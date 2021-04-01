
!> The SSASolver edited for use with MPM (GIMPM or sMPM). The essential difference
!! for the MPM version is that particles (material points) serve as the integration
!! points. Variables that update each iteration are mapped to material points in the
!! subroutine UpdateSSAParticleVals. Matrix assembly using the particles occurs in
!! subroutine LocalMatrixUVSSAMPM. Optionally, matrix assembly for some parts of the domain may
!! be performed using the usual FEM routines (LocalMatrixUVSSAFEM, and when using damage,
!! LocalMatrixUVSSAFEMDamage) when MPM isn't necessary for the entire ice domain. The FEM subroutines
!! are also called for bulk element matrix assembly at the ice front (see Huth et al 2020, Part I).
!! The usual FEM routine (LocalMatrixBCSSA) is used to enforce the ice front boundary condition.

!Alex Huth, 2020
!ahuth@princeton.edu


SUBROUTINE MPM_SSA( Model,Solver,dt,TransientSimulation )
  !--------------------------------------------------------------------------
  !**************************************************************************
  !
  !  Solve the 1D or 2D horizontal velocity with the SSA formulated for MPM
  !  NOTE: -Assumes density and gravity are constants
  !
  !  ARGUMENTS:
  !
  !  TYPE(Model_t) :: Model,
  !     INPUT: All model information (mesh, materials, BCs, etc...)
  !
  !  TYPE(Solver_t) :: Solver
  !     INPUT: Linear & nonlinear equation solver options
  !
  !  REAL(KIND=dp) :: dt,
  !     INPUT: Timestep size for time dependent simulations
  !
  !  LOGICAL :: TransientSimulation
  !     INPUT: Steady state or transient simulation
  !
  !**************************************************************************
  USE DefUtils
  USE SolverUtils
  USE MeshUtils
  USE Types
  USE MPMUtils

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Valuelist_t), POINTER :: Params

  TYPE(Nodes_t)   :: ElementNodes,Nodes
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, BoundaryElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material, BC
  TYPE(Variable_t), POINTER :: PointerToVariable,  VeloSol, Hsol, PVsol,ZsSol,PM,NSol,&
       MeshDamageVar

  LOGICAL :: AllocationsDone = .FALSE., Found, GotIt, CalvingFront, UnFoundFatal=.TRUE., Visited=.FALSE.
  LOGICAL :: Newton, RevisitMesh=.FALSE., firsttime = .FALSE.,&
       reducedam,musl,UseBfScale,DefGradIters

  INTEGER :: i, n, m, t, istat, DIM, p, STDOFs, iFriction, b, bb, bbb,xpicm,pind,countdam
  INTEGER :: NonlinearIter, NewtonIter, iter, other_body_id,np,newtonts


  INTEGER, POINTER :: Permutation(:), HPerm(:), NodeIndexes(:), PVPerm(:), &
       PMPerm(:),ZsPerm(:),NPerm(:),MeshDamagePerm(:)
  INTEGER, ALLOCATABLE :: LElem(:), RElem(:)
  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), H(:), Zb(:), Nval(:), PV(:), PMVal(:),Zs(:),MeshDamage(:)

  REAL(KIND=dp) :: UNorm, cn, dd, NonlinearTol, NewtonTol, MinSRInv, MinH, rhow, sealevel, usf, &
       PrevUNorm, relativeChange, minv, fm, PostPeak, MinN, gravity, rhoi,&
       CriticalDav,gridres,scale,stressscale,LinVelo,FricMaxVal,maxunorm,hf,Pval

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:),NodalBTrack(:),NodalBM(:), &
       NodalGravity(:), NodalViscosity(:), NodalDensity(:), NodalU(:),NodalV(:),NodalH(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, Friction
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, at0, starttime, endtime
#else
  REAL(KIND=dp) :: at, at0, starttime, endtime, CPUTime, RealTime
#endif
  INTEGER :: ei,j, ElemFirst, ElemLast, sz, IsAllocated, NoNewElements, count,No
  TYPE(Mesh_t), POINTER :: Mesh, OrigMesh => NULL()
  TYPE(Element_t), POINTER :: Telems(:), P1, P2, NewElements(:),OrigElems(:)
  TYPE(Particle_t), POINTER :: Particles
  !TYPE(MPforSSA_t), POINTER :: MP
  LOGICAL :: SEP,UseFEdit,applyzerostress
  LOGICAL :: UseZeroStressDamage,UseZeroStressDamageFirstOnly,shelfonly,reweightmpm,first,fpgroundonly,restarted

  LOGICAL :: UseFemMinMax,UseFem
  LOGICAL :: AdaptRelaxation,VisitedAdaptRelax,RelaxationAdapted
  REAL(KIND=dp) :: FemMinX,FemMaxX,SaveRelax,NewRelax,StartNewtonNorm,NormMultThres

  INTEGER :: ind, ind2, xpicinterval,zssurf,numberofparticles
  INTEGER :: GLnIP ! number of Integ. Points for GL Sub-element parametrization
  TYPE(Variable_t), POINTER :: GMSol,BedrockSol,Vstore
  INTEGER, POINTER :: VstorePerm(:)
  REAL(KIND=dp), POINTER :: VstoreVal(:)

  REAL(KIND=dp), ALLOCATABLE :: NodalZs(:), NodalZb(:),NodalBeta(:), NodalLinVelo(:),&
       NodalGM(:),NodalBed(:),NodalD(:,:)

  REAL(KIND=dp), POINTER :: height(:),depth(:)
  INTEGER :: ii,numlayers,NoP
  REAL(KIND=dp) :: oneovernumlayersminus1,melfrac
  REAL(KIND=dp) :: invvisc,btzav,dav,newdamval,zero=0.0_dp,one=1.0_dp,diff,db,ds
  REAL(KIND=dp) :: EigenVec(2,2),EigVals(2),strainrate(2,2),ww,xx,yy,zz,inc
  LOGICAL ::assigntolayers,convbefzs,usedcompress,usemelfrac,alwaysisodir
  INTEGER :: bnode,inode,nn,en,ni(2),jj
  REAL(KIND=dp) :: s1,s2,norm
  INTEGER, POINTER :: ENodeIndexes(:)
  TYPE(Variable_t), POINTER :: BTrack,mfvar,bm,opvar
  INTEGER, POINTER :: BTrackPerm(:),mfperm(:),bmperm(:),opperm(:)
  REAL(KIND=dp), POINTER :: BtrackVal(:),mf(:) ,bmval(:), op(:)

  REAL(KIND=dp) :: a_cm,a_secondsperyear,a_H0,a_v0,a_Q0,a_B0,a_A,a_C,a_EeExp,a_Acm,a_m1,a_m2

  SAVE rhow,sealevel, gravity, rhoi, gridres, NodalDensity,zssurf,convbefzs
  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName, ElementNodes
  SAVE NodalU, NodalV, NodalD, NodeIndexes, NodalH,CriticalDav,Visited, count
  SAVE Nodes,xpicm
  SAVE reweightmpm, NodalZs, NodalZb, NodalBeta,UseZeroStressDamage,&
       UseZeroStressDamageFirstOnly
  SAVE NodalLinVelo, NodalGM, NodalBed,NodalGravity,NodalViscosity,xpicinterval,fpgroundonly
  SAVE UseFemMinMax,FemMinX,FemMaxX,SaveRelax,AdaptRelaxation,NormMultThres,newtonts


  !****************
  !INITIALIZATION
  !****************
  restarted = .FALSE.
  starttime = RealTime()

  WRITE(SolverName, '(A)') 'MPM_SSA'



  IF (.NOT. Visited) THEN
     count = -1
     Visited = .TRUE.
  END IF

  count = count + 1


  Particles => GlobalParticles
  !MP => GlobalMPforSSA
  NoP = Particles % NumberOfParticles


  !> Each iteration, Update particle values needed to assemble stiffness and force matrices.
  !! The structure (type) MP is used for storage, which repurposes the memory otherwise
  !! allocated for the xpic routine (Particles % xpic) and damage layer increments
  !! (Particles % dD). For efficiency, a vectorized approach is used to update the
  !! particle values, whenever possible

  ! MP is a structure (type) used for storage when updating the particle values needed
  ! to assemble the stiffness and force matrices each iteration (see Subroutine
  ! UpdateSSAParticleVals(-----). It makes use of the
  ! memory already allocated for Particles % xpic and Particles % dD.

  Particles % xpic = 0.0_dp
  Particles % dD = 0.0_dp

  MP % DSRxx => Particles % dD(1:NoP,1,1)
  MP % DSRyy => Particles % dD(1:NoP,1,2)
  MP % DSRxy => Particles % dD(1:NoP,1,3)
  MP % eta   => Particles % dD(1:NoP,1,4)

  MP % muder      => Particles % dD(1:NoP,2,1)
  MP % slip       => Particles % dD(1:NoP,2,2)
  MP % driveforce => Particles % dD(1:NoP,2,3:4)

  MP % GradVel      => Particles % dD(:,3,1:4)
  MP % GridVelocity => Particles % dD(:,4,1:2)

  MP % Ezz => Particles % dD(1:NoP,4,3)
  MP % Exy => Particles % dD(1:NoP,4,4)

  MP % Hf => Particles % dD(1:NoP,5,1)
  MP % fN => Particles % dD(1:NoP,5,2)
  MP % fB => Particles % dD(1:NoP,5,3)
  MP % Ee => Particles % dD(1:NoP,5,4)

  MP % Dxx => Particles % dD(1:NoP,6,1)
  MP % Dyy => Particles % dD(1:NoP,6,2)
  MP % Dzz => Particles % dD(1:NoP,6,3)
  MP % Dxy => Particles % dD(1:NoP,6,4)

  MP % exxd1m1 => Particles % dD(1:NoP,7,1)
  MP % eyyd2m1 => Particles % dD(1:NoP,7,2)
  MP % ezzd3m1 => Particles % dD(1:NoP,7,3)
  MP % exyd4   => Particles % dD(1:NoP,7,4)

  !MP % falpha => Particles % xpic(1:NoP,1:2)  !not used in this version
  MP % slip2  => Particles % xpic(1:NoP,3)
  MP % ub     => Particles % xpic(1:NoP,4)
  MP % Velo   => Particles % xpic(1:NoP,5:6)


  ! Get some constants
  gridres = Particles % gridres
  rhow = Particles % rhow
  sealevel = Particles % sealevel
  CriticalDav = Particles % criticaldav
  MinH = 1.0_dp
  rhoi = Particles % rhoi
  Gravity = ABS(Particles % gravity)
  SEP = Particles % SEP

  VisitedAdaptRelax = .FALSE.
  RelaxationAdapted = .FALSE.

  IF (count == 0) THEN

     !---------------------------------------------------------------------------------!
     !everthing in this commented section was done in SSAPrepMeshToParticles!!!
     !so only need to do first time

     IF (.NOT. Particles % moveGL) THEN
        CALL Info('MPM_Initialize','update mesh zs',Level=4)
        CALL UpdateMeshZs( Particles,Model)
     END IF


     CALL Info(SolverName,'MPM interpolation of gradvel to particles...',Level=4)
     !3 is a dummy
     CALL MPMMeshVectorToParticle(Particles, Model, 1,3 )
     CALL Info(SolverName,'interpolation done',Level=4)

     !gradzs and falpha to particles
     CALL Info(SolverName,'MPM interpolation of gradzs to particles...',Level=4)
     CALL MPMMeshScalarToParticle(Particles, Model, 4)
     CALL Info(SolverName,'interpolation done',Level=4)


     IF (Particles % movegl) THEN
        !gmask, mb,ef,fp. Also bed if movegl
        CALL Info(SolverName,'MeshToParticles: Gmask',Level=4)
        CALL MPMMeshScalarToParticle( Particles, Model, 2)
     END IF

     IF (SEP .AND. Particles % movegl) THEN

        DO No = 1,NoP
           IF (Particles % Gmask(No) > 0.99_dp) THEN
              Particles % Gmask(No) = 1.0_dp
           ELSE IF (Particles % Gmask(No) < -0.99_dp) THEN
              Particles % Gmask(No) = -1.0_dp
           ELSE
              Hf = rhow * (sealevel-Particles % bedrock(No))/rhoi
              IF (Particles % H(No) .LT. Hf) THEN
                 Particles % Gmask(No) = 1.0_dp
              ELSE
                 Particles % Gmask(No) = -1.0_dp
              END IF
           END IF
        END DO
     END IF
  END IF



  Particles % NextCoordinate = 0.0_dp


  !**************************************************************************
  !--------------------------------------------------------------------------
  !                              Mesh Edits                                 !
  !                                                                         !
  ! set some bulk passive elements at the ice front to NULL so that ice     !
  ! front neumann boundary condition can be applied                         !
  !--------------------------------------------------------------------------

  CALL Info( SolverName, ' ', Level=4 )
  CALL Info( SolverName, 'Editing Mesh for Passive Neumann at P/A bounds', Level=4 )

  Mesh => Solver % Mesh

  ElemFirst = Mesh % NumberOfBulkElements+1
  ElemLast = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
  b = ElemLast-ElemFirst+1

  ALLOCATE(LElem(b))
  ALLOCATE(RElem(b))

  j = 0
  DO i = ElemFirst,ElemLast
     j = j+1
     Element => Mesh % Elements( i )
     Model % CurrentElement => Element

     P1 => Element % BoundaryInfo % Left
     P2 => Element % BoundaryInfo % Right

     IF ( ASSOCIATED(P1)) THEN
        LElem( j ) = P1 % ElementIndex
     ELSE
        LElem( j ) = -1
     END IF
     IF ( ASSOCIATED(P2)) THEN
        RElem( j ) = P2 % ElementIndex
     ELSE
        RElem(j) = -1
     END IF
  END DO

  j = 0
  DO i = ElemFirst,ElemLast
     j=j+1
     Element => Mesh % Elements( i )
     Model % CurrentElement => Element
     P1 => Element % BoundaryInfo % Left
     P2 => Element % BoundaryInfo % Right

     !THIS WAS <= before
     ! IF( i < ( ElemLast - Mesh % PassBCcnt )) THEN
     !    CYCLE
     ! ELSE
     IF ( CheckPassiveElement(P1) ) THEN
        NULLIFY ( Element % BoundaryInfo % Left )
     END IF

     IF ( CheckPassiveElement(P2) ) THEN
        NULLIFY ( Element % BoundaryInfo % Right )
     END IF
     ! END IF
  END DO

  CALL Info( SolverName, 'Done with Mesh Edits', Level=4 )
  CALL Info( SolverName, ' ', Level=4 )

  !**************************************************************************


  !--------------------------------------------------------------------------
  !                           Initialization
  !--------------------------------------------------------------------------

  Params => GetSolverParams()
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs


  DIM = CoordinateSystemDimension()

  ! IF DIM = STDOFs+1  Normal-Tangential can not be used => trick temporary set Model Dimension to STDOFs
  IF (DIM.eq.(STDOFs+1)) CurrentModel % Dimension = STDOFs

  HSol => VariableGet( Solver % Mesh % Variables, 'H',UnFoundFatal=UnFoundFatal)
  H => HSol % Values
  HPerm => HSol % Perm


  ZsSol => VariableGet( Solver % Mesh % Variables, 'Zs',UnFoundFatal=UnFoundFatal)
  Zs => ZsSol % Values
  ZsPerm => ZsSol % Perm


  ! Vplus is just used for storage purposes
  Vstore => VariableGet(Model % Mesh % Variables, 'Vplus' )
  IF (.NOT. ASSOCIATED(Vstore)) CALL Fatal(SolverName,'Vplus does not exist ')
  VstorePerm => Vstore % Perm
  VstoreVal => Vstore % Values

  PVSol => VariableGet( Solver % Mesh % Variables, 'PrevVel',UnFoundFatal=UnFoundFatal)
  PV => PVSol % Values
  PVPerm => PVSol % Perm

  PV(STDOFs*(PVPerm(:)-1)+1) = VariableValues(STDOFs*(Permutation(:)-1)+1)
  PV(STDOFs*(PVPerm(:)-1)+2) = VariableValues(STDOFs*(Permutation(:)-1)+2)

  IF (Particles % usedamage) THEN
     MeshDamageVar => VariableGet( Solver % Mesh % Variables, 'Mesh Damage',UnFoundFatal=UnFoundFatal)
     MeshDamage => MeshDamageVar % Values
     MeshDamagePerm => MeshDamageVar % Perm
  END IF

  PM => VariableGet( Mesh % Variables, 'surface')
  PMVal => PM % Values
  PMPerm => PM % Perm

  IF (SEP) THEN
     GLnIP=ListGetInteger( Solver % Values, &
          'GL integration points number',UnFoundFatal=.TRUE. )
  END IF

  GMSol => VariableGet( Solver % Mesh % Variables, 'Mask',UnFoundFatal=.TRUE. )
  BedrockSol => VariableGet( Solver % Mesh % Variables, 'bed',UnFoundFatal=.TRUE. )


  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN


     !Portions of the domain can be processed using FEM. Largely used for
     !testing purposes at this point, but with Use FEM MinMax = true,
     !the bounds in the x direction in which FEM is to be used can be
     !specified with FEM Min x an FEM Max x
     UseFemMinMax =  GetLogical( Params, 'Use FEM MinMax', GotIt)
     IF (.NOT. GotIt) THEN
        UseFemMinMax = .FALSE.
        CALL WARN(Solvername,'UseFemMinMax not specified -- assuming false')
     END IF

     IF (UseFemMinMax) THEN
        FemMinX =  GetCReal( Params, 'FEM Min x', GotIt)
        IF (.NOT. GotIt) THEN
           CALL FATAL(Solvername,'Must Define FEM Min x and FEM Max x if using FEMMinMax')
        END IF

        FemMaxX =  GetCReal( Params, 'FEM Max x', GotIt)
        IF (.NOT. GotIt) THEN
           CALL FATAL(Solvername,'Must Define FEM Min x and FEM Max x if using FEMMinMax')
        END IF
     END IF



     !is this useful?
     fpgroundonly =  GetLogical( Params, 'Friction only on grounded', GotIt)
     IF (.NOT. GotIt) THEN
        fpgroundonly = .FALSE.
        CALL WARN(Solvername,'Friction on grounded not specified -- assuming false')
     END IF



     !the reweighting scheme for MPM
     reweightmpm = GetLogical(Params, 'Particle Reweighting',GotIt)
     IF (.NOT. GotIt) THEN
        IF (Particles % ShapeFunctions == 'gimpm') THEN
           reweightmpm = .FALSE.
           CALL WARN(Solvername,'Particle Reweighting not specified -- assuming false')
        ELSE
           reweightmpm = .TRUE.
           CALL WARN(Solvername,'Particle Reweighting not specified -- assuming true')
        END IF
     END IF


     !interval in which xpic is used to update particle velocities/positions.
     !Otherwise, FLIP is used
     xpicinterval = GetInteger(Params, 'xpic interval',GotIt)
     IF (.NOT. GotIt) THEN
        xpicinterval = 1
        CALL WARN(Solvername,'xpic interval not specified. Assuming 1.')
     END IF


     UseZeroStressDamage =  GetLogical( Params, 'Use Zero Stress Damage', GotIt)
     IF (.NOT. GotIt) THEN
        UseZeroStressDamage = .FALSE.
        CALL WARN(Solvername,'UseZeroStressDamage not specified -- assuming false')
     END IF


     !experimental: for initializing a 2-D damage field from zero-stress damage on the first timestep
     UseZeroStressDamageFirstOnly =  GetLogical( Params, 'Use Zero Stress Damage First Timestep Only', GotIt)
     IF (.NOT. GotIt) THEN
        UseZeroStressDamageFirstOnly = .FALSE.
        CALL WARN(Solvername,'UseZeroStressDamageFirstOnly not specified -- assuming false')
     END IF

     !experiemntal: for turning the 2-D damage field from UseZeroStressDamageFirstOnly into a 3-D field
     !for later creep damage evolution
     assigntolayers = GetLogical( Params, 'Assign Zero Stress Damage To Layers', GotIt)
     IF (.NOT. GotIt) THEN
        assigntolayers = .FALSE.
        CALL WARN(Solvername,'Assign Zero Stress Damage To Layers not specified -- assuming false')
     END IF


     zssurf = 0
     IF (UseZeroStressDamage) THEN

        !Experimental: If using the (implicit) zero-stress damage model,
        !better convergence may be possible if the SSA solution is allowed to temporarily converge
        !without the zero-stress damage evolution first, and then turning on the zero-stress model
        convbefzs =  GetLogical( Params, 'Converge before zs', GotIt)
        IF (.NOT. GotIt) THEN
           convbefzs = .FALSE.
           CALL WARN(Solvername,'converge before zsn not specified -- assuming false')
        END IF

        zssurf = GetInteger(Params,'zero stress surface',GotIt)
        IF (.NOT. GotIt) THEN
           CALL WARN(Solvername,'zero stress surface = Integer $ not specified, assuming 0 (both surfs)')
           zssurf = 0
        ELSE
           IF (zssurf == 1) THEN
              CALL INFO(Solvername,'using zero stress damage for surface crevasses only!',Level=1)
           ELSEIF (zssurf == -1) THEN
              CALL INFO(Solvername,'using zero stress damage for basal crevasses only!',Level=1)
           ELSE
              CALL INFO(Solvername,'using zero stress damage for both surface and basal crevasses!',Level=1)
           END IF
        END IF
     END IF

     SaveRelax = GetCReal( Params, &
          'Nonlinear System Relaxation Factor',GotIt )
     IF ( .NOT. GotIt ) SaveRelax = 1.0_dp

     AdaptRelaxation =  GetLogical( Params, 'Use Adaptive Relaxation', GotIt)
     IF (.NOT. GotIt) THEN
        AdaptRelaxation = .FALSE.
        CALL Info(Solvername,'Use Adaptive Relaxation not specified -- assuming false',Level=1)
     END IF

     IF (AdaptRelaxation) THEN
        NormMultThres =  GetCReal( Params, 'Adaptive Norm Mult Threshold', GotIt)
        IF (.NOT. GotIt) THEN
           NormMultThres = 1.5_dp
        END IF
     END IF


     ! Allocate
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes


     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, &
          NodalH, NodalGravity,NodalViscosity,& !NodalBTrack,NodalBM,&
          NodalZb, NodalZs, NodalBeta, NodalLinVelo,  NodalGM,&
          NodalBed,NodalU, NodalV, NodalDensity, ElementNodes % x, &
          ElementNodes % y, ElementNodes % z,&
          Nodes % x, Nodes % y, Nodes % z )

     IF (Particles % usedamage .AND. AllocationsDone) THEN
        IF (ALLOCATED(NodalD)) DEALLOCATE(NodalD)
     END IF


     ALLOCATE( FORCE(STDOFs*N), LOAD(N), STIFF(STDOFs*N,STDOFs*N), &
          NodalH(N), NodalGravity(N), NodalViscosity(N), & !NodalBTrack(N),NodalBM(N),&
          NodalZb(N), NodalZs(N),NodalBeta(N), NodalLinVelo(N),  &
          NodalGM(N),NodalBed(N),NodalU(N), NodalV(N), NodalDensity(N), &
          ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N), &
          Nodes % x(N), Nodes % y(N), Nodes % z(N),STAT=istat )

     IF (Particles % usedamage) THEN
        ALLOCATE(NodalD(N,4))
     END IF

     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  NonlinearTol = GetConstReal( Solver % Values, &
       'Nonlinear System Convergence Tolerance' )

  NonlinearIter = GetInteger( Solver % Values, &
       'Nonlinear System Max Iterations',GotIt )

  IF ( .NOT.GotIt ) NonlinearIter = 1

  NewtonTol = ListGetConstReal( Solver % Values, &
       'Nonlinear System Newton After Tolerance', minv=0.0d0 )

  NewtonIter = ListGetInteger( Solver % Values, &
       'Nonlinear System Newton After Iterations', GotIt )
  if (.NOT.Gotit) NewtonIter = NonlinearIter + 1

  MaxUNorm = ListGetConstReal( Solver % Values, &
       'Max Norm', GotIt )
  IF (.NOT. GotIt) MaxUNorm = 1.0E99_dp

  Newton=.False.



  !order of xpic
  xpicm = ListGetInteger( Solver % Values, &
       'xpic m', GotIt )
  IF (.NOT.Gotit) THEN
     xpicm = 0
     CALL WARN(Solvername,'xpic not specified, assuming 0 (using FLIP)')
  END IF

  NewtonTS = ListGetInteger( Solver % Values, &
       'Timestep to start newton',GotIt)
  IF (.NOT. GotIt) NewtonTS = 0

  IF (count < NewtonTS) THEN
     NewtonIter = 1000
     NewtonTol = 0.0_dp
  END IF

  ! Read the Viscosity exponent m in MMaterial Section
  ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n)
  cn = GetConstReal( Model % Constants, 'Viscosity Exponent', GotIt )
  IF (.NOT. GotIt) CALL Fatal(Solvername,'Need to define "Viscosity Exponent = Real $1/n" in constants')
  MinSRInv = GetConstReal( Model % Constants, 'Critical Shear Rate', GotIt )
  IF (.NOT. GotIt) CALL Fatal(Solvername,'Need to define "Critical Shear Rate = Real $1crit" in constants')

  !Friction
  DO t = 1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     EXIT
  END DO

  Material => GetMaterial(Element)
  Friction = GetString(Material, 'SSA Friction Law', Found)
  IF (.NOT.Found) CALL FATAL(SolverName,'Could not find Material keyword >SSA Friction Law<')
  SELECT CASE(Friction)
  CASE('linear')
     iFriction = 1
     fm = 1.0_dp
  CASE('weertman')
     iFriction = 2
  CASE('coulomb')
     iFriction = 3
  CASE DEFAULT
     CALL FATAL(SolverName,'Friction should be linear, Weertman or Coulomb')
  END SELECT

  ! for Weertman and Coulomb friction
  IF (iFriction > 1) THEN
     fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found , UnFoundFatal=UnFoundFatal)
     MinN = ListGetConstReal( Material, 'SSA Min Effective Pressure', Found, UnFoundFatal=UnFoundFatal)
     !Previous default value: MinN = 1.0e-6_dp
     LinVelo = 0.0_dp
     LinVelo = ListGetConstReal( Material, 'SSA Friction Linear Velocity', Found , UnFoundFatal=UnFoundFatal)
  END IF
  ! only for Coulomb friction
  IF (iFriction > 2) THEN
     PostPeak = ListGetConstReal( Material, 'SSA Friction Post-Peak', Found, UnFoundFatal=UnFoundFatal )
     FricMaxVal = ListGetConstReal( Material, 'SSA Friction Maximum Value', Found , UnFoundFatal=UnFoundFatal)
  END IF
  !End Friction

  !--------------------------------------------------------------------------
  !                           End Initialization
  !--------------------------------------------------------------------------


  !**************************************************************************
  !                             START Iterations
  !**************************************************************************

  IF (UseZeroStressDamageFirstOnly .AND. count > 0) UseZeroStressDamage = .FALSE.

  applyzerostress = .FALSE.

  IF (UseZeroStressDamage) THEN
     IF (count > 0 .AND. (.NOT. convbefzs) ) applyzerostress = .TRUE.
     Particles % damage(:,3,3) = 0.0_dp
     Particles % damage(:,2,1:4) = Particles % Dav(:,1:4)
  END IF



  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS


100 Newton = .FALSE.



  DO iter=1,NonlinearIter

300  Newton = Newton

     at  = CPUTime()
     at0 = RealTime()

     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, &
          '-------------------------------------',Level=4 )
     WRITE( Message, * ) 'SSA BASAL VELOCITY NON-LINEAR ITERATION', iter
     CALL Info( SolverName, Message, Level=4 )
     If (Newton) Then
        WRITE( Message, * ) 'Newton linearisation is used'
        CALL Info( SolverName, Message, Level=4 )
     Endif
     CALL Info( SolverName, ' ', Level=4 )
     CALL Info( SolverName, &
          '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' ', Level=4 )

     !Initialize the system and do the assembly:
     !------------------------------------------
     CALL DefaultInitialize()



     !THIS IS WHERE VARIABLES ARE NOW UPDATED ON EACH PARTICLE:
     CALL UpdateSSAParticleVals( NoP, Gravity,rhoi,rhow,&
          iFriction,fm,cn,MinSRInv,minH,model,&
          minN,LinVelo,PostPeak,FricMaxVal,sealevel,&
          newton,usezerostressdamage,count,applyzerostress,UseFemMinMax,&
          FemMinX,FemMaxX,gridres,zssurf)



     !-------------------------------------------------------------------------------
     !                         Bulk Assembly
     !-------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t)
        IF (ParEnv % myPe .NE. Element % partIndex) CYCLE

        n = GetElementNOFNodes()

        NodeIndexes => Element % NodeIndexes

        ei = Element % ElementIndex


        ! set coords of highest occuring dimension to zero (to get correct path element)
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (STDOFs == 1) THEN !1D SSA
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (STDOFs == 2) THEN !2D SSA
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                STDOFs, ' . Aborting'
           CALL Fatal( SolverName, Message)
           STOP
        END IF

        !Previous Velocity
        NodalU(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+1)
        NodalV = 0.0_dp
        IF (STDOFs.EQ.2) NodalV(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+2)

        Material => GetMaterial(Element)

        IF (UseFemMinMax) THEN
           IF (ALL(ElementNodes % x > FemMinX) .AND. ALL(ElementNodes % x < FemMaxX)) THEN
              UseFem = .TRUE.
           ELSE
              UseFem = .FALSE.
           END IF
        ELSE
           UseFem = .FALSE.
        END IF

        !ei = t
        !Element % ElementIndex

        !----------------------------------- FEM STUFF -------------------------------------!
        IF (ElemTrack(ei) % Status == FEM .OR. UseFem) THEN

           Friction = GetString(Material, 'SSA Friction Law', Found)
           IF (.NOT.Found) &
                CALL FATAL(SolverName,'Could not find Material keyword >SSA Friction Law<')
           SELECT CASE(Friction)
           CASE('linear')
              iFriction = 1
              fm = 1.0_dp
           CASE('weertman')
              iFriction = 2
           CASE('coulomb')
              iFriction = 3
           CASE DEFAULT
              CALL FATAL(SolverName,'Friction should be linear, Weertman or Coulomb')
           END SELECT
           ! for Weertman and Coulomb friction
           IF (iFriction > 1) THEN
              fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found , UnFoundFatal=UnFoundFatal)
              MinN = ListGetConstReal( Material, 'SSA Min Effective Pressure', Found, UnFoundFatal=UnFoundFatal)
              !Previous default value: MinN = 1.0e-6_dp
              LinVelo = 0.0_dp
              LinVelo = ListGetConstReal( Material, 'SSA Friction Linear Velocity', Found , UnFoundFatal=UnFoundFatal)
           END IF
           ! only for Coulomb friction
           IF (iFriction > 2) THEN
              PostPeak = ListGetConstReal( Material, 'SSA Friction Post-Peak', Found, UnFoundFatal=UnFoundFatal )
              FricMaxVal = ListGetConstReal( Material, 'SSA Friction Maximum Value', Found , UnFoundFatal=UnFoundFatal)
           END IF

           IF (Particles % usedamage) THEN
              NodalD(1:n,1) = MeshDamage(4*(Permutation(NodeIndexes(1:n))-1)+1)
              NodalD(1:n,2) = MeshDamage(4*(Permutation(NodeIndexes(1:n))-1)+2)
              NodalD(1:n,3) = MeshDamage(4*(Permutation(NodeIndexes(1:n))-1)+3)
              NodalD(1:n,4) = MeshDamage(4*(Permutation(NodeIndexes(1:n))-1)+4)
           END IF

           ! Get the Nodal value of Zb and Zs
           NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))
           NodalH(1:n) = H(HPerm(NodeIndexes(1:n)))
           NodalZb(1:n) = NodalZs(1:n)-NodalH(1:n)
           NodalGravity(1:n) = Gravity
           NodalDensity(1:n) = rhoi
           NodalViscosity=0.0_dp
           NodalViscosity(1:n) = ListGetReal( Material, 'SSA Mean Viscosity',n, NodeIndexes,Found,&
                UnFoundFatal=UnFoundFatal)
           NodalBeta = 0.0_dp
           NodalBeta(1:n) = ListGetReal( &
                Material, 'SSA Friction Parameter', n, NodeIndexes(1:n), Found,&
                UnFoundFatal=UnFoundFatal)
           NodalLinVelo(1:n) = LinVelo
           IF (SEP) THEN
              NodalGM(1:n)=GMSol%Values(GMSol%Perm(NodeIndexes(1:n)))
              NodalBed(1:n)=BedrockSol%Values(BedrockSol%Perm(NodeIndexes(1:n)))
           ENDIF

           IF (ElemTrack(ei) % Status == FEM .OR. UseFem) THEN

              IF (Particles % usedamage) THEN


                 CALL LocalMatrixUVSSAFEMDamage (  STIFF, FORCE, Element, n, ElementNodes, NodalGravity, &
                      NodalDensity, NodalViscosity, NodalZb, NodalZs, NodalU, NodalV, &
                      NodalD,iFriction, NodalBeta, fm, &
                      NodalLinVelo, PostPeak, fricmaxval, NodalGM,NodalBed,SEP,rhow, &
                      cn, MinSRInv, MinH , STDOFs, Newton,minn,reweightmpm)
              ELSE
                 CALL LocalMatrixUVSSAFEM (  STIFF, FORCE, Element, n, ElementNodes, NodalGravity, &
                      NodalDensity, NodalViscosity, NodalZb, NodalZs, NodalU, NodalV, &
                      iFriction, NodalBeta, fm, NodalLinVelo, PostPeak, fricmaxval, &
                      NodalGM,NodalBed,SEP,rhow, &
                      cn, MinSRInv, MinH , STDOFs, Newton,minn)
              END IF
           END IF
           !----------------------------------- END FEM STUFF -------------------------------------!
        END IF


        IF  ( (ElemTrack(ei) % Status .NE. FEM) .AND. (.NOT. UseFem) ) THEN

           ! NodalH(1:n) = H(HPerm(NodeIndexes(1:n)))

           CALL LocalMatrixUVSSAMPM (  STIFF, FORCE, Element, n, ElementNodes, Gravity, &
                rhoi, rhow, NodalU, NodalV, iFriction, fm,cn, &
                STDOFs, Newton,model,PointerToVariable,HSol,reweightmpm,&
                count,usezerostressdamage)

        END IF

        !------------------------------------------------------------------------------

        CALL DefaultUpdateEquations( STIFF, FORCE )

        !-------------------------------------------------------------------------------
     END DO


     CALL DefaultFinishBulkAssembly()


     !-------------------------------------------------------------------------------
     !                         Neumann Boundary
     !-------------------------------------------------------------------------------
     DO t=1,GetNOFBoundaryElements()
        BoundaryElement => GetBoundaryElement(t)

        IF (STDOFS.NE.1) then
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE
        END IF
        IF ( GetElementFamily() == 1 ) CYCLE

        NodeIndexes => BoundaryElement % NodeIndexes

        !! The passive-active BC is put in the -1 partition, so turned off the following:
        ! IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE
        ! IF (BoundaryElement % partIndex > 0) CYCLE
        !! will need to fix this for parallel implementation

        n = GetElementNOFNodes()


        FORCE = 0.0_dp
        STIFF = 0.0_dp

        ! set coords of highest occuring dimension to zero (to get correct path element)
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (STDOFs == 1) THEN
           ElementNodes % y(1:n) = 0.0_dp
           ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (STDOFs == 2) THEN
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA with SSA var DOFs=',&
                STDOFs, '. Aborting'
           CALL Fatal( SolverName, Message)
           STOP
        END IF

        BC => GetBC()
        IF (.NOT.ASSOCIATED( BC ) ) CYCLE

        CalvingFront=.False.
        CalvingFront = ListGetLogical( BC, 'Calving Front', GotIt )


        P1 => BoundaryElement % BoundaryInfo % Left
        P2 => BoundaryElement % BoundaryInfo % Right

        pind = 0


        IF (ASSOCIATED(P1)) THEN
           pind = P1 % ElementIndex
        END IF
        IF (pind<1) THEN
           IF (ASSOCIATED(P2)) THEN
              pind = P2 % ElementIndex
           END IF
           IF (pind < 1) THEN
              CalvingFront = .False.
           END IF
        END IF


        IF (CalvingFront) THEN

           NodalH(1:n) = H(HPerm(NodeIndexes(1:n)))
           MinH = 1.0_dp
           NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))
           NodalZb(1:n) = NodalZs(1:n)-NodalH(1:n)

           CALL LocalMatrixBCSSA(  STIFF, FORCE, BoundaryElement, n, ElementNodes, rhoi, &
                gravity, NodalZb, NodalZs, rhow, sealevel,MinH )
        END IF

        CALL DefaultUpdateEquations( STIFF, FORCE )
        !------------------------------------------------------------------------------
     END DO

     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()

     !------------------------------------------------------------------------------
     !      Solve the system, update particle GradVel, check for convergence
     !------------------------------------------------------------------------------
     PrevUNorm = UNorm
     UNorm = DefaultSolve()


     IF (Unorm /= Unorm) THEN
        CALL WARN(Solvername,'nan detected, saving particles and exiting')
        CALL ParticleOutputVtu( Particles,Model )
        CALL FATAL(Solvername,'Exited')
     END IF


     !---
     !If solution is clearly not converging:
     !try setting velocity to 0, and restart solver with picard iters only
     IF (UNorm > MaxUNorm) THEN

        IF (restarted) THEN
           CALL Warn(Solvername,&
                'will not converge after restart, saving particles and exiting')
           CALL ParticleOutputVtu( Particles,Model )
           CALL FATAL(Solvername,'Exited')
        END IF

        CALL Warn(Solvername,&
             'SOLUTION IS DIVERGING, RESTARTING WITH VELOCITY=0')

        VariableValues = 0.0_dp
        restarted = .TRUE.

        CALL Info(Solvername,'interpolation of gradvel and gridvel to particles',Level=4)
        CALL MPMMeshVectorToParticle(Particles, Model, 1, 3 )
        CALL Info(Solvername,'interpolation done',Level=4)

        NewtonTol = 0.0_dp
        NewtonIter = NonlinearIter + 1
        Newton = .FALSE.

        IF (usezerostressdamage) THEN
           applyzerostress = .FALSE.
           Particles % Dav(:,1:4) = Particles % damage(:,2,1:4)

           Call Info(Solvername,&
                'Resetting zerostress damage. Will allow to converge before reapplying')

           NewRelax = SaveRelax*0.85_dp
           CALL ListAddConstReal( Solver % Values,  &
                'Nonlinear System Relaxation Factor', NewRelax )

           WRITE( Message, * ) 'New Relaxation Factor : ', NewRelax
           CALL Info(SolverName, Message, Level=1 )
           RelaxationAdapted = .TRUE.
        END IF
        GOTO 100
     END IF
     !----

     RelativeChange = Solver % Variable % NonlinChange

     WRITE( Message, * ) 'Result Norm   : ', UNorm, PrevUNorm
     CALL Info(SolverName, Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ', RelativeChange
     CALL Info(SolverName, Message, Level=4 )


     CALL Info(Solvername,&
          'interpolation of gradvel and gridvel to particles',Level=4)
     IF (count>0 .AND. Particles % flipvelfric) THEN
        !-1 does it flip style
        CALL MPMMeshVectorToParticle(Particles, Model, 1, -1 )
        CALL Info(Solvername,'interpolation done',Level=4)
     ELSE
        CALL Info(Solvername,&
             'interpolation of gradvel and gridvel to particles',Level=4)
        CALL MPMMeshVectorToParticle(Particles, Model, 1, 3 )
        CALL Info(Solvername,'interpolation done',Level=4)
     END IF


     !If using adaptrelaxation, and newton has been activated, and the
     !solution has diverged 2x its norm when newton was first activated,
     !then revert to the solution when newton was first activated.
     !Halve the relaxation factor.  Hopefully, now it will converge.
     !Otherwise, will try to convert to purely picard iterations
     !until convergence...
     IF (VisitedAdaptRelax .AND. UNorm> StartNewtonNorm*NormMultThres) THEN

        !Here, vstoreval is  used for storage
        VariableValues(STDOFs*(Permutation(:)-1)+1) = VstoreVal(2*(VstorePerm(:)-1)+1)
        VariableValues(STDOFs*(Permutation(:)-1)+2) = VstoreVal(2*(VstorePerm(:)-1)+2)

        Particles % GradVel(:,1:4) = MP % GradVel(:,1:4)
        Particles % GridVelocity(:,1:2) = MP % GridVelocity(:,1:2)

        IF (.NOT. RelaxationAdapted) THEN
           CALL Info(Solvername,'SOL DIVERGING, HALVING RELAXATION AND RESTARTING NEWTON ITERS',Level=1)
           NewRelax = SaveRelax*0.5_dp
           CALL ListAddConstReal( Solver % Values,  &
                'Nonlinear System Relaxation Factor', NewRelax )

           WRITE( Message, * ) 'New Relaxation Factor : ', NewRelax
           CALL Info(SolverName, Message, Level=1 )

           RelaxationAdapted = .TRUE.
           Newton = .TRUE.
           GOTO 300
        ELSE

           CALL Info(Solvername,'SOL DIVERGING AGAIN, RESTARTING FROM VEL = 0',Level=1)
           VariableValues = 0.0_dp
           restarted = .TRUE.

           CALL Info(Solvername,'interpolation of gradvel and gridvel to particles',Level=4)
           CALL MPMMeshVectorToParticle(Particles, Model, 1, 3 )
           CALL Info(Solvername,'interpolation done',Level=4)

           VisitedAdaptRelax = .FALSE.
           NewtonTol = 0.0_dp
           NewtonIter = NonlinearIter + 1
           Newton = .FALSE.
           GOTO 100
        END IF
     END IF


     !START NEWTON!!
     IF ( RelativeChange < NewtonTol .OR. &
          iter > NewtonIter ) THEN
        Newton = .TRUE.

        !if newton is being activated for the first time and
        !using adapt relaxation, save current velocity values on the
        !mesh and particles.  And gradvel on particles.
        IF (AdaptRelaxation .AND. (.NOT. VisitedAdaptRelax) ) THEN

           MP % GradVel(:,1:4) = Particles % GradVel(:,1:4)
           MP % GridVelocity(:,1:2) = Particles % GridVelocity(:,1:2)

           !Vstore just used for storage
           VstoreVal(2*(VstorePerm(:)-1)+1) = VariableValues(STDOFs*(Permutation(:)-1)+1)
           VstoreVal(2*(VstorePerm(:)-1)+2) = VariableValues(STDOFs*(Permutation(:)-1)+2)

           StartNewtonNorm = UNorm
           VisitedAdaptRelax = .TRUE.
        END IF
     END IF

     IF ( RelativeChange < NonLinearTol ) THEN
        IF (iter > 1) THEN

           !reset relaxation if it was adapted
           IF (RelaxationAdapted) THEN
              CALL ListAddConstReal( Solver % Values,  &
                   'Nonlinear System Relaxation Factor', SaveRelax )
           END IF

           EXIT
        END IF
     END IF

  END DO

  !solution converged and iters are complete, unless applying zero-stress after
  !initial convergence:
  IF (usezerostressdamage) THEN

     IF (.NOT. applyzerostress) THEN
        applyzerostress = .TRUE.
        CALL INFO(Solvername,&
             'INITIAL SSA CONVERGENCE ACHIEVED. NOW APPLY ZERO STRESS DAMAGE', Level = 1)
        GOTO 100
     END IF

     Particles % zsmaxdd = MAXVAL(Particles % damage(:,3,2)-Particles % damage(:,3,1))

     WRITE( Message, * ) 'Maximum zero stress damage increase  : ', Particles % zsmaxdd
     CALL Info(SolverName, Message, Level=1 )

  END IF

  !----------------------------------------------------------------------------!
  !  SOLUTION CONVERGED: END LOOP NON-LINEAR ITERATIONS
  !----------------------------------------------------------------------------!

  CALL Info(Solvername,'MPM interpolation of converged grid velocity to particles...',Level=4)
  CALL MPMMeshVectorToParticle(Particles, Model, 3, count)
  CALL Info(Solvername,'interpolation done',Level=4)

  Particles % dd = 0.0_dp
  !Particles % xpic = 0.0_dp

  !--------some additional options for zero-stress damage simulations-------
  IF (UseZeroStressDamageFirstOnly .AND. count == 0) THEN
     DO No = 1,NoP
        IF (Particles % damage(No,3,2) >= Particles % CriticalDav) Particles % Damstatus(No) = 1
        IF (Particles % Damstatus(No) == 1) THEN
           Particles % Dav(No,1:3) = Particles % riftdmax
           Particles % Dav(No,4) = 0.0_dp
        END IF
     END DO
  END IF

  !---- A routine to assign a 3-D damage distribution based on 2-D zero-stress damage
  IF (UseZeroStressDamageFirstOnly  .AND. count == 0 .AND. assigntolayers) THEN

     numlayers = Particles % numberofparticlelayers
     oneovernumlayersminus1 = 1.0_dp/DBLE(numlayers-1)
     ALLOCATE(height(numlayers))
     ALLOCATE(depth(numlayers))

     DO No = 1,NoP
        IF (Particles % damstatus(No) == 1) THEN
           Particles % Dav(No,1:3) = Particles % riftdmax
           Particles % Dav(No,4) = 0.0_dp
           Particles % damage(No,:,1:3) = Particles % riftdmax
           Particles % damage(No,:,4) = 0.0_dp
           CYCLE
        END IF

        ds = Particles % Damage(No,1,1)
        db = Particles % Damage(No,1,2)
        Particles % Damage(No,:,:) = 0.0_dp

        inc = oneovernumlayersminus1*Particles % H(No)

        DO ii = 1,numlayers
           height(ii) = (DBLE(ii-1))*inc
           depth(ii) = Particles % H(No)-height(ii)
        END DO

        !basal crev
        DO ii = 1,numlayers

           diff = (height(ii) - db)/inc
           IF (ii == 1 .AND. diff == 0.0_dp) EXIT
           IF (diff < 0.0_dp) THEN
              Particles % damage(No,ii,1) = Particles % DmaxI !1.0
           ELSE
              IF (diff <= 0.5_dp) THEN
                 Particles % damage(No,ii,1) = 0.5_dp-diff
              ELSE
                 IF (ii-1 .NE. 1) THEN
                    Particles % damage(No,ii-1,1) = 1.5_dp-diff
                 ELSE
                    Particles % damage(No,ii-1,1) = 2.0_dp*(1.0_dp-diff)
                 END IF
              END IF

              EXIT
           END IF
        END DO

        !surface crevs
        DO ii = numlayers,1,-1
           diff = (depth(ii) - ds)/inc
           IF (ii == numlayers .AND. diff == 0.0_dp) EXIT
           IF (diff < 0.0_dp) THEN
              Particles % damage(No,ii,1) = Particles % DmaxI
           ELSE
              IF (diff <= 0.5_dp) THEN
                 Particles % damage(No,ii,1) = Particles % damage(No,ii,1) + 0.5_dp-diff
              ELSE
                 IF (ii+1 .NE. numlayers) THEN
                    Particles % damage(No,ii+1,1) = (Particles % damage(No,ii+1,1)-1.0_dp) + 1.5_dp-diff
                 ELSE
                    Particles % damage(No,ii+1,1) = (Particles % damage(No,ii+1,1)-1.0_dp) + 2.0_dp*(1.0_dp-diff)
                 END IF
              END IF

              EXIT

           END IF
        END DO

        DO ii = 2,3
           Particles % Damage(No,:,ii) = Particles % Damage(No,:,1)
        END DO

        IF (Particles % gamma > 0 .AND. Particles % damstatus(No) .NE. 1) THEN
           strainrate(1,1) = Particles % Dav(No,1)
           strainrate(2,2) = Particles % Dav(No,2)
           strainrate(1,2) = Particles % Dav(No,4)
           strainrate(2,1) = strainrate(1,2)
           CALL Eigen2DSym_TryGenFirst(strainrate,EigVals,EigenVec)

           DO ii = 1,Particles % numberofparticlelayers
              newdamval = Particles % damage(No,ii,1)

              EigVals(1) = newdamval * (1.0_dp-Particles % Gamma)
              EigVals(2) = newdamval

              ww = EigVals(1)*EigenVec(1,1)
              xx = EigVals(2)*EigenVec(1,2)
              yy = EigVals(1)*EigenVec(2,1)
              zz = EigVals(2)*EigenVec(2,2)

              Particles % damage(No,ii,1) = EigenVec(1,1)*ww + EigenVec(1,2)*xx
              Particles % damage(No,ii,2) = EigenVec(2,1)*yy + EigenVec(2,2)*zz
              Particles % damage(No,ii,3) = EigVals(1)
              Particles % damage(No,ii,4) = EigenVec(2,1)*ww + EigenVec(2,2)*xx
           END DO
        END IF
        CALL VertIntDamFromVisc(Particles, No, numlayers, Model)
     END DO

     DEALLOCATE(depth,height)
  END IF
  !---- end routine to assign a 3-D damage distribution based on 2-D zero-stress damage -----


  !------------------------------------------------------------------------------
  !         Solution Converged: now do FLIP,PIC, or XPIC
  !------------------------------------------------------------------------------

  IF (count == 0 .AND. Particles % firsttimestepzero) THEN
     !if dt is zero on the first timestep (useful for initialization), use PIC
     CALL Info(Solvername,'PIC particle update on first timestep',Level=3)

     Particles % Velocity = Particles % GridVelocity
     Particles % NextCoordinate = Particles % GridVelocity

     IF (Particles % analytictest) THEN
        ! Pval = rhoi * (rhow-rhoi)/rhow
        ! Particles % Velocity(:,1) = ((Pval * Gravity * Particles % H(:)/ &
        !      (4.0_dp * Particles % Binit(:)))**(3.0_dp)) * Particles % Coordinate(:,1)
        ! Particles % Velocity(:,2) = 0.0_dp
        ! Particles % GridVelocity = Particles % Velocity

        !Analytical 1D solution on first timestep...

        PRINT *,'ANALYTICAL 1D SOLUTION ON FIRST TIMESTEP'

        a_H0 = GetConstReal( Model % Constants,'H0',GotIt )
        IF (.NOT. GotIt) CALL Fatal('USF_1dtest:', &
             'initH: Need to define "H0 = Real $" in constants')

        a_v0 = GetConstReal( Model % Constants,'v0',GotIt )
        IF (.NOT. GotIt) CALL Fatal('USF_1dtest:', &
             'initH: Need to define "H0 = Real $" in constants')

        a_cm = 1.0_dp/3.0_dp
        a_secondsperyear = 31556926.0_dp
        !a_H0 = 600.0_dp
        !a_v0 = 300.0_dp
        a_Q0 = a_H0*a_v0
        a_B0 = 1.9E8_dp
        a_A = ((a_B0*1.0E-6_dp)**(-3.0_dp))*a_secondsperyear !Mpa^(-3) a^(-1)
        a_C = (((910.0_dp*1.0e-6_dp*9.81_dp)/&
             (4.0_dp*(a_A**(-a_cm))))*(1.0_dp-910.0_dp/1028.0_dp))**3.0_dp
        !C is the weertman constant !C =2.45E-18; !m?3 s?1

        a_EeExp = (a_cm-1.0_dp)/2.0_dp
        a_Acm = a_A**(-a_cm)
        a_m1 = 4.0_dp*a_C/a_Q0
        a_m2 = 1.0_dp/(a_H0*a_H0*a_H0*a_H0)

        DO No = 1,Particles % NumberOfParticles
           IF (Particles % Coordinate(No,1)<0.0_dp) THEN
              Particles % H(No) = a_H0
              Particles % Velocity(No,1)=a_v0
              Particles % GradVel(No,1) = 0.0_dp
           ELSE
              Particles % H(No) = (a_m1*Particles % Coordinate(No,1) + a_m2)**(-0.25_dp)
              Particles % Velocity(No,1) = a_Q0/Particles % H(No)
              Particles % GradVel(No,1) = a_C * Particles % H(No)**(3.0_dp)
           ENDIF
        END DO

        Particles % Gradvel(:,2) = 0.0_dp
        Particles % Velocity(:,2) = 0.0_dp
        Particles % GridVelocity = Particles % Velocity
        Particles % NextCoordinate = Particles % Velocity

     END IF
  ELSE
     !if timestep > 1 and/or dt > 0

     IF (xpicm < 1) THEN
        CALL Info(Solvername,'FLIP particle update',Level=1)

        !FOR 2nd ORDER: old velocity (V^(k)) on xpic(No,1:2)
        Particles % xpic(:,1:2) = Particles % Velocity(:,1:2)

        !use FLIP
        CALL Info(Solvername,'MPM interpolation of velocity to particles...',Level=1)

        !acceleration based velocity update [ V^(k+1) = V^(k) + S(v^(k+)-v^(k))  ]
        CALL MPMMeshVectorToParticle(Particles, Model, 2, count)
        CALL Info(Solvername,'interpolation done',Level=4)

        !2nd ORDER
        !now xpic(1 and 2 ) are acceleration*dt, which equals S(v^(k+)-v^(k)) on grid
        !here, we get S(v^(k+)-v^(k)) = V^(k+1) - V^(k)
        Particles % xpic(:,1:2) = Particles % Velocity(:,1:2)-Particles % xpic(:,1:2)

        !NextCoordinate (below) will be multiplied with dt and added to previous position (X^(k))
        !for the position update: X^(k+1) = X^(k) + dt * Particles % NextCoordinate
        !where Particles % NextCoordinate = S(v^(k+)) - 0.5*(S(v^(k+)-v^(k)) )
        Particles % NextCoordinate(:,1:2) = Particles % GridVelocity(:,1:2) - &
             0.5_dp*Particles % xpic(:,1:2)

        Particles % xpic = 0.0_dp
     END IF

     IF (xpicm == 1) THEN
        Particles % Velocity = Particles % GridVelocity
        Particles % NextCoordinate = Particles % GridVelocity
     END IF

     IF (xpicm > 1) THEN

        IF (MOD(count,xpicinterval) == 0) THEN
           CALL Info(Solvername,'Calling XPIC',Level=3)
           CALL XPIC(Particles,Model,xpicm)
           CALL Info(Solvername,'XPIC done',Level=3)
        ELSE
           CALL Info(Solvername,'FLIP particle update',Level=1)

           !FOR 2nd ORDER: old velocity (V^(k)) on xpic(No,1:2)
           Particles % xpic(:,1:2) = Particles % Velocity(:,1:2)

           !use FLIP
           CALL Info(Solvername,'MPM interpolation of velocity to particles...',Level=1)

           !acceleration based velocity update [ V^(k+1) = V^(k) + S(v^(k+)-v^(k))  ]
           CALL MPMMeshVectorToParticle(Particles, Model, 2, count)
           CALL Info(Solvername,'interpolation done',Level=4)

           !2nd ORDER
           !now xpic(1 and 2 ) are acceleration*dt, which equals S(v^(k+)-v^(k)) on grid
           !here, we get S(v^(k+)-v^(k)) = V^(k+1) - V^(k)
           Particles % xpic(:,1:2) = Particles % Velocity(:,1:2)-Particles % xpic(:,1:2)

           !NextCoordinate (below) will be multiplied with
           !dt and added to previous position (X^(k))
           !for the position update: X^(k+1) = X^(k) + dt * Particles % NextCoordinate
           !where Particles % NextCoordinate = S(v^(k+)) - 0.5*(S(v^(k+)-v^(k)) )
           Particles % NextCoordinate(:,1:2) = Particles % GridVelocity(:,1:2) - &
                0.5_dp*Particles % xpic(:,1:2)

           Particles % xpic = 0.0_dp
        END IF
     END IF
  END IF


  IF (Particles % usestaticparticles) THEN
     DO No = 1,Particles % NumberOfParticles
        IF (Particles % Static(No)) THEN
           Particles % GridVelocity(No,:) = 0.0_dp
           Particles % Velocity(No,:) = 0.0_dp
           Particles % NextCoordinate(No,:) = 0.0_dp
        END IF
     END DO
  END IF


  !----------------------------------------------------------------------------------------
  ! Solution Converged: Undo the mesh edits needed for the calving front boundary condition
  !----------------------------------------------------------------------------------------
  j = 0
  DO i = ElemFirst,ElemLast
     j = j+1
     Element => Mesh % Elements( i )
     Model % CurrentElement => Element
     IF (LElem(j) .NE. (-1)) Element % BoundaryInfo % Left => Mesh % Elements(LElem( j ))
     IF (RElem(j) .NE. (-1)) Element % BoundaryInfo % Right=> Mesh % Elements(RElem( j ))
  END DO

!!! reset Model Dimension to dim
  IF (DIM.eq.(STDOFs+1)) CurrentModel % Dimension = DIM

  DEALLOCATE(LElem,RElem)


  !nullify the MP pointers


  NULLIFY(MP % Velo, MP % ub, MP % slip2, MP % exyd4, MP % ezzd3m1, MP % eyyd2m1, &
       MP % exxd1m1, MP % Dxy, MP % Dzz, MP % Dyy, MP % Dxx, MP % Ee, MP % fB,&
       MP % fN, MP % Hf, MP % Exy, MP % Ezz, MP % GridVelocity, MP % GradVel,&
       MP % driveforce, MP % slip, MP % muder, MP % eta, MP % DSRxy, MP % DSRyy, &
       MP % DSRxx)


  endtime = RealTime()

  WRITE( Message, * ) 'SOLUTION TIME  : ', endtime-starttime
  CALL Info(SolverName, Message, Level=1 )

  !**************************************************************************
  !                             END Iterations
  !**************************************************************************


CONTAINS


  !> Matrix assembly using material points (particles)
  !------------------------------------------------------------------------------

  SUBROUTINE LocalMatrixUVSSAMPM(  STIFF, FORCE, Element, n, Nodes, g, &
       rho, rhow,  LocalU, LocalV,Friction, fm,  &
       cm, STDOFs , Newton,model,&
       GridVel,GridH,reweightmpm,count,usezerostressdamage)


    !------------------------------------------------------------------------------
    USE TYPES
    USE DefUtils
    USE SolverUtils
    USE MeshUtils
    USE MPMUtils

    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LocalU(:), LocalV(:)
    INTEGER :: n, cp , STDOFs, Friction, numberofparticles,count
    REAL(KIND=dp) :: cm, fm, fq, TempParams(18), OrigParams(18),Normal(2),norm,d1prev
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: GridVel,GridH
    REAL(KIND=dp) :: Zs,GradZs(3),Vel(3),GradVel(3,3),noforce,fulldamf(STDOFs)
    LOGICAL :: Newton,reducedam,SEP,alwaysisodir
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ,detj2,&
         Basis2(n), dBasisdx2(n,3),ddBasisddx2(n,3,3)
    REAL(KIND=dp) :: gtotdbasis(n,3), mtotdbasis(n,3)
    REAL(KIND=dp) :: II(3,3),ID(3,3),Coord(3),LocalCoord(3),&
         dF(2,2),F(2,2),detF,invdF(2,2),oldgradvel(4),modifiedD(4)
    REAL(KIND=dp) :: g, rho, dhdx, dhdy , gravity, density,gridres,Hnew,rhoscale
    REAL(KIND=dp) :: beta, LinVelo, fC, fN, alpha, fB, velgradscale
    REAL(KIND=dp) :: gradS(2), A(2,2), StrainA(2,2), StrainB(2,2)
    REAL(KIND=dp) :: Exx, Eyy, Exy, Ezz, Ee, MinSRInv ,MinH, h2
    REAL(KIND=dp) :: rhow,minN,PostPeak,FricMaxVal,sealevel,bedrock,hf,diff
    REAL(KIND=dp) :: Jac(2*n,2*n), SOL(2*n), scale,LocalDistance,&
         SqrtElementMetric,txx,tyy,tzz,txy,ds,hab,db,d1,t0,fedir(3,3),zero=0.0_dp
    REAL(KIND=dp) :: angle, angle1,angle2, s1,s2,s1a,s1b,h_im,mpmweight,area,volinelem
    LOGICAL :: Stat, NewtonLin, fNewtonLin,UseBfScale,&
         delpstress,UseFEdit,&
         editforzerostress,UseZeroStressDamage,shelfonly,&
         reweightmpm,fpgroundonly,visited,applyzerostress
    INTEGER :: i, j, t, p, q , dim, No,ind,mm,sf,count2
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: bfscale(2)
    REAL(KIND=dp),POINTER :: h,Dxx,Dyy,Dzz,Dxy,&
         eta,muder,slip,slip2,&
         Velo(:),DSRxx,DSRyy,&
         DSRxy,driveforce(:),ub,falpha(:)
    REAL(KIND=dp) :: North,South,East,West,newpg(2),nxmax,nxmin,nymax,nymin
    REAL(KIND=dp) :: EV(3,3),xc,yc,xf,yf,slope,bf,nxx,nyy,xi,yi,xi_diff,xc_diff,&
         yi_diff,yc_diff,fafin(n,2),xis,xcs,yis,ycs

#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: t1,t2,tp,tm
#else
    REAL(KIND=dp) :: t1,t2,tp,tm
#endif


    STIFF = 0.0_dp
    FORCE = 0.0_dp
    Jac=0.0_dp

    ! Use Newton Linearisation
    NewtonLin = (Newton.AND.(cm.NE.1.0_dp))
    fNewtonLin = (Newton.AND.(fm.NE.1.0_dp))



    ind = Element % ElementIndex
    numberofparticles = ElemParticles(ind) % NumberOfParticles

    IF (reweightmpm) THEN
       Area = ElementArea(Model % Mesh,Element,n)
       mpmweight = Area/ElemTrack(ind)%Volume
    END IF

    IF (Particles % uplag) THEN
       IP = GaussPoints( Element , np=INT(Particles % elementfraction) )
    END IF




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO t = 1,numberofparticles

       No = ElemParticles(ind) % p(t)
       IF (Particles % Status(No) == PARTICLE_LOST) CYCLE
       IF (Particles % ShapeFunctions == 'gimpm') THEN

          !implicit gimpm element-wise basis
          stat = GIMPMElementInfo( t,Particles, Model,Element, Nodes, No, &
               detJ, scale, .False., Basis,dBasisdx)

          stat = sMPMElementInfo( Element,Particles, Model, Nodes, No, &
               Particles % gridres, Basis2,dBasisdx2)

          detJ = Particles % PVolume(No)
          scale = 1.0_dp
       ELSE
          IF (Particles % uplag) THEN
             !stat = ParticleElementInfo( Element, Particles % Coordinate(No,:), &
             !   detJ, Basis2, dBasisdx2 )

             stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t),  detJ, Basis2, dBasisdx2, ddBasisddx2, .FALSE. )

             detJ = detJ * IP%S(t)
          ELSE
             stat = sMPMElementInfo( Element,Particles, Model, Nodes, No, &
                  Particles % gridres, Basis2,dBasisdx2)
             detJ = Particles % PVolume(No)
          END IF

          Basis = Basis2
          dBasisdx = dBasisdx2
          scale = 1.0_dp
       END IF


       IF (reweightmpm) THEN
          detJ = detJ * mpmweight
       END IF

       h => Particles % H(No)

       Dxx     => MP % Dxx(No)
       Dyy     => MP % Dyy(No)
       Dzz     => MP % Dzz(No)
       Dxy     => MP % Dxy(No)

       DSRxx      =>  MP % DSRxx(No)
       DSRyy      =>  MP % DSRyy(No)
       DSRxy      =>  MP % DSRxy(No)
       eta        =>  MP % eta(No)
       muder      =>  MP % muder(No)
       slip       =>  MP % slip(No)
       driveforce =>  MP % driveforce(No,1:2)


       IF (iFriction>1) THEN
          slip2   => MP % slip2(No)
          ub      => MP % ub(No)
          Velo    => MP % Velo(No,1:2)
       END IF

       IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1
       IF (iFriction==1) fNewtonLin = .FALSE.

       StrainA=0.0_dp
       StrainB=0.0_dp

       IF (NewtonLin) THEN
          StrainA(1,1)=2.0_dp*Particles % GradVel(No,1) !Exx

          IF (STDOFs.EQ.2) THEN
             StrainB(1,1)= 0.5_dp * Particles % GradVel(No,3)

             StrainA(1,2)=Particles % GradVel(No,2) !Eyy
             StrainB(1,2)=0.5_dp * Particles % GradVel(No,4)

             StrainA(2,1)=Particles % GradVel(No,1) !Exx
             StrainB(2,1)=0.5_dp * Particles % GradVel(No,3)

             StrainA(2,2)=2.0_dp*Particles % GradVel(No,2) !Eyy
             StrainB(2,2)=0.5_dp * Particles % GradVel(No,4)
          END IF
       END IF


       A = 0.0_dp
       DO p=1,n
          DO q=1,n
             !!IMPLEMENTATION OF ANISOTROPIC DAMAGE IN STIFFNESS MATRIX:

             A(1,1) = (dBasisdx2(q,1)*(2.0_dp-Dxx-Dzz) - &
                  0.5_dp*dBasisdx2(q,2)*Dxy) * dBasisdx(p,1)

             IF (STDOFs.EQ.2) THEN

                A(1,1) = A(1,1) + (0.25_dp*dBasisdx2(q,2)*(2.0_dp-Dxx-Dyy) - &
                     0.5_dp*dBasisdx2(q,1)*Dxy) * dBasisdx(p,2)

                A(1,2) = (dBasisdx2(q,2)*(1.0_dp-Dzz) - 0.5_dp*dBasisdx2(q,1)*Dxy) * dBasisdx(p,1) + &
                     (-0.5_dp*dBasisdx2(q,2)*Dxy + 0.25_dp*dBasisdx2(q,1)*(2.0_dp-Dxx-Dyy)) * dBasisdx(p,2)

                A(2,1) = (dBasisdx2(q,1)*(1.0_dp-Dzz) - 0.5_dp*dBasisdx2(q,2)*Dxy)*dBasisdx(p,2) + &
                     (-0.5_dp*dBasisdx2(q,1)*Dxy + 0.25_dp*dBasisdx2(q,2)*(2.0_dp-Dxx-Dyy))*dBasisdx(p,1)

                A(2,2) = (dBasisdx2(q,2)*(2.0_dp-Dyy-Dzz)-0.5_dp*dBasisdx2(q,1)*Dxy)*dBasisdx(p,2) + &
                     (0.25_dp*dBasisdx2(q,1)*(2.0_dp-Dxx-Dyy) - 0.5_dp*dBasisdx2(q,2)*Dxy ) *dBasisdx(p,1)

             END IF

             A = 2.0_dp * h * eta * A

             DO i=1,STDOFs

                STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) +&
                     Slip * Basis2(q) * Basis(p)  * detJ

                DO j=1,STDOFs
                   STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                        A(i,j) * detJ
                END DO
             END DO


             IF ((fNewtonLin).AND.(iFriction > 1)) THEN
                DO i=1,STDOFs
                   Do j=1,STDOFs

                      STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                           Slip2 * Velo(i) * Velo(j) * Basis2(q) * Basis(p) * detJ
                   END DO
                END DO
             END IF

             IF (NewtonLin) then
                IF (STDOFs.EQ.1) THEN
                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                        detJ * 2.0_dp * h * StrainA(1,1) *dBasisdx(p,1) * &
                        muder * 2.0_dp * DSRxx * dBasisdx2(q,1)
                ELSE IF (STDOFs.EQ.2) THEN
                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                        detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                        (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*DSRxx+DSRyy)*&
                        dBasisdx2(q,1)+DSRxy*dBasisdx2(q,2))

                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) +&
                        detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                        (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*DSRyy+DSRxx)*&
                        dBasisdx2(q,2)+DSRxy*dBasisdx2(q,1))

                   Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) +&
                        detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                        (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*DSRxx+DSRyy)*&
                        dBasisdx2(q,1)+DSRxy*dBasisdx2(q,2))

                   Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) +&
                        detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                        (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*DSRyy+DSRxx)*&
                        dBasisdx2(q,2)+DSRxy*dBasisdx2(q,1))
                END IF
             END IF

          END DO

          DO i=1,STDOFs

             FORCE((STDOFs)*(p-1)+i) = FORCE((STDOFs)*(p-1)+i) - &
                  driveforce(i) * detJ * Basis(p)


          END DO

          IF ((fNewtonLin).AND.(iFriction>1)) THEN
             DO i=1,STDOFs
                FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) + &
                     Slip2 * Velo(i) * ub * ub  * detJ * Basis(p)
             END DO
          END IF
       END DO
    END DO



    IF (NewtonLin) THEN
       SOL(1:STDOFs*n:STDOFs)=LocalU(1:n)
       IF (STDOFs.EQ.2) SOL(2:STDOFs*n:STDOFs)=LocalV(1:n)

       STIFF(1:STDOFs*n,1:STDOFs*n) = STIFF(1:STDOFs*n,1:STDOFs*n) + &
            Jac(1:STDOFs*n,1:STDOFs*n)
       FORCE(1:STDOFs*n) = FORCE(1:STDOFs*n) + &
            MATMUL(Jac(1:STDOFs*n,1:STDOFs*n),SOL(1:STDOFs*n))
    END IF

    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUVSSAMPM

  !------------------------------------------------------------------------------


  !> Essentially the regular FEM solver
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixUVSSAFEM(  STIFF, FORCE, Element, n, Nodes, gravity, &
       Density, Viscosity, LocalZb, LocalZs, LocalU, LocalV, &
       Friction, LocalBeta, fm, LocalLinVelo, fq, fricmaxval, &
       NodalGM,NodalBed,SEP,rhow, &
       cm, MinSRInv, MinH, STDOFs , Newton,minn )
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), gravity(:), Density(:), &
         Viscosity(:), LocalZb(:), LocalZs(:), &
         LocalU(:), LocalV(:) , LocalBeta(:), &
         LocalLinVelo(:)
    REAL(KIND=dp) :: NodalGM(:),NodalBed(:)
    LOGICAL :: SEP,PartlyGroundedElement
    REAL(KIND=dp) :: rhow,fricmaxval,minn
    REAL(KIND=dp) :: Bedrock,Hf
    INTEGER :: n, cp , STDOFs, Friction
    REAL(KIND=dp) :: cm, fm, fq
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Newton
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ
    REAL(KIND=dp) :: g, rho, eta, h, dhdx, dhdy , muder
    REAL(KIND=dp) :: beta, LinVelo, fC, fN, Velo(2), ub, alpha, fB
    REAL(KIND=dp) :: gradS(2), A(2,2), StrainA(2,2), StrainB(2,2), Exx, Eyy, Exy, Ezz, Ee, MinSRInv ,MinH
    REAL(KIND=dp) :: Jac(2*n,2*n), SOL(2*n), Slip, Slip2
    LOGICAL :: Stat, NewtonLin, fNewtonLIn
    INTEGER :: i, j, t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    !------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    Jac=0.0_dp

    ! Use Newton Linearisation
    NewtonLin = (Newton.AND.(cm.NE.1.0_dp))
    fNewtonLin = (Newton.AND.(fm.NE.1.0_dp))

    IF (SEP) THEN
       PartlyGroundedElement=(ANY(NodalGM(1:n).GE.0._dp).AND.ANY(NodalGM(1:n).LT.0._dp))
       IF (PartlyGroundedElement) THEN
          IP = GaussPoints( Element , np=GLnIP )
       ELSE
          IP = GaussPoints( Element )
       ENDIF
    ELSE
       IP = GaussPoints( Element )
    ENDIF

    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       ! Needed Integration Point value

       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )
       eta = SUM( Viscosity(1:n) * Basis(1:n) )
       gradS = 0.0_dp
       gradS(1) = SUM( LocalZs(1:n) * dBasisdx(1:n,1) )
       IF (STDOFs == 2) gradS(2) = SUM( LocalZs(1:n) * dBasisdx(1:n,2) )
       h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n) )
       h=max(h,MinH)

       beta = SUM( LocalBeta(1:n) * Basis(1:n) )


       IF (SEP) THEN
          IF (ALL(NodalGM(1:n).LT.0._dp)) THEN
             beta=0._dp
          ELSE IF (PartlyGroundedElement) THEN
             bedrock = SUM( NodalBed(1:n) * Basis(1:n) )
             Hf= rhow * (sealevel-bedrock) / rho
             if (h.lt.Hf) beta=0._dp
          END IF
       END IF
       IF (iFriction > 1) THEN
          LinVelo = SUM( LocalLinVelo(1:n) * Basis(1:n) )
          IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1
          Velo = 0.0_dp
          Velo(1) = SUM(LocalU(1:n) * Basis(1:n))
          IF (STDOFs == 2) Velo(2) = SUM(LocalV(1:n) * Basis(1:n))
          ub = SQRT(Velo(1)*Velo(1)+Velo(2)*Velo(2))
          Slip2=1.0_dp
          IF (ub < LinVelo) then
             ub = LinVelo
             Slip2=0.0_dp
          ENDIF
       END IF

       ! IF (iFriction==3) THEN
       !    fC = SUM( LocalC(1:n) * Basis(1:n) )
       !    fN = SUM( LocalN(1:n) * Basis(1:n) )

       IF (iFriction==3) THEN
          ! fC = SUM( LocalC(1:n) * Basis(1:n) )
          ! fN = SUM( LocalN(1:n) * Basis(1:n) )
          fC = FricMaxVal

          bedrock = SUM( NodalBed(1:n) * Basis(1:n) )
          Hf= rhow * (sealevel-bedrock) / rho

          Hf = MAX(0.0_dp,Hf)
          fN = rho*g*(h-Hf)

          ! Effective pressure should be >0 (for the friction law)
          fN = MAX(fN, MinN)
       END IF

       IF (iFriction==1) THEN
          Slip = beta
          fNewtonLin = .FALSE.
       ELSE IF (iFriction==2) THEN
          Slip = beta * ub**(fm-1.0_dp)
          Slip2 = Slip2*Slip*(fm-1.0_dp)/(ub*ub)
       ELSE IF (iFriction==3) THEN
          IF (PostPeak.NE.1.0_dp) THEN
             alpha = (PostPeak-1.0_dp)**(PostPeak-1.0_dp) / PostPeak**PostPeak
          ELSE
             alpha = 1.0_dp
          END IF
          fB = alpha * (beta / (fC*fN))**(PostPeak/fm)
          Slip = beta * ub**(fm-1.0_dp) / (1.0_dp + fB * ub**PostPeak)**fm
          Slip2 = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
               fm*PostPeak*fB*ub**(PostPeak-2.0_dp)/(1.0_dp+fB*ub**PostPeak))
       END IF

       !------------------------------------------------------------------------------
       ! In the non-linear case, effective viscosity
       IF (cm.NE.1.0_dp) THEN
          Exx = SUM(LocalU(1:n)*dBasisdx(1:n,1))
          Eyy = 0.0_dp
          Exy = 0.0_dp
          IF (STDOFs.EQ.2) THEN
             Eyy = SUM(LocalV(1:n)*dBasisdx(1:n,2))
             Ezz = -Exx - Eyy
             Exy = SUM(LocalU(1:n)*dBasisdx(1:n,2))
             Exy = 0.5_dp*(Exy + SUM(LocalV(1:n)*dBasisdx(1:n,1)))
             Ee = 0.5_dp*(Exx*Exx + Eyy*Eyy + Ezz*Ezz) + Exy*Exy
          ELSE
             Ee = Exx*Exx
          END IF
          muder = eta * 0.5_dp * (2.0_dp**cm) * ((cm-1.0_dp)/2.0_dp) *  Ee**((cm-1.0_dp)/2.0_dp - 1.0_dp)
          IF (sqrt(Ee) < MinSRInv) THEN
             Ee = MinSRInv*MinSRInv
             muder = 0.0_dp
          END IF
          eta = eta * 0.5_dp * (2.0_dp**cm) * Ee**((cm-1.0_dp)/2.0_dp)
       END IF


       StrainA=0.0_dp
       StrainB=0.0_dp
       IF (NewtonLin) THEN
          StrainA(1,1)=SUM(2.0_dp*dBasisdx(1:n,1)*LocalU(1:n))

          IF (STDOFs.EQ.2) THEN
             StrainB(1,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

             StrainA(1,2)=SUM(dBasisdx(1:n,2)*LocalV(1:n))
             StrainB(1,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))

             StrainA(2,1)=SUM(dBasisdx(1:n,1)*LocalU(1:n))
             StrainB(2,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

             StrainA(2,2)=SUM(2.0_dp*dBasisdx(1:n,2)*LocalV(1:n))
             StrainB(2,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))
          END IF
       END IF

       A = 0.0_dp
       DO p=1,n
          DO q=1,n
             A(1,1) = 2.0_dp*dBasisdx(q,1)*dBasisdx(p,1)
             IF (STDOFs.EQ.2) THEN
                A(1,1) = A(1,1) + 0.5_dp*dBasisdx(q,2)*dBasisdx(p,2)

                A(1,2) = dBasisdx(q,2)*dBasisdx(p,1) + &
                     0.5_dp*dBasisdx(q,1)*dBasisdx(p,2)

                A(2,1) = dBasisdx(q,1)*dBasisdx(p,2) + &
                     0.5_dp*dBasisdx(q,2)*dBasisdx(p,1)

                A(2,2) = 2.0*dBasisdx(q,2)*dBasisdx(p,2) +&
                     0.5_dp*dBasisdx(q,1)*dBasisdx(p,1)
             END IF
             A = 2.0_dp * h * eta * A
             DO i=1,STDOFs
                STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) +&
                     Slip * Basis(q) * Basis(p) * IP % S(t) * detJ
                DO j=1,STDOFs
                   STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                        A(i,j) * IP % S(t) * detJ
                END DO
             END DO

             IF ((fNewtonLin).AND.(iFriction > 1)) THEN
                DO i=1,STDOFs
                   Do j=1,STDOFs
                      STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                           Slip2 * Velo(i) * Velo(j) * Basis(q) * Basis(p) * IP % S(t) * detJ
                   End do
                END DO
             END IF

             IF (NewtonLin) then
                IF (STDOFs.EQ.1) THEN
                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                        IP % S(t) * detJ * 2.0_dp * h * StrainA(1,1)*dBasisdx(p,1) * &
                        muder * 2.0_dp * Exx*dBasisdx(q,1)
                ELSE IF (STDOFs.EQ.2) THEN
                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                        IP % S(t) * detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                        (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2))

                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) +&
                        IP % S(t) * detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                        (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1))

                   Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) +&
                        IP % S(t) * detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                        (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2))

                   Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) +&
                        IP % S(t) * detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                        (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1))
                END IF
             END IF

          END DO

          DO i=1,STDOFs
             FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) - &
                  rho*g*h*gradS(i) * IP % s(t) * detJ * Basis(p)
             !was basis
          END DO

          IF ((fNewtonLin).AND.(iFriction>1)) THEN
             DO i=1,STDOFs
                FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) + &
                     Slip2 * Velo(i) * ub * ub * IP % s(t) * detJ * Basis(p)
                !was basis
             END DO
          END IF

       END DO
    END DO

    IF (NewtonLin) THEN
       SOL(1:STDOFs*n:STDOFs)=LocalU(1:n)
       IF (STDOFs.EQ.2) SOL(2:STDOFs*n:STDOFs)=LocalV(1:n)

       STIFF(1:STDOFs*n,1:STDOFs*n) = STIFF(1:STDOFs*n,1:STDOFs*n) + &
            Jac(1:STDOFs*n,1:STDOFs*n)
       FORCE(1:STDOFs*n) = FORCE(1:STDOFs*n) + &
            MATMUL(Jac(1:STDOFs*n,1:STDOFs*n),SOL(1:STDOFs*n))
    END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUVSSAFEM
  !------------------------------------------------------------------------------



  !> The regular FEM solver, edited for damage
  SUBROUTINE LocalMatrixUVSSAFEMDamage(  STIFF, FORCE, Element, n, Nodes, gravity, &
       Density, Viscosity, LocalZb, LocalZs, LocalU, LocalV,LocalD, &
       Friction, LocalBeta, fm, LocalLinVelo, fq, fricmaxval, &
       NodalGM,NodalBed,SEP,rhow, &
       cm, MinSRInv, MinH, STDOFs , Newton,minn,reweightmpm)

    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), gravity(:), Density(:), &
         Viscosity(:), LocalZb(:), LocalZs(:), &
         LocalU(:), LocalV(:) , LocalBeta(:), &
         LocalLinVelo(:),LocalD(:,:)
    REAL(KIND=dp) :: NodalGM(:),NodalBed(:)
    LOGICAL :: SEP,PartlyGroundedElement
    REAL(KIND=dp) :: rhow,fricmaxval,minn
    REAL(KIND=dp) :: Bedrock,Hf
    INTEGER :: n, cp , STDOFs, Friction
    REAL(KIND=dp) :: cm, fm, fq
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Newton
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ
    REAL(KIND=dp) :: g, rho, eta, h, dhdx, dhdy , muder
    REAL(KIND=dp) :: beta, LinVelo, fC, fN, Velo(2), ub, alpha, fB
    REAL(KIND=dp) :: gradS(2), A(2,2), StrainA(2,2), StrainB(2,2), Exx, Eyy, Exy, Ezz, Ee, MinSRInv ,MinH
    REAL(KIND=dp) :: Jac(2*n,2*n), SOL(2*n), Slip, Slip2
    REAL(KIND=dp) :: Dxx,Dyy,Dzz,Dxy,exxd1m1,eyyd2m1,ezzd3m1,exyd4,DSRxx,DSRyy,DSRxy
    LOGICAL :: Stat, NewtonLin, fNewtonLIn
    INTEGER :: i, j, t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    INTEGER :: ind,No
    LOGICAL :: reweightmpm
    REAL(KIND=dp) :: Area,scale,mpmweight
    !------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    Jac=0.0_dp

    ! Use Newton Linearisation
    NewtonLin = (Newton.AND.(cm.NE.1.0_dp))
    fNewtonLin = (Newton.AND.(fm.NE.1.0_dp))

    IF (SEP) THEN
       PartlyGroundedElement=(ANY(NodalGM(1:n).GE.0._dp).AND.ANY(NodalGM(1:n).LT.0._dp))
       IF (PartlyGroundedElement) THEN
          IP = GaussPoints( Element , np=GLnIP )
       ELSE
          IP = GaussPoints( Element )
       ENDIF
    ELSE
       IP = GaussPoints( Element )
    ENDIF

    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       ! Needed Integration Point value

       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )
       eta = SUM( Viscosity(1:n) * Basis(1:n) )
       gradS = 0.0_dp
       gradS(1) = SUM( LocalZs(1:n) * dBasisdx(1:n,1) )
       IF (STDOFs == 2) gradS(2) = SUM( LocalZs(1:n) * dBasisdx(1:n,2) )
       h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n) )
       h=max(h,MinH)

       beta = SUM( LocalBeta(1:n) * Basis(1:n) )


       IF (SEP) THEN
          IF (ALL(NodalGM(1:n).LT.0._dp)) THEN
             beta=0._dp
          ELSE IF (PartlyGroundedElement) THEN
             bedrock = SUM( NodalBed(1:n) * Basis(1:n) )
             Hf= rhow * (sealevel-bedrock) / rho
             if (h.lt.Hf) beta=0._dp
          END IF
       END IF
       IF (iFriction > 1) THEN
          LinVelo = SUM( LocalLinVelo(1:n) * Basis(1:n) )
          IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1
          Velo = 0.0_dp
          Velo(1) = SUM(LocalU(1:n) * Basis(1:n))
          IF (STDOFs == 2) Velo(2) = SUM(LocalV(1:n) * Basis(1:n))
          ub = SQRT(Velo(1)*Velo(1)+Velo(2)*Velo(2))
          Slip2=1.0_dp
          IF (ub < LinVelo) then
             ub = LinVelo
             Slip2=0.0_dp
          ENDIF
       END IF

       ! IF (iFriction==3) THEN
       !    fC = SUM( LocalC(1:n) * Basis(1:n) )
       !    fN = SUM( LocalN(1:n) * Basis(1:n) )

       IF (iFriction==3) THEN
          ! fC = SUM( LocalC(1:n) * Basis(1:n) )
          ! fN = SUM( LocalN(1:n) * Basis(1:n) )
          fC = FricMaxVal

          bedrock = SUM( NodalBed(1:n) * Basis(1:n) )
          Hf= rhow * (sealevel-bedrock) / rho

          Hf = MAX(0.0_dp,Hf)
          fN = rho*g*(h-Hf)

          ! Effective pressure should be >0 (for the friction law)
          fN = MAX(fN, MinN)
       END IF

       IF (iFriction==1) THEN
          Slip = beta
          fNewtonLin = .FALSE.
       ELSE IF (iFriction==2) THEN
          Slip = beta * ub**(fm-1.0_dp)
          Slip2 = Slip2*Slip*(fm-1.0_dp)/(ub*ub)
       ELSE IF (iFriction==3) THEN
          IF (PostPeak.NE.1.0_dp) THEN
             alpha = (PostPeak-1.0_dp)**(PostPeak-1.0_dp) / PostPeak**PostPeak
          ELSE
             alpha = 1.0_dp
          END IF
          fB = alpha * (beta / (fC*fN))**(PostPeak/fm)
          Slip = beta * ub**(fm-1.0_dp) / (1.0_dp + fB * ub**PostPeak)**fm
          Slip2 = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
               fm*PostPeak*fB*ub**(PostPeak-2.0_dp)/(1.0_dp+fB*ub**PostPeak))
       END IF

       !------------------------------------------------------------------------------
       ! In the non-linear case, effective viscosity
       IF (cm.NE.1.0_dp) THEN
          Exx = SUM(LocalU(1:n)*dBasisdx(1:n,1))
          Eyy = 0.0_dp
          Exy = 0.0_dp
          IF (STDOFs.EQ.2) THEN
             Eyy = SUM(LocalV(1:n)*dBasisdx(1:n,2))
             Ezz = -Exx - Eyy
             Exy = SUM(LocalU(1:n)*dBasisdx(1:n,2))
             Exy = 0.5_dp*(Exy + SUM(LocalV(1:n)*dBasisdx(1:n,1)))
             Ee = 0.5_dp*(Exx*Exx + Eyy*Eyy + Ezz*Ezz) + Exy*Exy
          ELSE
             Ee = Exx*Exx
          END IF
          muder = eta * 0.5_dp * (2.0_dp**cm) * ((cm-1.0_dp)/2.0_dp) *  Ee**((cm-1.0_dp)/2.0_dp - 1.0_dp)
          IF (sqrt(Ee) < MinSRInv) THEN
             Ee = MinSRInv*MinSRInv
             muder = 0.0_dp
          END IF
          eta = eta * 0.5_dp * (2.0_dp**cm) * Ee**((cm-1.0_dp)/2.0_dp)
       END IF



       Dxx = SUM( LocalD(1:n,1) * Basis(1:n) )
       Dyy = SUM( LocalD(1:n,2) * Basis(1:n) )
       Dzz = SUM( LocalD(1:n,3) * Basis(1:n) )
       Dxy = SUM( LocalD(1:n,4) * Basis(1:n) )


       IF (NewtonLin) THEN
          exxd1m1 = Exx*(Dxx-1.0_dp)
          eyyd2m1 = Eyy*(Dyy-1.0_dp)
          ezzd3m1 = Ezz*(Dzz-1.0_dp)
          exyd4 = Exy*Dxy

          DSRxx = (1.0_dp/3.0_dp) * (-2.0_dp*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
          DSRyy = (1.0_dp/3.0_dp) * (exxd1m1-2.0_dp*eyyd2m1+ezzd3m1-exyd4)
          DSRxy = -0.5_dp * (Exy*(Dxx + Dyy - 2.0_dp) + &
               Dxy*(Exx+Eyy))
       END IF


       StrainA=0.0_dp
       StrainB=0.0_dp
       IF (NewtonLin) THEN
          StrainA(1,1)=SUM(2.0_dp*dBasisdx(1:n,1)*LocalU(1:n))

          IF (STDOFs.EQ.2) THEN
             StrainB(1,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

             StrainA(1,2)=SUM(dBasisdx(1:n,2)*LocalV(1:n))
             StrainB(1,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))

             StrainA(2,1)=SUM(dBasisdx(1:n,1)*LocalU(1:n))
             StrainB(2,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

             StrainA(2,2)=SUM(2.0_dp*dBasisdx(1:n,2)*LocalV(1:n))
             StrainB(2,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))
          END IF
       END IF

       A = 0.0_dp
       DO p=1,n
          DO q=1,n

             A(1,1) = (dBasisdx(q,1)*(2.0_dp-Dxx-Dzz) - &
                  0.5_dp*dBasisdx(q,2)*Dxy) * dBasisdx(p,1)

             IF (STDOFs.EQ.2) THEN

                A(1,1) = A(1,1) + (0.25_dp*dBasisdx(q,2)*(2.0_dp-Dxx-Dyy) - &
                     0.5_dp*dBasisdx(q,1)*Dxy) * dBasisdx(p,2)

                A(1,2) = (dBasisdx(q,2)*(1.0_dp-Dzz) - 0.5_dp*dBasisdx(q,1)*Dxy) * dBasisdx(p,1) + &
                     (-0.5_dp*dBasisdx(q,2)*Dxy + 0.25_dp*dBasisdx(q,1)*(2.0_dp-Dxx-Dyy)) * dBasisdx(p,2)

                A(2,1) = (dBasisdx(q,1)*(1.0_dp-Dzz) - 0.5_dp*dBasisdx(q,2)*Dxy)*dBasisdx(p,2) + &
                     (-0.5_dp*dBasisdx(q,1)*Dxy + 0.25_dp*dBasisdx(q,2)*(2.0_dp-Dxx-Dyy))*dBasisdx(p,1)

                A(2,2) = (dBasisdx(q,2)*(2.0_dp-Dyy-Dzz)-0.5_dp*dBasisdx(q,1)*Dxy)*dBasisdx(p,2) + &
                     (0.25_dp*dBasisdx(q,1)*(2.0_dp-Dxx-Dyy) - 0.5_dp*dBasisdx(q,2)*Dxy ) *dBasisdx(p,1)
             END IF

             A = 2.0_dp * h * eta * A


             DO i=1,STDOFs
                STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) +&
                     Slip * Basis(q) * Basis(p) * IP % S(t) * detJ
                DO j=1,STDOFs
                   STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                        A(i,j) * IP % S(t) * detJ
                END DO
             END DO

             IF ((fNewtonLin).AND.(iFriction > 1)) THEN
                DO i=1,STDOFs
                   Do j=1,STDOFs
                      STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                           Slip2 * Velo(i) * Velo(j) * Basis(q) * Basis(p) * IP % S(t) * detJ
                   End do
                END DO
             END IF

             IF (NewtonLin) then
                IF (STDOFs.EQ.1) THEN
                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                        detJ * 2.0_dp * h * StrainA(1,1) *dBasisdx(p,1) * &
                        muder * 2.0_dp * DSRxx * dBasisdx(q,1)
                ELSE IF (STDOFs.EQ.2) THEN
                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                        detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                        (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*DSRxx+DSRyy)*&
                        dBasisdx(q,1)+DSRxy*dBasisdx(q,2))

                   Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) +&
                        detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                        (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*DSRyy+DSRxx)*&
                        dBasisdx(q,2)+DSRxy*dBasisdx(q,1))

                   Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) +&
                        detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                        (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*DSRxx+DSRyy)*&
                        dBasisdx(q,1)+DSRxy*dBasisdx(q,2))

                   Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) +&
                        detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                        (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*DSRyy+DSRxx)*&
                        dBasisdx(q,2)+DSRxy*dBasisdx(q,1))
                END IF
             END IF

          END DO

          DO i=1,STDOFs
             FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) - &
                  rho*g*h*gradS(i) * IP % s(t) * detJ * Basis(p)
             !was basis
          END DO

          IF ((fNewtonLin).AND.(iFriction>1)) THEN
             DO i=1,STDOFs
                FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) + &
                     Slip2 * Velo(i) * ub * ub * IP % s(t) * detJ * Basis(p)
                !was basis
             END DO
          END IF

       END DO
    END DO


    IF (NewtonLin) THEN
       SOL(1:STDOFs*n:STDOFs)=LocalU(1:n)
       IF (STDOFs.EQ.2) SOL(2:STDOFs*n:STDOFs)=LocalV(1:n)

       STIFF(1:STDOFs*n,1:STDOFs*n) = STIFF(1:STDOFs*n,1:STDOFs*n) + &
            Jac(1:STDOFs*n,1:STDOFs*n)
       FORCE(1:STDOFs*n) = FORCE(1:STDOFs*n) + &
            MATMUL(Jac(1:STDOFs*n,1:STDOFs*n),SOL(1:STDOFs*n))
    END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUVSSAFEMDamage


  !------------------------------------------------------------------------------

  !> Here, matrix assembly for boundaries is the same for MPM as FEM
  SUBROUTINE LocalMatrixBCSSA(  STIFF, FORCE, Element, n, ENodes, rhoi, &
       g, LocalZb, LocalZs, rhow, sealevel,MinH )
    !------------------------------------------------------------------------------
    USE TYPES
    USE DefUtils
    USE SolverUtils
    USE MeshUtils
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) ::  ENodes
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:),rhow, sealevel, MinH, zb,LocalZs(:),LocalZb(:),gravity
    !  REAL(KIND=dp) :: LocalH(:)
    INTEGER :: n
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
         DetJ,Normal(3), rhoi, g, alpha, h, h_im,norm,scale
    LOGICAL :: Stat
    INTEGER :: t, i
    TYPE(GaussIntegrationPoints_t) :: IP

    !------------------------------------------------------------------------------
    STIFF = 0.0_dp
    FORCE = 0.0_dp


    ! 1D-SSA Case : concentrated force at each nodes
    IF (STDOFs==1) THEN  !1D SSA but should be 2D problem (does elmer work in 1D?)
       DO i = 1, n

          ! h = SUM( (LocalH(1:n)) * Basis(1:n))
          ! h = max(h,MinH)
          ! h_im = max(0.0_dp ,  h * (rhoi/rhow) )


          h = LocalZs(i)-LocalZb(i)
          h = max(h,MinH)
          h_im = max(0.0_dp,sealevel-LocalZb(i))

          alpha=0.5_dp * g * (rhoi * h*h - rhow * h_im*h_im)

          FORCE(i) = FORCE(i) + alpha
       END DO

       ! 2D-SSA Case : force distributed along the line
       ! This will work in DIM=3D only if working with Extruded Mesh and Preserve
       ! Baseline as been set to True to keep the 1D-BC
    ELSE IF (STDOFs==2) THEN

       IP = GaussPoints( Element )
       DO t=1,IP % n

          stat = ElementInfo( Element, ENodes, IP % U(t), IP % V(t), &
               IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

          g = ABS(g)

          !  h = SUM( (LocalH(1:n)) * Basis(1:n))
          !  h = max(h,MinH)

          h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n))
          h_im = max(0.0_dp , SUM( (sealevel-LocalZb(1:n)) * Basis(1:n)) )

          ! h_im = max(0.0_dp ,  h * (rhoi/rhow) )

          alpha=0.5_dp * g * (rhoi * h*h - rhow * h_im*h_im)


          ! Normal in the (x,y) plane
          Normal = NormalVector( Element, ENodes, IP % U(t), IP % V(t), .TRUE.)
          norm=SQRT(normal(1)*normal(1) +normal(2)*normal(2))
          Normal(1) = Normal(1)/norm
          Normal(2) = Normal(2)/norm


          DO p=1,n
             DO i=1,STDOFs

                FORCE(STDOFs*(p-1)+i) =   FORCE(STDOFs*(p-1)+i) +&
                     alpha * Normal(i) * IP % s(t) * detJ * Basis(p)
             END DO
          END DO
       END DO
    ELSE
       CALL FATAL('SSASolver-SSABasalSolver','Do not work for STDOFs <> 1 or 2')
    END IF

  END SUBROUTINE LocalMatrixBCSSA


  !------------------------------------------------------------------------------


  !> Each iteration, Update particle values needed to assemble stiffness and force matrices.
  !! The structure (type) MP is used for storage, which repurposes the memory otherwise
  !! allocated for the xpic routine (Particles % xpic) and damage layer increments
  !! (Particles % dD). For efficiency, a vectorized approach is used to update the
  !! particle values, whenever possible


  SUBROUTINE UpdateSSAParticleVals( NoP, g, rho, rhow,&
       Friction, fm, cm, MinSRInv, MinH, model,&
       minN,LinVelo,PostPeak,FricMaxVal,sealevel,&
       newton,usezerostressdamage,count,applyzerostress,UseFemMinMax,&
       FemMinX,FemMaxX,gridres,zssurf)


    USE TYPES
    USE DefUtils
    USE SolverUtils
    USE MeshUtils
    USE MPMUtils

    INTEGER :: n, cp , STDOFs, Friction, numberofparticles,count,zssurf
    REAL(KIND=dp) :: cm, fm, fq, TempParams(18), OrigParams(18),Normal(2),norm,d1prev
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: GridVel,GridH
    REAL(KIND=dp) :: Zs,GradZs(3),Vel(3),GradVel(3,3),noforce
    LOGICAL :: Newton,reducedam,SEP
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: II(3,3),ID(3,3),DSR(3,3),Coord(3),LocalCoord(3),&
         dF(2,2),F(2,2),detF,invdF(2,2),oldgradvel(4),modifiedD(4)
    REAL(KIND=dp) :: g, rho, eta, h, dhdx, dhdy , muder, gravity, density,gridres,Hnew,rhoscale
    REAL(KIND=dp) :: beta, LinVelo, fC, Velo(2), ub, alpha, velgradscale
    REAL(KIND=dp) :: gradS(2), A(2,2), StrainA(2,2), StrainB(2,2),DDD(2,2)
    REAL(KIND=dp) :: Dxx,Dyy,Dzz,Dxy,Exx, Eyy, Exy, Ezz, Ee, MinSRInv ,MinH, h2, zb
    REAL(KIND=dp) :: rhow,minN,PostPeak,FricMaxVal,sealevel,bedrock,diff,TT,DD
    REAL(KIND=dp) :: tmf,difference

    REAL(KIND=dp) :: Slip, Slip2,scale,bfscale(2),dtot,fax,fay,&
         SqrtElementMetric,txx,tyy,tzz,txy,tp,ds,hab,db,d1,t0,fedir(3,3),q,qw,dw,tau(2,2)
    REAL(KIND=dp) :: angle, angle1,angle2, s1,s2,s1a,s1b,h_im,falpha,mpmweight,area,volinelem
    LOGICAL :: Stat, NewtonLin, fNewtonLIn,UseBfScale,&
         delpstress,UseFEdit,&
         editforzerostress,UseZeroStressDamage,shelfonly,applyzerostress,&
         reweightmpm,fpgroundonly,Visited=.FALSE.,usedcompress,usedcompress2
    INTEGER :: i, j, t, p, dim, No,ind,iter,mm,sf,count2,elind,NoP
    TYPE(Nodes_t) :: Nodes

    REAL(KIND=dp) :: one=1.0_dp,zero=0.0_dp,half=0.5_dp,onethird=1.0_dp/3.0_dp,two=2.0_dp
    REAL(KIND=dp) :: exxd1m1,eyyd2m1,ezzd3m1,exyd4

    LOGICAL :: UseFemMinMax,usemelfrac,alwaysisodir
    REAL(KIND=dp) :: FemMinX,FemMaxX,FemMinXEdit,FemMaxXEdit,eps1,melfrac,melfracmod
    REAL(KIND=dp) :: EigValues(2),EigenVec(2,2),EigenVec3d(3,3),xx,yy,ww,zz,tau2
    REAL(KIND=dp) :: denom,t1d2m1,t2d1m1,t3od3m1,d4t4,ETau(2,2),SR(2,2)


    SAVE :: II,Visited,eps1,FemMinXEdit,FemMaxXEdit


    IF (.NOT. Visited) THEN

       !cn is visc exponent
       II = 0.0_dp
       II(1,1) = 1.0_dp
       II(2,2) = 1.0_dp
       II(3,3) = 1.0_dp

       eps1 = TINY(1.0_dp)

       IF (UseFemMinMax) THEN
          FemMinXEdit = FemMinX+gridres
          FemMaxXEdit = FemMaxX-gridres
       END IF

       Visited = .TRUE.
    END IF


    ! The particle values needed to assemble stiffness and force matrices are
    ! saved using the type 'MP', a set of pointers that repurpose the
    ! memory that was allocated for the xpic routine (Particles % xpic) and
    ! layer damage increments (Particles % dD). For efficiency, a vectorized
    ! approach is used wherever possible to update the values of MP. This is
    ! done in the subroutine

    Particles % dD(:,1:2,:) = 0.0_dp
    Particles % xpic = 0.0_dp
    NewtonLin = (Newton.AND.(cm.NE.1.0_dp))
    DSR = 0.0_dp

    !IF (iter == 1) THEN
    !assign damage to MP, recalling that:
    !MP % Dxx => Particles % dD(1:NoP,6,1)
    !MP % Dyy => Particles % dD(1:NoP,6,2)
    !MP % Dzz => Particles % dD(1:NoP,6,3)
    !MP % Dxy => Particles % dD(1:NoP,6,4)
    Particles % dD(1:NoP,6,1:4) = Particles % Dav(1:NoP,1:4)
    !END IF

    !----------- VECTORIZABLE PROCESSING -----------

    !friction Part 1/2
    !MP % slip(1:NoP) = max(Particles % fp(1:NoP),0.0_dp)
    MP % slip(1:NoP) = Particles % fp(1:NoP)
    WHERE (MP % slip(1:NoP) < 0.0_dp) MP % slip(1:NoP) = 0.0_dp


    IF (Particles % movegl .OR. fpgroundonly) THEN
       WHERE (Particles % Gmask(1:NoP) > 0.0_dp) MP % slip(1:NoP) = 0.0_dp
    END IF

    IF (iFriction > 1) THEN
       IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1

       IF ( Particles % VelocityDependentFriction) THEN
          !use current velocity interpolated from grid for friction
          MP % Velo(1:NoP,1) = Particles % GridVelocity(1:NoP,1)
          MP % Velo(1:NoP,2) = Particles % GridVelocity(1:NoP,2)
       ELSE
          !use particle velocity from previous timestep for friction
          MP % Velo(1:NoP,1) = Particles % Velocity(1:NoP,1)
          MP % Velo(1:NoP,2) = Particles % Velocity(1:NoP,2)
       END IF

       MP % ub(1:NoP) = SQRT( MP % Velo(1:NoP,1) * MP % Velo(1:NoP,1) &
            + MP % Velo(1:NoP,2) * MP % Velo(1:NoP,2) )

       MP % Slip2(1:NoP) = 1.0_dp

       WHERE (MP % ub < LinVelo)
          MP % ub = LinVelo
          MP % Slip2 = 0.0_dp
       END WHERE

    END IF

    IF (iFriction == 1) THEN
       fNewtonLin = .FALSE.
    END IF



    !------------------------------------------------------------------------------
    ! In the non-linear case, effective viscosity

    !your binit is B = A**(-1/n).
    MP % eta(1:NoP) = Particles % Binit(1:NoP) * (Particles % EF(1:NoP) * 2.0_dp)**(-1.0_dp * cm)


    !Exx = Particles % GradVel(No,1)
    !Eyy = Particles % GradVel(No,2)
    !Ezz = -Exx - Eyy

    MP % Ezz(1:NoP) = -Particles % GradVel(1:NoP,1) - Particles % GradVel(1:NoP,2)
    MP % Exy(1:NoP) = 0.5_dp*(Particles % GradVel(1:NoP,3)  + Particles % GradVel(1:NoP,4) )

    !Ee = 0.5_dp*(Exx*Exx + Eyy*Eyy + Ezz*Ezz) + Exy*Exy
    MP % Ee(1:NoP) = 0.5_dp*(Particles % GradVel(1:NoP,1)*Particles % GradVel(1:NoP,1) + &
         Particles % GradVel(1:NoP,2)*Particles % GradVel(1:NoP,2) + &
         MP % Ezz(1:NoP) * MP % Ezz(1:NoP)) + MP % Exy(1:NoP) * MP % Exy(1:NoP)

    !muder = eta * 0.5_dp * (2.0_dp**cm) * ((cm-1.0_dp)/2.0_dp) &
    !     *  Ee**((cm-1.0_dp)/2.0_dp - 1.0_dp)

    MP % muder = MP % eta * 0.5_dp * (2.0_dp**cm) * ((cm-1.0_dp)/2.0_dp) &
         *  (MP % Ee)**((cm-1.0_dp)/2.0_dp - 1.0_dp)

    ! IF (Ee < MinSRInv*MinSRInv) THEN
    !    Ee = MinSRInv*MinSRInv
    !    muder = 0.0_dp
    ! END IF
    WHERE (MP % Ee < MinSRInv*MinSRInv)
       MP % muder = 0.0_dp
       MP % Ee = MinSRInv*MinSRInv
    END WHERE

    !eta = eta * 0.5_dp * (2.0_dp**cm) * Ee**((cm-1.0_dp)/2.0_dp)

    MP % eta = MP % eta * 0.5_dp * (2.0_dp**cm) * (MP % Ee**((cm-1.0_dp)/2.0_dp))


    !----------- LOOP PROCESSING ----------
    !basically only contains the zero stress damage model in this version

    DO No = 1,NoP

       IF (UseFemMinMax) THEN
          IF ((Particles % Coordinate(No,1)+0.5*Particles % Length(No,1)) &
               < FemMaxXEdit .AND. &
               (Particles % Coordinate(No,1)-0.5_dp * Particles % Length(No,1))&
               > FemMinXEdit) CYCLE
       END IF


       !------------
       IF (usezerostressdamage ) THEN

          !not vectorized yet
          Exx = Particles % GradVel(No,1)
          Eyy = Particles % GradVel(No,2)
          Exy = MP % Exy(No)
          Ezz = MP % Ezz(No)
          eta = MP % eta(No)
          h = Particles % H(No)


          IF (ApplyZeroStress .AND. Particles % Damstatus(No) .NE. 1) THEN
             IF (Particles % gamma == 0.0_dp) THEN
                Tau(1,1) = 2.0_dp * Exx + Eyy
                Tau(2,2) = 2.0_dp * Eyy + Exx
                Tau(1,2) = Exy
                Tau(2,1) = Exy

                !effective resistive stress
                Tau = Tau * 2.0_dp * eta

                TT = Tau(1,1) + Tau(2,2)
                DD = Tau(1,1) * Tau(2,2) - Tau(1,2) * Tau(1,2)

                !max principal deviatoric stress
                tp = 0.5_dp * TT + sqrt(MAX(0.25_dp*TT*TT-DD,0.0_dp))
                tp = tp*h
             ELSE

                exxd1m1 = Exx*(Particles % Dav(No,1)-one)
                eyyd2m1 = Eyy*(Particles % Dav(No,2)-one)
                ezzd3m1 = Ezz*(Particles % Dav(No,3)-one)
                exyd4 = Exy*Particles % Dav(No,4)

                !tau here is actually regular tau, not effective
                Tau(1,1) = two*eta*onethird * (-two*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
                Tau(2,2) = two*eta*onethird * (exxd1m1-two*eyyd2m1+ezzd3m1-exyd4)
                !Tau(3,3) = onethird * (exxd1m1+eyyd2m1-two*ezzd3m1+two*exyd4)
                Tau(1,2) = -two*eta*half * (Exy*(Particles % Dav(No,1) + Particles % Dav(No,2) - two) + &
                     Particles % Dav(No,4)*(Exx+Eyy))
                Tau(2,1) = Tau(1,2)

                Dxx = Particles % Dav(No,1)
                Dyy = Particles % Dav(No,2)
                Dzz = Particles % Dav(No,3)
                Dxy = Particles % Dav(No,4)

                denom = one/( Dxy*Dxy + Dxx + Dyy -Dxx*Dyy - one)
                t1d2m1 = Tau(1,1)*(Dyy-one)*denom
                t2d1m1 = Tau(2,2)*(Dxx-one)*denom
                t3od3m1 = two*eta*onethird * (exxd1m1+eyyd2m1-two*ezzd3m1+two*exyd4)/(Dzz-one)
                d4t4 = Tau(1,2)*Dxy*denom

                Etau(1,1) = onethird*(two*t1d2m1 - t2d1m1 +t3od3m1 -d4t4)
                Etau(2,2) = onethird*(-t1d2m1 + two*t2d1m1 +t3od3m1 -d4t4)
                !  Etau(3,3) = onethird*(-t1d2m1 - t2d1m1 - two*t3od3m1 + two*d4t4)
                Etau(1,2) = half*denom*(Tau(1,2)*(Dxx+Dyy-two) - Dxy*(Tau(1,1)+Tau(2,2)))
                Etau(2,1) = Etau(1,2)

                Tau(1,1) = 2.0_dp * ETau(1,1) + ETau(2,2)
                Tau(2,2) = 2.0_dp * ETau(2,2) + ETau(1,1)

                Tau(1,2) = Etau(1,2)
                Tau(2,1) = Tau(1,2)


                IF (count == 0) THEN

                   IF (Particles % damage(No,3,3) == 0.0_dp) THEN
                      CALL Eigen2DSym_TryGenFirst(Tau,EigValues,EigenVec)

                      Particles % damage(No,3,3) = 1.0_dp
                      Particles % damage(No,4,1) = EigenVec(1,1)
                      Particles % damage(No,4,2) = EigenVec(2,2)
                      Particles % damage(No,4,3) = EigenVec(1,2)
                      Particles % damage(No,4,4) = EigenVec(2,1)

                      Tau = MATMUL(MATMUL(TRANSPOSE(EigenVec),Tau),EigenVec)

                      EigValues(2) = MAXVAL(Tau)

                      ETau = MATMUL(MATMUL(TRANSPOSE(EigenVec),ETau),EigenVec)

                      IF (EigValues(2) == Tau(1,1)) THEN
                         tau2 = ETau(2,2)
                      ELSE
                         tau2 = ETau(1,1)
                      END IF

                   ELSE
                      EigenVec(1,1) = Particles % damage(No,4,1)
                      EigenVec(2,2) = Particles % damage(No,4,2)
                      EigenVec(1,2) = Particles % damage(No,4,3)
                      EigenVec(2,1) = Particles % damage(No,4,4)

                      Tau = MATMUL(MATMUL(TRANSPOSE(EigenVec),Tau),EigenVec)
                      EigValues(2) = MAXVAL(Tau)
                      ETau = MATMUL(MATMUL(TRANSPOSE(EigenVec),ETau),EigenVec)

                      IF (EigValues(2) == Tau(1,1)) THEN
                         tau2 = ETau(2,2)
                      ELSE
                         tau2 = ETau(1,1)
                      END IF
                      ! EigValues(1) = 0.0_dp
                   END IF

                ELSE
                   DDD(1,1) = Particles % Dav(No,1)
                   DDD(2,2) = Particles % Dav(No,2)
                   DDD(1,2) = Particles % Dav(No,4)
                   DDD(2,1) = DDD(1,2)
                   CALL Eigen2DSym_TryGenFirst(DDD,EigValues,EigenVec)

                   DDD = MATMUL(MATMUL(TRANSPOSE(EigenVec),DDD),EigenVec)
                   Tau = MATMUL(MATMUL(TRANSPOSE(EigenVec),Tau),EigenVec)
                   ETau = MATMUL(MATMUL(TRANSPOSE(EigenVec),ETau),EigenVec)

                   IF (DDD(1,1) > DDD(2,2)) THEN
                      EigValues(2) = Tau(1,1)
                      tau2 = ETau(2,2)
                   ELSE
                      EigValues(2) = Tau(2,2)
                      tau2 = ETau(1,1)
                   END IF
                END IF

                tp = EigValues(2)*h
             END IF

             IF (tp<0.0_dp) tp = 0.0_dp
             IF (tp .NE. tp) tp = 0.0_dp

             Q = rho * g * h * h
             Qw = rhow * g * h * h
             dw = 0.0_dp !surface crev water pressure. for now, set to zero

             IF (Particles % Gmask(No) < 0.0_dp) THEN
                hab = h - (rhow/rho) * (sealevel-Particles % bedrock(No))
                hab = MAX(0.0_dp,hab)
             ELSE
                hab = 0.0_dp
             END IF


             IF (Particles % Gamma > 0.0_dp) THEN
                !tau2 = h * (2.0_dp/3.0_dp) * (EigValues(1) - 0.5_dp*EigValues(2))
                tau2 = tau2*h


                !both surface and basal crevasses
                dtot  = ( -(dw * rhow + hab*rho)*Q &
                     + rhow * (Qw * dw + h * tp) ) &
                     / ( (rhow-rho)*Q + rhow*(tp-tau2) )
                IF (dtot .NE. dtot) THEN
                   dtot  = ( -(dw * rhow + hab*rho)*Q &
                        + rhow * (Qw * dw + h * tp) ) &
                        / ( (rhow-rho)*Q + rhow*(tp-tau2) + eps1)
                END IF

                !only surface crevasses
                ds = (h * tp)/(Q + (tp-tau2))
                IF (ds .NE. ds) ds = (h * tp)/(Q + (tp-tau2) +eps1)

                !only basal crevasses
                db = rho*(H*tp - hab*Q)/(Q*(rhow-rho)+(tp-tau2)*rho)
                IF (db .NE. db) db = rho*(H+tp - hab*Q)/(Q*(rhow-rho)+(tp-tau2)*rho + eps1)

             ELSE
                !both surface and basal crevasses
                dtot  = ( -(dw * rhow + hab*rho)*Q &
                     + rhow * (Qw * dw + h * tp) ) &
                     / ( (rhow-rho)*Q + tp*rhow )

                !only surface crevasses
                ds = (h * tp)/(Q + tp)

                IF ((ds .NE. ds) .OR. (dtot .NE. dtot)) THEN
                   dtot  = ( -(dw * rhow + hab*rho)*Q &
                        + rhow * (Qw * dw + h * tp) ) &
                        / ( (rhow-rho)*Q + tp*rhow + eps1 )
                   ds = (h * tp)/(Q + tp + eps1)
                   IF (dtot .NE. dtot) dtot = 0.0_dp
                   IF (ds .NE. ds) ds = 0.0_dp
                END IF

                !only basal crevasses
                db = rho*(H*tp - hab*Q)/(Q*(rhow-rho)+tp*rho)
                IF (db .NE. db) THEN
                   db = rho*(H+tp - hab*Q)/(Q*(rhow-rho)+tp*rho + eps1)
                END IF
             END IF

             IF (dtot<zero) dtot = zero
             IF (db<zero) db = zero

             IF (zssurf == 1) THEN
                !only want surface crevasses
                dtot = ds
                db = ds
             ELSE IF (zssurf == -1) THEN
                !only want basal crevasses
                dtot = db
                ds = db
             END IF

             !IF (Particles % Gmask(No) < 0.0_dp) THEN
             dtot = MAX(ds,dtot)
             !END IF

             IF (dtot .NE. dtot) dtot = 0.0_dp

             dtot = MAX(0.0_dp,MIN(dtot,h*Particles % riftdmax))

             !damage(No,3,1) is the max principal damage from the last step
             IF (dtot/h >= Particles % damage(No,3,1)) THEN
                Particles % damage(No,3,2) = dtot/H

                IF (Particles % gamma > 0.0_dp) THEN
                   EigValues(2) = dtot/H
                   EigValues(1) = 0.0_dp

                   ww = EigValues(1)*EigenVec(1,1)
                   xx = EigValues(2)*EigenVec(1,2)
                   yy = EigValues(1)*EigenVec(2,1)
                   zz = EigValues(2)*EigenVec(2,2)

                   Particles % Dav(No,1) = EigenVec(1,1)*ww + EigenVec(1,2)*xx
                   Particles % Dav(No,2) = EigenVec(2,1)*yy + EigenVec(2,2)*zz
                   Particles % Dav(No,4) = EigenVec(2,1)*ww + EigenVec(2,2)*xx
                   Particles % Dav(No,3) = 0.0_dp !dtot/H
                ELSE
                   Particles % Dav(No,1:3) = dtot/h
                   Particles % Dav(No,4) = 0.0_dp
                END IF

             ELSE
                Particles % damage(No,3,2) = Particles % damage(No,3,1)
                Particles % Dav(No,1:4) = Particles % damage(No,2,1:4)
             END IF

             IF (UseZeroStressDamageFirstOnly) THEN
                !calculate actual ds and db and save them
                !on particles % damage (No,1:2,1)
                !These calculated heights add up to
                !dtot above.

                tp = tp*(1.0_dp-Particles % damage(No,3,2))/H
                ds  = tp/(rhoi * g)
                ds = MAX(ds,zero)
                db = (rhoi/(rhow-rhoi)) * ( ds -hab)
                db = MAX(db,zero)
                Particles % damage(No,1,1) = ds
                Particles % damage(No,1,2) = db

             END IF
          END IF

          Particles % dD(No,6,1:4) = Particles % Dav(No,1:4)
       END IF
    END DO
    !end of loop

    !-----------MORE VECTORIZABLE PROCESSING----------

    !------------ Friction Part 2/2

    IF (iFriction > 1) THEN

       IF (iFriction ==2 ) THEN
          MP % Slip = MP % Slip * MP % ub**(fm-1.0_dp)
          MP % Slip2 = MP % Slip2 * MP % Slip * &
               (fm-1.0_dp)/(MP % ub * MP % ub)

       ELSEIF (iFriction == 3) THEN

          !Effective pressure. Really just need to do this once.
          MP % Hf(1:NoP) = rhow * (sealevel-Particles % bedrock(1:NoP))/rho
          WHERE (MP % Hf < 0.0_dp) MP % Hf = 0.0_dp
          MP % fN(1:NoP) = rho*g*(Particles % H(1:NoP) - MP % Hf(1:NoP))
          ! Effective pressure should be >0 (for the friction law)
          !MP % fN = MAX(MP % fN, MinN)
          WHERE(MP % fN < MinN) MP % fN = MinN

          IF (PostPeak.NE.1.0_dp) THEN
             alpha = (PostPeak-1.0_dp)**(PostPeak-1.0_dp) / PostPeak**PostPeak
          ELSE
             alpha = 1.0_dp
          END IF
          MP % fB = alpha * (MP % slip / (FricMaxVal * MP % fN))**(PostPeak/fm)

          MP % Slip = (MP % slip) * (MP % ub)**(fm-1.0_dp) / (1.0_dp + (MP % fB) * (MP % ub)**PostPeak)**fm
          MP % Slip2 = (MP % Slip2) * (MP % Slip) * ((fm-1.0_dp) / ((MP % ub) * (MP % ub)) - &
               fm*PostPeak*(MP % fB)*(MP % ub)**(PostPeak-2.0_dp)/(1.0_dp+(MP % fB)*(MP % ub)**PostPeak))
       END IF

    END IF



    IF (NewtonLin) THEN

       !slow way
       ! ID = 0.0_dp
       ! ID(1,1) = Particles % Dav(No,1) !Dxx
       ! ID(2,1) = Particles % Dav(No,4); ID(1,2) = Particles % Dav(No,4) !Dxy
       ! ID(2,2) = Particles % Dav(No,2) !Dyy
       ! ID(3,3) = Particles % Dav(No,3) !Dzz
       ! ID = II - ID

       ! !DSR = deviatoric damage eff strain rate tensor
       ! DSR = 0.0_dp
       ! DSR(1,1) = Exx; DSR(2,2) = Eyy; DSR(3,3) = Ezz
       ! DSR(2,1) = Exy; DSR(1,2) = Exy
       ! DSR = (MATMUL(ID,DSR) + MATMUL(DSR,ID))
       ! DSR = 0.5_dp*(DSR - (II * ( ( DSR(1,1)+DSR(2,2)+DSR(3,3) )/3.0_dp) ) )

       !fast way
       MP % exxd1m1(1:NoP) = Particles % GradVel(1:NoP,1) * (MP % Dxx(1:NoP) - one)
       MP % eyyd2m1(1:NoP) = Particles % GradVel(1:NoP,2) * (MP % Dyy(1:NoP) - one)
       MP % ezzd3m1(1:NoP) = MP % Ezz(1:NoP) * (MP % Dzz(1:NoP)-one)
       MP % exyd4(1:NoP) = MP % Exy(1:NoP) * MP % Dxy(1:NoP)

       MP % DSRxx = onethird * (-two* MP % exxd1m1 + MP % eyyd2m1 + MP % ezzd3m1 - MP % exyd4)
       MP % DSRyy = onethird * (MP % exxd1m1-two* MP % eyyd2m1 + MP % ezzd3m1 - MP % exyd4)
       !DSRzz = onethird * (exxd1m1+eyyd2m1-two*ezzd3m1+two*exyd4) !doesn't need saving
       MP % DSRxy(1:NoP) = -half * (MP % Exy(1:NoP) *(MP % Dxx(1:NoP) + MP % Dyy(1:NoP) - two) + &
            MP % Dxy(1:NoP) * (Particles % GradVel(1:NoP,1) + Particles % GradVel(1:NoP,2) ))
       !DSRyx = DSRxy !doesn't need saving

    ELSE
       MP % DSRxx = zero
       MP % DSRyy = zero
       MP % DSRxy = zero
    END IF

    !driveforce
    MP % driveforce(1:NoP,1) = rho*g*Particles % H(1:NoP)*Particles % GradZs(1:NoP,1)
    MP % driveforce(1:NoP,2) = rho*g*Particles % H(1:NoP)*Particles % GradZs(1:NoP,2)

  END SUBROUTINE UpdateSSAParticleVals

!------------------------------------------------------------------------------
END SUBROUTINE MPM_SSA
!------------------------------------------------------------------------------

!!! Procedures for SSA-MPM (sMPM or GIMPM), to be called from the .sif
!!! -Alex Huth, 2020 (ahuth@princeton.edu)


!> Initialize the MPM Simulation
SUBROUTINE MPM_Initialize( Model,Solver,dt,TransientSimulation )

  USE MPMUtils

  IMPLICIT NONE
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: dim
  TYPE(Particle_t), POINTER :: Particles
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName  

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: starttime
#else
  REAL(KIND=dp) :: starttime
#endif

  WRITE(SolverName, '(A)') 'MPM_Initialize'  

  Particles => GlobalParticles
  dim = CoordinateSystemDimension()

  starttime = RealTime()
  Particles % simtime = starttime
  Particles % firstsimtime = Particles % simtime

  CALL Info(SolverName,'sMPM/GIMPM Initialization',Level=3)

  CALL Info(SolverName,'Setting Particle Preliminaries',Level=4)
  CALL SetMPMParticlePreliminaries( Particles, Model,dim )

  CALL Info(SolverName,'Initializing Particles',Level=4)
  CALL InitializeParticles( Particles, Model )

  CALL Info(SolverName,'Initializing InterpLayers Map',Level=4)
  CALL MakeInterpLayers( Particles, Model)

  IF (Particles % ShapeFunctions == 'gimpm') THEN
     CALL Info(SolverName,&
          'Assigning Particles to Elements and assembling Passive Mask',Level=4)
     CALL GetElemParticles_GIMPM(Particles, Model )
  ELSE
     CALL Info(SolverName,&
          'Assigning Particles to Elements and assembling Passive Mask',Level=4)
     CALL GetElemParticles_sMPM(Particles, Model )
  END IF

  IF (Particles % shapefunctions == 'smpm' .AND. Particles % hoop) THEN
     CALL Info(SolverName,'sMPM hoop: skipping initparticlevars',Level=1)
  ELSE
     CALL Info(SolverName,'Initializing Particle Variables',Level=4)
     CALL InitParticleVars(Particles, Model, Solver)
  END IF

  CALL Info(SolverName,'-------- Done --------',Level=3)

END SUBROUTINE MPM_Initialize

!**************************************************************************

!> Update material point def grad, gimpm shape, area, position, splitting, etc
SUBROUTINE ParticleUpdates( Model, Solver, dt, TransientSimulation)

  USE MPMUtils 
  USE DefUtils
  USE Lists

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Mesh_t), POINTER :: Mesh  
  REAL(KIND=dp) :: dt,rhoi, davsplitthres
  LOGICAL :: TransientSimulation

  INTEGER :: No, ii, infor, numoflayers,ElementIndex,VisitedTimes
  REAL(KIND=dp) :: F(2,2),G(2,2),Q(2,2),U(2)
  REAL(KIND=dp) :: gridres,maxlength, maxcalclength, &
       mincalclength,maxDPlength,mlength
  LOGICAL :: Visited=.FALSE.,SplitP,GotIt,constlintemp
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  REAL(KIND=dp) :: detF,maskav
  INTEGER :: maxcplength,maxstretchlength,maxstrain
  TYPE(Element_t), POINTER :: BulkElement

  TYPE(Variable_t), POINTER :: mask  
  INTEGER, POINTER :: maskPerm(:)
  REAL(KIND=dp), POINTER :: maskVal(:)

  !For use corner update
  TYPE(Element_t),Pointer :: Element
  INTEGER, POINTER :: NodeIndexes(:)  
  INTEGER :: ind,nn,selem(4)
  REAL(KIND=dp) :: xmin,xmax,ymin,ymax,newy,newx,Coord(3),midpcoords(4,3),SqrtElementMetric,&
       midpvel(2),tl1,tl2,sterm,Basis(4),dBasisdx(4,3),midpcoords2(4,2)
  Logical :: stat,UseCorners
  TYPE(Variable_t), POINTER :: GridVel=>NULL()  


  SAVE :: gridres,maxlength,maxDPlength,numoflayers,&
       rhoi,davsplitthres,Visited,VisitedTimes,&
       mask, maskperm, maskval,UseCorners,ElementNodes

  Particles => GlobalParticles
  Params => ListGetSolverParams()
  WRITE(SolverName, '(A)') 'ParticleUpdates'

  IF (.NOT. Visited) THEN

     nn = Model % Mesh % MaxElementNodes

     ALLOCATE(ElementNodes % x(nn),ElementNodes % y(nn),ElementNodes % z(nn))  
     mask => VariableGet(Model % Mesh % Variables, 'Mask' )
     maskPerm => mask % Perm
     maskVal => mask % Values     

     gridres = Particles % gridres
     maxlength = Particles % maxlength
     maxDPlength = Particles % maxdplength !max length of damaged particles
     numoflayers = Particles % numberofparticlelayers
     rhoi = Particles % rhoi
     davsplitthres = Particles % davsplitthres

     IF (Particles % ShapeFunctions == 'gimpm') THEN

        UseCorners = ListGetLogical( Params, 'Update GIMPM with Corners', GotIt )
        IF (.NOT. GotIt) UseCorners = .TRUE.

        IF (UseCorners) THEN
           CALL Info( SolverName, 'USING GIMPM CORNER UPDATE', Level=1 )           
        END IF

     END IF

     VisitedTimes = 0
     Visited = .TRUE.

     IF (Particles % firsttimestepzero) RETURN
  END IF


  Mesh=>GetMesh()

  IF (usecorners) THEN
     midpcoords = 0.0_dp
     GridVel => VariableGet( Mesh % Variables, 'SSAVelocity' )
     IF(.NOT. ASSOCIATED( GridVel ) ) THEN
        CALL Fatal(SolverName, 'SSAVelocity field variable does not exist: SSAVelocity')
     END IF
  END IF

  CALL Info(SolverName,'Updating F, Vol, Lengths, Splitting, and Coords',Level=3)  

  VisitedTimes = VisitedTimes + 1

  SplitP = .FALSE. 

  maxcalclength = 0.0_dp 
  mincalclength = maxlength

  dt = Particles % dtime

  Particles % xpic = 0.0_dp

  DO No = 1,Particles % NumberOfParticles

     IF (Particles % Static(No)) CYCLE

     IF (MAXVAL(ABS(Particles % Dav(No,:))) < davsplitthres) THEN
        mlength = maxlength
     ELSE
        mlength = maxDPlength
     END IF

     !for particle splitting around grounding line
     ElementIndex = GetParticleElement( Particles, No )
     IF( ElementIndex .NE. 0 ) THEN
        BulkElement => Model % Mesh % Elements( ElementIndex )
        NodeIndexes => BulkElement % NodeIndexes
        nn = BulkElement % TYPE % NumberOfNodes
        maskav = SUM(maskval(maskperm(nodeindexes)))/nn

        IF (ABS(maskav) < 0.9_dp) THEN
           !we have an ungrounding particle
           mlength = MIN(mlength,Particles % maxGLlength)

           !mark it on xpic for the splitting routine.
           Particles % xpic(No,1) = 1.0_dp
        END IF
     END IF

     !UPDATE PARTICLE DEFORMATION GRADIENT
     F = 0.0_dp
     F(1,1) = Particles % F(No,1)
     F(2,2) = Particles % F(No,2)
     F(1,2) = Particles % F(No,3)
     F(2,1) = Particles % F(No,4)
     G = 0.0_dp
     G(1,1) = 1.0_dp + dt * Particles % GradVel(No,1) 
     G(2,2) = 1.0_dp + dt * Particles % GradVel(No,2) 
     G(1,2) = dt * Particles % GradVel(No,3) 
     G(2,1) = dt * Particles % GradVel(No,4)

     !F = MATMUL(G,F)
     F(1,1) = F(1,1)*G(1,1) + F(2,1)*G(1,2)
     F(1,2) = F(1,2)*G(1,1) + F(2,2)*G(1,2)
     F(2,1) = F(1,1)*G(2,1) + F(2,1)*G(2,2)
     F(2,2) = F(1,2)*G(2,1) + F(2,2)*G(2,2)

     Particles % F(No,1) = F(1,1)
     Particles % F(No,2) = F(2,2)
     Particles % F(No,3) = F(1,2)
100  Particles % F(No,4) = F(2,1)


     IF (Particles % ShapeFunctions == 'smpm') THEN
        !
        ! sMPM
        !

        detF = F(1,1)*F(2,2) - F(1,2)*F(2,1)

        Particles % pvolume(No) = detF * Particles % GVolume(No)

        Particles % Strain(No,1:2) = Particles % Strain(No,1:2) + dt * Particles % GradVel(No,1:2)

        Particles % Length(No,1:2) = Particles % OrigLength(No,1:2) * &
             (1.0_dp + Particles % Strain(No,1:2))

        DO ii = 1,2
           maxcalclength = MAX(maxcalclength,Particles % Length(No,ii))
           mincalclength = MIN(mincalclength,Particles % Length(No,ii))
        END DO

        IF ( ANY(Particles % Length(No,:) > mlength) ) SplitP = .TRUE.

     ELSE
        !
        ! GIMPM
        !

        ind = Particles % ElementIndex(No)

        IF (ind>0 .AND. UseCorners .AND. Particles % Status(No) == PARTICLE_ACTIVE) THEN
           !
           ! "CORNERS" LENGTH UPDATE (here, using midpoints shortcut)
           !

           Element => Mesh % Elements(ind)           
           nn = Element % TYPE % NumberOfNodes
           NodeIndexes => Element % NodeIndexes
           xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
           xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
           ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
           ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))

           Coord = GetParticleCoord(Particles,No)

           !if using mid point coords
           midpcoords(1,1) = Coord(1) - 0.5_dp * Particles % Length(No,1)
           midpcoords(1,2) = Coord(2)
           midpcoords(2,1) = Coord(1)
           midpcoords(2,2) = Coord(2) + 0.5_dp * Particles % Length(No,2)
           midpcoords(3,1) = Coord(1) + 0.5_dp * Particles % Length(No,1)
           midpcoords(3,2) = Coord(2)
           midpcoords(4,1) = Coord(1)
           midpcoords(4,2) = Coord(2) - 0.5_dp * Particles % Length(No,2)

           IF (midpcoords(1,1) < xmin) THEN
              selem(1) = ElemParticles(ind) % SurroundingElems(4)
           ELSE
              selem(1) = ind
           END IF

           IF (midpcoords(2,2) > ymax) THEN
              selem(2) = ElemParticles(ind) % SurroundingElems(2)
           ELSE
              selem(2) = ind
           END IF

           IF (midpcoords(3,1) > xmax) THEN
              selem(3) = ElemParticles(ind) % SurroundingElems(6)
           ELSE
              selem(3) = ind
           END IF

           IF (midpcoords(4,2) < ymin) THEN
              selem(4) = ElemParticles(ind) % SurroundingElems(8)
           ELSE
              selem(4) = ind
           END IF

           IF (ANY(selem==0)) THEN
              Particles % Status(No) = PARTICLE_LEAVING
              goto 100
           END IF

           DO ii = 1,4
              IF (selem(ii)>0) THEN

                 Element => Mesh % Elements(selem(ii))
                 CALL GetElementNodes(ElementNodes,Element)

                 stat = sMPMElementInfoFromCoords( Particles, Model, ElementNodes, &
                      Particles % gridres,midpcoords(ii,1:2), Basis,dBasisdx)

                 CALL GetVectorFieldInMesh(GridVel,Element,Basis,midpvel)

                 midpcoords(ii,1:2) = midpcoords(ii,1:2) + dt * midpvel
              END IF
           END DO

           tl1 = MAXVAL(midpcoords(:,1))-MINVAL(midpcoords(:,1))
           tl2 = MAXVAL(midpcoords(:,2))-MINVAL(midpcoords(:,2))

           detF = F(1,1)*F(2,2) - F(1,2)*F(2,1)           
           sterm = sqrt( (detF * Particles % OrigLength(No,1)*Particles % OrigLength(No,2))/&
                (tl1*tl2) )

           IF (tl1 > 0.0_dp .AND. tl2 > 0.0_dp .AND. sterm > 0.0_dp) THEN
              Particles % Length(No,1) = tl1*sterm
              Particles % Length(No,2) = tl2*sterm
           ELSE

              !for debug, mostly. Something is probably wrong if this occurs
              WRITE( Message, * ) 'ERROR: fixing misshapen particle no...', no
              CALL Warn( SolverName, Message )                 

              tl1 = MAX(tl1,1.0_dp)
              tl2 = MAX(tl2,1.0_dp)
              sterm = sqrt( (detF * Particles % OrigLength(No,1)*Particles % OrigLength(No,2))/&
                   (tl1*tl2) )              
           END IF

           !END "CORNERS" LENGTH UPDATE

        ELSE
           ! "STRETCH TENSOR" LENGTH UPDATE

           !UPDATE PARTICLE LENGTHS
           ! F = MATMUL( TRANSPOSE(F) , F )
           G(1,1) = F(1,1)*F(1,1) + F(2,1)*F(2,1)
           G(1,2) = F(1,1)*F(1,2) + F(2,1)*F(2,2)
           G(2,1) = G(1,2)
           G(2,2) = F(1,2)*F(1,2) + F(2,2)*F(2,2)



           !U is eigenvalues, Q is eigenvectors (columns)
           !ascending order (returns smallest eigvalues first)
           U = 0.0_dp           
           CALL Eigen2D(G,U,Q)         
           U=SQRT(U)


           ! Update lengths
           ! Slow version:
           ! The final form of G is the stretch tensor
           ! G = 0.0_dp; G(1,1) = U(1); G(2,2) = U(2)           
           ! G = MATMUL(Q,G); G = MATMUL(G,TRANSPOSE(Q))
           ! Particles % Length(No,1) = Particles % OrigLength(No,1) * G(1,1)
           ! Particles % Length(No,2) = Particles % OrigLength(No,2) * G(2,2)           

           ! Fast version:
           Particles % Length(No,1) = Particles % OrigLength(No,1) * &
                ( U(1)*Q(1,1)*Q(1,1) + U(2)*Q(1,2)*Q(1,2) )

           Particles % Length(No,2) = Particles % OrigLength(No,2) * &
                (   U(1)*Q(2,1)*Q(2,1) + U(2)*Q(2,2)*Q(2,2) )

           !END "STRECH TENSOR" LENGTH UPDATE
        END IF

        !Prevents particles from getting really small (is this necessary?):
        !Currently min length of particles is set to 0.5 m
        DO ii = 1,2
           IF (Particles % Length(No,ii) < 0.5_dp) THEN
              Particles % Length(No,ii) = 0.5_dp
           END IF
        END DO

        DO ii = 1,2
           maxcalclength = MAX(maxcalclength,Particles % Length(No,ii))
           mincalclength = MIN(mincalclength,Particles % Length(No,ii))
        END DO

        Particles % pvolume(No) = Particles % Length(No,1) * Particles % Length(No,2)

        IF ( ANY(Particles % Length(No,:) > mlength) ) SplitP = .TRUE.        

        !end sMPM or GIMPM
     END IF

     Particles % mass(No) = Particles % pvolume(No) * Particles % H(No) * Particles % rhoi
  END DO

  !Update particle positions
  DO No = 1,Particles % NumberOfParticles
     IF (Particles % Static(No)) CYCLE
     Particles % Coordinate(No,1:2) = &
          Particles % Coordinate(No,1:2) + dt * Particles % NextCoordinate(No,1:2)
  END DO

  WRITE( Message, * ) 'Maximum length of particles', maxcalclength
  CALL Info( SolverName, Message, Level=1 )
  WRITE( Message, * ) 'Minimum length of particles', mincalclength
  CALL Info( SolverName, Message, Level=1 )       

  CALL Info(SolverName,'Checking for particle splitting',Level=4)  
  IF ( SplitP ) THEN
     CALL ParticleSplitting( Particles, Model, numoflayers )
     CALL Info(SolverName,'Particle Splitting Done',Level=4)  
  END IF

  Particles % mass(:) = Particles % pvolume(:) * Particles % H(:) * Particles % rhoi   

  !Assign particles to elements
  IF (Particles % ShapeFunctions == 'gimpm') THEN
     CALL Info(SolverName, &
          'Assigning GIMPM Particles to Elements and Assembling Passive Mask',Level=1)
     CALL GetElemParticles_GIMPM( Particles, Model )
  ELSE
     CALL Info(SolverName, &
          'Assigning MPM Particles to Elements and Assembling Passive Mask',Level=1)
     CALL GetElemParticles_sMPM( Particles, Model )
  END IF

  !intepolation of mass balance particles is done separately from other
  !variable interpolations
  IF (.NOT. Particles % ConstMB) THEN
     CALL MPMMeshScalarToParticle(Particles, Model, 7)
  END IF


END SUBROUTINE ParticleUpdates

!**************************************************************************

!> Updated Lagrangian mode: Update material point locs to Gauss point locs
!! Needed because mesh nodes are allowed to move in updated lagrangian mode
SUBROUTINE UpdateLagParticleUpdates( Model, Solver, dt, TransientSimulation)

  USE MPMUtils

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation,stat
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: No,nn,j,t,np,ind
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Basis(4), dBasisdx(4,3), ddBasisddx(4,3,3), detJ 
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t) :: Nodes
  INTEGER, POINTER :: NodeIndexes(:)

  Particles => GlobalParticles
  Mesh => GetMesh()
  Particles % uplag = .TRUE.

  ALLOCATE( Nodes % x(4), Nodes % y(4), Nodes % z(4))

  Do t=1,Solver % Mesh % NumberOfBulkElements
     Element => Solver % Mesh % Elements(t)
     ind = Element % ElementIndex
     IF ( Element % TYPE % NumberOfNodes .NE. 4) CYCLE     
     NodeIndexes => Element % NodeIndexes
     Nodes % x(1:4) = Solver % Mesh % Nodes % x(NodeIndexes)
     Nodes % y(1:4) = Solver % Mesh % Nodes % y(NodeIndexes)
     Nodes % z(1:4) = 0.0_dp       

     IP = GaussPoints( Element , np=INT(Particles % elementfraction) )

     DO j=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(j), IP % V(j), &
             IP % W(j),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

        No = ElemParticles(ind) % p(j)

        Particles % Coordinate(No,1) = SUM( Basis(1:4) * Nodes % x(1:4) ) !IP % U(j)
        Particles % Coordinate(No,2) = SUM( Basis(1:4) * Nodes % y(1:4) ) !IP % V(j)        
        Particles % PVolume(No) = IP % S(j) * detJ
        Particles % GVolume(No) = IP % S(j) * detJ

     END DO
  END DO

  DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)  

END SUBROUTINE UpdateLagParticleUpdates

!**************************************************************************

!> Updated Lagrangian mode: Update mesh (advect nodes)
SUBROUTINE UpdateLagMeshUpdates( Model, Solver, dt, TransientSimulation)

  USE MPMUtils
  USE DefUtils
  USE Lists
  IMPLICIT NONE

  TYPE(Solver_t), TARGET :: Solver
  TYPE(Element_t), POINTER :: Element  
  TYPE(Model_t) :: Model  
  REAL(KIND=dp) :: dt,x,y,dist
  LOGICAL :: TransientSimulation,UnFoundFatal=.TRUE.
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: Velvar,heightvar,prevvar
  REAL(KIND=dp), POINTER :: Velval(:),prevval(:)
  INTEGER, POINTER :: VelPerm(:),NodeIndexes(:),prevperm(:)
  INTEGER :: ii,dim,jj,t,nn
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName


  WRITE(SolverName, '(A)') 'Updated Lagrangian Mesh Deformation'

  Mesh => GetMesh()       

  dim = CoordinateSystemDimension()

  velvar => VariableGet( Mesh % Variables, 'SSAVelocity' )
  velval => velvar % values
  velperm => velvar % perm

  prevvar => VariableGet( Mesh % Variables, 'prevvel' )
  prevval => prevvar % values
  prevperm => prevvar % perm

  DO t=1,Solver % Mesh % NumberOfBulkElements
     Element => Solver % Mesh % Elements(t)
     Model % CurrentElement => Solver % Mesh % Elements(t)
     nn = GetElementNOFNodes(Element)
     NodeIndexes => Element % NodeIndexes        

     DO jj = 1,nn
        ii = NodeIndexes(jj)
        prevval(2*(prevperm(ii)-1)+1) = Mesh % Nodes % x(ii)
        prevval(2*(prevperm(ii)-1)+2) = Mesh % Nodes % y(ii)
     END DO
  END DO


  Do t=1,Solver % Mesh % NumberOfBulkElements
     Element => Solver % Mesh % Elements(t)
     Model % CurrentElement => Solver % Mesh % Elements(t)
     nn = GetElementNOFNodes(Element)
     NodeIndexes => Element % NodeIndexes

     DO jj = 1,nn
        ii = NodeIndexes(jj)

        Mesh % Nodes % x(ii) = prevval(2*(prevperm(ii)-1)+1) + &
             (velval(2*(velperm(ii)-1) + 1) ) * dt

        Mesh % Nodes % y(ii) = prevval(2*(prevperm(ii)-1)+2) + &
             (velval(2*(velperm(ii)-1) + 2)) * dt
     END DO
  END DO

END SUBROUTINE UpdateLagMeshUpdates


!**************************************************************************

!> Update particle thickness and mass
SUBROUTINE UpdateParticleHandMass( Model, Solver, dt, TransientSimulation)

  USE MPMUtils
  USE DefUtils

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(ValueList_t), POINTER :: Params  
  TYPE(Model_t) :: Model
  LOGICAL :: TransientSimulation  
  TYPE(Mesh_t), POINTER :: Mesh  
  REAL(KIND=dp) :: dt,divu,lc
  INTEGER :: No, ii,smoothiters
  LOGICAL :: Visited=.FALSE.,GotIt,nohupdate,smoothdam
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE :: Visited, nohupdate,lc,smoothdam,smoothiters

  Particles => GlobalParticles
  Mesh => GetMesh()
  WRITE(SolverName, '(A)') 'Update_Particle_H_and_Mass'

  CALL Info(SolverName,'Updating Particle H and Mass',Level=3) 

  Params => GetSolverParams()

  IF( .NOT. Visited ) THEN

     smoothdam =  GetLogical( Params, 'Smooth Full Dam Particles H', GotIt)
     IF (.NOT. GotIt) THEN
        smoothdam = .FALSE.
        CALL Info(SolverName,'Use smooth full dam particles not specified -- assuming false',Level=4)         
     END IF

     IF (smoothdam) THEN
        lc = GetCReal( Params, 'smooth dam lc', GotIt)
        IF (.NOT. GotIt) CALL Fatal(SolverName,&
             'Need to define "smooth dam lc=Real $lc"')

        smoothiters = GetInteger( Params,'smoothing iters',GotIt)
        IF (.NOT. GotIt) CALL Fatal(SolverName,&
             'Need to define "smooth iters = Integer"')        
     END IF

     nohupdate = ListGetLogical(Params,'no h update',GotIt)
     IF (.NOT. GotIt) THEN
        nohupdate = .FALSE.
     END IF

     Visited = .TRUE.

     IF (Particles % firsttimestepzero) RETURN
  END IF

  dt = Particles % dtime 

  IF( .NOT. nohupdate) THEN

     !regular update, no additional basal melt
     DO No = 1,Particles % NumberOfParticles

        divu = Particles % GradVel(No,1) + Particles % GradVel(No,2)

        Particles % H(No) = Particles % H(No) * (-1.0_dp) * &
             divu * dt  + Particles % H(No) + Particles % MB(No) * dt

        IF (Particles % H(No) < 1.0_dp) Particles % H(No) = 1.0_dp
     END DO
  END IF

  IF (Particles % ShapeFunctions == 'gimpm') THEN
     Particles % pvolume(:) = Particles % Length(:,1)*Particles % Length(:,2)
  END IF

  IF (smoothdam) THEN
     CALL Info(SolverName,'Smoothing thickness for fully damaged material points',Level=4)      
     CALL smoothrupth(Particles, Mesh, lc, smoothiters)
  END IF


  IF (Particles % hoop) THEN

     !just using Particles % xpic as temporary storage...
     Particles % xpic(:,1) = Particles % Coordinate(:,1)*Particles % Coordinate(:,1)
     Particles % xpic(:,2) = Particles % Coordinate(:,2)*Particles % Coordinate(:,2)
     Particles % xpic(:,1) = sqrt(Particles % xpic(:,1)+Particles % xpic(:,2))
     WHERE (Particles % xpic(:,1)<=70000.0_dp) Particles % H(:) = 400.0_dp

     Particles % xpic = 0.0_dp
  END IF


  Particles % mass(:) = Particles % pvolume(:) * Particles % H(:) * Particles % rhoi  

END SUBROUTINE UpdateParticleHandMass

!**************************************************************************

!> (Old?) solver to update particle H, Zs, grounded mask, Mass according to mesh
!! Where is this used? 
SUBROUTINE Update_Particle_H_GZs_GM_Mass_From_Mesh( Model, Solver, dt, TransientSimulation)

  USE MPMUtils

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(ValueList_t), POINTER :: Params  
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt,Hf
  LOGICAL :: TransientSimulation
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: No, ii
  LOGICAL :: Visited=.FALSE.,GotIt,nohupdate
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE :: visited, nohupdate

  Particles => GlobalParticles
  WRITE(SolverName, '(A)') 'Update_Particle_H_GZs_GM_and_Mass'


  !h and gmsk from mesh to particle
  CALL MPMMeshScalarToParticle(Particles, Model, 8)


  IF (Particles % SEP) THEN
     DO No = 1, Particles % NumberOfParticles

        IF (Particles % Gmask(No) > 0.99_dp) THEN
           Particles % Gmask(No) = 1.0_dp
        ELSE IF (Particles % Gmask(No) < -0.99_dp) THEN
           Particles % Gmask(No) = -1.0_dp
        ELSE
           Hf = Particles % rhow * &
                (Particles % sealevel-Particles % bedrock(No)) &
                /Particles % rhoi

           IF (Particles % H(No) .LT. Hf) THEN
              Particles % Gmask(No) = 1.0_dp
           ELSE
              Particles % GMask(No) = -1.0_dp
           END IF
        END IF

     END DO
  END IF

  Particles % mass(:) = Particles % pvolume(:) * Particles % H(:) * Particles % rhoi 

END SUBROUTINE Update_Particle_H_GZs_GM_Mass_From_Mesh


!**************************************************************************

!> Update mesh H, Binit, Velocity, etc from particle values
SUBROUTINE ParticlesToMesh( Model, Solver, dt, TransientSimulation)

  USE MPMUtils

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt,dirichletmax
  LOGICAL :: TransientSimulation,Visited=.FALSE.,GotIt,Test1D=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  TYPE(Variable_t), POINTER :: H,Hi,Vel1,Vel1i
  REAL(KIND=dp), POINTER :: HVal(:),HiVal(:),Vel1Val(:),Vel1iVal(:)
  INTEGER, POINTER :: HPerm(:),HiPerm(:),Vel1Perm(:),Vel1iPerm(:)
  INTEGER :: i

  SAVE :: Visited,SolverName,Test1D,H,HPerm,HVal,Hi,HiPerm,HiVal,Vel1,Vel1Perm,&
       Vel1Val,Vel1i,Vel1iPerm,Vel1iVal,dirichletmax

  Particles => GlobalParticles

  IF (.NOT. Visited) THEN

     WRITE(SolverName, '(A)') 'ParticleToMesh'

     Particles % weighth = ListGetLogical(Solver % Values,'weight h',GotIt)
     IF (.NOT. GotIt) Particles % weighth = .TRUE.

     PRINT *,''
     IF (Particles % weighth) THEN
        CALL Info(SolverName,'Weighting H',Level=3) 
     ELSE
        CALL Info(SolverName,'Not Weighting H',Level=3) 
     END IF
     PRINT *,''

     Test1D = ListGetLogical(Solver % Values,'Test1D',GotIt)
     IF (.NOT. GotIt) Test1D = .FALSE.     

     IF (Test1D) THEN

        H => VariableGet(Model % Mesh % Variables, 'H' )
        HPerm => H % Perm
        HVal => H % Values

        Hi => VariableGet(Model % Mesh % Variables, 'Hinit' )
        HiPerm => Hi % Perm
        HiVal => Hi % Values

        Vel1 => VariableGet(Model % Mesh % Variables, 'SSAVelocity 1' )
        Vel1Perm => Vel1 % Perm
        Vel1Val => Vel1 % Values

        Vel1i => VariableGet(Model % Mesh % Variables, 'InitVel 1' )
        Vel1iPerm => Vel1i % Perm
        Vel1iVal => Vel1i % Values

        dirichletmax = GetConstReal(Solver % Values,'dirichlet max x',GotIt)

     END IF

     Visited = .TRUE.
  END IF

  CALL Info(SolverName,'Updating Mesh H, Binit, Velocity',Level=3) 
  CALL MPMParticlesToNodes( Particles, Model, 2)
  CALL Info(SolverName,'Done Updating Mesh H, Binit, Velocity',Level=3)



  IF (Test1D) THEN

     Do i = 1,Model % Mesh % NumberOfNodes
        IF (Model % Mesh % Nodes % x(i) <= dirichletmax) THEN
           HVal(HPerm(i)) = HiVal(HiPerm(i))
           Vel1Val(Vel1Perm(i)) = Vel1iVal(Vel1iPerm(i))
        END IF

     END DO
  END IF

  !update mesh zs if static grounding line.  Otherwise, this is taken care of
  !by floatation solver
  !TODO: combine update mesh zs and floatation solver into one file
  IF (.NOT. Particles % movegl) THEN
     CALL Info(SolverName,'update mesh zs',Level=4)
     CALL UpdateMeshZs( Particles,Model)
  END IF

END SUBROUTINE ParticlesToMesh

!**************************************************************************

!> Update particle values from mesh values in preparation for SSA solution
SUBROUTINE MeshToParticles( Model, Solver, dt, TransientSimulation)

  USE MPMUtils

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt  
  LOGICAL :: TransientSimulation
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  Particles => GlobalParticles
  Params => ListGetSolverParams()
  WRITE(SolverName, '(A)') 'MeshToParticles'

  CALL SSAPrepMeshToParticles( Particles, Model)

  ! 3D viscosity from temperature, if needed:
  IF ((.NOT.  Particles % constlintemp) .AND. (.NOT. Particles % useconsttemp)) THEN
     CALL Info(SolverName,'interpolation of bz to ungrounding particles...',Level=1)
     CALL MPMMeshVectorToParticle(Particles, Model, 5,2 )     
     CALL Info(SolverName,'interpolation done',Level=1)
  END IF

END SUBROUTINE MeshToParticles

!**************************************************************************

!> Interpolates damage defined on mesh to the particle.
!! Not typically used.
SUBROUTINE MeshDamageToParticle(Model, Solver, dt, TransientSimulation)

  USE MPMUtils
  USE Lists

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt,divu
  TYPE(Element_t), POINTER :: Element  
  LOGICAL :: TransientSimulation
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: No, ii, elem
  LOGICAL :: Visited=.FALSE.,GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  TYPE(Variable_t), POINTER :: DVar  
  INTEGER, POINTER :: DPerm(:)
  REAL(KIND=dp), POINTER :: DValues(:)
  REAL(KIND=dp) :: Basis(4),dBasisdx(4,3),Coord(3)
  LOGICAL :: Stat
  REAL(KIND=dp) :: D,SqrtElementMetric

  Particles => GlobalParticles
  WRITE(SolverName, '(A)') 'MeshDamageToParticle'

  Particles % Dav = 0.0_dp
  DVar => VariableGet(Model % Mesh % Variables, 'damage' )
  DValues => DVar % Values
  DPerm => DVar % Perm  

  DO No = 1,Particles % NumberOfParticles
     Coord = GetParticleCoord(Particles,No)
     elem = GetParticleElement(Particles,No)
     Element => Model % Mesh % Elements(elem)
     stat = ParticleElementInfo( Element, Coord, &
          SqrtElementMetric, Basis, dBasisdx )
     CALL GetScalarFieldInMesh(DVar,Element,Basis,D)
     Particles % Dav(No,:) = D
  END DO

END SUBROUTINE MeshDamageToParticle

!**************************************************************************

!>
SUBROUTINE SaveParticleData( Model,Solver,dt,TransientSimulation )

  USE MPMUtils
  USE DefUtils
  USE MeshUtils
  USE SolverUtils
  USE Lists
  USE GeneralUtils

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles => NULL()
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt, SaveInterval,SaveStep,EndTime,AlwaysSaveTime,ymax
  LOGICAL :: TransientSimulation
  TYPE(ValueList_t), POINTER :: Params     
  LOGICAL ::  Found,UseInterval,GotIt,UseAlwaysSaveTime,savetimestep,&
       UseMismipFinalDamSave,mismipfinalsave
  INTEGER :: VisitedTimes = 0, OutputInterval,No

  CHARACTER(LEN=MAX_NAME_LEN) :: OutputDirectory,FileNamePrefix,FileName

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: currenttime
#else
  REAL(KIND=dp) :: currenttime
#endif

  SAVE :: OutputInterval,VisitedTimes,SaveStep,SaveInterval,UseInterval,&
       OutputDirectory,FileNamePrefix,FileName,EndTime,&
       UseAlwaysSaveTime,AlwaysSaveTime,ymax,UseMismipFinalDamSave

  Particles => GlobalParticles
  Params => ListGetSolverParams()

  IF (VisitedTimes == 0) THEN

     CALL Info('Save Particle Data',&
          '----------Setting Save Params----------',Level=4)
     OutputInterval = GetInteger( Params,'Output Interval',Found)
     SaveInterval = GetCReal( Params, 'Save Interval', Found)
     IF (.NOT. Found) CALL Fatal('SaveParticleData',&
          'Need to define "Save Interval k=Real $k"')


     OutputDirectory = GetString( Solver % Values,'Filename Directory')
     FileNamePrefix = GetString( Solver % Values,'Filename Prefix')

     FileName = TRIM(OutputDirectory)// '' //TRIM(FilenamePrefix)// '.result'     

     SaveStep = 0.0_dp

     UseInterval = GetLogical( Params, 'Use Output Interval', GotIt )
     IF (.NOT. GotIt) UseInterval = .FALSE.

     UseAlwaysSaveTime = GetLogical( Params, 'Use Always Save Time', GotIt )
     IF (.NOT. GotIt) THEN
        Call Warn('SaveParticleData',&
             'Did not specify "Use Always Save Time" in Params so setting to false')        
        UseAlwaysSaveTime = .FALSE.
     END IF

     !after a specified time, saving will occur every time step
     IF (UseAlwaysSaveTime) THEN
        AlwaysSaveTime = GetCReal( Params,'Always Save Time', GotIt)
        IF (.NOT. GotIt) THEN
           Call Warn('SaveParticleData',&
                'Did not specify "Always Save Time = Real" in Params so will not use!')        
           UseAlwaysSaveTime = .FALSE.
        END IF
     END IF

     !For MISMIP+ damage simulation: stop simulation and save when rift propagates across domain
     UseMismipFinalDamSave = GetLogical( Params, 'Use MISMIP Final Damage Save', GotIt )
     IF (.NOT. GotIt) THEN
        Call Warn('SaveParticleData',&
             'Did not specify "Use MISMIP Final Damage Save" in Params so setting to false')         
        UseMismipFinalDamSave = .FALSE.
     END IF

     IF (UseMismipFinalDamSave) THEN

        ymax = 40000.0_dp - (Particles % gridres/sqrt(Particles % elementfraction))

        PRINT *,''
        PRINT *,'Using ymax = ',ymax
        PRINT *,''
     END IF

     EndTime = GetCReal( Params,'StopTime', Found)
     IF (.NOT. Found) THEN
        Call Warn('SaveParticleData',&
             'Did not specify "StopTime= Real" in Params so setting to HUGE')
        EndTime = HUGE(1.0_dp)
     END IF

  END IF

  currenttime = RealTime()
  Particles % simtime = currenttime - Particles % firstsimtime

  PRINT *,''
  PRINT *,'Filename: ',FileName
  PRINT *,'Simulation Time (minutes): ',Particles % simtime/60.0_dp
  PRINT *,''

  mismipfinalsave = .FALSE.

  IF (UseMISMIPfinaldamsave) THEN
     DO No = 1,Particles % NumberOfParticles
        IF (Particles % Coordinate(No,2) < ymax) CYCLE
        IF (Particles % damstatus(No) >= 1) THEN
           mismipfinalsave = .TRUE.
           EXIT
        END IF
     END DO
  END IF

  IF (UseInterval) THEN

     IF (OutputInterval > 0 .AND. MOD(VisitedTimes,OutputInterval) == 0) THEN
        savetimestep = .TRUE.
     ELSEIF (UseAlwaysSaveTime .AND. Particles % time >= AlwaysSaveTime) THEN
        savetimestep = .TRUE.
     ELSEIF (mismipfinalsave) THEN
        savetimestep = .TRUE.
     ELSE
        savetimestep = .FALSE.
     END IF

     IF (savetimestep) THEN
        CALL Info('Save Particle Data',&
             '----------Saving Particle Data----------',Level=3)

        IF(.NOT. ListCheckPresent( Params,'Filename Prefix') ) THEN
           CALL ListAddString( Params,'Filename Prefix','particles')
        END IF
        CALL ParticleOutputVtu( Particles,Model )

     END IF

  ELSE

     IF (Particles % time >= SaveStep) THEN
        savetimestep = .TRUE.
     ELSEIF (UseAlwaysSaveTime .AND. Particles % time >= AlwaysSaveTime) THEN
        savetimestep = .TRUE.
     ELSEIF (mismipfinalsave) THEN
        savetimestep = .TRUE.        
     ELSE
        savetimestep = .FALSE.
     END IF

     IF (savetimestep) THEN
        CALL Info('Save Particle Data',&
             '----------Saving Particle Data----------',Level=3)

        IF(.NOT. ListCheckPresent( Params,'Filename Prefix') ) THEN
           CALL ListAddString( Params,'Filename Prefix','particles')
        END IF
        CALL ParticleOutputVtu( Particles,Model )

        IF (Particles % time >= SaveStep) THEN
           SaveStep = SaveStep+SaveInterval
        END IF
     END IF
  END IF

  IF (mismipfinalsave) THEN
     CALL FATAL ('final mismip damage saved','EXITING')
  END IF

  VisitedTimes = VisitedTimes + 1

  IF (Particles % time > EndTime) THEN

     !save last timestep if you didn't already
     IF (Particles % time < SaveStep) THEN
        IF(.NOT. ListCheckPresent( Params,'Filename Prefix') ) THEN
           CALL ListAddString( Params,'Filename Prefix','particles')
        END IF
        CALL ParticleOutputVtu( Particles,Model )
     END IF

     CALL FATAL('END TIME REACHED','EXITING')
  END IF

END SUBROUTINE SaveParticleData

!**************************************************************************

!> Option 1 for setting timestepping at start of computational cycle
!! Used when another solver (e.g. damage) computes the timestep during the previous cycle
FUNCTION MPMTimestep( Model ) RESULT( dt )
  USE Types
  USE Lists
  USE MPMUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: vtime
  REAL(KIND=dp) :: dt
  INTEGER :: VisitedTimes = 0,TimeIntervals
  TYPE(Solver_t), POINTER :: Solver
  REAL(KIND=dp), POINTER :: TimestepSizes(:,:)
  REAL(KIND=dp) :: dt0,prevdt
  LOGICAL :: GotIt
  TYPE(Particle_t), POINTER :: Particles  => NULL()   

  SAVE :: VisitedTimes,Mesh

  Particles => GlobalParticles

  IF (VisitedTimes == 0) THEN
     Mesh => GetMesh()
     vtime = VariableGet(Mesh % Variables, 'time')
     vtime % Values(1) = 0.0_dp

     Particles % time = 0.0_dp
     Particles % dtime = 0.0_dp
  END IF

  dt = Particles % dtime

  WRITE(Message,'(a,ES12.3)') 'MPM timestep before damage',dt
  CALL Info( 'MPM timestep',Message, Level=1 )
  CALL ListAddConstReal(Model % Simulation,'res: MPM timestep',dt)

  VisitedTimes = VisitedTimes + 1

END FUNCTION MPMTimestep

!**************************************************************************

!> Option 2 for setting timestepping at start of computational cycle
!! Used when no other solver computes the timestep
FUNCTION SSATimestep( Model ) RESULT( dt )

  USE Types
  USE Lists
  USE DefUtils
  USE MPMUtils  

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: vmax,vtime,VVar
  REAL(KIND=dp) :: dt,gridres,mult,velmax,steadytimestep,cflconstant
  INTEGER :: VisitedTimes = 0,TimeIntervals
  TYPE(Solver_t), POINTER :: Solver
  REAL(KIND=dp), POINTER :: TimestepSizes(:,:),VValues(:)
  INTEGER, POINTER :: VPerm(:)
  LOGICAL :: Found,usesteady
  TYPE(Particle_t), POINTER :: Particles  => NULL()   


  SAVE :: VisitedTimes,mult,usesteady,steadytimestep,gridres,cflconstant 

  Particles => GlobalParticles 

  IF (VisitedTimes == 0) THEN

     usesteady = GetLogical( Model % Constants, 'Use Steady Timestep', Found )
     IF (.NOT. Found) THEN
        usesteady = .FALSE.
     END IF

     IF (usesteady) THEN
        steadytimestep = GetConstReal( Model % Constants,'Steady Timestep',Found )
        IF (.NOT. Found) CALL Fatal('SSATimestep:', &
             'Need to define "gridres = Real $" in constants')
     END IF

     IF (.NOT. usesteady) THEN
        !  vtime % Values(1) = 0.0_dp
        gridres = GetConstReal( Model % Constants,'Grid Resolution',Found )
        IF (.NOT. Found) CALL Fatal('SSATimestep:', &
             'Need to define "gridres = Real $" in constants')

        mult = gridres/sqrt(2.0_dp)

        velmax = 0.0_dp

        cflconstant = GetConstReal( Model % Constants,'cfl constant',Found)
        IF (.NOT. Found) THEN
           cflconstant = 0.9_dp
           CALL Warn('SSATimestep:',&
                'Did not define "cfl constant = Real $" in constants so set to 0.9!')
        END IF

     END IF
  END IF

  IF (usesteady) THEN
     dt = steadytimestep
  ELSE

     !CFL
     dt = MAXVAL(ABS(Particles % Velocity(:,1)) + ABS(Particles % Velocity(:,2))) / gridres
     dt = cflconstant/dt
  END IF

  IF (VisitedTimes == 0 .AND. Particles % firsttimestepzero) THEN
     dt = 0.0_dp
     Particles % time = 0.0_dp
  END IF


  VisitedTimes = VisitedTimes + 1

  Particles % dtime = dt
  Particles % time = Particles % time + dt


  CALL Info('','',Level=1)
  CALL Info('','',Level=1)
  WRITE(Message,'(a,I10)') 'Timestep',VisitedTimes  
  CALL Info('SSATimestep',Message,Level=1)
  
  WRITE(Message,'(a,ES12.3)') 'dt',dt  
  CALL Info('SSATimestep',Message,Level=1)
  
  WRITE(Message,'(a,ES12.3)') 'Time',Particles % Time  
  CALL Info('SSATimestep',Message,Level=1)  
  CALL Info('','',Level=1)


  ! WRITE(Message,'(a,ES12.3)') 'ssa timestep',dt
  !CALL Info( 'ssa timestep',Message, Level=1 )
  CALL ListAddConstReal(Model % Simulation,'res: ssa timestep',dt)

  CALL ListAddConstReal(Model % Simulation,'res: MPM timestep',dt)  

END FUNCTION SSATimestep

!**************************************************************************

!> The 3-D SSA-MPM creep damage solver.
!! CAUTION: There are several experimental functions in this solver.
!! The primary damage parameters/controls can be set from the .sif, but
!! usually, some editing of this routine will be needed to get all the details worked out

SUBROUTINE UpdateCreepDamage( Model,Solver,dt,TransientSimulation )

  USE MPMUtils
  USE DefUtils
  USE ElementUtils
  USE SolverUtils
  USE Lists
  USE GeneralUtils

  IMPLICIT NONE

  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  TYPE(Element_t), POINTER :: BulkElement
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Valuelist_t), POINTER :: Params
  REAL(KIND=dp),allocatable :: dD(:),D(:),zref(:),ddscalevec(:)
  LOGICAL, allocatable :: layerdone(:)

  INTEGER :: Status,ElementIndex, numlayers, &
       NoParticles,dim
  INTEGER :: No, ii, ifail,jj,slowp,slowl,fastp,fastl,zhighi,zlowi
  INTEGER :: slowrkmsteps,fastrkmsteps,avrkmsteps,rkmcount
  INTEGER :: whichsurf, start,finish,step,ruptind,&
       layercount,ruptlayercount,subcritlayercount
  REAL(KIND=dp) :: Coord(3), GradVel(4),stress(4),&
       pSR(3,3),count,count2,TTT,DDD,maxeig,dtvec(3)
  REAL(KIND=dp) :: maxtimestep, maxgiventimestep, maxalloweddD,a,b, &
       gravity, tol, u, v, Dmax, n, MinSRInv, rhow, rhoi, Ap, &
       Peff,maxdD, stressthres,srthres,kparam,maxdDlayer, &
       EffH, MaxPDamage, criticaldamage,gaussk,lc, gridres, &
       av, newbf, ddmax, zs, z, normheight, H, Arr, velmag, &
       maxvel,divu,pddmax,Ee,RHS,depth,pw,pressure1,RHS2,savegamma,&
       maxdam,mpd,rupttime,subcrittime,savebf,Dlay(2,2),rupttol,currenttol,&
       slowrkmtime,fastrkmtime,totaltime,sthresmod,nbrPi,sigmath_var,&
       vertlc,loceps=1.0e-20_dp,ellipsesthres,defaultcriticaldamage,diff,maxp,&
       eigdmax,rcount,zhigh,zlow,EigValues(2),EigenVec(2,2)

  INTEGER :: NoVec(15),damnancount

  !efficiency variables
  REAL(KIND=dp) :: Eeexp,EFexp,zsRHS,MinSRInvSquared,rhowtimesgravity,rhoitimesgravity
  REAL(KIND=dp) :: nldDthres
  REAL(KIND=dp) :: oneovernumlayersminus1,preptime,rkmtime,zero=0.0_dp,one=1.0_dp,dist,dist2

  LOGICAL :: Visited = .FALSE., GotIt, RKM,&
       allowgrounded,useeffp,isotropic,skipisodam,rupt,restart,&
       defaultmpdonly,groundbasalwaterp,defaultnodzz,tuneellipse,&
       waterpforbasalonly,useellipse,justusegaussian,&
       useborstad,troubleshoot,usewp,usenonlocskip
  TYPE(Particle_t),  POINTER :: Particles
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,nonlocal
  LOGICAL :: zerostresscompare,suture
  REAL(KIND=dp) :: maxdDsave,targetdd,ddscale,nonlocskip,cflconstant

  REAL(KIND=dp) :: DS(2,2)
  LOGICAL :: larcenfmel  

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: starttime, endtime,time2s,time2e
#else
  REAL(KIND=dp) :: starttime, endtime,time2s,time2e
#endif  


  SAVE :: Visited, RKM, maxgiventimestep, maxalloweddD,&
       Gravity,numlayers,tol,n,MinSRInv,rhow,rhoi,D,dD,layerdone,&
       dim,stressthres,kparam,criticaldamage,gaussk,lc, gridres,&
       SolverName,nonlocal,allowgrounded,isotropic,skipisodam,&
       Eeexp,EFexp, zsRHS,MinSRInvSquared,rhowtimesgravity,rhoitimesgravity,&
       oneovernumlayersminus1,defaultnodzz,rupttol,groundbasalwaterp,&
       waterpforbasalonly,sthresmod,nbrPi,vertlc,useellipse,&
       zref,ellipsesthres,justusegaussian,tuneellipse,useborstad, &
       maxddsave,defaultcriticaldamage,troubleshoot,ddscalevec,&
       targetdd,nldDthres, zerostresscompare, &
       nonlocskip,usenonlocskip,cflconstant,suture,&
       larcenfmel 


  Params => GetSolverParams()
  Mesh => GetMesh()

  Particles => GlobalParticles

!!!!------------------------------------ FIRST TIME ONLY -------------------------------!!!!!!

  IF( .NOT. Visited ) THEN

     WRITE(SolverName, '(A)') 'UpdateCreepDamage'

     dim = Mesh % Meshdim

     cflconstant = GetCReal( Params, 'cfl constant', GotIt)
     IF (.NOT. GotIt) THEN
        cflconstant = 0.9_dp
     END IF


     nonlocal = ListGetString(Params,'nonlocal scheme',GotIt)
     IF (.NOT. GotIt) THEN
        CALL Fatal(SolverName,&
             'Need to define "nonlocal scheme=string "none","integral",or "gradient"')
     END IF

     !actually, only integral is available now...
     IF ((nonlocal == 'integral') .OR. (nonlocal == 'gradient')) THEN

        nldDthres = GetCReal( Params, 'Nonlocal dD rate threshold', GotIt)
        IF (.NOT. GotIt) THEN
           CALL Warn(SolverName,&
                'did not define Nonlocal dD rate threshold, setting to zero!')
           nldDthres = 0.0_dp
        ELSE
           PRINT *,''
           PRINT *,'nldDthres',nldDthres
           PRINT *,''
        END IF

        lc = GetCReal( Params, 'Nonlocal Regularization lc', GotIt)
        IF (.NOT. GotIt) CALL Fatal(SolverName,&
             'Need to define "Nonlocal Regularization lc=Real $lc"')

        gaussk = GetCReal( Params, 'Nonlocal Regularization k', GotIt)
        IF (.NOT. GotIt) CALL Fatal(SolverName,&
             'Need to define "Nonlocal Regularization k=Real $k"')

        IF (lc <= 0.0_dp) THEN
           lc = 0.0_dp
           gaussk = 0.0_dp
           nonlocal = 'none'
        ELSE
           CALL nonlocalsurroundelems(Particles,Mesh,lc)
        END IF

     ELSE
        lc = 0.0_dp
        gaussk = 0.0_dp
        nonlocal = 'none'
     END IF


     !feature not available in this version
     larcenfmel = .FALSE.

     !skip the nonlocal scheme if x-coord < some threshold
     !not recommended, generally, but can be useful to speed things up when debugging
     nonlocskip = GetCReal( Params, 'Nonlocal skip lower than x', GotIt )
     IF (.NOT. GotIt) CALL Warn(SolverName,&
          'didnt define Nonlocal skip lower than x=Real $"')
     IF (GotIt) THEN
        usenonlocskip = .TRUE.
     ELSE
        usenonlocskip = .FALSE.
     END IF

     maxgiventimestep = Particles % maximumtimestep
     maxalloweddD = GetCReal( Params, 'Maximum Allowed dDdt', GotIt )
     IF (.NOT. GotIt) CALL Fatal(SolverName,&
          'Need to define "Maximum Allowed dDdt=Real $maxdDdt"')

     IF (maxalloweddD == -9999.0_dp) THEN
        maxalloweddD = HUGE(1.0_dp)
     END IF

     targetdD = GetCReal( Params, 'Target dDdt', GotIt )
     IF (.NOT. GotIt) CALL Fatal(SolverName,&
          'Need to define "Target dDdt=Real $"')

     vertlc = GetCReal( Params, 'Vertical Regularization lc', GotIt )
     IF (.NOT. GotIt) THEN
        CALL WARN(SolverName,&
             'Vertical Regularization lc Not Defined, setting to 0!!')
        vertlc = 0
     END IF


     !Borstad 2016 damage model. Might have bugs.
     useborstad = GetLogical( Params,'Use Borstad',GotIt)

     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "use borstad = Logical" in Params so assuming false')
        useborstad = .FALSE.
     END IF

     IF (useborstad) THEN
        stressthres = GetConstReal( Model % Constants, 'sthresborstad', GotIt )
        IF (.NOT. GotIt) CALL Fatal(SolverName,&
             'Need to define "sthresborstad=Real $sthresborstad" in constants')
        kparam = GetConstReal( Model % Constants, 'kparam', GotIt )
        IF (.NOT. GotIt) CALL Fatal(SolverName,&
             'Need to define "kparam=Real $kparam" in constants')    
     END IF


     !Ellipse nonlocal scheme, similar to Giry 2011. 
     useellipse = GetLogical( Params,'Use Ellipse',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "use ellipse = Logical" in Params so assuming true')
        useellipse = .TRUE.
     END IF

     IF (useellipse) THEN

        !not used if 'Just use gaussian for ellipse = true'
        ellipsesthres = GetCReal( Params,'Ellipse Stress Threshold',GotIt)
        IF (.NOT. GotIt) THEN
           Call FATAL(SolverName,&
                'Need to specify "Ellipse Stress Threshold = Real" to use nonlocal ellipse ')
        END IF

        !skip the ellipse, just use a gaussian nonlocal kernal
        justusegaussian = GetLogical( Params,'Just Use Gaussian for ellipse',GotIt)
        IF (.NOT. GotIt) THEN
           Call Warn(SolverName,&
                'Didnt specify "Just Use Gaussian for ellipse = Logical", assuming true')
           justusegaussian = .TRUE.
        END IF
     END IF

     !NOTE: even if not useellipse is not specified, it
     !now overrides as the default nonlocal scheme, but with a
     !Gaussian kernal and no threshold. It is superior to the old scheme.
     !TODO: eliminate the old nonlocal scheme entirely
     IF (nonlocal == 'integral' .AND. (.NOT. useellipse)) THEN
        CALL WARN(SolverName,&
             'Ellipse nonloc scheme is always used now by default, but with a Gaussian kernal')
        useellipse = .TRUE.
        justusegaussian = .TRUE.
        ellipsesthres = 0.0_dp
     END IF


     !if true, no more damage accumulation on damage layers that have fully ruptured in one principal direction
     !However, spin contributions will still occur
     Particles % noevolveruptlayers = GetLogical( Params,'no evolve ruptured layers',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Didnt specify "no evolve ruptured layers = Logical", assuming false')
        Particles % noevolveruptlayers = .FALSE.
     END IF


     troubleshoot = GetLogical( Params,'troubleshoot',GotIt)
     IF (.NOT. GotIt) THEN
        troubleshoot = .FALSE.
     END IF

     !Runge-Kutta-Merson or Forward Euler?
     RKM = GetLogical( Params,'RKM',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "RKM = Logical" in Params so assuming RKM = TRUE')
        RKM = .TRUE.
     END IF

     waterpforbasalonly = GetLogical( Params,'Water Pressure for Basal Only',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "Water Pressure for Basal Only= Logical" in Params. Assuming false')
        waterpforbasalonly = .FALSE.
     END IF

     allowgrounded = GetLogical( Params,'Allow Grounded Damage',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "Allow Grounded Damage= Logical" in Params so assuming shelf damage only')
        allowgrounded = .FALSE.
     END IF

     groundbasalwaterp = .FALSE.

     IF (allowgrounded) THEN
        groundbasalwaterp = GetLogical( Params,'Allow Grounded Basal Water Pressure',GotIt)
        IF (.NOT. GotIt) THEN
           Call Warn(SolverName,&
                'Did not specify "Allow Grounded Basal Water Pressure= Logical" in Params so assuming false')
           groundbasalwaterp = .FALSE.
        END IF
     END IF

     skipisodam = GetLogical( Params,'Skip IsoDam Particles',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "Skip IsoDam Particles= Logical" in Params so false')
        skipisodam = .FALSE.
     END IF



     CriticalDamage = Particles % CriticalDamage
     MinSRInv = Particles % criticalshearrate
     rhow = Particles % rhow
     rhoi = Particles % rhoi
     !particles % viscosityexponent is typically 1/3
     n = Particles % ViscosityExponent
     n = 1.0_dp/n
     gridres = Particles % gridres
     Gravity = ABS(Particles % Gravity)
     numlayers = Particles % numberofparticlelayers

     !used to set random noise on damage threshold as in Krug 2014
     sthresmod = GetConstReal( Params, 'Stress Threshold Modifier', GotIt )
     IF (.NOT. GotIt) CALL Fatal(SolverName,&
          'Need to define "Stress Threshold Modifier=Real $" in constants')

     nbrPi = 3.141592_dp

     IF (RKM) THEN
        dt = GetConstReal( Model % Constants, 'First rkm dt', GotIt )
        IF (.NOT. GotIt) CALL Fatal(SolverName,&
             'Need to define "first rkm dt=Real $" in constants')
        Particles % dtime = dt
        Particles % nextdtime = dt
     END IF

     tol = GetCReal( Params,'Damage convergence tolerance',GotIt)
     IF (.NOT. GotIt) CALL Fatal(SolverName, &
          'Need to define "Damage convergence tolerance = Real $Dtol"')

     rupttol = GetCReal( Params,'Ruptured Damage convergence tolerance',GotIt)
     IF (.NOT. GotIt) CALL Fatal(SolverName, &
          'Need to define "Ruptured Damage convergence tolerance = Real $Dtol"')     

     ALLOCATE( layerdone(numlayers) )
     ALLOCATE(zref(numlayers))
     ALLOCATE( dD(4) )
     ALLOCATE( D(4) )
     ALLOCATE( dDscalevec(numlayers) )

     a=0.0_dp

     IF (Particles % Gamma == 0.0_dp) THEN
        isotropic = .TRUE.
     ELSE
        isotropic = .FALSE.
     END IF

     Particles % stresshigh = -HUGE(1.0_dp)
     Particles % stresslow = HUGE(1.0_dp)

     !various variable shortcuts for damage calculations
     EeExp =  (1.0_dp-n)/(2.0_dp * n)
     EFexp = -1.0_dp/n
     zsRHS = (1.0_dp - rhoi/rhow)
     MinSRInvSquared = MinSRInv*MinSRInv
     rhowtimesgravity = rhow*gravity
     rhoitimesgravity = rhoi*gravity
     oneovernumlayersminus1 = 1.0_dp/DBLE(numlayers-1)

     DO ii = 1,numlayers
        zref(ii) = -one + DBLE(ii-1)*oneovernumlayersminus1
     END DO

     defaultcriticaldamage = Particles % criticaldamage
     defaultnodzz = Particles % nodzz


     Visited = .TRUE.

  END IF

!!!!------------------------------------ END FIRST TIME ONLY -------------------------------!!!!!

  CALL Info(SolverName,'Updating Particle Damage',Level=3)  

  starttime = RealTime()

  maxtimestep = cflconstant * (1.0_dp/ MAXVAL( ABS(Particles % Velocity(:,1))/Particles % gridres &
       + ABS(Particles % Velocity(:,2))/Particles % gridres))

  IF (RKM) THEN
     !if using Runge-Kutta-Merson, dt for the current timestep (Particles % dtime)
     !was determined during the previous timestep (Particles % nextdtime).
     !If this dt results in too much change in damage
     !over the current timestep, we restart the damage procedure
     !using 0.5*dt
     Particles % dtime = MIN(Particles % nextdtime,maxtimestep)
  END IF

  layercount = 0

  dt = Particles % dtime
  PRINT *,'Damage dt',dt

  !for Runge-kutta-merson
  a = 0.0_dp
  b = dt

  Particles % xpic(:,1) = Particles % Velocity(:,1)*Particles % Velocity(:,1) &
       + Particles % Velocity(:,2) * Particles % Velocity(:,2)

  maxvel = MAXVAL(Particles % xpic(:,1))

  Particles % xpic = zero


  !----------SKIPPING--------
  !note whether a particle should be skipped or not
  !on xpic(No,6)

  NoParticles = Particles % NumberOfParticles


  IF (maxalloweddD > zero) THEN

     DO No = 1,NoParticles
        !----------- can this particle be skipped? -------------!

        IF (Particles % damstatus(No) == 1 .AND. Particles % useriftdmax) CYCLE        

        IF (.NOT. allowgrounded) THEN
           IF (Particles % Gmask(No) < 0.0_dp) CYCLE
        END IF

        IF (Particles % nodamregion) THEN
           IF ( (Particles % Coordinate(No,1) < Particles % ndxmax) .AND. & 
                (Particles % Coordinate(No,1) > Particles % ndxmin) ) THEN
              IF ( (Particles % Coordinate(No,2) < Particles % ndymax) .AND. & 
                   (Particles % Coordinate(No,2) > Particles % ndymin) ) THEN
                 Particles % Damage(No,:,:) = 0.0_dp
                 Particles % Dav(No,:) = 0.0_dp
                 CYCLE
              END IF
           END IF
        END IF

        IF (Particles % restrictdam) THEN
           IF (Particles % Coordinate(No,1) < Particles % rdxmin) CYCLE
           IF (Particles % Coordinate(No,1) > Particles % rdxmax) CYCLE       
           IF (Particles % Coordinate(No,2) < Particles % rdymin) CYCLE
           IF (Particles % Coordinate(No,2) > Particles % rdymax) CYCLE  
        END IF


        IF ((Particles % Status(No) == PARTICLE_LOST) &
             .OR. (Particles % Status(No) == PARTICLE_ALLOCATED)) CYCLE

        IF (skipisodam) THEN
           IF (Particles % damstatus(No) == -1) CYCLE
        END IF

        Particles % xpic(No,6) = one

        !------------------- end skip particles ------------!
     END DO
  END IF



100 NoParticles = Particles % NumberOfParticles
  Particles % dD = zero
  maxdD = zero
  maxvel = zero
  maxdDlayer = zero
  Particles % currentgamma = Particles % gamma


  !Particles % xpic is used for the XPIC node to particle
  !mapping scheme, but here it is used to track which particles
  !are damaged. This speeds up the nonlocal scheme.
  Particles % xpic(:,1:5) = zero

  count = zero
  count2 = zero

  !particle strain rate
  pSR = 0.0_dp

  IF (UseBorstad) THEN
     DO No = 1,NoParticles
        IF (Particles % xpic(No,6) == zero) CYCLE        
        srthres = (stressthres/Particles % Binit(No))**n
        CALL BorstadDamage(Particles, No, srthres, Particles % GradVel(No,:),Particles % dD(No,1,1),n,kparam)

        maxdD = MAX(Particles % dD(No,1,1),maxdD)
     END DO

     Particles % dav(:,2) = Particles % dav(:,1)
     Particles % dav(:,3) = Particles % dav(:,1)
     Particles % dav(:,4) = 0.0_dp
     Particles % dD = 0.0_dp
     Particles % xpic = 0.0_dp

     PRINT *,'Borstad maxdD',maxdD

     RETURN
  END IF


  Particles % xpic(:,3) = 0.0_dp
  Particles % xpic(:,4) = numlayers

  !debugging 
  Particles % equaleigcount = 0

  IF (maxalloweddD > zero) THEN 

     DO No = 1, NoParticles

        IF (Particles % xpic(No,6) == zero) CYCLE


        !Particles % damstatus(No) == -1 occurs when
        !"Use Isotropic Damage for Initially Damaged Particles = Logical True"
        !is denoted in the SIF (under Constants)
        !This is for particles that are initially (isotropically) damaged at the
        !start of the simulation

        !There are 2 options on how to treat these particles:
        !1) do not allow any additional accumulation of damage.
        !   This option is activated using "Skip IsoDam Particles = Logical"
        !   in the damage solver section, and is denoted "skipisodam" below
        !2) default: allow additional damage accumulation, but only isotropically

        IF (Particles % useisodam) THEN
           IF (skipisodam) THEN
              IF (Particles % damstatus(No) == -1) CYCLE
           ELSE
              IF (Particles % damstatus(No) == -1) THEN
                 Particles % currentgamma = zero
                 Particles % nodzz = .FALSE.
              ELSE
                 Particles % currentgamma = Particles % gamma
                 Particles % nodzz = defaultnodzz
              END IF
           END IF
        END IF

        !if a particle is fully damaged and Particles % dmaxII > Particles % dmaxI,
        !we simply evolve the damage using the spin tensor.  The particle will also be
        !exempt from the nonlocal scheme and any restriction on the magnitude of
        !dD allowed over a timestep
        ! TODO: check that this still works, given recent changes??
        IF (Particles % dmaxII_dom .AND. Particles % damstatus(No)==1) THEN
           D(1:4) = Particles % Dav(No,1:4)

           Particles % dvdxmdudy = Particles % GradVel(No,4) - Particles % Gradvel(No,3)

           CALL runge_kutta_merson_fulldam_dmaxIIdom(Particles,D,dD,a,b,rupttol,ifail)

           IF (ifail == 1) THEN
              CALL Fatal(SolverName,'(dmaxII_dom) check rkm input params')
           ELSEIF (ifail == 2) THEN
              PRINT *,''
              CALL Warn(SolverName,'(dmaxII_dom) RKM: VERY SMALL STEP TAKEN!')
              PRINT *,'No',No
              PRINT *,'Layer from bottom',ii
              PRINT *,'Dav',Particles % Dav(No,:)
           END IF

           Particles % dD(No,:,1) = dD(1)
           Particles % dD(No,:,2) = dD(2)
           Particles % dD(No,:,3) = dD(3)
           Particles % dD(No,:,4) = dD(4)
           CYCLE
        END IF


        !----------------- Partially damaged particles ---------------!

        !particle strain rates
        pSR(1,1) = Particles % GradVel(No,1)
        pSR(2,1) = 0.5_dp*(Particles % Gradvel(No,3) + Particles % GradVel(No,4))
        pSR(1,2) = pSR(2,1)
        pSR(2,2) = Particles % GradVel(No,2)
        pSR(3,3) = -pSR(1,1)-pSR(2,2)
        Particles % dvdxmdudy = Particles % GradVel(No,4) - Particles % Gradvel(No,3)

        !2nd invariant of strain rate tensor
        Ee = 0.5_dp*(pSR(1,1)*pSR(1,1) + pSR(2,2)*pSR(2,2) + &
             pSR(3,3)*pSR(3,3)) + pSR(1,2)*pSR(1,2)


        !instead of passing these a bunch of times for RKM, we just define them
        !globally
        Particles % psr = pSR

        IF (Ee < MinSRInvSquared) THEN
           Ee = MinSRInvSquared
        END IF

        !RHS of Glen's flow law except (effective) strain-rate tensor (SR)
        !and temperature flow factor (Bt)
        !dev_stress = Bt*RHS*SR

        !EeExp =  (1.0_dp-n)/(2.0_dp * n)
        !EFexp = -1.0_dp/n

        RHS = (Ee**EeExp) * (Particles % EF(No)**EFexp )         


        IF (Particles % Gmask(No) < 0.0_dp) THEN
           zs = Particles % H(No)+Particles % Bedrock(No)
        ELSE
           !particle surface height
           !zsRHS = (1.0_dp - rhoi/rhow)             
           zs = Particles % H(No)*zsRHS
        END IF


        !----------------------------------------------------------------------!
        !---------         calculate change in damage                ----------!
        !----------------------------------------------------------------------!        

        !These calculations are done for each vertical layer of the particle.
        !This code currently assumes damage should nucleate at the surface or base
        !of the shelf and propagate towards the center.  A layer will not damage
        !unless all other layers above or below the layer are damaged.

        !Calculations are performed on each layer starting from the basal layer 
        !and looping towards the surface layer.  However, if one of these layers
        !accumulates no damage, the loop is restarted, but now working from
        !surface layer towards the bottom layer.  This surface-to-base loop ends
        !when a layer is reached that either:
        !-accumulates no damage 
        !-has already been processed from the basal-to-surface loop
        !This way, we don't have to process the middle layers that won't accumulate
        !damage, resulting in significantly faster code (especially if most of the
        !domain is undamaged).

        IF (Particles % damstatus(No) == 1) THEN
           currenttol = rupttol
        ELSE
           currenttol = tol
        END IF

200     layercount = 0

        !the vector layerdone = TRUE where D has been calculated, FALSE where not
        layerdone(:) = .FALSE.


        Particles % currentno = No


        DO whichsurf = 1,2

           IF (whichsurf == 1) THEN
              !"basal crevasses"

              !use basal water pressure term
              usewp = .TRUE.

              IF (Particles % Gmask(No) < zero .AND. (.NOT. groundbasalwaterp)) CYCLE

              start = 1
              finish= numlayers
              step = 1
           ELSE

              !"surface crevasses"

              IF (waterpforbasalonly) THEN
                 usewp = .FALSE.
              END IF

              !starting from top
              start = numlayers
              finish = 2
              step = -1
           END IF


           ddscalevec = 1.0_dp

           DO ii = start,finish,step

              IF ( layerdone(ii) ) EXIT
              layercount = layercount+1
              layerdone(ii) = .TRUE.

              dD = zero
              D(1:4) = Particles % Damage(No,ii,1:4)

              !---- get water pressure and depth ----
              !numlayersminus1 = numlayers-1
              !  z = zs - Particles % H(No) + (ii-1)*(Particles % H(No)/numlayersminus1)
              ! z = zs + Particles % H(No) * (-one + DBLE(ii-1)*oneovernumlayersminus1)
              z = zs + Particles % H(No) * zref(ii)


              !this is converted to Peff within dDdt subroutine
              !Peff = pressure1 - (Tau(1,1) + Tau(2,2))
              !rhoitimesgravity = rhoi*gravity
              !pressure1 = (rhoitimesgravity * (zs-z)) - Pw

              IF (z>=0.0_dp .OR. (Particles % GMask(No)<zero .AND. (.NOT. groundbasalwaterp))) THEN
                 !NO BASAL WATER PRESSURE
                 !  Pw = 0.0_dp
                 Particles % pressure1 = rhoitimesgravity * (zs-z)
              ELSE
                 !BASAL WATER PRESSURE
                 !Pw = rhowtimesgravity * (0.0_dp-z)
                 !pressure1 = (rhoitimesgravity * (zs-z)) - Pw

                 IF (usewp) THEN
                    Particles % pressure1 = rhoitimesgravity*(zs-z) + rhowtimesgravity*z
                 ELSE
                    Particles % pressure1 = rhoitimesgravity * (zs-z)
                 END IF
              END IF


              !Tau = RHS2*(effective strain rate)
              Particles % RHS = RHS*Particles % Bz(No,ii)

              Particles % currentlayer = ii

              !introduce some noise on the stress threshold:
              u = EvenRandom()
              v = EvenRandom()     
              sigmath_var = ABS(0+sthresmod*SQRT(-2.0_dp*LOG((one-u)))*COS(2.0_dp*nbrPi*v))
              Particles % sthresmod = Particles % sthres * (one+sigmath_var)



              IF (RKM) THEN

                 IF (Particles % noevolveruptlayers) THEN
                    !check if any principal damage component is > criticaldamage. Then it is rupt.
                    !and we just evolve via spin.
                    TTT = D(1)+D(2)
                    DDD = D(1)*D(2)-D(4)*D(4)
                    maxeig = 0.5_dp*TTT+sqrt(0.25_dp*TTT*TTT-DDD)

                    IF (maxeig > Particles % criticaldamage) THEN
                       Particles % prupt = .TRUE.
                    ELSE
                       Particles % prupt = .FALSE.
                    END IF

                    CALL runge_kutta_merson(Particles,D,dD,a,b,currenttol,ifail,ddscale)
                    Particles % prupt = .FALSE.

                 ELSE
                    CALL runge_kutta_merson(Particles,D,dD,a,b,currenttol,ifail,ddscale)
                    Particles % prupt = .FALSE.
                 END IF

                 Particles % tempgamma = Particles % gamma
                 Particles % currentgamma = Particles % gamma

                 IF (ifail == 1) THEN
                    CALL Fatal(SolverName,'check rkm input params')
                 ELSEIF (ifail == 2) THEN
                    PRINT *,''
                    CALL Warn(SolverName,'RKM: VERY SMALL STEP TAKEN!')
                    PRINT *,'No',No
                    PRINT *,'Layer from bottom',ii
                    PRINT *,'Dav',Particles % Dav(No,:)
                    PRINT *,'Dold',Particles % Damage(No,ii,:)
                    PRINT *,'Dnew',D+dD
                    PRINT *,'dd',dd
                 END IF

                 IF (waterpforbasalonly .AND. whichsurf == 1) THEN

                    TTT = D(1)+dD(1)+D(2)+dD(2)
                    DDD = (D(1)+dD(1))*(D(2)+dD(2))-(D(4)+dD(4))*(D(4)+dD(4))
                    eigdmax = 0.5_dp*TTT+sqrt(0.25_dp*TTT*TTT-DDD)                    

                    usewp = .TRUE.
                 END IF

              ELSE
                 !explicit forward euler.  Works fine for isotropic, though
                 !RKM may still be preferred. Doesn't work will with
                 !anisotropic, due to spin contributions

                 CALL dDdtfast(D,dD)

                 IF (.NOT. Isotropic) THEN
                    dD = D+dD
                    !dmax is now the max value of dmaxI,dmaxII,dmaxIII
                    CALL FixPrincipalDamage(dD,Particles % Dmax)
                    dD = dD-D
                    Particles % prupt = .FALSE.
                 END IF
              END IF

              IF (ALL(dD .EQ. zero) .AND. ALL(D .EQ. zero) ) EXIT

              Particles % dD(No,ii,1:4) = dD(1:4)
           END DO
           !done with layer loop from a surface

           IF (useellipse) THEN
              !save on xpic which layers from the bottom and top the local damage
              !propagates to.  This allows the nonlocal scheme to skip
              !nondamaged layers.
              IF (whichsurf == 1) THEN
                 Particles % xpic(No,3) = ii
              ELSE
                 Particles % xpic(No,4) = ii
              END IF
           END IF
        END DO
        !done with surface loop


        !either restrict the timestep based on the change in vertically integrated damage
        !or the change in damage of any layer EDIT: always dvert int  now

        IF (Particles % damstatus(No)<1 ) THEN 

           CALL GetMaxdDPrincipalDamageVert(Particles,No,numlayers,pddmax)

           IF (pddmax > maxdD) THEN
              maxp = No
              maxdD = pddmax
           END IF

        END IF


        IF (RKM) THEN

           IF (maxdD > maxalloweddd) THEN
              maxdD = 0.0_dp
              Particles % dtime = Particles % dtime/1.5_dp
              restart=.TRUE.

              dt = Particles % dtime 
              PRINT *,'Damage dt',dt
              PRINT *,'Restart due to No',No

              IF (dt < 1.0e-20_dp) THEN
                 !really small timestep...something has probably gone wrong.
                 PRINT *,'No',No
                 PRINT *,'dD 1',Particles % dD(No,:,1)
                 PRINT *,'dD 2',Particles % dD(No,:,2)
                 PRINT *,'dD 3',Particles % dD(No,:,3)
                 PRINT *,'dD 4',Particles % dD(No,:,4)
                 PRINT *,'D 1',Particles % Damage(No,:,1)
                 PRINT *,'D 2',Particles % Damage(No,:,2)
                 PRINT *,'D 3',Particles % Damage(No,:,3)
                 PRINT *,'D 4',Particles % Damage(No,:,4)
                 PRINT *,'Dav',Particles % Dav(No,:)
                 PRINT *,'damstatus',Particles % damstatus(No)
                 CALL FATAL('end','end')
              END IF

              a = zero
              b = dt
              Particles % dD(No,:,:) = zero
              go to 200
           END IF

           IF (restart) THEN
              restart = .FALSE.
              go to 100
           END IF
        END IF

        !track which particles are damaged
        IF (ANY(Particles % dD(No,:,:)>zero)) THEN

           count2 = count2 + one

           IF (ANY(Particles % dD(No,:,:)>(nldDthres*dt))) THEN
              count = count + one
              Particles % xpic(No,1) = count
           END IF
        ELSE
           IF (useellipse) THEN
              Particles % xpic(No,3) = 0.0_dp
              Particles % xpic(No,4) = numlayers
           END IF
        END IF


     END DO

  END IF


  Particles % currentgamma = Particles % gamma

  !------Set new timestep and (if RKM) an initial guess for the next timestep ------- !


  maxtimestep = MIN(maxtimestep,maxgiventimestep)

  IF (RKM) THEN

     !set current timestep and set a possible timestep for next time

     Particles % Time = Particles % Time + Particles % dtime


     dtvec(1) = Particles % dtime * 1.8_dp
     dtvec(2) = Particles % dtime * targetdD/maxdD
     dtvec(3) = maxtimestep
     Particles % nextdtime = MINVAL(dtvec)


     IF (Particles % nextdtime == maxtimestep) THEN
        PRINT *,'timestep restricted by velocity'
     END IF

  ELSE

     maxtimestep = (1.0_dp/(MAXVAL(ABS(Particles % Velocity(:,1)))/gridres + &
          MAXVAL(ABS(Particles % Velocity(:,2)))/gridres))

     maxtimestep = MIN(maxtimestep,maxgiventimestep)
     !forward euler can set the timestep on the fly

     IF (maxdD > 0.0_dp) THEN
        dt = targetdD/maxdD
     ELSE
        dt = maxtimestep
     END IF

     dt = MIN(dt,maxtimestep)

     Particles % dTime = dt
     Particles % Time = Particles % Time + dt
  END IF

  maxdDlayer = MAXVAL(Particles % dD)

  PRINT *,''
  PRINT *,'MAX dD layer over current timestep: ',maxdDlayer
  PRINT *,''
  PRINT *,'MAX drate layer over current timestep: ',maxdDlayer/dt
  PRINT *,''
  PRINT *,'MAX dD over current timestep: ',maxdD
  PRINT *,''
  PRINT *,'MAX particle: ',maxp
  PRINT *,' '
  PRINT *,'new dt:',dt
  PRINT *,' '
  PRINT *,'new Time:', Particles % time

  IF (RKM) THEN
     PRINT *,''
     PRINT *,'next dt:', Particles % nextdtime
  END IF


  !--------------- nonlocal scheme --------------!

  IF (lc>0.0_dp ) THEN
     IF (Nonlocal == 'integral') THEN

        IF (useellipse) THEN
           IF (usenonlocskip) THEN
              count = zero
              DO No = 1,NoParticles
                 IF (Particles % coordinate(No,1)<nonlocskip) THEN
                    Particles % xpic(No,1) = 0.0_dp
                 ELSE
                    IF (Particles % xpic(No,1) > 0.0_dp) THEN
                       count = count + 1
                    END IF
                 END IF
              END DO
           END IF
           !5 and 6 will be the furthest layers that basal/surface
           !damage propagates to of a particle and its surrounding
           !particles.  This is used to restrict processing to only
           !the damaged layers.

           !now in subroutine
           !Particles % xpic(:,5:6) = Particles % xpic(:,3:4)

           CALL Info(SolverName,'Start Nonlocal Damage',Level=3)
           time2s = RealTime()

           PRINT *,''
           PRINT *,'num particles with dD > 0',count2
           PRINT *,'num particles with dD > dDthres',count
           PRINT *,'diff',count2-count
           PRINT *,''

           CALL nonlocalintegraldDellipseRobust(Particles, numlayers, &
                INT(count), lc, gaussk, gridres,vertlc,groundbasalwaterp,&
                ellipsesthres,justusegaussian,Mesh)

           time2e = RealTime()-time2s
           WRITE( Message, '(a,f8.1,a)') 'Nonlocal Damage Finished in ', time2e, ' seconds'
           CALL Info(SolverName,Message,Level=3)
        ELSE
           !This is the old nonlocal scheme, which should no longer be called
           !by default
           !TODO: clean this up
           
           CALL Info(SolverName,'Start Nonlocal Damage',Level=3)   
           CALL nonlocalintegraldD(Particles, numlayers, &
                INT(count), lc, gaussk, gridres,vertlc)
           CALL Info(SolverName,'Nonlocal Damage Done',Level=3)
        END IF

     ELSE IF (nonlocal == 'gradient') THEN
        CALL Fatal(SolverName,'nonlocal implicit gradient no longer available.')
     END IF
  END IF



!!! --------------  UPDATE DAMAGE --------------- !!

  time2s = RealTime()
  damnancount = 0
  Particles % prupt = .FALSE.

  DO No = 1, NoParticles

     Particles % currentno = No

     IF (Particles % damstatus(No) == 1 .AND. Particles % useriftdmax) CYCLE


     IF (ANY(Particles % dD(No,:,:) .NE. 0.0_dp)) THEN

        IF (skipisodam) THEN
           IF (Particles % damstatus(No) == -1) CYCLE
        END IF

        IF (Particles % Static(No)) THEN
           Particles % Dav(No,:) = 0.0_dp
           Particles % damage(No,:,:) = 0.0_dp
        END IF


        DO ii = 1,numlayers

           Particles % currentlayer = ii

           IF (ANY(Particles % dD(No,:,:) .NE. 0.0_dp)) THEN

              IF (isotropic .AND. Particles % damstatus(No) .NE. -1) THEN

                 !ISOTROPIC (gamma == 0)

                 Particles % Damage(No,ii,1:3) = Particles % Damage(No,ii,1) + &
                      Particles % dD(No,ii,1)

                 IF (Particles % Damage(No,ii,1) < Particles % mindam ) THEN
                    Particles % Damage(No,ii,:) = 0.0_dp
                 END IF

                 IF (Particles % Damage(No,ii,1) > Particles % CriticalDamage) THEN
                    Particles % Damage(No,ii,1:2) = Particles % DmaxI

                    IF (.NOT. Particles % nodzz) THEN
                       Particles % Damage(No,ii,3) = Particles % DmaxIII
                    END IF


                    !this shouldn't happen
                    IF (Particles % Damage(No,ii,1) < 0.0_dp) THEN
                       Particles % Damage(No,ii,1:3) = 0.0_dp
                    END IF
                 END IF
              ELSE
                 !ANISOTROPIC (gamma > 0)

                 D(:) = Particles % Damage(No,ii,:) + Particles % dD(No,ii,:)

                 !makes sure principal damage components of the layer are bounded
                 !by 0 and DMax, and will set a component to DMax if
                 !it exceeds CriticalDamage.

                 IF (Particles % DamStatus(No) == -1) THEN
                    CALL FixPrincipalDamage(D,Particles % isodamcritdam)
                 ELSE
                    CALL FixPrincipalDamage(D,CriticalDamage)
                 END IF

                 Particles % prupt = .FALSE.

                 IF (ALL(D<Particles % mindam)) D = 0.0_dp

                 IF ( D(1) .NE. D(1) ) THEN
                    PRINT *,''
                    PRINT *,'No',No
                    PRINT *,'Particles % dd(No,ii,:)',Particles % dd(No,ii,:)
                    PRINT *,'Particles % damage(No,ii,:)',Particles % damage(No,ii,:)
                    CALL FATAL(SolverName,'nan detected on a layer')
                 END IF

                 Particles % Damage(No,ii,:) = D
              END IF
           END IF
        END DO


        Particles % currentlayer = -1

        !vertically integrated damage for DAv
        CALL VertIntDamFromVisc(Particles, No, numlayers,Model)

        IF (Particles % dmaxII_dom .AND. Particles % damstatus(No) == 1) CYCLE

        IF (Particles % damStatus(No) == -1) THEN
           !material points where damage is assigned, but not evolved

           IF (ANY(Particles % Dav(No,:)>Particles % isodamcritdav)) THEN

              Particles % Dav(No,1:2) = Particles % initdmax !DMaxI
              Particles % Dav(No,3) = Particles % initdmax !DmaxIII
              Particles % Dav(No,4) = 0.0_dp

              Particles % Damage(No,:,1:2) = Particles % DMaxI
              Particles % Damage(No,:,3) = Particles % DmaxIII
              Particles % Damage(No,:,4) = 0.0_dp
           END IF


        ELSE IF (Particles % Gamma == 0.0_dp) THEN
           !isotropic damage
           !Determine if Dav should rupture

           IF (ANY(Particles % Dav(No,:)>Particles % criticaldav)) THEN

              Particles % Dav(No,1:2) = Particles % riftdmax
              !Particles % DMaxI
              Particles % Dav(No,3) = Particles % riftdmax
              !Particles % DmaxIII
              Particles % Dav(No,4) = 0.0_dp

              Particles % Damage(No,:,1:2) = Particles % DMaxI
              Particles % Damage(No,:,3) = Particles % DmaxIII
              Particles % Damage(No,:,4) = 0.0_dp
              Particles % damstatus(No) = 1
           END IF

        ELSE
           !1. get the principal damage directions for DAv
           !2. get the principal damage values for each layer, and rotate back to cartesian
           !   reference frame based on the principal damage eigenvector of DAv.
           !   This guarantees that the principal damage values of each particle layer
           !   are pointed exactly in the same direction.  This eliminates the very slight
           !   errors (e.g. ~1.0E-8 per timestep?) in principal directions between layers that
           !   eventually adds up to affect the solution. For example, in the fully anisotropic
           !   case, no damage should develop in any principal direction except the direction of the
           !   max principal damage, which this scheme guarantees.
           !3. vertically integrate damage again.

           IF (troubleshoot) THEN
              IF (Particles % origno(No) == 111694) THEN
                 PRINT *,'===FIRST==='
                 CALL CheckPrincipalDamage(No)
              END IF
           END IF

           IF (Particles % fixdavdirections) THEN
              ! Experimental...reorients damage tensors on all layers of a material point so that
              ! they are aligned with their average orientation

              !also calls fixprincipaldamageverint if ruptured
              CALL FixPrincipalDavAndLayers(Particles,No,numlayers,Model)

              IF (troubleshoot) THEN
                 IF (Particles % origno(No) == 111694) THEN
                    PRINT *,'===SECOND==='
                    CALL CheckPrincipalDamage(No)
                 END IF
              END IF

           ELSE
              !DEFAULT ROUTINE:
              !processes ruptured components according to the settings in the sif
              !if any rupt, sets Particles % damstatus(No) = 1
              CALL FixPrincipalDamageVertInt(No,Model)

           END IF
        END IF
     END IF

     IF (ANY(Particles % Dav(No,:) .NE. Particles % Dav(No,:))) THEN

        damnancount = damnancount + 1
        IF (damnancount == 1) THEN
           PRINT *,'No',No
           PRINT *,'Dav(No,:)',Particles % Dav(No,:)
           PRINT *,'damage(No,:,1)',Particles % damage(No,:,1)
           PRINT *,'damage(No,:,2)',Particles % damage(No,:,2)
           PRINT *,'damage(No,:,3)',Particles % damage(No,:,3)
           PRINT *,'damage(No,:,4)',Particles % damage(No,:,4)
        END IF
        Particles % Dav(No,:) = 0.0_dp
        WHERE (Particles % Damage(No,:,:) .NE. Particles % Damage(No,:,:)) &
             Particles % Damage(No,:,:)  = 0.0_dp

        CALL VertIntDamFromVisc(Particles, No, numlayers,Model)     
     END IF

  END DO


  time2e = RealTime()

  PRINT *,'damage update time',time2e-time2s


  IF (damnancount > 0) THEN
     PRINT *,''
     PRINT *,'WARNING: DAMAGE NAN COUNT',damnancount
     PRINT *,''
  END IF


  ! Particles % mpdonly = defaultmpdonly
  Particles % currentgamma = Particles % gamma

  IF (lc>0.0_dp ) THEN
     ddmax = MAXVAL(ABS(Particles % dD(:,:,:)))     
     PRINT *,'Max particle layer dD over current timestep after nonlocal: ', ddmax
     PRINT *,' '
  END IF

  Particles % xpic = 0.0_dp

  PRINT *,''
  PRINT *,'Stress high',Particles % stresshigh
  PRINT *,'Stress low', Particles % stresslow

  !debugging
  PRINT *,'equal eig count',Particles % equaleigcount


  endtime = RealTime()
  PRINT *,''
  PRINT *,'DAMAGE SOLUTION TIME',endtime-starttime

  CALL Info(SolverName,'----- Damage Update Done -----',Level=3)

END SUBROUTINE UpdateCreepDamage

!**************************************************************************

!> The necking/mass balance modification (Bassis and Ma, 2015) for the
!! zero stress damage model (Sun et al., 2017).
!! In addition, if 'zero stress only = Logical True', then this solver just serves as
!! timestep control for the non-modified zero-stress damage model

SUBROUTINE UpdateDamageModifiedZeroStress( Model,Solver,dt,TransientSimulation )

  USE MPMUtils
  USE DefUtils
  USE ElementUtils
  USE SolverUtils
  USE Lists
  USE GeneralUtils

  IMPLICIT NONE

  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  TYPE(Element_t), POINTER :: BulkElement
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Valuelist_t), POINTER :: Params
  REAL(KIND=dp),allocatable :: zref(:)
  REAL(KIND=dp) :: dD(4),D(4)
  LOGICAL, allocatable :: layerdone(:)

  INTEGER :: Status,ElementIndex, numlayers, &
       NoParticles,dim
  INTEGER :: No, ii, ifail,jj,slowp,slowl,fastp,fastl
  INTEGER :: slowrkmsteps,fastrkmsteps,avrkmsteps,rkmcount
  INTEGER :: whichsurf, start,finish,step,ruptind,&
       layercount,ruptlayercount,subcritlayercount
  REAL(KIND=dp) :: Coord(3), GradVel(4),stress(4),&
       pSR(3,3),psrrot(2,2),count,TTT,DDD,maxeig,dtvec(3),rr(4),DNew(4)
  REAL(KIND=dp) :: maxtimestep, maxgiventimestep, maxalloweddD,a,b, &
       gravity, tol, u, v, Dmax, n, MinSRInv, rhow, rhoi, Ap, &
       Peff,maxdD, stressthres,srthres,kparam,maxdDlayer, &
       EffH, MaxPDamage, criticaldamage,gaussk,lc, gridres, &
       av, newbf, ddmax, zs, z, normheight, H, Arr, velmag, &
       maxvel,divu,pddmax,pddmax1,pddmin,Ee,RHS,depth,pw,pressure1,RHS2,savegamma,&
       maxdam,mpd,rupttime,subcrittime,savebf,Dlay(2,2),rupttol,currenttol,&
       slowrkmtime,fastrkmtime,totaltime,sthresmod,nbrPi,sigmath_var,&
       vertlc,loceps=1.0e-20_dp,defaultcriticaldamage,diff,maxp,&
       cflconst,maxvelx,maxvely,cfltimestep,temph,mb

  REAL(KIND=dp) w,x,y,Df(2,2)

  INTEGER :: NoVec(15),damnancount

  !efficiency variables
  REAL(KIND=dp) :: Eeexp,EFexp,zsRHS,MinSRInvSquared,rhowtimesgravity,rhoitimesgravity
  REAL(KIND=dp) :: oneovernumlayersminus1,preptime,rkmtime,zero=0.0_dp,one=1.0_dp


  LOGICAL :: Visited = .FALSE., GotIt, RKM,&
       allowgrounded,useeffp,isotropic,skipisodam,rupt,restart,&
       defaultmpdonly,groundbasalwaterp,defaultnodzz,&
       waterpforbasalonly,&
       useborstad,troubleshoot,usemelt
  TYPE(Particle_t),  POINTER :: Particles
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  LOGICAL :: linkwithzerostress,cflonly,zsonly
  REAL(KIND=dp) :: maxdDsave,targetdd,crevmelt,bmb

  REAL(KIND=dp) :: mu0,S0,EigVal(2),EigVec(2,2),beta,E1,E2,nstarnum,&
       nstardenom,nstar,rhos,S0gen

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: starttime, endtime,time2s,time2e
#else
  REAL(KIND=dp) :: starttime, endtime,time2s,time2e
#endif  


  SAVE :: Visited, RKM, maxgiventimestep, maxalloweddD,&
       Gravity,numlayers,tol,n,MinSRInv,rhow,rhoi,D,dD,layerdone,&
       dim,stressthres,kparam,criticaldamage,gridres,&
       SolverName,allowgrounded,skipisodam,&
       Eeexp,EFexp, zsRHS,MinSRInvSquared,rhowtimesgravity,rhoitimesgravity,&
       oneovernumlayersminus1,defaultnodzz,rupttol,&
       nbrPi,vertlc,zref,maxddsave,defaultcriticaldamage,&
       targetdd,S0gen,rhos,linkwithzerostress,cflconst,cflonly,usemelt,crevmelt,zsonly


  Params => GetSolverParams()
  Mesh => GetMesh()

  Particles => GlobalParticles

!!!!------------------------------------ FIRST TIME ONLY -------------------------------!!!!!!

  IF( .NOT. Visited ) THEN


     WRITE(SolverName, '(A)') 'UpdateDamageModifiedZeroStress'

     dim = Mesh % Meshdim

     CFLconst = GetCReal( Params, 'CFL Constant', GotIt )
     IF (.NOT. GotIt) THEN
        CALL WARN(SolverName,&
             'didnt define "CFL Constant=Real", setting to 0.9')
        cflconst = 0.9_dp
     END IF

     maxgiventimestep = Particles % maximumtimestep

     maxalloweddD = GetCReal( Params, 'Maximum Allowed dDdt', GotIt )
     IF (.NOT. GotIt) CALL Fatal(SolverName,&
          'Need to define "Maximum Allowed dDdt=Real $maxdDdt"')

     IF (maxalloweddD == -9999.0_dp) THEN
        maxalloweddD = HUGE(1.0_dp)
     END IF

     targetdD = GetCReal( Params, 'Target dDdt', GotIt )
     IF (.NOT. GotIt) CALL Fatal(SolverName,&
          'Need to define "Target dDdt=Real $"')


     cflonly = GetLogical( Params, 'cflonly',GotIt)
     IF (.NOT. GotIt) cflonly=.FALSE.

     zsonly = GetLogical( Params, 'zero stress only',GotIt)
     IF (.NOT. GotIt) zsonly=.FALSE.
     IF (zsonly) PRINT *,'USING ZERO STRESS ONLY!'

     vertlc = GetCReal( Params, 'Vertical Regularization lc', GotIt )
     IF (.NOT. GotIt) THEN
        CALL WARN(SolverName,&
             'Vertical Regularization lc Not Defined, setting to 0!!')
        vertlc = 0
     END IF

     !add 1 m/yr melt for floating particles?
     usemelt = GetLogical( Params,'Use melt for floating',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'not using melt for floating')
        usemelt = .FALSE.
     END IF

     IF (usemelt) THEN
        PRINT *,'USING MELT FOR FLOATING'
        crevmelt = GetCReal( Params,'Crevasse melt rate',GotIt)
        IF (GotIt) THEN
           PRINT *,'crevasse melt rate',crevmelt
        ELSE
           Call Fatal(SolverName,&
                'Did not specify "crevasse melt rate = real $" in Params!!')
        END IF
     END IF


     linkwithzerostress = GetLogical( Params,'link with zero stress',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "link with zero stress = Logical", assume true')
        linkwithzerostress = .TRUE.
     END IF

     !Runge-Kutta-Merson or Forward Euler?
     RKM = GetLogical( Params,'RKM',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "RKM = Logical" in Params so assuming RKM = TRUE')
        RKM = .TRUE.
     END IF

     allowgrounded = GetLogical( Params,'Allow Grounded Damage',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "Allow Grounded Damage= Logical" in Params so assuming shelf damage only')
        allowgrounded = .FALSE.
     END IF

     skipisodam = GetLogical( Params,'Skip IsoDam Particles',GotIt)
     IF (.NOT. GotIt) THEN
        Call Warn(SolverName,&
             'Did not specify "Skip IsoDam Particles= Logical" in Params so false')
        skipisodam = .FALSE.
     END IF


     CriticalDamage = Particles % CriticalDav 

     MinSRInv = Particles % criticalshearrate
     rhow = Particles % rhow
     rhoi = Particles % rhoi

     !particles % viscosityexponent is 1/3
     n = Particles % ViscosityExponent
     n = 1.0_dp/n
     gridres = Particles % gridres
     Gravity = ABS(Particles % Gravity)
     numlayers = Particles % numberofparticlelayers

     nbrPi = 3.141592_dp

     dt = GetConstReal( Model % Constants, 'First rkm dt', GotIt )
     IF (.NOT. GotIt) CALL Fatal(SolverName,&
          'Need to define "first rkm dt=Real $" in constants')
     Particles % dtime = dt
     Particles % nextdtime = dt

     tol = GetCReal( Params,'Damage convergence tolerance',GotIt)
     IF (.NOT. GotIt) CALL Fatal(SolverName, &
          'Need to define "Damage convergence tolerance = Real $Dtol"')

     rhos = GetCReal( Params, 'Snow Density', GotIt )
     IF (.NOT. GotIt) CALL Fatal(SolverName,&
          'Need to define "Snow Density=Real $"')


     S0gen = (rhoi-rhos)*(rhow-rhoi)*gravity/(4.0_dp*(rhow-rhos))


     !various variable shortcuts for damage calculations
     EeExp =  (1.0_dp-n)/(2.0_dp * n)
     EFexp = -1.0_dp/n
     MinSRInvSquared = MinSRInv*MinSRInv

     defaultcriticaldamage = Particles % criticaldav !mage

     DO No = 1,Particles % NumberOfParticles
        IF (Particles % damstatus(No) == -1) THEN
           IF (ANY(Particles % Dav(No,:)>Particles % isodamcritdav)) THEN
              Particles % Dav(No,1:3) = Particles % Dmax
              Particles % Dav(No,4) = 0.0_dp
              Particles % damstatus(No) = 1
              CYCLE
           END IF
        ELSE
           IF (Particles % gamma > 0.0_dp) THEN
              !check if any principal damage component is > criticaldamage. Then it is rupt.
              !and we just evolve via spin.
              TTT = Particles % Dav(No,1)+Particles % Dav(No,2)
              DDD = Particles % Dav(No,1)*Particles % Dav(No,2)-&
                   Particles % Dav(No,4)*Particles % Dav(No,4)
              maxeig = 0.5_dp*TTT+sqrt(0.25_dp*TTT*TTT-DDD)
           ELSE
              maxeig = Particles % Dav(No,1)
           END IF

           IF (maxeig >= Particles % criticaldav) THEN
              ! Particles % prupt = .TRUE.
              Particles % damstatus(No) = 1 
              Particles % Dav(No,1:3) = Particles % riftdmax
              Particles % Dav(No,4) = 0.0_dp
              Particles % damstatus(No) = 1
           END IF
        END IF
     END DO

     Visited = .TRUE.     
  END IF

!!!!------------------------------------ END FIRST TIME ONLY -------------------------------!!!!!

  CALL Info(SolverName,'Updating Particle Damage Modified Zero Stress',Level=3)  

  starttime = RealTime()

  !dt for the current timestep (Particles % dtime)
  !was determined during the previous timestep (Particles % nextdtime).
  !If this dt results in too much change in damage
  !over the current timestep, we restart the damage procedure
  !using 0.5*dt

  IF (cflonly) THEN
     Particles % zsmaxdd = 0.0_dp
  END IF


  !you link with zero stress, so need to adjust dtime accordingly,
  !hoping to keep zero stress damage in check.
  IF (linkwithzerostress) THEN
     IF (Particles % zsmaxdd > maxalloweddD) THEN
        dt = MIN(Particles % nextdtime, Particles % dtime * maxalloweddD/Particles % zsmaxdd)
     ELSE
        dt = Particles % nextdtime
     END IF
  ELSE
     dt = Particles % nextdtime
  END IF

  IF (Particles % nextdtime > dt) THEN
     Particles % dtime = dt
  ELSE
     Particles % dtime = Particles % nextdtime
  END IF

  dt = Particles % dtime
  PRINT *,'Damage dt',dt

  !for Runge-kutta-merson
  a = 0.0_dp
  b = dt

  maxvel = MAXVAL(Particles % Velocity(:,1)*Particles % Velocity(:,1) &
       + Particles % Velocity(:,2) * Particles % Velocity(:,2))

  !used for storage, has nothing to do with xpic here
  !TODO: define an alias
  Particles % xpic = zero


  NoParticles = Particles % NumberOfParticles



  !----------- mark particles to skip -------------!
  IF (maxalloweddD > zero) THEN

     DO No = 1,NoParticles

        IF (.NOT. allowgrounded) THEN
           IF (Particles % Gmask(No) < 0.0_dp) CYCLE
        END IF

        IF (Particles % nodamregion) THEN
           IF ( (Particles % Coordinate(No,1) < Particles % ndxmax) .AND. & 
                (Particles % Coordinate(No,1) > Particles % ndxmin) ) THEN
              IF ( (Particles % Coordinate(No,2) < Particles % ndymax) .AND. & 
                   (Particles % Coordinate(No,2) > Particles % ndymin) ) THEN
                 Particles % Damage(No,:,:) = 0.0_dp
                 Particles % Dav(No,:) = 0.0_dp
                 CYCLE
              END IF
           END IF
        END IF

        IF (Particles % restrictdam) THEN
           IF (Particles % Coordinate(No,1) < Particles % rdxmin) CYCLE
           IF (Particles % Coordinate(No,1) > Particles % rdxmax) CYCLE       
           IF (Particles % Coordinate(No,2) < Particles % rdymin) CYCLE
           IF (Particles % Coordinate(No,2) > Particles % rdymax) CYCLE  
        END IF

        IF ((Particles % Status(No) == PARTICLE_LOST) &
             .OR. (Particles % Status(No) == PARTICLE_ALLOCATED)) CYCLE


        IF (Particles % DamStatus(No) == 1) CYCLE

        Particles % xpic(No,6) = one
     END DO
  END IF
  !------------------- end marking particles to skip ------------!


100 NoParticles = Particles % NumberOfParticles
  Particles % dD = zero
  maxdD = zero

  Particles % currentgamma = Particles % gamma

  !Particles % xpic is used for the XPIC node to particle
  !mapping scheme, but here it is used to track which particles
  !are damaged. This speeds up the nonlocal scheme.
  Particles % xpic(:,1:5) = zero

  count = zero

  !particle strain rate
  pSR = 0.0_dp


  IF (maxalloweddD > zero) THEN 

     DO No = 1, NoParticles

        ! IF (Particles % xpic(No,6) == zero) CYCLE

        IF (Particles % useisodam) THEN
           IF (skipisodam) THEN
              IF (Particles % damstatus(No) == -1) CYCLE
           ELSE
              IF (Particles % damstatus(No) == -1) THEN
                 Particles % currentgamma = zero
                 Particles % nodzz = .FALSE.
              ELSE
                 Particles % currentgamma = Particles % gamma
                 Particles % nodzz = defaultnodzz
              END IF
           END IF
        END IF

        !----------------------------------------------------------------------!
        !---------         calculate change in damage                ----------!
        !----------------------------------------------------------------------!        

        !particle strain rates
        pSR(1,1) = Particles % GradVel(No,1)
        pSR(2,1) = 0.5_dp*(Particles % Gradvel(No,3) + Particles % GradVel(No,4))
        pSR(1,2) = pSR(2,1)
        pSR(2,2) = Particles % GradVel(No,2)
        pSR(3,3) = -pSR(1,1)-pSR(2,2)

        divu = Particles % GradVel(No,1) + Particles % GradVel(No,2)
        temph = Particles % H(No)


        Particles % dvdxmdudy = Particles % GradVel(No,4) - Particles % Gradvel(No,3)

        !2nd invariant of strain rate tensor
        Ee = 0.5_dp*(pSR(1,1)*pSR(1,1) + pSR(2,2)*pSR(2,2) + &
             pSR(3,3)*pSR(3,3)) + pSR(1,2)*pSR(1,2)

        IF (Ee < MinSRInvSquared) THEN
           Ee = MinSRInvSquared
        END IF

        !the 0.5 is needed due to how jeremy defines B,mu, and tau
200     mu0 = 0.5_dp*Particles % Binit(No) * (Ee**EeExp) * (Particles % EF(No)**EFexp )


        IF (Particles % gamma < 1.0_dp) THEN
           !principal horizontal strain rates
           CALL Eigen2DSym_TryGenFirst(pSR(1:2,1:2),EigVal,EigVec)

           E1 = EigVal(2)
           E2 = EigVal(1)

           Particles % bmd = Particles % Dav(No,1)

           !  mu0 = mu0*( 1.0_dp-Particles % Dav(No,1) )
        ELSE
           !get the strain rates coaxial with the max principal damage
           Df(1,1) = Particles % Dav(No,1)
           Df(2,2) = Particles % Dav(No,2)
           Df(1,2) = Particles % Dav(No,4)
           CALL Eigen2DSym_TryGenFirst(Df,EigVal,EigVec)


           !mu0 = mu0 * (1.0_dp-EigVal(2))
           Particles % bmd = EigVal(2)

           Df = MATMUL(MATMUL(TRANSPOSE(EigVec),Df),EigVec)

           psrrot = MATMUL(MATMUL(TRANSPOSE(EigVec),psr(1:2,1:2)),EigVec)

           IF (Df(1,1)>Df(2,2)) THEN
              E1 = psrrot(1,1)
              E2 = psrrot(2,2)
           ELSE
              E2 = psrrot(1,1)
              E1 = psrrot(2,2)
           END IF
        END IF

        mb = Particles % mb(No)
        bmb = Particles % mb(No)

        !add 1 m/yr melt to floating particles?
        IF (usemelt .AND. Particles % Gmask(No) >= 0.0_dp) THEN
           !  Particles % MB(No) = Particles % MB(No)-5.0_dp !1.0_dp
           bmb = bmb-crevmelt
        END IF

        IF (E1 == 0.0_dp) THEN
           E1 = TINY(1.0_dp)
        END IF

        IF (mu0 == 0.0_dp) THEN
           mu0 = TINY(1.0_dp)
        END IF

        beta = E2/E1

        nstarnum = 4.0_dp*(1.0_dp + beta + beta*beta)
        nstardenom = nstarnum + 3.0_dp*(n-1.0_dp)*beta*beta
        nstar = nstarnum*n/nstardenom

        IF (.NOT. RKM) THEN
           !  mu0 = mu0*( 1.0_dp-Particles % Dav(No,1) )
           mu0 = mu0*( 1.0_dp-Particles % bmd)
           S0 = S0gen*Particles % H(No) / (mu0*E1)
           !the bassis and ma 2015 sign convention is negative mb is accumulation, so
           !swap signs to be consistent with our convention (melting is negative)
           RHS = nstar*(1.0_dp-S0)*E1 - bmb/Particles % H(No)
        ELSE
           !bmd
           S0 = S0gen / (mu0*E1)
           RHS = nstar*E1
        END IF


        IF (zsonly) RHS = 0.0_dp        

        dD = zero
        D(1:4) = Particles % Dav(No,1:4)

        Particles % currentno = No


        IF (Particles % gamma > 0.0_dp) THEN
           !check if any principal damage component is > criticaldamage. Then it is rupt.
           !and we just evolve via spin.
           TTT = D(1)+D(2)
           DDD = D(1)*D(2)-D(4)*D(4)
           maxeig = 0.5_dp*TTT+sqrt(MAX(0.25_dp*TTT*TTT-DDD,0.0_dp))

           IF (maxeig >= Particles % criticaldav) THEN
              Particles % prupt = .TRUE.
           ELSE
              Particles % prupt = .FALSE.
           END IF
        ELSE
           IF (D(1) >= Particles % criticaldav) THEN
              Particles % prupt = .TRUE.
           ELSE
              Particles % prupt = .FALSE.
           END IF
        END IF


        IF (RKM) THEN
           CALL runge_kutta_merson_bassis(Particles,D,dD,a,b,tol,RHS,S0,ifail,temph,divu,bmb,bmb)

           IF (ifail == 1) THEN
              CALL Fatal(SolverName,'check rkm input params')
           ELSEIF (ifail == 2) THEN
              PRINT *,''
              CALL Warn(SolverName,'RKM: VERY SMALL STEP TAKEN!')
              PRINT *,'No',No
              PRINT *,'Layer from bottom',ii
              PRINT *,'Dav',Particles % Dav(No,:)
              PRINT *,'Dold',Particles % Damage(No,ii,:)
              PRINT *,'Dnew',D+dD
              PRINT *,'dd',dd
           END IF

           !h update
           Particles % dD(No,2,1) = tempH
        ELSE
           CALL bassisinc(Particles,D,RHS,rr,dd)
           dD=dD*dt
           rr=rr*dt
           Dnew = D+dD !*dt
           CALL FixPrincipalDamageInc(Dnew,dD,rr)
           dD = Dnew-D

           !h update
           Particles % dD(No,2,1) = temph*(1.0_dp-divu*dt) +bmb*dt
           IF (Particles % dD(No,2,1) < 1.0_dp) Particles % dD(No,2,1) = 1.0_dp
        END IF

        Particles % prupt = .FALSE.        

        Particles % dD(No,1,1:4) = dD(1:4)

        IF (Particles % damstatus(No)<1 ) THEN
           IF (Particles % Gamma > 0.0_dp) THEN
              TTT = dD(1)+dD(2)
              DDD = dD(1)*dD(2)-dD(4)*dD(4)
              pddmax = 0.5_dp*TTT+sqrt(0.25_dp*TTT*TTT-DDD)
              pddmin = 0.5_dp*TTT-sqrt(0.25_dp*TTT*TTT-DDD)

              pddmax1 = MAX(ABS(pddmax),ABS(pddmin))

              IF (ABS(pddmax) >= ABS(pddmin)) THEN
                 pddmax1 = pddmax
              ELSE
                 pddmax1 = pddmin
              END IF
           ELSE
              pddmax1 = dD(1)
           END IF

           IF (ABS(pddmax1) > ABS(maxdD)) THEN
              maxp = No
              maxdD = pddmax1
           END IF
        END IF

        !if using melt parameterization, reset for particle
        ! IF (usemelt .AND. Particles % Gmask(No) >= 0.0_dp) THEN
        !    Particles % MB(No) = Particles % MB(No)+meltmb
        !5          !.0_dp
        !  END IF        

        IF (ABS(maxdD) > maxalloweddd) THEN
           maxdD = 0.0_dp
           Particles % dtime = Particles % dtime/1.5_dp
           restart=.TRUE.

           dt = Particles % dtime 
           PRINT *,'Damage dt',dt
           PRINT *,'Restart due to No',No

           IF (dt < 1.0e-20_dp) THEN
              PRINT *,'No',No
              PRINT *,'dD 1',Particles % dD(No,1,1)
              PRINT *,'dD 2',Particles % dD(No,1,2)
              PRINT *,'dD 3',Particles % dD(No,1,3)
              PRINT *,'dD 4',Particles % dD(No,1,4)
              PRINT *,'D 1',Particles % Damage(No,:,1)
              PRINT *,'D 2',Particles % Damage(No,:,2)
              PRINT *,'D 3',Particles % Damage(No,:,3)
              PRINT *,'D 4',Particles % Damage(No,:,4)
              PRINT *,'Dav',Particles % Dav(No,:)
              PRINT *,'damstatus',Particles % damstatus(No)
              CALL FATAL('end','end')
           END IF

           a = zero
           b = dt
           Particles % dD(No,:,:) = zero
           go to 200
        END IF

        IF (restart) THEN
           restart = .FALSE.
           go to 100
        END IF

        !track which particles are damaged
        IF (ANY(Particles % dD(No,1,:)>zero)) THEN
           !  count2 = count2 + one
           ! IF (ANY(Particles % dD(No,:.:)>dDthres)) THEN
           count = count + one
           Particles % xpic(No,1) = count
           ! END IF
        END IF
     END DO
  END IF

  !------Set new timestep and an initial guess for the next timestep ------- !

  !maxtimestep = MIN(maxgiventimestep,((gridres/sqrt(2.0_dp))/(2.0_dp * sqrt(maxvel))))

  ! maxvelx = MAXVAL(Particles % Velocity(:,1))
  !  maxvely = MAXVAL(Particles % Velocity(:,2))

  !cfltimestep = cflconst * (1.0_dp/ (maxvelx/gridres + maxvely/gridres))

  cfltimestep = cflconst * (1.0_dp/ MAXVAL( ABS(Particles % Velocity(:,1))/gridres &
       + ABS(Particles % Velocity(:,2))/gridres))

  maxtimestep = MIN(maxgiventimestep,cfltimestep)


  Particles % Time = Particles % Time + Particles % dtime

  dtvec(1) = Particles % dtime * 1.8_dp
  dtvec(2) = Particles % dtime * targetdD/ABS(maxdD)
  dtvec(3) = maxtimestep

  Particles % nextdtime = MINVAL(dtvec)

  IF (cflonly) THEN
     Particles % nextdtime = MIN(dtvec(2),cfltimestep)
     ! Particles % nextdtime = cfltimestep
  END IF


  IF (Particles % nextdtime == maxtimestep) THEN
     PRINT *,'timestep restricted by velocity'
  END IF

  PRINT *,'MAX dD over current timestep: ',maxdD
  PRINT *,''
  PRINT *,'MAX particle: ',maxp
  PRINT *,' '
  PRINT *,'new dt:',dt
  PRINT *,' '
  PRINT *,'new Time:', Particles % time
  PRINT *,''
  PRINT *,'next dt:', Particles % nextdtime


!!! --------------  UPDATE DAMAGE --------------- !!

  DO No = 1, NoParticles

     Particles % H(No) = Particles % dd(No,2,1)

     IF (Particles % damstatus(No) == 1) THEN
        Particles % damage(No,3,1) = Particles % Dav(No,1)
        Particles % dbassis(No) = 0.0_dp
        CYCLE
     END IF

     IF (skipisodam) THEN
        IF (Particles % damstatus(No) == -1) CYCLE
     END IF

     IF (Particles % Static(No)) THEN
        Particles % Dav(No,:) = 0.0_dp
        Particles % damage(No,:,:) = 0.0_dp
     END IF

     IF (Particles % damstatus(No) < 1) THEN

        Particles % Dav(No,:) = Particles % Dav(No,:) + &
             Particles % dD(No,1,:)

        IF (Particles % Gamma < 1.0_dp) THEN

           IF (Particles % Dav(No,1) >= Particles % criticaldav) THEN
              Particles % Dav(No,1:3) = Particles % riftDmax
              Particles % Damstatus(No) = 1
           ELSE
              Particles % Damstatus(No) = 0
           END IF

           IF (Particles % Dav(No,1) < 0.0_dp) THEN
              Particles % Dav(No,1:3) = 0.0_dp
           END IF

           Particles % Dav(No,4) = 0.0_dp

           Particles % damage(No,3,1) = Particles % Dav(No,1)

           IF (Particles % outputdbassis) THEN
              Particles % dbassis(No) = Particles % dD(No,1,1)/dt
           END IF

        ELSE

           IF (Particles % outputdbassis) THEN
              Df(1,1) = Particles % dD(No,1,1)
              Df(2,2) = Particles % dD(No,1,2)
              Df(1,2) = Particles % dD(No,1,4)
              Df(2,1) = Df(1,2)
              CALL Eigen2DSym_TryGenFirst(Df,EigVal,EigVec) 
              Particles % dbassis(No) = EigVal(2)/dt
           END IF

           Df(1,1) = Particles % Dav(No,1)
           Df(2,2) = Particles % Dav(No,2)
           Df(1,2) = Particles % Dav(No,4)
           Df(2,1) = Df(1,2)

           CALL Eigen2DSym_TryGenFirst(Df,EigVal,EigVec)        

           IF (EigVal(2) > Particles % riftDmax) EigVal(2) = Particles % Riftdmax
           IF (EigVal(2) < 0.0_dp) EigVal(2) = 0.0_dp
           EigVal(1) = 0.0_dp


           !rotate back
           w = EigVal(1)*EigVec(1,1)
           x = EigVal(2)*EigVec(1,2)
           y = EigVal(1)*EigVec(2,1) 
           z = EigVal(2)*EigVec(2,2)

           Particles % Dav(No,1) = EigVec(1,1)*w + EigVec(1,2)*x
           Particles % Dav(No,2) = EigVec(2,1)*y + EigVec(2,2)*z
           Particles % Dav(No,4) = EigVec(2,1)*w + EigVec(2,2)*x
           Particles % Dav(No,3) = 0.0_dp 

           IF (EigVal(2)>=Particles % criticaldav) THEN
              EigVal(2) = Particles % RiftDmax
              Particles % Dav(No,1:3) = Particles % Riftdmax
              Particles % Dav(No,4) = 0.0_dp
              Particles % damstatus(No) = 1
           ELSE
              Particles % Damstatus(No) = 0
           END IF

           Particles % damage(No,3,1) = EigVal(2)

        END IF
     ELSE
        Particles % damage(No,3,1) = Particles % Dav(No,1)
     END IF
     !  END IF
  END DO


  IF (lc>0.0_dp ) THEN
     ddmax = MAXVAL(ABS(Particles % dD(:,:,:)))     
     PRINT *,'Max particle  dD over current timestep after nonlocal: ', ddmax
     PRINT *,' '
  END IF

  Particles % xpic = 0.0_dp

  endtime = RealTime()
  PRINT *,''
  PRINT *,'DAMAGE SOLUTION TIME',endtime-starttime

  CALL Info(SolverName,'----- Damage Update Done -----',Level=3)

END SUBROUTINE UpdateDamageModifiedZeroStress



SUBROUTINE SaveMaxFrontX_1Dtest( Model, Solver, dt, TransientSimulation)

  USE MPMUtils

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt,maxx,fx,H,vel,gradvel
  LOGICAL :: TransientSimulation
  INTEGER :: No
  LOGICAL :: Visited=.FALSE.,saveh=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: OutputDirectory,FileNamePrefix,FileName,SolverName

  SAVE :: Visited, OutputDirectory,FileNamePrefix,FileName,saveh

  Particles => GlobalParticles
  WRITE(SolverName, '(A)') 'SaveMaxFrontX'

  CALL Info(SolverName,'Saving Max Front X',Level=3) 


  IF( .NOT. Visited ) THEN

     saveh = GetLogical( Solver % Values,'Save thickness')

     OutputDirectory = GetString( Solver % Values,'Filename Directory')
     FileNamePrefix = GetString( Solver % Values,'Filename Prefix')

     FileName = TRIM(OutputDirectory)// '/' //TRIM(FilenamePrefix)// '.result'

     Visited = .TRUE.
  END IF

  maxx = -HUGE(1.0_dp)

  DO No = 1,Particles % NumberOfParticles

     IF (Particles % shapefunctions == 'gimpm') THEN
        fx = Particles % Coordinate(No,1) + 0.5_dp * (Particles % Length(No,1))
     ELSE
        fx = Particles % Coordinate(No,1) + 0.5_dp * (Particles % pvolume(No)/Particles % Length(No,2))
     END IF


     IF (fx>maxx) THEN
        maxx = fx

        IF (saveh) THEN
           H = Particles % H(No)
           vel = Particles % Velocity(No,1)
           gradvel = Particles % GradVel(No,1)
        END IF
     END IF
  END DO

  IF (saveh) THEN
     OPEN (10, FILE=FileName,POSITION='APPEND')
     WRITE (10,'(5ES22.12E3)') Particles % time, maxx,H,vel,gradvel
     CLOSE(10)
  ELSE
     OPEN (10, FILE=FileName,POSITION='APPEND')
     WRITE (10,'(2ES22.12E3)') Particles % time, maxx
     CLOSE(10)
  END IF

  PRINT *,'Time',Particles % time
  PRINT *,'maxx',maxx 

END SUBROUTINE SaveMaxFrontX_1Dtest


SUBROUTINE SaveParticleLoc_1Dtest( Model, Solver, dt, TransientSimulation)

  USE MPMUtils
  USE DefUtils
  USE MeshUtils
  USE SolverUtils
  USE Lists
  USE GeneralUtils  

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt,count,fx,h,vel,gradvel,nearest,miny,closestx,targetinitloc
  REAL(KIND=dp) :: scale, minx,maxx,KE,M,rr,rl
  REAL(KIND=dp) :: pxmin,pxmax,totvol,vol1d
  LOGICAL :: TransientSimulation
  INTEGER :: No,whichno
  LOGICAL :: Visited=.FALSE.,savehvelgradveltime
  CHARACTER(LEN=MAX_NAME_LEN) :: OutputDirectory,FileNamePrefix,FileName,SolverName

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: currenttime
#else
  REAL(KIND=dp) :: currenttime
#endif    

  SAVE :: Visited, OutputDirectory,FileNamePrefix,FileName,whichno,savehvelgradveltime,&
       minx,maxx

  Particles => GlobalParticles
  WRITE(SolverName, '(A)') 'SaveParticleLoc'

  CALL Info(SolverName,'Saving Particle Location',Level=3) 


  IF( .NOT. Visited ) THEN    

     OutputDirectory = GetString( Solver % Values,'Filename Directory')
     FileNamePrefix = GetString( Solver % Values,'Filename Prefix')

     FileName = TRIM(OutputDirectory)// '/' //TRIM(FilenamePrefix)// '.result'

     TargetInitLoc = GetCReal(Solver % Values,'Target Initial Location')

     miny = MINVAL(Particles % Coordinate(1:Particles % NumberOfParticles,2))

     miny = miny+10.0_dp

     closestx = HUGE(1.0_dp)
     whichno = 0

     minx = 0.0_dp
     maxx = 250000.0_dp

     DO No = 1,Particles % NumberOfParticles

        IF (Particles % Coordinate(No,2)>miny) CYCLE

        IF (Particles % Coordinate(No,1)<TargetInitLoc) CYCLE

        IF (Particles % Coordinate(No,1)<closestx) THEN
           closestx = Particles % Coordinate(No,1)
           whichno = No
        END IF
     END DO

     PRINT *,''
     PRINT *,'TargetInitLoc',TargetInitLoc
     PRINT *,'miny',miny
     PRINT *,'closestx',closestx
     PRINT *,'WHICHNO',whichno
     PRINT *,''

     savehvelgradveltime = GetLogical( Solver % Values,'Save hvelgradveltime')

     Visited = .TRUE.
  END IF

  count = 0.0_dp
  fx = 0.0_dp
  h = 0.0_dp
  vel = 0.0_dp
  gradvel = 0.0_dp
  vol1d = 0.0_dp

  rr = -HUGE(1.0_dp)
  rl = HUGE(1.0_dp)

  KE = 0.0_dp
  M = 0.0_dp

  DO No = 1,Particles % NumberOfParticles
     IF (Particles % origno(No) == whichno) THEN
        fx = fx + Particles % Coordinate(No,1)
        h = h + Particles % H(No)
        vel = vel + Particles % Velocity(No,1)
        gradvel = gradvel + Particles % gradvel(No,1)

        IF (Particles % shapefunctions == 'gimpm') THEN
           vol1d = vol1d + Particles % gvolume(No)
        ELSE
           vol1d = vol1d + Particles % pvolume(No)
        END IF

        count = count + 1

        rr = MAX(rr,Particles % Coordinate(No,1)+Particles % Length(No,1)/2.0_dp)
        rl = MIN(rl,Particles % Coordinate(No,1)-Particles % Length(No,1)/2.0_dp)
     END IF

     IF (Particles % shapefunctions == 'gimpm') THEN
        pxmin = Particles % Coordinate(No,1) - Particles % Length(No,1)/2.0_dp        
        pxmax = Particles % Coordinate(No,1) + Particles % Length(No,1)/2.0_dp
        totvol = pxmax-pxmin

        IF (pxmax<minx) CYCLE
        IF (pxmin>maxx) CYCLE

        IF (pxmax > maxx) THEN
           !overlaps R bound
           scale = (maxx-pxmin)/Particles % Length(No,1)
        ELSEIF (pxmin < minx) THEN
           !overlaps L bound
           scale = (pxmax-minx)/Particles % Length(No,1)
        ELSE
           scale = 1.0_dp
        END IF

     ELSE
        IF (Particles % Coordinate(No,1) > maxx) CYCLE
        IF (Particles % Coordinate(No,1) < minx) CYCLE
        scale = 1.0_dp
     END IF

     M = M + Particles % Mass(No)*scale
     KE = KE + 0.5_dp*Particles % Mass(No)*scale*Particles % Velocity(No,1)*Particles % Velocity(No,1)     
  END DO

  IF (count>0) THEN
     fx = fx/count
     h = h/count
     vel = vel/count
     gradvel = gradvel/count
  ELSE
     fx = -9999
     h =-9999
     vel = -9999
     gradvel = -9999
     vol1d = -9999
  END IF


  IF (savehvelgradveltime) THEN

     currenttime = RealTime()
     Particles % simtime = currenttime - Particles % firstsimtime

     OPEN (10, FILE=FileName,POSITION='APPEND')
     WRITE (10,'(2ES27.20,F4.1,9ES27.20)') &
          rr,rl,count,Particles % time, fx,H,vel,gradvel,Particles % simtime,M,KE,vol1d
     CLOSE(10)
  ELSE
     OPEN (10, FILE=FileName,POSITION='APPEND')
     WRITE (10,'(2ES22.12E3)') Particles % time, fx
     CLOSE(10)
  END IF

END SUBROUTINE SaveParticleLoc_1Dtest


!> Several different measures of stress error are saved here. The one used in the
!! paper is W_T_E
SUBROUTINE SteadyStressError_1Dtest( Model, Solver, dt, TransientSimulation)

  USE MPMUtils
  USE DefUtils
  USE MeshUtils
  USE SolverUtils
  USE Lists
  USE GeneralUtils  

  IMPLICIT NONE
  TYPE(Particle_t), POINTER :: Particles  
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  LOGICAL :: TransientSimulation
  INTEGER :: No
  LOGICAL :: Visited=.FALSE.
  REAL(KIND=dp) :: dt
  CHARACTER(LEN=MAX_NAME_LEN) :: OutputDirectory,FileNamePrefix,FileName,SolverName
  REAL(KIND=dp) :: miny,minx,maxx,cm,secondsperyear,H0,v0,Q0,B0,A,C,EeExp,Axm,m1,m2,acm
  REAL(KIND=dp) :: count,RMSnum_Tau,RMSnum_RTau,tot_P_Tau,tot_P_Rtau
  REAL(KIND=dp) ::  Exx,Ee,P_Tau,P_RTau,vol1d,Ha,A_Tau,A_RTau,totvolT,totvolRT
  REAL(KIND=dp) :: RMS_Tau,RMS_RTau
  REAL(KIND=dp) :: mean_P_Tau,mean_P_RTau,T_SS_tot,RT_SS_tot,T_SS_RES,RT_SS_RES
  REAL(KIND=dp) :: RR_T,RR_RT
  REAL(KIND=dp) :: NW_T_E,NW_RT_E,W_T_E,W_RT_E,W_T_E_b,W_RT_E_b
  REAL(KIND=dp) :: RSE_T,RSE_RT,nonanT,nonanRT,minsrinvsquared,nancount


  SAVE :: Visited, OutputDirectory,FileNamePrefix,FileName,miny,minx,maxx,&
       cm,secondsperyear,H0,v0,Q0,B0,A,C,EeExp,Acm,m1,m2,minsrinvsquared

  Particles => GlobalParticles
  WRITE(SolverName, '(A)') 'SaveParticleLoc'

  CALL Info(SolverName,'Saving Particle Location',Level=3) 


  IF( .NOT. Visited ) THEN    

     OutputDirectory = GetString( Solver % Values,'Filename Directory')
     FileNamePrefix = GetString( Solver % Values,'Filename Prefix')

     FileName = TRIM(OutputDirectory)// '/' //TRIM(FilenamePrefix)// '.result'

     miny = MINVAL(Particles % Coordinate(1:Particles%NumberOfParticles,2))

     miny = miny+5.0_dp

     minx = 0.0_dp
     maxx = 250000.0_dp

     MinSRInvSquared = Particles % criticalshearrate * Particles % criticalshearrate


     !ANALYTICAL STUFF
     cm = 1.0_dp/3.0_dp
     secondsperyear = 31556926.0_dp
     H0 = 600.0_dp
     v0 = 300.0_dp
     Q0 = H0*v0
     B0 = 1.9E8_dp
     A = ((B0*1.0E-6_dp)**(-3.0_dp))*secondsperyear !Mpa^(-3) a^(-1)
     C = (((910.0_dp*1.0e-6_dp*9.81_dp)/&
          (4.0_dp*(A**(-cm))))*(1.0_dp-910.0_dp/1028.0_dp))**3.0_dp
     !C is the weertman constant !C =2.45E-18; !m?3 s?1

     EeExp = (cm-1.0_dp)/2.0_dp
     Acm = A**(-cm)
     m1 = 4.0_dp*C/Q0     
     m2 = 1.0_dp/(H0*H0*H0*H0)


     Visited = .TRUE.
  END IF


  count = 0.0_dp
  RMSnum_Tau = 0.0_dp
  RMSnum_RTau = 0.0_dp
  tot_P_Tau = 0.0_dp
  tot_P_RTau = 0.0_dp
  nancount = 0.0_dp


  !XPIC AND dD are simply used for storage here, and have nothing to do with the
  !xpic routine nor the damage solver in this context
  Particles % xpic = 0.0_dp
  Particles % dD = 0.0_dp

  DO No = 1,Particles % NumberOfParticles

     IF (Particles % Coordinate(No,1) < minx) CYCLE
     IF (Particles % Coordinate(No,1) > maxx) CYCLE
     IF (Particles % Coordinate(No,2) > miny) CYCLE

     IF (Particles % ShapeFunctions == 'gimpm') THEN
        IF (Particles % Coordinate(No,1) + 0.5_dp*Particles % Length(No,1) > maxx) CYCLE
        IF (Particles % Coordinate(No,1) - 0.5_dp*Particles % Length(No,1) < minx) CYCLE
     END IF

     !----------particle stresses-----
     Exx = Particles % GradVel(No,1)
     Ee = Exx*Exx

     IF (Ee < MinSRInvSquared) THEN
        Ee = MinSRInvSquared
     END IF

     !dev stress
     P_Tau = Acm*(Ee**EeExp)*Exx
     !resistive stress
     P_RTau = 2.0_dp*P_Tau*Particles % H(No)

     vol1d = Particles % pvolume(No)/Particles % Length(No,2)


     Particles % xpic(No,2) = P_Tau
     Particles % xpic(No,3) = P_RTau
     Particles % dD(No,1,2) = vol1d

     !------analytical stresses-------
     Ha = (m1*Particles % Coordinate(No,1) + m2)**(-0.25_dp)
     Exx = C*Ha*Ha*Ha
     Ee = Exx*Exx
     !dev stress
     A_Tau = Acm*(Ee**EeExp)*Exx
     !resistive stress
     A_RTau = 2.0_dp*A_Tau*Ha


     IF ((A_Tau .NE. A_Tau) .OR. (P_Tau .NE. P_Tau)) THEN
        nancount = nancount + 1
        CYCLE
     END IF



     !For RMS:
     RMSnum_Tau = RMSnum_Tau + (A_Tau-P_Tau)*(A_Tau-P_Tau)
     RMSnum_RTau = RMSnum_RTau + (A_RTau-P_RTau)*(A_RTau-P_RTau)

     !For R Squared:
     tot_P_Tau = tot_P_Tau + P_Tau
     tot_P_RTau = tot_P_RTau + P_RTau

     !Miscellaneous
     Particles % XPIC(No,4) = A_Tau-P_Tau
     Particles % XPIC(No,5) = A_RTau-P_RTau
     Particles % XPIC(No,6) = A_Tau
     Particles % dD(No,1,1) = A_RTau

     count = count+1.0_dp
     Particles % xpic(No,1) = 1.0_dp        

  END DO

  !-----RMS-----
  RMS_Tau = sqrt(RMSnum_Tau/count)
  RMS_RTau = sqrt(RMSnum_RTau/count)

  !-----R Squared-----
  mean_P_Tau = tot_P_Tau/count
  mean_P_RTau = tot_P_RTau/count
  WHERE (Particles % xpic(:,1) == 1.0_dp)
     Particles % xpic(:,2) = Particles % xpic(:,2)-mean_P_Tau
     Particles % xpic(:,3) = Particles % xpic(:,3)-mean_P_RTau
  END WHERE
  
  T_SS_tot = SUM(Particles % xpic(:,2)*Particles % xpic(:,2))
  RT_SS_tot = SUM(Particles % xpic(:,3)*Particles % xpic(:,3))
  T_SS_RES = SUM(Particles % xpic(:,4)*Particles % xpic(:,4))
  RT_SS_RES = SUM(Particles % xpic(:,5)*Particles % xpic(:,5))

  RR_T = 1.0_dp - T_SS_RES/T_SS_tot
  RR_RT = 1.0_dp - RT_SS_RES/RT_SS_tot


  !----Stress Error, Weighted, like Bing et al 2019 BC paper ----

  W_T_E_b = 0.0_dp
  W_RT_E_b = 0.0_dp
  totvolT = 0.0_dp
  totvolRT = 0.0_dp
  nonanT = 0.0_dp
  nonanRT = 0.0_dp

  DO No = 1,Particles % NumberOfParticles

     IF (Particles % xpic(No,1) .NE. 1.0_dp) CYCLE

     !relative stress error at material point:     
     IF (Particles % xpic(No,6) .NE. 0.0_dp) THEN
        RSE_T = (ABS(Particles % xpic(No,4))/ABS(Particles % xpic(No,6)))*Particles % dD(No,1,2)
        IF (RSE_T == RSE_T) THEN
           W_T_E_b = W_T_E_b + RSE_T
           totvolT = totvolT + Particles % dD(No,1,2)
        ELSE
           nonanT = nonanT + 1
        END IF
     END IF

     IF (Particles % dD(No,1,1) .NE. 0.0_dp) THEN
        RSE_RT = (ABS(Particles % xpic(No,5))/ABS(Particles % dD(No,1,1)))*Particles % dD(No,1,2)
        IF (RSE_RT == RSE_RT) THEN
           W_RT_E_b = W_RT_E_b + RSE_RT
           totvolRT = totvolRT + Particles % dD(No,1,2)
        ELSE
           nonanRT = nonanRT + 1
        END IF
     END IF
  END DO

  W_T_E_b = W_T_E_b/totvolT
  W_RT_E_b = W_RT_E_b/totvolRT



  !---Stress Error, Not Weighted ----
  !num = sum(sqrt((data-actual)**2)*vol)
  !denom = sum(sqrt(actual**2)*vol)

  !ABS( actual tau - particle tau)
  Particles % xpic(:,4) = ABS(Particles % xpic(:,4))
  !ABS( actual Rtau - particle Rtau)
  Particles % xpic(:,5) = ABS(Particles % xpic(:,5))
  !ABS( actual tau)
  Particles % xpic(:,6) = ABS(Particles % xpic(:,6))
  !ABS( actual RTau)
  Particles % dD(:,1,1) = ABS(Particles % dD(:,1,1))

  NW_T_E = SUM(Particles % xpic(:,4))/SUM(Particles % xpic(:,6))
  NW_RT_E = SUM(Particles % xpic(:,6))/SUM(Particles % dD(:,1,1))



  !----Stress Error,  Weighted ----
  !ABS( actual tau - particle tau) * vol
  Particles % xpic(:,4) = Particles % xpic(:,4) * Particles % dD(:,1,2)
  !ABS( actual Rtau - particle Rtau) * vol
  Particles % xpic(:,5) = Particles % xpic(:,5) * Particles % dD(:,1,2)
  !ABS( actual tau) * vol
  Particles % xpic(:,6) = Particles % xpic(:,6) * Particles % dD(:,1,2)
  !ABS( actual RTau) * vol
  Particles % dD(:,1,1) = Particles % dD(:,1,1) * Particles % dD(:,1,2)

  W_T_E = SUM(Particles % xpic(:,4))/SUM(Particles % xpic(:,6))
  W_RT_E = SUM(Particles % xpic(:,6))/SUM(Particles % dD(:,1,1))


  Particles % xpic = 0.0_dp
  Particles % dD = 0.0_dp


  OPEN (10, FILE=FileName,POSITION='APPEND') 
  WRITE (10,'(11ES19.12,3F4.1)') &
       Particles % time, RMS_Tau, RMS_RTau, RR_T, RR_RT, NW_T_E, &
       NW_RT_E, W_T_E, W_RT_E, W_T_E_b, W_RT_E_b,nonanT,nonanRT,nancount
  CLOSE(10)

END SUBROUTINE SteadyStressError_1Dtest


SUBROUTINE AdvectRiftsTest( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE MPMUtils
  USE SolverUtils
  USE Lists
  USE GeneralUtils

  IMPLICIT NONE
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  TYPE(Element_t), POINTER :: BulkElement  
  REAL(KIND=dp) :: dt, PVal, Basis(4),xt(4),SqrtElementMetric,Coord(3)
  REAL(KIND=dp) :: xnodes(4),ynodes(4),scale,xmin
  LOGICAL :: TransientSimulation
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: dim,No,ii
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp), POINTER :: Val(:)
  TYPE(Valuelist_t), POINTER :: Params
  LOGICAL :: GotIt,InterpToParticles,Stat
  TYPE(Particle_t), POINTER :: Particles
  INTEGER, POINTER :: NodeIndexes(:)    

  Particles => GlobalParticles  
  Params => GetSolverParams()  
  Mesh => Solver % Mesh
  dim = CoordinateSystemDimension()

  InterpToParticles = GetLogical( Params,'Interp To Particles',GotIt)
  IF (.NOT. GotIt) THEN
     InterpToParticles = .FALSE.
  END IF


  Var => VariableGet(Model % Mesh % Variables, 'TempVar' )
  Perm => Var % Perm
  Val => Var % Values
  Val = 0.0_dp


  scale = 1.0e3_dp

  xt =(/40.0_dp,45.0_dp,50.0_dp,55.0_dp/)*scale*10.0_dp


  DO ii = 1,Mesh % NumberOfBulkElements
     BulkElement => Model % Mesh % Elements( ii )
     NodeIndexes => BulkElement % NodeIndexes
     xnodes = Mesh % Nodes % x(NodeIndexes(1:4))
     ynodes = Mesh % Nodes % y(NodeIndexes(1:4))

     xmin = MINVAL(xnodes)


     IF (ANY(xt==xmin) ) THEN !.OR. ANY( (xt+Particles%gridres)==xmin)) THEN

        IF (ALL(ynodes >= 37.5_dp*scale) .OR. &
             ALL(ynodes <= 2.5_dp*scale) .OR. &
             (ALL(ynodes >= 7.5_dp*scale) .AND. ALL(ynodes<=12.5_dp*scale)) .OR. &
             (ALL(ynodes >= 17.5_dp*scale) .AND. ALL(ynodes<=22.5_dp*scale)) .OR. &
             (ALL(ynodes >= 27.5_dp*scale) .AND. ALL(ynodes<=32.5_dp*scale))) THEN


           Val(Perm(NodeIndexes(1:4))) = 1.0_dp

        END IF
     END IF
  END DO


  IF (InterpToParticles) THEN
     PRINT *,'Interpolating to Particles'

     Particles % usetracer = .TRUE.
     DO No = 1,Particles % NumberOfParticles

        BulkElement => Model % Mesh % Elements(Particles % ElementIndex(No) )        
        Coord = GetParticleCoord( Particles, No)
        stat = ParticleElementInfo( BulkElement, Coord, &
             SqrtElementMetric, Basis )

        CALL GetScalarFieldInMesh(Var, BulkElement, Basis, PVal )

        Particles % Tracer(No) = PVal
     END DO
  END IF
END SUBROUTINE AdvectRiftsTest

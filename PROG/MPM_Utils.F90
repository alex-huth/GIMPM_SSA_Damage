!-----------------------------------------------------------------------------
!> Module including generic utilities for the Material Point Method (MPM)
!! for the shallow shelf approximation (SSA) of ice flow.
!! Includes standard MPM (sMPM) and Generalized Interpolation MPM (GIMPM).
!!
!! Note this code is not yet parallelized.
!! However, some subroutines are similar to (or directly edited from)
!! the particle routines/utils already available in Elmer,
!! which will provide some direction for future parallelization...
!!
!! -Alex Huth, 2020 (ahuth@princeton.edu)
!-----------------------------------------------------------------------------

MODULE MPMUtils

  USE DefUtils
  USE Lists
  USE MeshUtils
  USE GeneralUtils
  USE Interpolation

  IMPLICIT NONE

  TYPE Particle_t
     INTEGER :: Dim, NumberOfParticles=0, MaxNumberOfParticles=0,rkmsteps

     REAL(KIND=dp) :: time, dtime, nextdtime

     !damage stuff
     REAL(KIND=dp) :: psr(3,3),pressure1,RHS,dvdxmdudy !, P_i, P_w
     REAL(KIND=dp) :: CurrentEVal(3),CurrentEVect(2,2)

     REAL(KIND=dp), POINTER :: Coordinate(:,:) => NULL()
     REAL(KIND=dp), POINTER :: NextCoordinate(:,:) => NULL()

     INTEGER, POINTER :: Status(:) => NULL()
     INTEGER, POINTER :: ElementIndex(:) => NULL()
     INTEGER, POINTER :: InterpElem(:) => NULL()

     ! Mark the internal elements without any interface nodes
     LOGICAL, POINTER :: InternalElements(:) => NULL()

     ! YOUR ADDITIONS
     CHARACTER(LEN=MAX_NAME_LEN) :: DamageModel,ShapeFunctions

     REAL(KIND=dp) :: CriticalDamage,CriticalDav,Dmax,&
          ah,bf,bh,k1,k2,gamma,rf,sthres,kf,bfmgamma,Gravity,gridres,&
          criticalshearrate,maximumtimestep,viscosityexponent,maxlength,&
          maxdplength,maxGLlength,davsplitthres,dinittolerance,warningtol,&
          rhoi,rhow,sealevel,elementfraction,consttemp,giveneta,&
          rdxmin,rdxmax,rdymin,rdymax,initdmax,zsmaxdd,&
          stresshigh,stresslow,ndxmin,ndxmax,ndymin,ndymax,melangefraction,&
          currentgamma,tempgamma,fullrupt,fillpercent,DMaxI,DMaxII,DMaxIII,&
          beginisothres,isodamcritdam,&
          isodamcritdav,mindam,sthresmod,dDscale,EFscale,OnePPCXEnd,&
          riftdmax,bmd,bmdsave,tonguelfp,damstrainmax,rmelfrac,DavDMaxI,&
          DavDmaxII,DavDmaxIII,contactfrac

     INTEGER :: dsyevcountpd,generalcountpd,dsyevcountet,&
          generalcountet,returnedpd,equaleigcount

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: simtime,firstsimtime,&
          dsyevtimepd,dsyevtimeet,generaltimepd,generaltimeet,&
          ttime1,ttime2,fpdtime,dddttime,times,timee
#else
     REAL(KIND=dp) :: simtime,firstsimtime,&
          dsyevtimepd,dsyevtimeet,generaltimepd,generaltimeet,&
          ttime1,ttime2,fpdtime,dddttime,times,timee
#endif


     INTEGER :: NumberOfTemperatureLayers,buffer,&
          dinitnumiters,numberofparticlelayers,&
          currentno,currentlayer,bcracklayers,dbcracklayers

     LOGICAL :: ConstLinTemp,nodamfront,useconsttemp,&
          movegl,constmb,constef,constfric,SEP,simpleadvecttest,&
          usegiveneta,setyvelzero,usebcforprevvel,usestaticparticles,UseHBC,&
          usesavedbasis,uplag,velocitydependentfriction,coulombfriction,&
          useisodam,peig,femforallnotfull,lmod,modifiedmurakami,&
          restrictdam,rkmcritdam,nodzz,forcedzz,usedavtensor,&
          fixdavdirections,nodamregion,frontmelange,noxpiconbounds,analytictest,&
          firsttimestepzero,isorift,dmaxII_dom,&
          flipvelfric,weighth,trackstrain,&
          alwayssplitfour,savedDscale,usedamscale,FEMifGrounded,nodaminflow,&
          binitlowerboundfromfirsttimestep,alwaysfemfront,prupt,noevolveruptlayers,&
          outputdbassis,mix49,Usedamage,mismipmelt2,IncNumPAtStart,UseOnePPC,simpleinitdam,&
          MISMIP,useriftdmax,nospin,trackdamstrain,initriftdam,&
          hoop,unstructuredmesh,&
          initcrack,noriftbounduntilfullsep,efboundsfromfirsttimestep,&
          useriftmelange,gradualdyz,&
          larcmelfracmod,LarCFixDavTest,LarCDamTraj,usetracer,usetruecauchydamage

     REAL(KIND=dp) :: icrackwidth,icrackx1,icrackx2,icracky1,icracky2

     LOGICAL, POINTER :: UseInterpElem(:) => NULL()
     LOGICAL, POINTER :: Static(:) => NULL()

     REAL(KIND=dp), POINTER :: Damage(:,:,:) => NULL()
     REAL(KIND=dp), POINTER :: dD(:,:,:) => NULL()
     REAL(KIND=dp), POINTER :: damstrain(:) => NULL()
     REAL(KIND=dp), POINTER :: F(:,:) => NULL()
     REAL(KIND=dp), POINTER :: Dav(:,:) => NULL()
     REAL(KIND=dp), POINTER :: GradVel(:,:) => NULL()
     REAL(KIND=dp), POINTER :: GradZs(:,:) => NULL()
     REAL(KIND=dp), POINTER :: GradH(:,:) => NULL()
     REAL(KIND=dp), POINTER :: Velocity(:,:) => NULL()
     REAL(KIND=dp), POINTER :: GridVelocity(:,:) => NULL()
     REAL(KIND=dp), POINTER :: Length(:,:) => NULL()
     REAL(KIND=dp), POINTER :: OrigLength(:,:) => NULL()
     REAL(KIND=dp), POINTER :: Bz(:,:) => NULL()
     REAL(KIND=dp), POINTER :: xpic(:,:) => NULL()
     REAL(KIND=dp), POINTER :: EF(:) => NULL()
     REAL(KIND=dp), POINTER :: GMask(:) => NULL()
     REAL(KIND=dp), POINTER :: Bedrock(:) => NULL()
     REAL(KIND=dp), POINTER :: MB(:) => NULL()
     REAL(KIND=dp), POINTER :: Binit(:) => NULL()
     REAL(KIND=dp), POINTER :: FP(:) => NULL()
     REAL(KIND=dp), POINTER :: H(:) => NULL()
     REAL(KIND=dp), POINTER :: GVolume(:) => NULL()
     REAL(KIND=dp), POINTER :: PVolume(:) => NULL()
     REAL(KIND=dp), POINTER :: dbassis(:) => NULL()
     REAL(KIND=dp), POINTER :: Mass(:) => NULL()
     REAL(KIND=dp), POINTER :: Strain(:,:) => NULL()
     REAL(KIND=dp), POINTER :: Tracer(:) => NULL()

     INTEGER, POINTER :: OrigNo(:) => NULL()
     INTEGER, POINTER :: DamStatus(:) => NULL()
  END TYPE Particle_t

  TYPE MPforSSA_t
     REAL(KIND=dp), POINTER :: DSRxx(:) => NULL()
     REAL(KIND=dp), POINTER :: DSRyy(:) => NULL()
     REAL(KIND=dp), POINTER :: DSRxy(:) => NULL()
     REAL(KIND=dp), POINTER :: eta(:)   => NULL()
     REAL(KIND=dp), POINTER :: muder(:) => NULL()
     REAL(KIND=dp), POINTER :: slip(:) => NULL()
     REAL(KIND=dp), POINTER :: Ezz(:) => NULL()
     REAL(KIND=dp), POINTER :: Exy(:) => NULL()
     REAL(KIND=dp), POINTER :: Hf(:) => NULL()
     REAL(KIND=dp), POINTER :: fN(:) => NULL()
     REAL(KIND=dp), POINTER :: fB(:) => NULL()
     REAL(KIND=dp), POINTER :: Ee(:) => NULL()
     REAL(KIND=dp), POINTER :: Dxx(:) => NULL()
     REAL(KIND=dp), POINTER :: Dyy(:) => NULL()
     REAL(KIND=dp), POINTER :: Dzz(:) => NULL()
     REAL(KIND=dp), POINTER :: Dxy(:) => NULL()
     REAL(KIND=dp), POINTER :: exxd1m1(:) => NULL()
     REAL(KIND=dp), POINTER :: eyyd2m1(:) => NULL()
     REAL(KIND=dp), POINTER :: ezzd3m1(:) => NULL()
     REAL(KIND=dp), POINTER :: exyd4(:)   => NULL()
     REAL(KIND=dp), POINTER :: slip2(:)  => NULL()
     REAL(KIND=dp), POINTER :: ub(:)     => NULL()
     REAL(KIND=dp), POINTER :: driveforce(:,:) => NULL()
     REAL(KIND=dp), POINTER :: GradVel(:,:) => NULL()
     REAL(KIND=dp), POINTER :: GridVelocity(:,:) => NULL()
     REAL(KIND=dp), POINTER :: Velo(:,:)   => NULL()
     ! REAL(KIND=dp), POINTER :: falpha => NULL !not used in this version

  END type MPforSSA_t


  TYPE :: ElementParticles_t
     INTEGER :: NumberOfParticles
     INTEGER, allocatable :: p(:)
     INTEGER, allocatable :: SurroundingElems(:)
     INTEGER, allocatable :: NLSelems(:)
     INTEGER :: numNLSelems

     REAL(KIND=dp), allocatable :: Basis(:,:)
     REAL(KIND=dp), allocatable :: dBasisdx(:,:,:)
  END TYPE ElementParticles_t

  TYPE :: ElementTrack_t
     REAL(KIND=dp) :: Volume
     INTEGER :: InitialSurface
     INTEGER :: Status
     INTEGER :: GStatus
     INTEGER :: ddlayfrombottom1
     INTEGER :: ddlayfromtop1
     INTEGER :: ddlayfrombottom2
     INTEGER :: ddlayfromtop2
  END TYPE ElementTrack_t

  TYPE :: SSAHInterpolationNodes_t
     INTEGER :: NumberOfNodes
     INTEGER, allocatable :: Node(:)
     INTEGER, allocatable :: ClosestParticle(:)
     REAL(KIND=dp), allocatable :: Distance(:)
  END TYPE SSAHInterpolationNodes_t


  TYPE :: InterpLayers_t
     INTEGER :: NumOfTLayers
     INTEGER :: NumOfDLayers
     INTEGER, allocatable :: Map(:,:)
     REAL(KIND=dp), allocatable :: InterpFun(:,:)
  END type InterpLayers_t


  INTEGER, PARAMETER :: &

       PARTICLE_ALLOCATED = 1, &
       PARTICLE_ACTIVE = 2, &
       PARTICLE_FIXEDCOORD = 3, &
       PARTICLE_LEAVING = 4, &
       PARTICLE_LOST = 5, &

       EMPTY = 1, &
       IGNORE = 2, &
       NOTFULL = 3, &
       FEM = 4, &
       FULL = 5

  TYPE(Particle_t), TARGET, SAVE :: GlobalParticles
 ! TYPE(MPforSSA_t), TARGET, SAVE :: GlobalMPforSSA
  TYPE(MPforSSA_t), TARGET, SAVE :: MP
  TYPE(ElementParticles_t), public, SAVE, Allocatable,target :: ElemParticles(:)
  TYPE(ElementTrack_t), public, SAVE, Allocatable,target :: ElemTrack(:)
  TYPE(InterpLayers_t), public, SAVE :: InterpLayers

  !**************************************************************************
  !**************************************************************************

CONTAINS

  !> Gets the elements where the particle is located
  FUNCTION GetParticleElement(Particles, No) RESULT ( Index )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Index

    Index = Particles % ElementIndex(No)
  END FUNCTION GetParticleElement

  !> Gets the Cartesian coordinates for a particle
  FUNCTION GetParticleCoord(Particles, No) RESULT ( Coord )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Coord(3)
    INTEGER :: dim

    Coord(3) = 0.0_dp
    dim = Particles % dim
    Coord(1:dim) = Particles % Coordinate(no,:)
    Coord(1) = Particles % Coordinate(No,1)
    Coord(2) = Particles % Coordinate(No,2)
  END FUNCTION GetParticleCoord

  !> Gets the Status for a particle
  !! 1. PARTICLE_ALLOCATED; 2. PARTICLE_ACTIVE; 3. PARTICLE_FIXEDCOORD
  !! 4. PARTICLE_LEAVING; 5. PARTICLE_LOST
  FUNCTION GetParticleStatus(Particles,No) RESULT ( Status )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Status

    Status = Particles % Status(No)
  END FUNCTION GetParticleStatus

  !> Determines whether a point is within a square grid cell
  FUNCTION coordsinelement( Mesh, ElementInd, Coord) RESULT (inelement)

    REAL(KIND=dp) :: xmin,xmax, ymin, ymax,Coord(3)
    INTEGER :: nn,ElementInd
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: inelement

    Element => Mesh % Elements(ElementInd)
    nn = Element % TYPE % NumberOfNodes
    NodeIndexes => Element % NodeIndexes
    xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
    xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
    ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
    ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
    IF ((Coord(1) <= xmax) .AND. (Coord(1)>=xmin) &
         .AND. (Coord(2)>=ymin) .AND. (Coord(2)<=ymax)) THEN
       inelement = .TRUE.
    ELSE
       inelement = .FALSE.
    END IF

  END FUNCTION coordsinelement

  !> Determines what volume of a GIMPM domain is within a square grid cell
  FUNCTION GetParticleVolumeInElement(Particles,No,Element,Model) RESULT (Volume)
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Element_t), POINTER :: Element
    TYPE(Model_t) :: Model
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: No,nn
    REAL(KIND=dp) :: Volume,Coord(3),xmin,xmax,ymin,ymax,LX,LY,N,S,E,W

    nn = Element % TYPE % NumberOfNodes
    NodeIndexes => Element % NodeIndexes

    xmin = MINVAL(Model % Mesh % Nodes % x(NodeIndexes(1:nn)))
    xmax = MAXVAL(Model % Mesh % Nodes % x(NodeIndexes(1:nn)))
    ymin = MINVAL(Model % Mesh % Nodes % y(NodeIndexes(1:nn)))
    ymax = MAXVAL(Model % Mesh % Nodes % y(NodeIndexes(1:nn)))

    Coord = GetParticleCoord( Particles, No)

    LX = Particles % Length(No,1)
    LY = Particles % Length(No,2)

    W = Coord(1)-LX*0.5_dp
    E = Coord(1)+LX*0.5_dp
    S = Coord(2)-LY*0.5_dp
    N = Coord(2)+LY*0.5_dp

    Volume = 0.0_dp
    IF (xmax < Coord(1)) THEN
       Volume =  xmax-W
    ELSEIF (xmin > Coord(1)) THEN
       Volume =  E-xmin
    ELSE
       Volume = MIN(E,xmax)-MAX(W,xmin)
    END IF

    IF ( ymin > Coord(2) ) THEN
       Volume = (N-ymin) * Volume
    ELSEIF (ymax < Coord(2)) THEN
       Volume = (ymax-S) * Volume
    ELSE
       Volume = (MIN(N,ymax)-MAX(S,ymin)) * Volume
    END IF

    Volume = ABS(Volume)
  END FUNCTION GetParticleVolumeInElement

  !> sMPM shape function from Cartesian coordinates (not element-wise)
  FUNCTION sMPMElementInfoFromCoords( Particles, Model, Nodes, h, &
       coord, Basis, dBasisdx) RESULT (stat)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: Basis(4),dBasisdx(4,3),Sv(4,2),dSv(4,2),diff(4,2)
    REAL(KIND=dp) :: one=1.0_dp,two = 2.0_dp,coord(2)
    REAL(KIND=dp) :: h,oogr,e,n,half=0.5_dp,midx,midy
    LOGICAL :: Stat

    diff(:,1) = coord(1) - Nodes % x
    diff(:,2) = coord(2) - Nodes % y

    dSv = one/h
    dSv = -sign( dSv,diff);

    midx = MINVAL(Nodes % x) + half*h
    midy = MINVAL(Nodes % y) + half*h

    oogr = two/h
    e = oogr*(coord(1)-midx)
    n = oogr*(coord(2)-midy)

    WHERE(Nodes % x < midx)
       Sv(:,1) = one-e
    ELSEWHERE
       Sv(:,1) = one+e
    END WHERE

    WHERE(Nodes % y < midy)
       Sv(:,2) = one-n
    ELSEWHERE
       Sv(:,2) = one+n
    END WHERE

    Sv = Sv*half

    WHERE(ABS(Sv)>one) dSv = -dSv

    Basis = Sv(:,1) * Sv(:,2)
    dBasisdx(:,1) = Sv(:,2) * dSv(:,1)
    dBasisdx(:,2) = Sv(:,1) * dSv(:,2)

    Stat = .TRUE.
  END FUNCTION sMPMElementInfoFromCoords


  !> sMPM shape functions (element-wise)
  !! for rectangular elements, unless Particle % unstructuredmesh = .TRUE.
  FUNCTION sMPMElementInfo( Element,Particles, Model, Nodes, No, h, &
       Basis, dBasisdx) RESULT (stat)

    USE TYPES
    USE DefUtils
    USE SolverUtils
    USE MeshUtils

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Element_t), POINTER :: Element
    TYPE(Model_t) :: Model
    TYPE(Nodes_t) :: Nodes
    INTEGER :: No
    REAL(KIND=dp) :: Basis(4),dBasisdx(4,3),ddBasisddx(4,3,3),Sv(4,2),dSv(4,2),diff(4,2)
    REAL(KIND=dp) :: one=1.0_dp,two = 2.0_dp,half=0.5_dp
    REAL(KIND=dp) :: h,oogr,e,n,midx,midy,detj,Coord(3)
    LOGICAL :: Stat,Visited=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: np,t

    SAVE :: t,Visited

    !If you run the sMPM in Updated Lagrangian mode, then
    !both particle and nodal position are updated each timestep
    !The particle locations are always assumed to be at the usual FEM
    !Gauss points:
    IF (Particles % uplag) THEN
       IF (.NOT. Visited) THEN
          t = 1
          Visited = .TRUE.
       END IF

       IP = GaussPoints( Element , np=INT(Particles % elementfraction) )
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       t = t+1

       IF (t>Particles % elementfraction) t = 1
       RETURN
    END IF

    IF (Particles % unstructuredmesh) THEN
       !just call particleelementinfo for the hoop test,
       !because it uses a non-rectangular mesh. The rest of th
       Coord = GetParticleCoord(Particles,No)
       stat = ParticleElementInfo( Element, Coord, &
            detJ, Basis, dBasisdx )
       ddBasisddx = 0.0_dp

       RETURN
    END IF

    diff(:,1) = Particles % Coordinate(No,1) - Nodes % x
    diff(:,2) = Particles % Coordinate(No,2) - Nodes % y

    dSv = one/h
    dSv = -sign( dSv,diff)

    midx = MINVAL(Nodes % x) + half*h
    midy = MINVAL(Nodes % y) + half*h

    oogr = two/h
    e = oogr*(Particles % Coordinate(No,1)-midx)
    n = oogr*(Particles % Coordinate(No,2)-midy)

    WHERE(Nodes % x < midx)
       Sv(:,1) = one-e
    ELSEWHERE
       Sv(:,1) = one+e
    END WHERE

    WHERE(Nodes % y < midy)
       Sv(:,2) = one-n
    ELSEWHERE
       Sv(:,2) = one+n
    END WHERE

    Sv = Sv*half

    WHERE(ABS(Sv)>one) dSv = -dSv

    Basis = Sv(:,1) * Sv(:,2)
    dBasisdx(:,1) = Sv(:,2) * dSv(:,1)
    dBasisdx(:,2) = Sv(:,1) * dSv(:,2)

    Stat = .TRUE.
  END FUNCTION sMPMElementInfo


  !> GIMPM shape functions from Cartesian coordinates (not element-wise)
  !! not used currently
  FUNCTION GIMPMElementInfoFromCoords( Particles, Model, Nodes, No, h, &
       Basis, dBasisdx,iter) RESULT (stat)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    TYPE(Nodes_t) :: Nodes
    INTEGER :: No,j,k,iter
    REAL(KIND=dp) :: Basis(4),dBasisdx(4,3),Sv(4,2),dSv(4,2),diff(4,2)
    REAL(KIND=dp), DIMENSION (2) :: t1,t2,t3,t4,t5,t6,lp
    REAL(KIND=dp) :: h
    REAL(KIND=dp), POINTER :: pc(:)
    LOGICAL :: Stat

    Basis = 0.0_dp
    dBasisdx = 0.0_dp

    lp = 0.5_dp * Particles % Length(No,1:2)

    t1 = -h - lp
    t2 = -h + lp
    t3 = -lp
    t4 = lp
    t5 = h - lp
    t6 = h + lp

    Sv = 0.0_dp
    dSv = 0.0_dp

    diff(:,1) = Particles % Coordinate(No,1) - Nodes % x
    diff(:,2) = Particles % Coordinate(No,2) - Nodes % y


    DO j = 1,2
       DO k = 1,4

          !basis
          IF ( (t1(j)<diff(k,j)) .AND. (diff(k,j) <= t2(j)) ) THEN
             Sv(k,j) = ( (h+lp(j)+diff(k,j))*(h+lp(j)+diff(k,j)))/(4.0_dp*h*lp(j))
             dSv(k,j) = (h+lp(j)+diff(k,j))/(2.0_dp*h*lp(j))
          ELSEIF ( (t2(j)<diff(k,j)) .AND. (diff(k,j) <= t3(j)) ) THEN
             Sv(k,j) = 1.0_dp + diff(k,j)/h
             dSv(k,j) = 1.0_dp/h
          ELSEIF ( (t3(j)<diff(k,j)) .AND. (diff(k,j) <= t4(j)) ) THEN
             Sv(k,j) = 1.0_dp - (diff(k,j)*diff(k,j) + lp(j)*lp(j)) / (2.0_dp*h*lp(j))
             dSv(k,j) = -diff(k,j)/(h*lp(j))
          ELSEIF ((t4(j)<diff(k,j)) .AND. (diff(k,j) <= t5(j))) THEN
             Sv(k,j) = 1.0_dp - diff(k,j)/h
             dSv(k,j) = -1.0_dp/h
          ELSEIF ((t5(j)<diff(k,j)) .AND. (diff(k,j) <= t6(j))) THEN
             Sv(k,j) = ( (h+lp(j)-diff(k,j))*(h+lp(j)-diff(k,j)))/(4.0_dp*h*lp(j))
             dSv(k,j) = -(h+lp(j)-diff(k,j))/(2.0_dp*h*lp(j))
             !ELSE
             !   Sv(k,j) = 0.0_dp
             !   dSv(k,j) = 0.0_dp
          END IF

       END DO
    END DO

    Basis = Sv(:,1) * Sv(:,2)
    dBasisdx(:,1) = Sv(:,2) * dSv(:,1)
    dBasisdx(:,2) = Sv(:,1) * dSv(:,2)
  END FUNCTION GIMPMElementInfoFromCoords


  !> GIMPM shape functions (element-wise, following Charlton et al., 2017)
  FUNCTION GIMPMElementInfo( pperm, Particles, Model,Element, Nodes, No, &
       detj, scale, calcscale, Basis, dBasisdx) RESULT (stat)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    TYPE(Element_t), TARGET :: Element
    TYPE(Element_t), POINTER :: PointElement
    TYPE(Nodes_t), OPTIONAL   :: Nodes
    REAL(KIND=dp) :: detj
    REAL(KIND=dp) :: Basis(4),pg(2),nx(4),ny(4), &
         pl(2),Lg(2),Ll(2),E1(2),E2(2),Sa(2),Sb(2),dSa(2),dSb(2)
    REAL(KIND=dp) :: weight,nxmax,nxmin,nymax,nymin,Sx,dSx,Sy,dSy,&
         N,S,E,W,scale,newpg(2),newLg(2)
    REAL(KIND=dp), OPTIONAL :: dBasisdx(4,3)
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: calcscale
    INTEGER :: pperm
    INTEGER :: nn,k,No,xx,yy
    LOGICAL :: Stat,ordered=.FALSE.,actuallycalcscale

    IF ( Element % TYPE % ElementCode .NE. 404) THEN
       CALL FATAL('GIMPMElementInfo','Currently only works for ElemType 404!!')
    END IF

    !just to be consistent with elementinfo and solvers
    Stat = .TRUE.

    IF ( (Particles % Usesavedbasis) .AND. (pperm .NE. 0) .AND. &
         ((Particles % Status(No) < PARTICLE_LEAVING) &
         .OR. (.NOT. calcscale) ) ) THEN

       !use saved basis
       Basis(1:4) = ElemParticles(Element % ElementIndex) % Basis(pperm,1:4)
       dBasisdx(1:4,1:3) = ElemParticles(Element % ElementIndex) % dBasisdx(pperm,1:4,1:3)

       detj = (Particles % Length(No,1)) * (Particles % Length(No,2))
       scale = 1.0_dp

    ELSE

       Basis(:) = 0.0_dp
       dBasisdx(:,:) = 0.0_dp

       nn = Element % Type % NumberOfNodes
       NodeIndexes => Element % NodeIndexes

       nx(1:nn) = Model % Mesh % Nodes % x(NodeIndexes)
       ny(1:nn) = Model % Mesh % Nodes % y(NodeIndexes)

       !global element max and min coords
       nxmax = MAXVAL(nx)
       nxmin = MINVAL(nx)

       nymax = MAXVAL(ny)
       nymin = MINVAL(ny)

       !Global particle coords
       pg(1) = Particles % Coordinate(No,1)
       pg(2) = Particles % Coordinate(No,2)

       !Global particle lengths
       Lg(1) = Particles % Length(No,1)
       Lg(2) = Particles % Length(No,2)

       scale = 1.0_dp

       IF (calcscale) THEN
          IF (Particles % Status(No) >= PARTICLE_LEAVING) THEN

             !we are overlapping a boundary.
             !we find where the current element is relative to the
             !particle element. Order of elements:
             ! 1-2-3
             ! 4-5-6
             ! 7-8-9
             ! where 5 is the particle element.
             !we do this using vector rc, where
             ! rc(1) is the row, rc(2) is the column

             N = Particles % Coordinate(No,2) + Particles % Length(No,2)/2.0_dp
             S = Particles % Coordinate(No,2) - Particles % Length(No,2)/2.0_dp
             E = Particles % Coordinate(No,1) + Particles % Length(No,1)/2.0_dp
             W = Particles % Coordinate(No,1) - Particles % Length(No,1)/2.0_dp

             IF (nxmax < pg(1)) THEN
                !  rc(2) = 1
                newLg(1) =  ABS(nxmax-W)
                newpg(1) = (nxmax+W)/2.0_dp
                xx=1
             ELSEIF (nxmin > pg(1)) THEN
                !  rc(2) = 3
                newLg(1) =  ABS(E-nxmin)
                newpg(1) = (E+nxmin)/2.0_dp
                xx=2
             ELSE
                !  rc(2) = 2
                newLg(1) = ABS(MIN(E,nxmax)-MAX(W,nxmin))
                newpg(1) = (MIN(E,nxmax)+MAX(W,nxmin))/2.0_dp
                xx=3
             END IF

             IF ( nymin > pg(2) ) THEN
                !  rc(1) = 1
                newLg(2) = ABS(N-nymin)
                newpg(2) = (N+nymin)/2.0_dp
                yy=1

             ELSEIF (nymax < pg(2)) THEN
                !  rc(1) = 3
                newLg(2) = ABS(nymax-S)
                newpg(2) = (nymax+S)/2.0_dp
                yy=2
             ELSE
                !  rc(1) = 2
                newLg(2) = ABS(MIN(N,nymax)-MAX(S,nymin))
                newpg(2) = (MIN(N,nymax)+MAX(S,nymin))/2.0_dp
                yy=3
             END IF

             scale = (newLg(1) * newLg(2))/Particles % GVolume(No)

             IF (scale < 0) PRINT *,'negative scale',scale

             pg(1) = newpg(1)
             pg(2) = newpg(2)

             Lg(1) = newLg(1)
             Lg(2) = newLg(2)

             IF (scale <= 0.0_dp) scale = 0.0_dp
          END IF
       END IF


       ! Local particle coords
       pl(1) = (2.0_dp * (pg(1) - nxmin) / (nxmax - nxmin) ) - 1.0_dp
       pl(2) = (2.0_dp * (pg(2) - nymin) / (nymax - nymin) ) - 1.0_dp

       !Local particle lengths
       !the 2 is necessary because local coords are -1 to 1,
       !not 0 to 1...
       Ll(1) = 2.0_dp*Lg(1)/(nxmax - nxmin)
       Ll(2) = 2.0_dp*Lg(2)/(nymax - nymin)

       DO k = 1,2
          E1(k) = MAX( -1.0_dp, pl(k) - ( Ll(k) / 2.0_dp ) )
          E2(k) = MIN(  1.0_dp, pl(k) + ( Ll(k) / 2.0_dp ) )

          Sa(k) = ( 1.0_dp / (4.0_dp * Ll(k) ) ) * &
               ( ( 2.0_dp * E2(k) ) - ( E2(k) * E2(k) ) - &
               ( 2.0_dp * E1(k) ) + ( E1(k) * E1(k) ) )

          Sb(k) = ( 1.0_dp / ( 4.0_dp * Ll(k) ) ) * &
               ( ( 2.0_dp * E2(k) ) + ( E2(k) * E2(k) ) - &
               ( 2.0_dp * E1(k) ) - ( E1(k) * E1(k) ) )

          IF ( PRESENT(dBasisdx) ) THEN
             dSa(k) = ( E1(k) - E2(k) ) / ( 2.0_dp * Lg(k) )
             dSb(k) = ( E2(k) - E1(k) ) / ( 2.0_dp * Lg(k) )
          END IF
       END DO

       DO k = 1,nn
          IF (nx(k) == nxmin) THEN
             Sx = Sa(1)
             IF ( PRESENT(dBasisdx) ) dSx = dSa(1)
          ELSE
             Sx = Sb(1)
             IF ( PRESENT(dBasisdx) ) dSx = dSb(1)
          END IF

          IF (ny(k) == nymin) THEN
             Sy = Sa(2)
             IF ( PRESENT(dBasisdx) ) dSy = dSa(2)
          ELSE
             Sy = Sb(2)
             IF ( PRESENT(dBasisdx) ) dSy = dSb(2)
          END IF

          Basis(k) = Sx * Sy

          IF ( PRESENT(dBasisdx) ) THEN
             dBasisdx(k,1) = dSx*Sy
             dBasisdx(k,2) = dSy*Sx
          END IF
       END DO

       ! detJ = Lg(1)*Lg(2)
       detJ = Particles % PVolume(No)

    END IF

  END FUNCTION GIMPMElementInfo

  !**************************************************************************

  !> linspace function (as in MATLAB)
  SUBROUTINE linspace(from, to, array)
    REAL(KIND=dp), intent(in) :: from, to
    REAL(KIND=dp), intent(out) :: array(:)
    REAL(KIND=dp) :: range
    INTEGER :: n, ii

    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
       array(1) = from
       return
    end if

    do ii=1, n
       array(ii) = from + range * (ii - 1) / (n - 1)
    end do
  END SUBROUTINE linspace

  !**************************************************************************

  !> Subroutine sets up some preliminary information needed for MPM
  !! Defines some global particle constants so that we can keep things clean and
  !! avoid calling them for multiple subroutines
  !! CAUTION: Some of these global parameters may no longer be actively used
  SUBROUTINE SetMPMParticlePreliminaries(Particles,Model,dim)

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: dim
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: MinCoord(3), MaxCoord(3), s(3)
    INTEGER :: ierr
    LOGICAL :: GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

    WRITE(SolverName, '(A)') 'InitializeParticles: SetParticlePreliminaries'

    Mesh => GetMesh()
    IF( .NOT. ASSOCIATED( Mesh ) ) THEN
       CALL Fatal('SetParticleDimensions','No Mesh associated')
    END IF

    IF( PRESENT( dim ) ) THEN
       IF( dim == 2 .OR. dim == 3 ) THEN
          Particles % dim = dim
       ELSE
          CALL Fatal('SetParticleDimensions','Invalid dimension')
       END IF
    ELSE
       Particles % dim = Mesh % Meshdim
    END IF

    ! Create list of faces / edges
    !Mesh => GetMesh()
    CALL FindMeshEdges( Mesh, .FALSE.)
    IF ( ParEnv % PEs > 1 ) THEN
       CALL SParEdgeNumbering(Mesh,Allmesh=.TRUE.)
       CALL SParFaceNumbering(Mesh,Allmesh=.TRUE.)
    END IF


    Particles % ShapeFunctions = ListGetString( Model % Constants,'Shape Functions',GotIt )
    IF (.NOT. GotIt) Particles % ShapeFunctions = 'gimpm'

    PRINT *,''
    PRINT *,''
    PRINT *,''
    PRINT *,''
    PRINT *,'SHAPE FUNCTIONS: ',Particles % ShapeFunctions
    PRINT *,''
    PRINT *,''
    PRINT *,''
    PRINT *,''


    Particles % DamageModel = ListGetString( Model % Constants,'Damage Model',GotIt )
    IF (.NOT. GotIt) Particles % DamageModel = 'creep'


    Particles % rhow = GetConstReal( Model % Constants, 'Water Density', GotIt )
    If (.NOT.GotIt) Then
       WRITE(Message,'(A)') 'Constant Water Density not found. &
            &Setting to 1.03225e-18'
       CALL INFO(SolverName, Message, level=20)
       Particles % rhow = 1.03225e-18_dp
    End if

    Particles % rhoi = GetConstReal( Model % Constants, 'Ice Density', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Ice Density = Real $" in constants section')

    Particles % SeaLevel = GetConstReal( Model % Constants, 'Sea Level', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Sea Level = Real $" in constants section')


    Particles % CriticalDamage = GetConstReal( Model % Constants, 'Critical Damage', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Critical Damage=Real $" in constants')

    Particles % CriticalDav = GetConstReal( Model % Constants, 'Critical DAv', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Critical Dav=Real $" in constants')


    Particles % mindam = GetConstReal( Model % Constants, 'Min Damage Threshold', GotIt )
    IF (.NOT. GotIt) THEN
       CALL WARN(SolverName,&
            'Min Damage Threshold=Real $" not defined in constants; setting to 1.0e-20')
       Particles % mindam = 1.0e-20_dp
    END IF


    Particles % useisodam = GetLogical( Model % Constants, &
         'Use Isotropic Damage for Initially Damaged Particles', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Isotropic Damage for Initially Damaged Particles not defined, assume false!')
       Particles % useisodam = .FALSE.
    END IF


    IF (Particles % ShapeFunctions == 'gimpm') THEN
       Particles % trackstrain = .FALSE.
    ELSE
       Particles % trackstrain = .TRUE.
    END IF
    !Particles % trackstrain = .TRUE.

    Particles % usetracer = GetLogical( Model % Constants, &
         'Use Tracer', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Tracer not defined, assume false!')
       Particles % usetracer = .FALSE.
    END IF


    Particles % alwayssplitfour = .FALSE.


    Particles % isodamcritdam = GetConstReal( Model % Constants, 'Iso Critical Damage', GotIt )
    IF (.NOT. GotIt) THEN
       CALL WARN(SolverName,&
            'Iso Critical Damage=Real $" not defined in constants; setting to criticaldamage')
       Particles % isodamcritdam = Particles % CriticalDamage
    END IF


    Particles % isodamcritdav = GetConstReal( Model % Constants, 'Iso Critical Dav', GotIt )
    IF (.NOT. GotIt) THEN
       CALL WARN(SolverName,&
            'Iso Critical Dav=Real $" not defined in constants; setting to criticaldav')
       Particles % isodamcritdav = Particles % CriticalDav
    END IF



    Particles % ah = GetConstReal( Model % Constants, 'ah', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "ah = Real $" in constants')

    Particles % Bf = GetConstReal( Model % Constants, 'Bf', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Bf = Real $" in constants')

    IF (Particles % Bf == -9999) THEN
       Particles % Bf = HUGE(1.0_dp)
    END IF



    Particles % Bh = GetConstReal( Model % Constants, 'Bh', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Bh = Real $" in constants')

    Particles % k1 = GetConstReal( Model % Constants, 'k1', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "k1 = Real $" in constants')

    Particles % k2 = GetConstReal( Model % Constants, 'k2', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "k2=Real $" in constants')

    Particles % gamma = GetConstReal( Model % Constants, 'gamma', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "gamma=Real $" in constants')

    Particles % currentgamma = Particles % gamma
    Particles % tempgamma = Particles % gamma



    Particles % DmaxI  = GetConstReal( Model % Constants, 'DMax I', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "DMax I=Real $" in constants')

    Particles % DmaxII  = GetConstReal( Model % Constants, 'DMax II', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "DMax II=Real $" in constants')

    Particles % DmaxIII  = GetConstReal( Model % Constants, 'DMax III', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "DMax III=Real $" in constants')

    Particles % Dmax = MAX(Particles % DmaxI,Particles % DmaxII)
    Particles % Dmax = MAX(Particles % DmaxIII,Particles % Dmax)


    IF (Particles % DmaxII>Particles % DmaxI) THEN
       Particles % DmaxII_Dom = .TRUE.
    ELSE
       Particles % DmaxII_Dom = .FALSE.
    END IF


    Particles % DavDmaxI  = GetConstReal( Model % Constants, 'Dav DMax I', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not define "Dav DMax I=Real $" in constants so setting to Dmax I')
       Particles % DavDmaxI = Particles % DmaxI
    END IF


    Particles % DavDmaxII  = GetConstReal( Model % Constants, 'Dav DMax II', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not define "Dav DMax II=Real $" in constants so setting to Dmax II')
       Particles % DavDmaxII = Particles % DmaxII
    END IF

    Particles % DavDmaxIII  = GetConstReal( Model % Constants, 'Dav DMax III', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not define "Dav DMax III=Real $" in constants so setting to Dmax III')
       Particles % DavDmaxIII = Particles % DmaxIII
    END IF

    Particles % gradualdyz =  GetLogical( Model % Constants, &
         'gradualdyz', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not define "gradualdyz = Logical" in constants so setting to False')
       Particles % gradualdyz = .FALSE.
    ELSE
       IF (Particles % gradualdyz) THEN
          PRINT *,'USING GRADUALDYZ'
       END IF
    END IF



    Particles % useriftdmax = GetLogical( Model % Constants, &
         'Use Rift Dmax', GotIt )
    IF (Particles % DavDmaxI == Particles % DavDmaxII) THEN
       PRINT *,'DavDmaxI == DavDMaxII, so automatically using rift dmax'
       Particles % useriftdmax = .TRUE.
    ELSE
       IF (.NOT. GotIt) THEN
          CALL Warn(SolverName, &
               'didnt define useriftdmax, assuming true!')
          Particles % useriftdmax = .TRUE.
       END IF
    END IF

    IF (Particles % DamageModel .NE. 'creep' .AND. .NOT. Particles % useriftdmax) THEN
       CALL Warn(SolverName, &
            'must use riftdmax with selected damage model. Setting useriftdmax = true!')
       Particles % useriftdmax = .TRUE.
    END IF


    Particles % riftdmax = GetConstReal(Model % Constants, &
         'Rift Dmax',GotIt)
    IF (.NOT. GotIt) THEN
       IF (Particles % useriftdmax) THEN
          CALL Warn(SolverName, &
               'no layer dmax defined, assume same as davdmaxI!')
       END IF
       Particles % riftdmax = Particles % DavDmaxI
    END IF

    IF (Particles % gamma == 0.0_dp) THEN
       IF (Particles % DmaxII>Particles % DmaxI) THEN
          CALL FATAL(SolverName,'DmaxII cannot currently exceed DmaxI when gamma==0!')
       END IF
    END IF

    Particles % noDzz = GetLogical( Model % Constants, &
         'No Dzz', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'No Dzz not defined, assume false!')
       Particles % noDzz = .FALSE.
    END IF

    Particles % forceDzz = GetLogical( Model % Constants, &
         'Force Dzz', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Force Dzz not defined, assume false!')
       Particles % forceDzz = .FALSE.
    END IF

    Particles % unstructuredmesh = GetLogical( Model % Constants, 'Unstructured mesh',GotIt)
    IF (.NOT. GotIt) Particles % unstructuredmesh = .FALSE.

    Particles % hoop = GetLogical( Model % Constants, 'hoop', GotIt )
    IF (.NOT. GotIt) Particles % hoop = .FALSE.

    IF (Particles % hoop) THEN
       Particles % unstructuredmesh = .TRUE.
    END IF


    Particles % nospin = GetLogical( Model % Constants, 'no damage spin', GotIt )
    IF (.NOT. GotIt) Particles % nospin = .FALSE.
    IF (Particles % nospin) PRINT *,'no spin for damage!!!'


    Particles % constlintemp = GetLogical( Model % Constants, &
         'Constant Linear Temperature', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Constant Linear Temperature not defined in constants, assuming false!')
       Particles % constlintemp = .FALSE.
    END IF

    Particles % useconsttemp = GetLogical( Model % Constants, &
         'Use Constant Temperature', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Constant Temperature not defined in constants, assuming false!')
       Particles % useconsttemp = .FALSE.
    END IF

    IF (Particles % useconsttemp) THEN
       Particles % consttemp = GetConstReal( Model % Constants, 'Constant Temperature', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "Constant Temperature=Real $" in constants')
    END IF

    Particles % setyvelzero = GetLogical( Model % Constants, &
         'SetYVelZero', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'SetYVelZero not defined in constants, assuming false!')
       Particles % SetYVelZero = .FALSE.
    END IF



    Particles % restrictdam = GetLogical( Model % Constants, &
         'Restrict Damage', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Restrict Damage not defined in constants, assuming false!')
       Particles % restrictdam = .FALSE.
    END IF


    Particles % rdxmin = GetConstReal( Model % Constants, 'Restrict Damage X Min', GotIt )
    IF (.NOT. GotIt)  Particles % rdxmin = -HUGE(1.0_dp)
    Particles % rdxmax = GetConstReal( Model % Constants, 'Restrict Damage X Max', GotIt )
    IF (.NOT. GotIt)  Particles % rdxmax = -HUGE(1.0_dp)
    Particles % rdymin = GetConstReal( Model % Constants, 'Restrict Damage Y Min', GotIt )
    IF (.NOT. GotIt)  Particles % rdymin = -HUGE(1.0_dp)
    Particles % rdymax = GetConstReal( Model % Constants, 'Restrict Damage Y Max', GotIt )
    IF (.NOT. GotIt)  Particles % rdymax = -HUGE(1.0_dp)


    Particles % nodamregion = GetLogical( Model % Constants, &
         'Use No Damage Region', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use No Damage Region not defined in constants, assuming false!')
       Particles % nodamregion = .FALSE.
    END IF

    IF (Particles % nodamregion) THEN

       Particles % ndxmin = GetConstReal( Model % Constants, 'No Damage Region X Min', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "No Damage Region X Min=Real $" in constants')

       Particles % ndxmax = GetConstReal( Model % Constants, 'No Damage Region X Max', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "No Damage Region X Max=Real $" in constants')

       Particles % ndymin = GetConstReal( Model % Constants, 'No Damage Region Y Min', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "No Damage Region Y Min=Real $" in constants')

       Particles % ndymax = GetConstReal( Model % Constants, 'No Damage Region Y Max', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "No Damage Region Y Max=Real $" in constants')
    END IF


    Particles % alwaysfemfront = GetLogical( Model % Constants, &
         'Always Use FEM at Front', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Did not specify Always Use FEM at Front in constants, assuming false!')
       Particles % alwaysfemfront = .FALSE.
    END IF

    Particles % prupt = .FALSE.


    Particles % analytictest = GetLogical( Model % Constants, &
         'Analytic Test', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'No Analytic Test = Logical in constants, assuming false!')
       Particles % analytictest = .FALSE.
    END IF

    Particles % firsttimestepzero = GetLogical( Model % Constants, &
         'First Timestep Zero', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'No First Timestep Zero = Logical in constants, assuming True!')
       Particles % firsttimestepzero = .TRUE.
    END IF


    Particles % flipvelfric = GetLogical( Model % Constants, &
         'FLIP updates for friction velocity', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'FLIP updates for friction velocity = Logical $ not found in constants, assuming False!')
       Particles % flipvelfric = .FALSE.
    END IF



    Particles % usesavedbasis = GetLogical( Model % Constants, &
         'Use Saved Basis', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Saved Basis not defined in constants, assuming false!')
       Particles % usesavedbasis = .FALSE.
    END IF

    IF (Particles % shapefunctions .NE. 'gimpm') THEN
       CALL Warn(SolverName,'Saved basis only works with gimpm, setting false!!')
       Particles % usesavedbasis = .FALSE.
    END IF

    Particles % uplag = .FALSE.

    Particles % VelocityDependentFriction = GetLogical( Model % Constants, &
         'Update Particle Velocities for Friction', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Update Particle Velocities for Friction not defined in constants, assuming false!')
       Particles % VelocityDependentFriction = .FALSE.
    END IF

    Particles % CoulombFriction = GetLogical( Model % Constants, &
         'Use Coulomb Friction', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Coulomb Friction not defined in constants, assuming false!')
       Particles % CoulombFriction = .FALSE.
    END IF



    Particles % isorift = GetLogical( Model % Constants, &
         'Rupture All Damage Components for Rift', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Rupture All Damage Components for Rift = Logical in constants, assuming True!')
       Particles % isorift = .TRUE.
    END IF

    IF (Particles % useisodam) THEN
       Particles % initDmax  = GetConstReal( Model % Constants, 'Iso Max Damage', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'If using Iso dam Need to define "Iso Max Damage=Real $" in constants')
    END IF


    Particles % rkmcritdam = GetLogical( Model % Constants, &
         'Use Critical Damage for RKM', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Critical Damage for RKM not defined, assume false!')
       Particles % rkmcritdam = .FALSE.
    END IF


    Particles % usedavtensor = GetLogical( Model % Constants, &
         'Use Dav Tensor for RKM', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Dav Tensor for RKM not defined, assume false!')
       Particles % usedavtensor = .FALSE.
    END IF


    Particles % fixdavdirections = GetLogical( Model % Constants, &
         'Fix Dav Principal Directions', GotIt )

    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Fix Dav Principal Directions not defined, assume false!')
       Particles % fixdavdirections = .FALSE.
    END IF



    Particles % modifiedmurakami = GetLogical( Model % Constants, &
         'Use Modified Murakami', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Modified Murakami not found, assume false!')
       Particles % modifiedmurakami = .False.
    END IF

    Particles % usetruecauchydamage = GetLogical( Model % Constants, &
         'Use True Cauchy Damage', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use True Cauchy Damage not found, assume false!')
       Particles % usetruecauchydamage = .False.
    END IF


    Particles % femforallnotfull = GetLogical( Model % Constants, &
         'Always fill not full elements', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Always fill not full elements not defined, assume false!')
       Particles % femforallnotfull = .FALSE.
    END IF


    Particles % usestaticparticles = GetLogical( Model % Constants, &
         'Use Static Particles', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Static Particles not defined in constants, assuming false!')
       Particles % usestaticparticles = .FALSE.
    END IF


    Particles % UseBCforPrevVel = GetLogical( Model % Constants, &
         'Use BC for PrevVel', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'UseBCforPrevVel not defined in constants, assuming true!')
       Particles % UseBCforPrevVel = .TRUE.
    END IF

    Particles % UseHBC = GetLogical( Model % Constants, &
         'Use Thickness BC', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Thickness BC not defined in constants, assuming false!')
       Particles % UseHBC = .TRUE.
    END IF


    Particles % FEMifGrounded = GetLogical( Model % Constants, &
         'Use FEM if Grounded', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use FEM if Grounded not defined in constants, assuming false!')
       Particles % FEMifGrounded = .FALSE.
    END IF

    ! --  not available in this version
    !these are features specific to Larsen C and
    !internal rift BCs...
    Particles % useriftmelange = .FALSE.
    Particles % rmelfrac = 0.0_dp
    Particles % larcmelfracmod = .FALSE.
    Particles % trackdamstrain = .FALSE.
    Particles % noriftbounduntilfullsep = .FALSE.
    Particles % frontmelange = .FALSE.
    ! ------


    Particles % initcrack = GetLogical( Model % Constants, &
         'initcrack', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'initcrack = $ not def in const, assumed false!')
       Particles % initcrack = .FALSE.
    END IF

    IF (Particles % initcrack) THEN
       Particles % icrackwidth = GetConstReal( Model % Constants, 'initcrack width', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "initcrack width=Real $" in constants')
       Particles % icrackwidth = Particles % icrackwidth * 0.5_dp

       Particles % icrackx1 = GetConstReal( Model % Constants, 'icrack x1', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "icrack x1=Real $" in constants')

       Particles % icrackx2 = GetConstReal( Model % Constants, 'icrack x2', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "icrack x2=Real $" in constants')

       Particles % icracky1 = GetConstReal( Model % Constants, 'icrack y1', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "icrack y1=Real $" in constants')

       Particles % icracky2 = GetConstReal( Model % Constants, 'icrack y2', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "icrack y2=Real $" in constants')

       Particles % bcracklayers = GetInteger( Model % Constants, &
            'basalcrack layers', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "Number of basalcrack Layers = Integer $" in constants')

       Particles % dbcracklayers = GetInteger( Model % Constants, &
            'diffuse basalcrack layers', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "Number of diffuse basalcrack Layers = Integer $" in constants')

    END IF



    Particles % mismip = GetLogical( Model % Constants, &
         'MISMIP Tests', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'MISMIP Tests = $ not def in const, assumed false!')
       Particles % mismip = .FALSE.
    END IF


    !dont allow Binit (viscosity from the mesh) on a particle decrease
    !below what existed at the first timestep
    Particles % binitlowerboundfromfirsttimestep = GetLogical( Model % Constants, &
         'Binit Lower Bound from First Timestep', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Binit Lower Bound from First Timestep = $ not def in const, assumed false!')
       Particles % binitlowerboundfromfirsttimestep = .FALSE.
    END IF

    !dont allow EF (the mesh-based enhancement factor) on a particle decrease
    !below what existed at the first timestep
    Particles % efboundsfromfirsttimestep = GetLogical( Model % Constants, &
         'EF bounds from First Timestep', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'EF bounds from First Timestep = $ not def in const, assumed false!')
       Particles % efboundsfromfirsttimestep = .FALSE.
    END IF

    Particles % usegiveneta = GetLogical( Model % Constants, &
         'Use Given Eta', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName, &
            'Use Given Eta not defined in constants, assuming false!')
       Particles % usegiveneta = .FALSE.
    ELSE
       Particles % giveneta = GetConstReal( Model % Constants, 'Given Eta', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName,&
            'Need to define "Given Eta=Real $" in constants')
    END IF

    !use 4 particles per cell (ppc) in part of the domain, 9 ppc elsewhere.
    !for some reason, this may be broken
    Particles % mix49 =  GetLogical( Model % constants,'mixfrac4and9',GotIt)
    IF (.NOT. GotIt) THEN
       Particles % mix49 = .FALSE.
    END IF


    Particles % rf = GetConstReal( Model % Constants, 'rf', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "rf=Real $" in constants')

    Particles % nodaminflow = GetLogical( Model % constants,'no evolve damage on inflow particles',GotIt)
    IF (.NOT. GotIt) THEN
       Call Warn(SolverName,&
            'Did not specify "no evolve damage on inflow particles = Logical" in constants. Assuming false')
       Particles % nodaminflow = .FALSE.
    END IF

    Particles % mismipmelt2 = GetLogical( Model % constants,'mismipmelt 2 test',GotIt)
    IF (.NOT. GotIt) THEN
       Call Warn(SolverName,&
            'Did not specify "mismipmelt2 test = Logical" in constants. Assuming false')
       Particles % mismipmelt2 = .FALSE.
    END IF

    Particles % outputdbassis = GetLogical( Model % constants,'output dbassis',GotIt)
    IF (.NOT. GotIt) THEN
       Call Warn(SolverName,&
            'Did not specify "output dbassis = Logical" in constants. Assuming false')
       Particles % outputdbassis = .FALSE.
    END IF

    Particles % sthres = GetConstReal( Model % Constants, 'sthres', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "sthres=Real $" in constants')

    Particles % Gravity = GetConstReal( Model % Constants, 'Gravity', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Gravity=Real $" in constants')


    Particles % NumberOfTemperatureLayers = GetInteger( Model % Constants, &
         'Number of Temperature Layers', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Number of Temperature Layers = Real $" in constants')


    Particles % NumberOfParticleLayers = GetInteger( Model % Constants, &
         'Number of Particle Layers', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Number of Particle Layers = Real $" in constants')

    IF (Particles % NumberOfParticleLayers < 8) THEN
       Particles % NumberOfParticleLayers = 8
       CALL Warn(Solvername,&
            'Automatically increasing number of particle layers to 8, as required')
    END IF


    IF (Particles % DamageModel == 'zero stress') THEN
       Particles % NumberOfParticleLayers = Particles % NumberOfTemperatureLayers
    END IF


    Particles % gridres = GetConstReal( Model % Constants,'Grid Resolution',gotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "gridres = Real $" in constants')


    Particles % fillpercent = GetConstReal( Model % Constants,'FEM fill element under percent',gotIt )
    IF (.NOT. GotIt) THEN
       Particles % fillpercent = 0.0_dp
       CALL Warn(SolverName,&
            'Did not define "FEM fill element under percent" = Real $, setting to 0.0')
    END IF


    Particles % CriticalShearRate = GetConstReal( Model % Constants, &
         'Critical Shear Rate', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Critical Shear Rate = Real $" in constants')

    Particles % MaximumTimeStep = GetConstReal( Model % Constants, 'Maximum Time Step', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Maxtimum Time Step=Real $" in constants')

    Particles % ViscosityExponent = GetConstReal( Model % Constants, &
         'Viscosity Exponent', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName,&
         'Need to define "Viscosity Exponent = Real $" in constants')

    Particles % maxlength = GetConstReal( Model % Constants, &
         'Maximum Particle Length', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Maximum Particle Length = Real $" in constants')

    Particles % maxDPlength = GetConstReal( Model % Constants, &
         'Maximum Damaged Particle Length', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Maximum Damaged Particle Length = Real $" in constants')

    Particles % maxGLlength = GetConstReal( Model % Constants, &
         'Maximum Grounding Line Particle Length', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Maximum Grounding Line Particle Length = Real $" in constants')

    Particles % davsplitthres = GetConstReal( Model % Constants, 'Dav Split Threshold', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Dav Split Threshold = Real $" in constants section')


    Particles % buffer = GetInteger( Model % Constants, 'Number of Buffer Particles', GotIt )
    IF (.NOT. GotIt) Particles % buffer = 0


    Particles % DInitTolerance = GetConstReal( Model % Constants, 'Dinit Tolerance', GotIt )
    If (.NOT.GotIt) Then
       WRITE(Message,'(A)') 'Constant >Dinit Tolerance< not found. &
            &Setting to 1.0e0-10_dp'
       CALL INFO(SolverName, Message, level=20)
       Particles % DInitTolerance=1.0e-10_dp
    END IF

    Particles % WarningTol = GetConstReal( Model % Constants, &
         'Dinit Warning Tolerance', GotIt )
    If (.NOT.GotIt) Then
       WRITE(Message,'(A)') 'Constant >Dinit Warning Tolerance< not found. &
            &Setting to 1.0e-04_dp'
       CALL INFO(SolverName, Message, level=20)
       Particles % WarningTol=1.0e-4_dp
    END IF

    Particles % DInitNumIters = GetInteger( Model % Constants, 'Dinit Iters', GotIt )
    If (.NOT.GotIt) Then
       WRITE(Message,'(A)') 'Constant >Dinit Iters< not found. &
            &Setting to 1.0e05_dp'
       CALL INFO(SolverName, Message, level=20)
       Particles % DInitNumIters = 100000
    END IF


    Particles % elementfraction = GetConstReal( Model % Constants, 'Particle Element Fraction', GotIt )
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Particle Element Fraction = Real $frac" in constants')

    Particles % nodamfront = GetLogical( Model % Constants,'Ignore Front Damage',GotIt)
    IF (.NOT. GotIt) Particles % nodamfront = .FALSE.


    Particles % usedamage = GetLogical( Model % Constants, 'Use Damage', GotIt)
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not specify "Use Damage = Logical" in constants, so assuming true!!')
       Particles % usedamage = .FALSE.
    END IF


    Particles % movegl = GetLogical( Model % Constants, 'Move GL', GotIt)
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not specify "Move GL = Logical" in constants, so assuming false!!')
       Particles % movegl = .FALSE.
    END IF

    Particles % constmb = GetLogical( Model % Constants, 'Constant MB Parameter', GotIt)
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not specify "Constant MB Parameter= Logical" in constants, so assuming false!!')
       Particles % constmb = .FALSE.
    END IF


    Particles % constef = GetLogical( Model % Constants, 'Constant EF Parameter', GotIt)
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not specify "Constant EF Parameter= Logical" in constants, so assuming false!!')
       Particles % constef = .FALSE.
    END IF


    Particles % constfric = GetLogical( Model % Constants, 'Constant Friction Parameter', GotIt)
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not specify "Constant Friction Parameter = Logical" in constants, so assuming false!!')
       Particles % constfric = .FALSE.
    END IF


    Particles % SEP = GetLogical( Model % Constants, 'Use SEP', GotIt)
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not specify "Use SEP = Logical" in constants, so assuming false!!')
       Particles % SEP = .FALSE.
    END IF

    Particles % SimpleAdvectTest = GetLogical( Model % Constants, 'SimpleAdvectTest', GotIt)
    IF (.NOT. GotIt) THEN
       CALL Warn(SolverName,&
            'Did not specify "SimmpleAdvectTest = Logical" in constants, so assuming false!!')
       Particles % SimpleAdvectTest = .FALSE.
    END IF

    Particles % IncNumPAtStart = GetLogical( Model % Constants, &
         'Increase Number Of Particles At Start', GotIt)
    IF (.NOT. GotIt) THEN
       Particles % IncNumPAtStart = .TRUE.
    END IF

    Particles % UseOnePPC = GetLogical( Model % Constants, &
         'Use One PPC Edit', GotIt)
    IF (.NOT. GotIt) THEN
       Particles % UseOnePPC = .False.
    END IF



    IF (Particles % UseOnePPC) THEN
       Particles % OnePPCXEnd = GetConstReal( Model % Constants, &
            'One PPC XEnd', GotIt )
       IF (.NOT. GotIt) THEN
          CALL Fatal(SolverName,'Need to specify One PPC X End in constants')
       END IF
       PRINT *,'USING ONE PPC BEHIND X = ',Particles % OnePPCXEnd
    END IF

  END SUBROUTINE SetMPMParticlePreliminaries

  !**************************************************************************

  !> Subroutine initializes particle locations and allocates space for particle vars
  !! Largely borrowed from ParticleUtils.F90
  !! Parallel does not work yet (for this subroutine nor for the greater MPM model in general)
  SUBROUTINE InitializeParticles( Particles, Model, InitParticles, AppendParticles )
    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: InitParticles
    LOGICAL, OPTIONAL :: AppendParticles
    TYPE(Model_t) :: Model
    TYPE(ValueList_t), POINTER :: Params, BodyForce
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t), POINTER :: CurrentElement,Element
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t) :: Nodes
    INTEGER :: NewParticles,LastParticle,NoElements
    INTEGER :: dim, ElementIndex, body_id, bf_id,ElementInd
    REAL(KIND=dp), POINTER :: rWork(:,:),Coordinate(:,:),Velocity(:,:)
    REAL(KIND=dp) :: Velo(3), Coord(3), Center(3), time0, dist,ymove,xmove,xc(2),yc(2),cm(4)
    INTEGER :: i,ii,j,k,l,n,vdofs,nonodes, InitStatus, TotParticles
    INTEGER, POINTER :: MaskPerm(:) => NULL(), InvPerm(:) => NULL(), NodeIndexes(:) => NULL()
    LOGICAL :: Found, GotIt, GotMask, RequirePositivity, ContinuousNumbering, GotWeight
    REAL(KIND=dp), POINTER :: InitialValues(:,:)
    REAL(KIND=dp) :: mass,boltz,temp,coeff,eps,frac,meanval ,loc,gridres,len
    REAL(KIND=dp) :: MinCoord(3), MaxCoord(3), Diam, DetJ, MinDetJ, MaxDetJ, &
         MinWeight, MaxWeight, MeanWeight
    REAL(KIND=dp), POINTER :: MaskVal(:) => NULL()
    REAL(KIND=dp), ALLOCATABLE :: Weight(:)
    INTEGER :: nx,ny,nz,nmax,ix,iy,iz,ind,No,np
    LOGICAL :: CheckForSize, Parallel,shiftxlocstoelementedges,Stat
    LOGICAL, POINTER :: DoneParticle(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, str
    REAL(KIND=dp) :: loc1,loc2,loc3,loc4,stag,stag2
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: Basis(4), dBasisdx(4,3), ddBasisddx(4,3,3)

    INTEGER :: ntheta
    REAL(KIND=dp) :: r,r1,r2,Area,thetaseg,H,numazimelems,pmaskval
    REAL(KIND=dp), ALLOCATABLE :: theta(:)
    TYPE(Variable_t), POINTER :: GridH=>NULL()

    SAVE :: Nodes,ntheta

    Mesh => GetMesh()
    Params => ListGetSolverParams()
    dim = Particles % Dim
    Parallel = ( ParEnv % PEs > 1 )

    GotMask = .FALSE.
    VariableName = ListGetString( Params,'Initialization Condition Variable',GotIt )

    IF(GotIt) THEN
       RequirePositivity = .TRUE.
       Var => VariableGet( Mesh % Variables, TRIM(VariableName) )
       IF( .NOT. ASSOCIATED( Var ) ) THEN
          CALL Fatal('InitializeParticles','Mask / Condition variable does not exist!')
       END IF

       MaskPerm => Var % Perm
       MaskVal => Var % Values

       IF(.NOT. ( ASSOCIATED( MaskPerm ) .AND. ASSOCIATED(MaskVal)) ) THEN
          CALL Warn('InitializeParticles','Initialization variable does not exist?')
       ELSE IF( MAXVAL( MaskPerm ) == 0 ) THEN
          CALL Warn('InitializeParticles','Initialization variable of size zero?')
          nonodes = 0
          noelements = 0
          InvPerm => NULL()
       ELSE
          GotMask = .TRUE.

          ALLOCATE (InvPerm( Mesh % NumberOfBulkElements ))
          InvPerm = 0

          j = 0
          DO i=1,Mesh % NumberOfBulkElements
             CurrentElement => Mesh % Elements(i)
             IF (CheckPassiveElement(CurrentElement)) CYCLE
             NodeIndexes =>  CurrentElement % NodeIndexes
             n = CurrentElement % TYPE % NumberOfNodes

             IF( ANY( MaskPerm( NodeIndexes )==0 ) ) CYCLE

             !   IF (ANY (MaskVal( MaskPerm(NodeIndexes) )<0.0_dp)) CYCLE
             IF (ALL(MaskVal( MaskPerm(NodeIndexes) )<0.0_dp)) CYCLE
             ! If some of bulk elements have been found avtive
             j = j + 1
             InvPerm(j) = i
          END DO
          noelements = j

          PRINT *,'Total elements vs. masked',Mesh % NumberOfBulkElements,noelements
       END IF
    ELSE
       nonodes = Mesh % NumberOfNodes
       noelements = Mesh % NumberOfBulkElements
    END IF

    ! Now decide on the number of particles.
    !-------------------------------------------------------------------------

    frac = Particles % elementfraction

    IF( GotIt ) THEN
       NewParticles = NINT( frac * noelements )
    ELSE
       CALL Fatal('InitializeParticles','Need to specify particle element fraction!')
    END IF

    IF( ParEnv% PEs == 1 ) THEN
       TotParticles = NewParticles
    ELSE
       TotParticles = NINT( ParallelReduction( 1.0_dp * NewParticles ) )
    END IF

    IF( TotParticles == 0 ) THEN
       CALL Fatal('InitializeParticles','No Particles to Initialize')
    ELSE
       WRITE( Message,'(A,I0)') 'Number of Particles: ',TotParticles
       CALL Info('InitializeParticles',Message,Level=6)
    END IF

    shiftxlocstoelementedges = GetLogical( Params,'Shift X Locs to Element Edges',GotIt )

    IF (shiftxlocstoelementedges) THEN
       NewParticles = NewParticles + NINT(sqrt(frac))
    END IF


    ! If there are no particles in this partition, nothing to do
    !-------------------------------------------------------------------------
    IF( NewParticles == 0 ) RETURN

    IF (Particles % IncNumPAtStart) THEN
       LastParticle = FLOOR(NewParticles * 1.5)
    ELSE
       LastParticle = NewParticles
    END IF

    IF (Particles % mix49) THEN
       k = 0
       DO i = 1,noelements
          CurrentElement => Mesh % Elements(i)
          NodeIndexes =>  CurrentElement % NodeIndexes
          IF (MAXVAL(Mesh % Nodes % x(NodeIndexes)) > 300000.0_dp) THEN
             k = k+9
          ELSE
             k = k+4
          END IF
       END DO

       LastParticle = INT(FLOOR(REAL(k) * 1.25_dp))
    END IF


    gridres = GetConstReal( Model % Constants,'Grid Resolution',gotIt )
    IF (.NOT. GotIt) CALL Fatal('InitializeParticles', &
         'Need to define "gridres = Real $" in constants')

  IF (Particles % hoop .AND. Particles % shapefunctions == 'smpm') THEN

       GridH => VariableGet( Mesh % Variables,'Hinit',GotIt)
       IF (.NOT. GotIt) THEN
          CALL FATAL('cant find hinit','fail')
       END IF


       No = 0
       CALL AllocateParticles( Particles,Model,LastParticle)
       Particles % NumberOfParticles = LastParticle

       len = Particles % gridres/sqrt(frac)
       r1 = 30000.0-0.5_dp*len
       r2 = 70000.0_dp !120000.0

       numazimelems = 41.0_dp
       ntheta = INT(sqrt(frac)*numazimelems)
       ALLOCATE(theta(ntheta))
       thetaseg = (10.0_dp*(pi/180.0))/ntheta
       theta(1) = 0.5_dp*thetaseg

       !for rotated mesh
       theta(1) = theta(1) - (5.0_dp*(pi/180.0))

       DO i = 2,ntheta
          theta(i)=theta(i-1)+thetaseg
       END DO

       r = r1
       i = 1
       DO No = 1,LastParticle
          r = r+len

          IF (r>r2) THEN
             r = r1+len
             i = i+1
             IF (i>ntheta) EXIT
          END IF

          Particles % Coordinate(No,1) = r*COS(theta(i))
          Particles % Coordinate(No,2) = r*SIN(theta(i))

          Coord = GetParticleCoord( Particles, No)

          ElementInd = 0
          CALL LocateParticleInMeshOctree( ElementInd, Coord(1:3) )

          IF (ElementInd == 0) THEN
             PRINT *,'No',No
             PRINT *,'coords',Coord(:)
             Particles % Status(No) = PARTICLE_LOST
             Particles % Coordinate(No,:) = -10000.00
             CYCLE
          END IF


          CurrentElement => Mesh % Elements(ElementInd)
          NodeIndexes =>  CurrentElement % NodeIndexes
          n = CurrentElement % TYPE % NumberOfNodes


          IF (No == 1) THEN
             ALLOCATE(Nodes % x(n),Nodes % y(n),Nodes % z(n))
          END IF

          Area = ElementArea(Model % Mesh,CurrentElement,n)

          Area = Area/frac

          Particles % PVolume(No) = Area
          Particles % GVolume(No) = Area

          Particles % Length(No,1) = (MAXVAL(Mesh % Nodes % x(NodeIndexes))&
               -MINVAL(Mesh % Nodes % x(NodeIndexes)))/(sqrt(frac))
          Particles % Length(No,2) = Area/Particles % Length(No,1)
          Particles % OrigLength(No,1:2) = Particles % Length(No,1:2)

          IF (sqrt(Coord(1)*Coord(1) + Coord(2)*Coord(2)) <= 70000.0_dp) THEN
             Particles % H(No) = 400.0_dp
          ELSE
             Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
             Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
             Nodes % z(1:n) = 0.0_dp

             stat = sMPMElementInfo( CurrentElement, Particles, Model, Nodes, No, &
                  Particles % gridres, Basis,dBasisdx)
             CALL GetScalarFieldInMesh(GridH, CurrentElement, Basis, H)

             Particles % H(No) = H
          END IF

          Particles % EF(No) = 1.0_dp
          Particles % Mass(No) = Particles % pvolume(No) * &
               Particles % H(No) * Particles % rhoi
          Particles % binit(No) = Particles % giveneta
          Particles % bz(No,:) = Particles % giveneta
          Particles % Gmask(No) = 1.0_dp
          Particles % origno(No) = No
          Particles % F(No,:) = 0.0_dp
          Particles % F(No,1) = 1.0_dp
          Particles % F(No,2) = 1.0_dp
          Particles % ElementIndex(No) = ElementInd
          Particles % InterpElem(No) = ElementInd

       END DO

       IF (No<LastParticle) THEN
          Particles % Status(No:LastParticle) = PARTICLE_LOST
          Particles % Coordinate(No:LastParticle,:) = -10000000.0_dp
       END IF

       PRINT *,'number of particles after start', No-1

       RETURN
    END IF

    ! Allocate particles
    !-------------------------------------------------------------------------

    CALL AllocateParticles( Particles, Model, LastParticle )

    Particles % NumberOfParticles = LastParticle

    Coordinate => Particles % Coordinate

    CALL Info('InitializeParticles',&
         'Initializing particles evenly among elements',Level=10)

    Particles % NumberOfParticles = NewParticles
    PRINT *,'Initializing particles in elements:',NewParticles,noelements

    IF (Particles % mix49) THEN

       k=0

       DO i = 1,noelements

          j = InvPerm(i)
          CurrentElement => Mesh % Elements(j)
          NodeIndexes =>  CurrentElement % NodeIndexes
          n = CurrentElement % TYPE % NumberOfNodes


          IF (MAXVAL(Mesh % Nodes % x(NodeIndexes)) <= 300000.0_dp) THEN
             loc = gridres/4.0_dp

             DO ii=1,4

                k = k + 1
                IF (ii == 1) THEN
                   ymove = loc
                   xmove = loc
                ELSE IF (ii==2)  THEN
                   ymove = loc
                   xmove = -loc
                ELSE IF (ii==3) THEN
                   ymove = -loc
                   xmove = loc
                ELSE IF (ii==4) THEN
                   ymove = -loc
                   xmove = -loc
                END IF

                IF( j > Mesh % NumberOfBulkElements ) THEN
                   PRINT *,'j too large',j,i,k,(NoElements-1)*(i-1)/(NewParticles-1)+1
                END IF

                Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove
                Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
                IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n
                ! Only a bulk element may own a particle
                IF( j <= Mesh % NumberOfBulkElements ) THEN
                   Particles % ElementIndex(k) = j
                END IF
             END DO

          ELSE

             loc = gridres/3.0_dp

             DO ii=1,9

                k = k + 1
                j = InvPerm(i)

                IF (ii == 1) THEN
                   ymove = loc
                   xmove = -loc
                ELSE IF (ii==2)  THEN
                   ymove = loc
                   xmove = 0.0_dp
                ELSE IF (ii==3) THEN
                   ymove = loc
                   xmove = loc
                ELSE IF (ii==4) THEN
                   ymove = 0.0_dp
                   xmove = -loc
                ELSE IF (ii == 5) THEN
                   ymove = 0.0_dp
                   xmove = 0.0_dp
                ELSE IF (ii==6)  THEN
                   ymove = 0.0_dp
                   xmove = loc
                ELSE IF (ii==7) THEN
                   ymove = -loc
                   xmove = -loc
                ELSE IF (ii==8) THEN
                   ymove = -loc
                   xmove = 0.0_dp
                ELSE IF (ii==9) THEN
                   ymove = -loc
                   xmove = loc
                END IF

                IF( j > Mesh % NumberOfBulkElements ) THEN
                   PRINT *,'j too large',j,i,k,(NoElements-1)*(i-1)/(NewParticles-1)+1
                END IF

                Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove
                Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
                IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n

                ! Only a bulk element may own a particle
                IF( j <= Mesh % NumberOfBulkElements ) THEN
                   Particles % ElementIndex(k) = j
                END IF
             END DO
          END IF
       END DO

    ELSE

       IF (frac == 1.0_dp) THEN
          loc = gridres/2.0_dp
          k=0
          DO i = 1,noelements

             k = k + 1
             j = InvPerm(i)

             ymove = 0.0_dp !loc
             xmove = 0.0_dp !loc

             IF( j > Mesh % NumberOfBulkElements ) THEN
                PRINT *,'j too large',j,i,k,(NoElements-1)*(i-1)/(NewParticles-1)+1
             END IF

             CurrentElement => Mesh % Elements(j)
             NodeIndexes =>  CurrentElement % NodeIndexes
             n = CurrentElement % TYPE % NumberOfNodes
             Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove
             Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
             IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n

             ! Only a bulk element may own a particle
             IF( j <= Mesh % NumberOfBulkElements ) THEN
                Particles % ElementIndex(k) = j
             END IF
          END DO

       ELSE IF (frac == 4.0_dp) THEN
          loc = gridres/4.0_dp
          !loc = 166.0_dp
          k=0
          DO i = 1,noelements


             IF (Particles % UseOnePPC) THEN
                j = InvPerm(i)
                CurrentElement => Mesh % Elements(j)
                NodeIndexes =>  CurrentElement % NodeIndexes
                n = CurrentElement % TYPE % NumberOfNodes

                IF (ANY(Mesh % Nodes % x(NodeIndexes(1:n))<Particles % OnePPCXEnd)) THEN
                   k = k + 1
                   Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n
                   Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n
                   IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n
                   IF( j <= Mesh % NumberOfBulkElements ) THEN
                      Particles % ElementIndex(k) = j
                   END IF

                   CYCLE
                END IF
             END IF

             DO ii=1,4

                j = InvPerm(i)
                k = k + 1
                IF (ii == 1) THEN
                   ymove = loc
                   xmove = loc
                ELSE IF (ii==2)  THEN
                   ymove = loc
                   xmove = -loc
                ELSE IF (ii==3) THEN
                   ymove = -loc
                   xmove = loc
                ELSE IF (ii==4) THEN
                   ymove = -loc
                   xmove = -loc
                END IF

                IF( j > Mesh % NumberOfBulkElements ) THEN
                   PRINT *,'j too large',j,i,k,(NoElements-1)*(i-1)/(NewParticles-1)+1
                END IF

                CurrentElement => Mesh % Elements(j)
                NodeIndexes =>  CurrentElement % NodeIndexes
                n = CurrentElement % TYPE % NumberOfNodes
                Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove
                Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
                IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n

                ! Only a bulk element may own a particle
                IF( j <= Mesh % NumberOfBulkElements ) THEN
                   Particles % ElementIndex(k) = j
                END IF
             END DO

             IF (shiftxlocstoelementedges .AND. i == noelements) THEN

                DO ii = 1,2
                   k = k+1

                   IF (ii == 1) THEN
                      xmove = -loc
                      ymove = loc

                   ELSEIF (ii==2) THEN
                      xmove = -loc
                      ymove = -loc
                   END IF

                   CurrentElement => Mesh % Elements(j)
                   NodeIndexes =>  CurrentElement % NodeIndexes
                   n = CurrentElement % TYPE % NumberOfNodes
                   Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove + gridres
                   Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
                   IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n
                END DO

                Coordinate(:,1) = Coordinate(:,1) - loc
             END IF
          END DO


       ELSE IF (frac == 9.0_dp) THEN
          loc = gridres/3.0_dp
          !loc = 166.0_dp
          k=0

          stag2 = 0.0_dp !loc/2.0_dp
          stag = 0.0_dp !loc/3.0_dp
          DO i = 1,noelements
             DO ii=1,9

                k = k + 1
                j = InvPerm(i)

                IF (ii == 1) THEN
                   ymove = loc
                   xmove = -loc
                ELSE IF (ii==2)  THEN
                   ymove = loc
                   xmove = 0.0_dp
                ELSE IF (ii==3) THEN
                   ymove = loc
                   xmove = loc
                ELSE IF (ii==4) THEN
                   ymove = 0.0_dp
                   xmove = -loc +stag2
                ELSE IF (ii == 5) THEN
                   ymove = 0.0_dp
                   xmove = 0.0_dp +stag2
                ELSE IF (ii==6)  THEN
                   ymove = 0.0_dp
                   xmove = loc +stag2
                ELSE IF (ii==7) THEN
                   ymove = -loc
                   xmove = -loc -stag
                ELSE IF (ii==8) THEN
                   ymove = -loc
                   xmove = 0.0_dp -stag
                ELSE IF (ii==9) THEN
                   ymove = -loc
                   xmove = loc -stag
                END IF

                IF( j > Mesh % NumberOfBulkElements ) THEN
                   PRINT *,'j too large',j,i,k,(NoElements-1)*(i-1)/(NewParticles-1)+1
                END IF

                CurrentElement => Mesh % Elements(j)
                NodeIndexes =>  CurrentElement % NodeIndexes
                n = CurrentElement % TYPE % NumberOfNodes
                Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove
                Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
                IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n

                ! Only a bulk element may own a particle
                IF( j <= Mesh % NumberOfBulkElements ) THEN
                   Particles % ElementIndex(k) = j
                END IF
             END DO


             IF (shiftxlocstoelementedges .AND. i == noelements) THEN

                DO ii = 1,3
                   k = k+1

                   IF (ii == 1) THEN
                      xmove = -loc
                      ymove = loc

                   ELSEIF (ii==2) THEN
                      xmove = -loc
                      ymove = 0.0_dp
                   ELSEIF (ii==3) THEN
                      xmove = -loc
                      ymove = -loc
                   END IF

                   CurrentElement => Mesh % Elements(j)
                   NodeIndexes =>  CurrentElement % NodeIndexes
                   n = CurrentElement % TYPE % NumberOfNodes
                   Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove + gridres
                   Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
                   IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n
                END DO

                Coordinate(:,1) = Coordinate(:,1) - loc
             END IF
          END DO


       ELSE IF (frac == 16.0_dp) THEN
          loc = gridres/8.0_dp
          !loc = 166.0_dp
          loc1 = loc
          loc2 = 3*loc
          loc3 = 5*loc
          loc4 = 7*loc
          k=0
          DO i = 1,noelements
             DO ii=1,16

                k = k + 1
                j = InvPerm(i)

                IF (ii == 1) THEN
                   ymove = loc1
                   xmove = loc1
                ELSE IF (ii==2)  THEN
                   ymove = loc2
                   xmove = loc1
                ELSE IF (ii==3) THEN
                   ymove = loc3
                   xmove = loc1
                ELSE IF (ii==4) THEN
                   ymove = loc4
                   xmove = loc1
                ELSEIF (ii == 5) THEN
                   ymove = loc1
                   xmove = loc2
                ELSE IF (ii==6)  THEN
                   ymove = loc2
                   xmove = loc2
                ELSE IF (ii==7) THEN
                   ymove = loc3
                   xmove = loc2
                ELSE IF (ii==8) THEN
                   ymove = loc4
                   xmove = loc2
                ELSEIF (ii == 9) THEN
                   ymove = loc1
                   xmove = loc3
                ELSE IF (ii==10)  THEN
                   ymove = loc2
                   xmove = loc3
                ELSE IF (ii==11) THEN
                   ymove = loc3
                   xmove = loc3
                ELSE IF (ii==12) THEN
                   ymove = loc4
                   xmove = loc3
                ELSEIF (ii == 13) THEN
                   ymove = loc1
                   xmove = loc4
                ELSE IF (ii==14)  THEN
                   ymove = loc2
                   xmove = loc4
                ELSE IF (ii==15) THEN
                   ymove = loc3
                   xmove = loc4
                ELSE IF (ii==16) THEN
                   ymove = loc4
                   xmove = loc4
                END IF

                IF( j > Mesh % NumberOfBulkElements ) THEN
                   PRINT *,'j too large',j,i,k,(NoElements-1)*(i-1)/(NewParticles-1)+1
                END IF

                CurrentElement => Mesh % Elements(j)
                NodeIndexes =>  CurrentElement % NodeIndexes
                n = CurrentElement % TYPE % NumberOfNodes
                Coordinate(k,1) = MINVAL(Mesh % Nodes % x(NodeIndexes)) + xmove
                Coordinate(k,2) = MINVAL(Mesh % Nodes % y(NodeIndexes)) + ymove
                !Coordinate(k,1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n + xmove
                !Coordinate(k,2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n + ymove
                IF( dim == 3 ) Coordinate(k,3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n

                ! Only a bulk element may own a particle
                IF( j <= Mesh % NumberOfBulkElements ) THEN
                   Particles % ElementIndex(k) = j
                END IF
             END DO
          END DO
       ELSE
          CALL Fatal('InitializeParticles', &
               'Particle element fraction can only be 1,4, 9,or 16 currently')
       END IF
    END IF


    DEALLOCATE (InvPerm)

    PRINT *,'maxcoord',MAXVAL(Particles % Coordinate(1:Particles % NumberOfParticles,1))


  END SUBROUTINE InitializeParticles

  !**************************************************************************

  !>MakeInterpLayers holds the interpolation map and functions to
  !!interpolate between the temperature layers
  !!(e.g. variable Temperature with dofs = # temperature layers)
  !!particle layers (# of particles layers is specified in constants).
  !!InterpLayers % Map(particlelayer,1) and InterpLayers % Map(particlelayer,2) give the
  !!two temperature layers that a particle layers lies between.
  !!InterpLayers % InterpFun(particlelayer,1) and InterpLayers % InterpFun(particlelayer,2)
  !!give the two weights to interpolate between the particlelayer and the two closest temperature
  !!layers specified in InterpLayers % Map(particlelayer,1:2)
  SUBROUTINE MakeInterpLayers(Particles, Model)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    INTEGER :: numoflayers1,numoflayers2,H1,H2,ii,jj,I
    REAL(KIND = dp), POINTER :: x1(:),x2(:),z1(:),z2(:),interp1(:),interp2(:)
    REAL(KIND = dp) :: H,intsize1,intsize2
    INTEGER, POINTER :: L1(:),L2(:)
    LOGICAL :: GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

    WRITE(SolverName, '(A)') 'MakeInterpLayers'

    numoflayers1 = Particles % NumberOfTemperatureLayers
    numoflayers2 = Particles % NumberOfParticleLayers

    H1 = numoflayers1-1

    H = DBLE(H1)

    !spacing 1 (temperature)
    ALLOCATE(x1(numoflayers1),z1(numoflayers1))

    x1 = (/ (I,I=1, numoflayers1) /)

    intsize1 = H/DBLE(H1)
    z1 = intsize1*(DBLE(x1)-1.0_dp)

    !spacing 2 (damage)
    ALLOCATE(x2(numoflayers2),z2(numoflayers2))

    H2 = numoflayers2-1

    x2 = (/ (I,I=1, numoflayers2) /)
    intsize2 = H/DBLE(H2)
    z2 = intsize2*(DBLE(x2)-1.0_dp)

    !interpolation vars
    ALLOCATE(interp1(numoflayers2),interp2(numoflayers2))
    ALLOCATE(L1(numoflayers2),L2(numoflayers2))

    DO ii = 1,numoflayers2
       DO jj = 2,numoflayers1

          IF (z2(ii) <= z1(jj)) THEN

             L1(ii) = jj-1
             L2(ii) = jj

             interp1(ii) = 1.0_dp-(z2(ii)-z1(jj-1))/intsize1
             interp2(ii) = 1.0_dp-(z1(jj) - z2(ii))/intsize1
             EXIT
          END IF
       END DO
    END DO

    InterpLayers % NumOfTLayers = numoflayers1
    InterpLayers % NumOfDLayers = numoflayers2

    ALLOCATE(InterpLayers % Map(numoflayers2,2), InterpLayers % InterpFun(numoflayers2,2))

    InterpLayers % Map(:,1) = L1
    InterpLayers % Map(:,2) = L2

    PRINT *,'min map1',MINVAL(InterpLayers % Map(:,1))
    PRINT *,'max map1',MAXVAL(InterpLayers % Map(:,1))

    PRINT *,'min map2',MINVAL(InterpLayers % Map(:,2))
    PRINT *,'max map2',MAXVAL(InterpLayers % Map(:,2))

    InterpLayers % InterpFun(:,1) = interp1
    InterpLayers % InterpFun(:,2) = interp2

    DEALLOCATE(x1,z1)
    DEALLOCATE(x2,z2)
    DEALLOCATE(interp1,interp2)
    DEALLOCATE(L1,L2)

  END SUBROUTINE MakeInterpLayers

  !**************************************************************************

  !> Initialize particle variable values,
  !! which as specified in the sif or interpolated from the mesh
  SUBROUTINE InitParticleVars(Particles, Model, Solver )

    USE ElementDescription
    USE Lists

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Nodes_t)   :: ElementNodes
    TYPE(Model_t) :: Model
    TYPE(Element_t), POINTER :: BulkElement
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: GridInvVisc=>NULL(),GridBtz=>NULL(),GridH=>NULL(), &
         GridVelo=>NULL(), GridExtraDamage=>NULL(), GridFP=>NULL(), GridZs=>NULL(), &
         GridDMask=>NULL(), GridEF=>NULL(),GridBed=>NULL(),GridStatic=>NULL(),&
         DamTrajMask,Slope,damtraj1,damtraj2
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: No, Status,ElementIndex,NoParticles,nn,n,layers,btzlayers,ii,jj
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:),newviscz(:),Dam(:)
    REAL(KIND=dp) :: Coord(3), dtime,StressThres,SqrtElementMetric, &
         anisoparam, gridres, rhow, rhoi, Dmax, btzav,g,newdamval
    REAL(KIND=dp) :: InvVisc,Btz,ExtraDamage,H,FP,Dmask,Zs,GradZs(3),&
         Velocity(3),GradVelocity(3,3),x(7),y(7),dist,prevdist,Bedrock, GradH(3),&
         pdamtrajmask,pslope,pdamtraj1,pdamtraj2
    REAL(KIND=dp) :: frac,D,EF,Binit,criticaldamage,criticaldav,maxpdamage,damh,&
         fricparam,efparam,mbparam,sealevel,Hf,static,dscale,LarCRiftWidth
    REAL(KIND=dp) :: EigenVec(2,2),EigVals(2),strainrate(2,2),ww,xx,yy,zz,en
    LOGICAL :: Stat,GotIt,LarC,NoInitDam,LarCDamTraj,noprevdamatrift,&
         ConstantMB,ConstantEF,ConstantFric,testfedit,testpass,bumptest,test1d
    REAL(KIND=dp) :: xbump,hbump,abump,bbump,cbump,ge
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName,SolverName
    REAL(KIND=dp) :: cm,eta,Exx,Eyy,Ezz,Exy,Ee,Tau(3,3),tp,TT,DD,Q,Qw,dw,hab,dtot,ds,db

    REAL(KIND=dp) :: H0 = 600.0_dp,v0 = 300.0_dp,Q0,termA,termB,C,tempx,temph,dh
    LOGICAL :: damrift,isodamatrift,pir
    REAL(KIND=dp) :: dsave,cslope,ang,DDD(2,2),rot(2,2)

    INTEGER :: i
    REAL(KIND=dp) :: dirichletmax
    TYPE(Variable_t), POINTER :: GH,GHi,GVel1,GVel1i
    REAL(KIND=dp), POINTER :: GHVal(:),GHiVal(:),GVel1Val(:),GVel1iVal(:)
    INTEGER, POINTER :: GHPerm(:),GHiPerm(:),GVel1Perm(:),GVel1iPerm(:)


    Params => GetSolverParams()
    Mesh => GetMesh()
    n = Mesh % MaxElementNodes

    WRITE(SolverName, '(A)') 'InitParticleVars'

    ALLOCATE(ElementNodes % x(n),ElementNodes % y(n),ElementNodes % z(n))

    IF (Particles % useisodam) THEN
       Dmax = Particles % initDmax
       CriticalDamage = Particles % isodamcritdam
       !initDmax
       CriticalDav = Particles % isodamcritdav
       !initDmax
    ELSE
       Dmax  = Particles % DmaxI
       CriticalDamage = Particles % CriticalDamage
       CriticalDav = Particles % CriticalDav
    END IF


    anisoparam = Particles % gamma
    gridres = Particles % gridres
    rhow = Particles % rhow
    rhoi = Particles % rhoi
    layers = Particles % numberofparticlelayers
    frac = Particles % elementfraction
    sealevel = Particles % sealevel

    g = ABS(Particles % gravity)

    ConstantMB = Particles % constmb
    ConstantEF = Particles % constef
    ConstantFric = Particles % constfric

    IF (ConstantMB) THEN
       mbparam = GetConstReal( Model % Constants, 'mbparam', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "mbparam = Real $mbparam" in constants')
       Particles % MB(:) = mbparam
    END IF

    IF (ConstantEF) THEN
       efparam = GetConstReal( Model % Constants, 'efparam', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "efparam = Real $efparam" in constants')

       Particles % EF(:) = efparam
    END IF

    IF (ConstantFric) THEN
       fricparam = GetConstReal( Model % Constants, 'fricparam', GotIt )
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "fricparam = Real $fricparam" in constants')

       Particles % FP(:) = fricparam
    END IF


    Particles % simpleinitdam = GetLogical(Params,'simple init dam',GotIt)
    IF (Particles % simpleinitdam) THEN
       PRINT *,'initializing damage evenly for all layers'
    END IF

    noprevdamatrift = GetLogical(Params,'no prev dam at rift',GotIt)
    IF (noprevdamatrift) THEN
       PRINT *,'not including previous damage at rift'
    END IF


    isodamatrift = GetLogical(Params,'isodam at rift',GotIt)
    IF (.NOT. GotIt) isodamatrift = .TRUE.
    IF (isodamatrift) THEN
       PRINT *,'isodam at rift'
       Particles % initriftdam = .TRUE.
    ELSE
       Particles % initriftdam = .FALSE.
    END IF

    test1d = GetLogical( Params,'Test1d',GotIt)
    IF (.NOT. GotIt) THEN
       Call Warn(SolverName,&
            'Did not specify "Test1d = Logical" in Params, so assuming false!!')
       test1d = .FALSE.
    END IF

    testpass = GetLogical( Params,'Test Pass',GotIt)
    IF (.NOT. GotIt) THEN
       Call Warn(SolverName,&
            'Did not specify "Test Pass = Logical" in Params, so assuming false!!')
       testpass = .FALSE.
    END IF

    testfedit = GetLogical( Params,'Test Fedit',GotIt)
    IF (.NOT. GotIt) THEN
       Call Warn(SolverName,&
            'Did not specify "Test Fedit = Logical" in Params, so assuming false!!')
       testfedit = .FALSE.
    END IF

    NoInitDam = GetLogical( Model % Constants,'No Init Dam',GotIt)
    IF (.NOT. GotIt) THEN
       Call Warn(SolverName,&
            'Did not specify "not Init Dam= Logical" in Constants, so assuming false!!')
       NoInitDam = .FALSE.
    END IF

    IF (Particles % usestaticparticles) THEN
       VariableName = ListGetString(Params,'Static Particle Variable Name',GotIt)
       IF( GotIt ) THEN
          GridStatic => VariableGet( Mesh % Variables, TRIM(VariableName) )
          IF(.NOT. ASSOCIATED( GridStatic ) ) THEN
             CALL Fatal(SolverName,'Static Particle variable does not exist: '//TRIM(VariableName))
          END IF
       END IF
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "Static Particle Variable Name = String ---- ')
    END IF


    VariableName = ListGetString(Params,'InvVisc Variable Name',GotIt)
    IF( GotIt ) THEN
       GridInvVisc => VariableGet( Mesh % Variables, TRIM(VariableName) )
       IF(.NOT. ASSOCIATED( GridInvVisc ) ) THEN
          CALL Fatal(SolverName,'InvVisc variable does not exist: '//TRIM(VariableName))
       END IF
    END IF
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "InvVisc Variable Name = String ---- ')

    VariableName = ListGetString(Params,'Velocity Variable Name',GotIt)
    IF( GotIt ) THEN
       GridVelo => VariableGet( Mesh % Variables, TRIM(VariableName) )
       IF(.NOT. ASSOCIATED( GridVelo ) ) THEN
          CALL Fatal(SolverName, &
               'Velocity field variable does not exist: '//TRIM(VariableName))
       END IF
    END IF
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Velocity Variable Name = String ---- ')


    VariableName = ListGetString(Params,'Additional Initial D Variable Name',GotIt)
    IF( GotIt ) THEN
       GridExtraDamage => VariableGet( Mesh % Variables, TRIM(VariableName) )
    END IF

    VariableName = ListGetString(Params,'Thickness Variable Name',GotIt)
    IF( GotIt ) THEN
       GridH => VariableGet( Mesh % Variables, TRIM(VariableName) )
       IF(.NOT. ASSOCIATED( GridH ) ) THEN
          CALL Fatal(SolverName, &
               'Thickness variable does not exist: '//TRIM(VariableName))
       END IF
    END IF
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Thickness Variable Name = String ---- ')


    IF (.NOT. ConstantFric) THEN
       VariableName = ListGetString(Params,'Friction Parameter Name',GotIt)
       IF( GotIt ) THEN
          GridFP => VariableGet( Mesh % Variables, TRIM(VariableName) )
          IF(.NOT. ASSOCIATED( GridFP ) ) THEN
             CALL Fatal(SolverName, &
                  'Friction Parameter does not exist: '//TRIM(VariableName))
          END IF
       END IF
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "Friction Parameter Name = String ---- ')
    END IF


    VariableName = ListGetString(Params,'Surface Height Variable Name',GotIt)
    IF( GotIt ) THEN
       GridZs => VariableGet( Mesh % Variables, TRIM(VariableName) )
       IF(.NOT. ASSOCIATED( GridZs ) ) THEN
          CALL Fatal(SolverName, &
               'Surface Height Variable does not exist: '//TRIM(VariableName))
       END IF
    END IF
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Surface Height Variable Name = String ---- ')


    VariableName = ListGetString(Params,'Damage Mask Variable Name',GotIt)
    IF( GotIt ) THEN
       GridDMask => VariableGet( Mesh % Variables, TRIM(VariableName) )
       IF(.NOT. ASSOCIATED( GridDMask ) ) THEN
          CALL Fatal(SolverName, &
               'Damage Mask Variable does not exist: '//TRIM(VariableName))
       END IF
    END IF
    IF (.NOT. GotIt) CALL Fatal(SolverName, &
         'Need to define "Damage Mask Variable Name = String ---- ')


    IF (.NOT. ConstantEF) THEN
       VariableName = ListGetString(Params,'EF Variable Name',GotIt)
       IF( GotIt ) THEN
          GridEF => VariableGet( Mesh % Variables, TRIM(VariableName) )
          IF(.NOT. ASSOCIATED( GridEF ) ) THEN
             CALL Fatal(SolverName, &
                  'EF Variable does not exist: '//TRIM(VariableName))
          END IF
       END IF
       IF (.NOT. GotIt) CALL Fatal(SolverName, &
            'Need to define "EF Variable Name = String ---- ')
    END IF

    VariableName = ListGetString(Params,'Bedrock Variable Name',GotIt)
    IF( GotIt ) THEN
       GridBed => VariableGet( Mesh % Variables, TRIM(VariableName) )
    END IF

    IF (test1d) THEN
        GHi => VariableGet(Model % Mesh % Variables, 'Hinit' )
        GHiPerm => GHi % Perm
        GHiVal => GHi % Values

        GVel1i => VariableGet(Model % Mesh % Variables, 'InitVel 1' )
        GVel1iPerm => GVel1i % Perm
        GVel1iVal => GVel1i % Values

        GVel1 => VariableGet(Model % Mesh % Variables, 'SSAVelocity 1' )
        GVel1Perm => GVel1 % Perm
        GVel1Val => GVel1 % Values

        GH => VariableGet(Model % Mesh % Variables, 'H' )
        GHPerm => GH % Perm
        GHVal => GH % Values

        dirichletmax = GetConstReal(Solver % Values,'dirichlet max x',GotIt)
    ENDIF



    ALLOCATE( Basis(n), dBasisdx(n, 3),newviscz(layers),Dam(layers) )


!!! ------ some features not available in this version .... !!!
    !ignore them for now
    Particles % LarCFixDavTest = .FALSE.
    damrift = .FALSE.
    pir = .FALSE.
!!! -------------------------------------------------------- !!!


!!!!------------------------------ Procedure -------------------------------!!!!!!

    NoParticles = Particles % NumberOfParticles

    PRINT *,'NoParticles',NoParticles

    Particles % Status(1:NoParticles) = PARTICLE_ACTIVE
    Particles % InterpElem(:) = Particles % ElementIndex(:)

    IF (Particles % mix49) THEN
       DO No = 1,NoParticles

          IF (Particles % Coordinate(No,1) > 300000.0_dp) THEN
             Particles % Length(No,:) = gridres/3.0_dp
          ELSE
             Particles % Length(No,:) = gridres/2.0_dp
          END IF
       END DO
    END IF


    IF (Particles % ShapeFunctions == 'gimpm') THEN

       IF (.NOT. Particles % mix49) THEN
          IF (frac == 1.0_dp) THEN
             Particles % Length(:,:) = gridres
          ELSE IF (frac == 4.0_dp) THEN
             Particles % Length(:,:) = gridres/2.0_dp
          ELSE IF (frac == 9.0_dp) THEN
             Particles % Length(:,:) = gridres/3.0_dp
          ELSE IF (frac == 16.0_dp) THEN
             Particles % Length(:,:) = gridres/4.0_dp
          ELSE
             CALL Fatal(SolverName, &
                  'Particle Element Fraction can currently only be 16,9,4, or 1')
          END IF

          IF (Particles % UseOnePPC) THEN
             WHERE (Particles % Coordinate(:,1)<Particles % OnePPCXEnd) Particles % Length(:,1) = gridres
             WHERE (Particles % Coordinate(:,1)<Particles % OnePPCXEnd) Particles % Length(:,2) = gridres
          END IF
       END IF


       Particles % OrigLength = Particles % Length

       Particles % GVolume(:) = Particles % Length(:,1) * &
            Particles % Length(:,2)
       Particles % PVolume(:) = Particles % GVolume(:)

    ELSE

       IF (frac == 1.0_dp) THEN
          Particles % PVolume = gridres*gridres
          Particles % GVolume = Particles % PVolume
       ELSE IF (frac == 4.0_dp) THEN
          Particles % PVolume = (gridres/2.0_dp)*(gridres/2.0_dp)
          Particles % GVolume = Particles % PVolume
       ELSE IF (frac == 9.0_dp) THEN
          Particles % PVolume = (gridres/3.0_dp)*(gridres/3.0_dp)
          Particles % GVolume = Particles % PVolume
       ELSE IF (frac == 16.0_dp) THEN
          Particles % PVolume = (gridres/4.0_dp)*(gridres/4.0_dp)
          Particles % GVolume = Particles % PVolume
       ELSE
          CALL Fatal(SolverName, &
               'Particle Element Fraction can currently only be 16,9,4, or 1')
       END IF
    END IF


    DO No = 1,NoParticles
       ElementIndex = GetParticleElement( Particles, No )
       IF( ElementIndex == 0 ) CYCLE
       BulkElement => Mesh % Elements( ElementIndex )

       Coord = GetParticleCoord( Particles, No )

       NodeIndexes => BulkElement % NodeIndexes
       nn = BulkElement % TYPE % NumberofNodes
       CALL GetElementNodes(ElementNodes,BulkElement)

       stat = sMPMElementInfo( BulkElement,Particles, Model, ElementNodes, No, &
            Particles % gridres, Basis,dBasisdx)

       CALL GetScalarFieldInMesh(GridDmask, BulkElement, Basis, Dmask )

       Particles % Gmask(No) = Dmask
    END DO

    IF (Particles % movegl) Particles % Gmask = -Particles % Gmask


    IF (particles % usegiveneta) THEN
       Particles % Bz(:,:) = Particles % giveneta
    ELSE
       !get temperature
       PRINT *,'interpolating temperature to particles...'
       CALL MPMMeshVectorToParticle( Particles, Model, 5, 1)
       PRINT *,'temperature interpolation complete.'
    END IF

    DO No = 1, NoParticles


       Particles % OrigNo(No) = No

       Dam = 0.0_dp
       ElementIndex = GetParticleElement( Particles, No )
       IF( ElementIndex == 0 ) CYCLE
       BulkElement => Mesh % Elements( ElementIndex )
       IF(.NOT. ASSOCIATED( BulkElement ) ) CYCLE

       Coord = GetParticleCoord( Particles, No )

       NodeIndexes => BulkElement % NodeIndexes
       nn = BulkElement % TYPE % NumberofNodes

       CALL GetElementNodes(ElementNodes,BulkElement)

       stat = sMPMElementInfo( BulkElement, Particles, Model, ElementNodes, No, &
            Particles % gridres, Basis,dBasisdx)

       IF( .NOT. stat ) THEN
          CALL Warn(SolverName,'Particle not in element')
          CYCLE
       END IF

       CALL GetScalarFieldInMesh(GridInvVisc, BulkElement, Basis, InvVisc )


       IF (ASSOCIATED(GridExtraDamage)) THEN
          CALL GetScalarFieldInMesh(GridExtraDamage, BulkElement, Basis, ExtraDamage )
       ELSE
          ExtraDamage = 0.0_dp
       END IF

       CALL GetScalarFieldInMesh(GridH, BulkElement, Basis, H, dBasisdx, GradH )
       CALL GetScalarFieldInMesh(GridDmask, BulkElement, Basis, Dmask )
       CALL GetScalarFieldInMesh(GridZs,BulkElement, Basis, Zs, dBasisdx, GradZs )
       CALL GetVectorFieldInMesh(GridVelo,BulkElement, Basis, Velocity, &
            dBasisdx, GradVelocity)


       IF (.NOT. ConstantFric) THEN
          CALL GetScalarFieldInMesh(GridFP, BulkElement, Basis, FP )
          Particles % FP(No) = FP
       END IF


       IF (.NOT. ConstantEF) THEN
          CALL GetScalarFieldInMesh(GridEF, BulkElement, Basis, EF )
          ! Particles % EF(No) = MIN(EF,1.0_dp)
          Particles % EF(No) = MAX(EF,0.0_dp)
       END IF


       IF (ASSOCIATED(GridBed)) THEN
          CALL GetScalarFieldInMesh(GridBed,BulkElement, Basis, Bedrock )
          Particles % Bedrock(No) = Bedrock
       ELSE
          IF (Particles % SEP) THEN
             CALL Fatal(SolverName,'Must Define Bedrock Var if using SEP!')
          END IF
       END IF


       Particles % Velocity(No,1:2) = Velocity(1:2)
       Particles % GridVelocity(No,1:2) = Velocity(1:2)
       Particles % GradVel(No,1) = GradVelocity(1,1)
       Particles % GradVel(No,2) = GradVelocity(2,2)
       Particles % GradVel(No,3) = GradVelocity(1,2)
       Particles % GradVel(No,4) = GradVelocity(2,1)
       Particles % GradZs(No,1:2) = GradZs(1:2)
       Particles % GradH(No,1:2) = GradH(1:2)

       !deformation gradient should be identity to start
       Particles % F(No,1) = 1.0_dp
       Particles % F(No,2) = 1.0_dp
       Particles % F(No,3) = 0.0_dp
       Particles % F(No,4) = 0.0_dp

       H = MAX(H,1.0_dp)

       ! IF (test1d) THEN
       !    cm = 1.0_dp/3.0_dp
       !    secondsperyear = 31556926.0_dp
       !    H0 = 600.0_dp
       !    v0 = 300.0_dp
       !    Q0 = H0*v0
       !    B0 = 1.9E8_dp
       !    A = ((B0*1.0E-6_dp)**(-3.0_dp))*secondsperyear !Mpa^(-3) a^(-1)
       !    C = (((910.0_dp*1.0e-6_dp*9.81_dp)/&
       !         (4.0_dp*(A**(-cm))))*(1.0_dp-910.0_dp/1028.0_dp))**3.0_dp
       !    !C is the weertman constant !C =2.45E-18; !m?3 s?1
       !    EeExp = (cm-1.0_dp)/2.0_dp
       !    Acm = A**(-cm)
       !    m1 = 4.0_dp*C/Q0
       !    m2 = 1.0_dp/(H0*H0*H0*H0)
       !    Ha = (m1*Particles % Coordinate(No,1) + m2)**(-0.25_dp)
       !    !Velocity
       !    Va = Q0/Ha
       !    Exx = C*Ha*Ha*Ha
       !    Particles % Velocity(No,1) = Va
       !    Particles % H(No) = Ha
       ! END IF

       Particles % H(No) = H
       Particles % Mass(No) = Particles % pvolume(No) * Particles % H(No) * rhoi


       ! VISCOSITY AND DAMAGE INITIALIZATION
       ! invvisc and btz are both in form B = A^(-1/3).

       newviscz = Particles % bz(No,:)

       IF (layers >1) THEN
          btzav = (((SUM(newviscz)-newviscz(1)/2.0_dp-newviscz(layers)/2.0_dp)*H/(DBLE(layers)-1.0_dp) ))/H
       ELSE
          btzav = Particles % bz(No,1)
       END IF


       Particles % Binit(No) = btzav

       btzav = btzav * Particles % EF(No)**(-1.0_dp/3.0_dp)

       Particles % Gmask(No) = Dmask


       IF (Particles % SEP) THEN
          IF (.NOT. Particles % movegl) THEN
             CALL Fatal(SolverName,'Must define movegl = true if using SEP!')
          END IF

          IF (Particles % Gmask(No) > 0.99_dp) THEN
             Particles % Gmask(No) = 1.0_dp
          ELSE IF (Particles % Gmask(No) < -0.99_dp) THEN
             Particles % Gmask(No) = -1.0_dp
          ELSE
             Hf = Particles % rhow * &
                  (Particles % sealevel-Particles % bedrock(No))/Particles % rhoi
             IF (Particles % H(No) .LT. Hf) THEN
                Particles % Gmask(No) = 1.0_dp
             ELSE
                Particles % Gmask(No) = -1.0_dp
             END IF
          END IF
       END IF



       IF (noinitdam) InvVisc = btzav

       IF (Dmask < 0.9_dp) THEN
          D = 0.0_dp
          Particles % Binit(No) = invvisc
          Particles % EF(No) = 1.0_dp
       ELSE

          D = 1.0_dp- (InvVisc/btzav)
       END IF


       IF (.NOT. Particles % SimpleAdvectTest) THEN
          IF (Particles % Gamma == 0.0_dp) THEN
             IF (ExtraDamage > 0.5_dp) THEN
                IF (D > 0.01_dp) THEN

                   D = D*2.0_dp
                   IF (D >= CriticalDav) D = Dmax
                END IF
             END IF
          END IF


          IF (ExtraDamage > 0.5_dp) THEN
             D = Dmax
          END IF
       ELSE

          D = ExtraDamage
       END IF



       IF (Particles % usestaticparticles) THEN

          CALL GetScalarFieldInMesh(GridStatic,BulkElement, Basis, Static )

          IF (Static > 0.95_dp) THEN
             Particles % Static(No) = .TRUE.
             Particles % Damage(No,:,:) = 0.0_dp
             Particles % dD(No,:,:) = 0.0_dp
             Particles % GradVel(No,:) = 0.0_dp
             Particles % Velocity(No,:) = 0.0_dp
             Particles % GridVelocity(No,:) = 0.0_dp
             Particles % xpic(No,:) = 0.0_dp
             Particles % EF(No) = 1.0_dp
             Particles % MB(No) = 0.0_dp
             Particles % Gmask(No) = 1.0_dp
          ELSE
             Particles % Static(No) = .FALSE.
          END IF
       ELSE
          Particles % Static(No) = .FALSE.
       END IF


       IF (damrift) Particles % DamStatus(No) = 1

       IF (Particles % UseIsoDam .AND. (.NOT. damrift)) THEN

          IF (ANY(Particles % Dav(No,:)>0.0_dp) .OR. Particles % Gmask(No) < 0.9_dp) THEN

             !isodam is now damstatus, where damstatus = -1 is the equivalent of
             !specifying isodam on the particle.
             !if DamStatus = 0, regular damage rules apply
             !is DamStatus = 1, the particle is fully damaged and can be evolved
             !using just the spin tensor

             IF ((.NOT. LarCDamTraj) .OR. Particles % Gmask(No)<0.9_dp) THEN
                Particles % DamStatus(No) = -1
             END IF
          END IF
       END IF


       damrift = .FALSE.
    END DO


    IF (Particles % initcrack) THEN

       !the D1 direction
       cslope = -(Particles % icrackx2 - Particles % icrackx1)/&
            (Particles % icracky2 - Particles % icracky1)

       ang = ATAN(cslope)

       rot(1,1) = COS(ang)
       rot(1,2) = -SIN(ang)
       rot(2,1) = SIN(ang)
       rot(2,2) = COS(ang)

       DDD = 0.0_dp
       DDD(1,1) = 1.0_dp

       DDD = MATMUL(MATMUL(rot,DDD),TRANSPOSE(rot))


       Particles % damage = 0.0_dp
       Particles % dD = 0.0_dp
       Particles % Dav = 0.0_dp

       DO No = 1,Particles % numberofparticles

          dist = DistToLineBetweenTwoPoints(Particles % Coordinate(No,1),&
               Particles % Coordinate(No,2),Particles % icrackx1,Particles % icracky1,&
               Particles % icrackx2,Particles % icracky2)

          IF (Particles % Coordinate(No,2)>Particles % icracky2) THEN
             dist = sqrt((Particles % Coordinate(No,1)-Particles % icrackx2)* &
                  (Particles % Coordinate(No,1)-Particles % icrackx2) + &
                  (Particles % Coordinate(No,2)-Particles % icracky2)* &
                  (Particles % Coordinate(No,2)-Particles % icracky2) )
          END IF

          IF (dist <= Particles % icrackwidth) THEN
             Particles % Damage(No,1:Particles % bcracklayers,1) = DDD(1,1)
             Particles % Damage(No,1:Particles % bcracklayers,2) = DDD(2,2)
             Particles % Damage(No,1:Particles % bcracklayers,4) = DDD(1,2)
          ELSE
             Particles % Damage(No,1:Particles % dbcracklayers,1) = DDD(1,1)
             Particles % Damage(No,1:Particles % dbcracklayers,2) = DDD(2,2)
             Particles % Damage(No,1:Particles % dbcracklayers,4) = DDD(1,2)
          END IF


          CALL VertIntDamFromVisc(Particles, No, Particles % NumberOfParticleLayers,Model)
          CALL FixPrincipalDamageVertInt(No,Model)
       END DO
    END IF


    CALL MPMParticlesToNodes( Particles, Model, 2)

  IF (Test1D) THEN
     print *,'TEST1d FIXING INITIAL MESH VALS'
     Do i = 1,Model % Mesh % NumberOfNodes
        IF (Model % Mesh % Nodes % x(i) <= dirichletmax) THEN
           GHVal(GHPerm(i)) = GHiVal(GHiPerm(i))
           GVel1Val(GVel1Perm(i)) = GVel1iVal(GVel1iPerm(i))
        END IF

     END DO
  END IF


    IF ( ASSOCIATED( Basis )    )  DEALLOCATE( Basis )
    IF ( ASSOCIATED( dBasisdx ) )  DEALLOCATE(dBasisdx)
    IF ( ASSOCIATED( newviscz ) ) DEALLOCATE( newviscz )
    IF ( ASSOCIATED( Dam ) )      DEALLOCATE ( Dam )
  END SUBROUTINE InitParticleVars


  !**************************************************************************

  !> Distance to a line specified between two points +/- 1 km
  !! used for Larsen C initialization
  FUNCTION DistToLineBetweenTwoPoints(x0,y0,x1,y1,x2,y2) RESULT(dist)

    REAL(KIND=dp) :: x0,y0,x1,x2,y1,y2
    REAL(KIND=dp) :: xmax,xmin,ymax,ymin,dist1,dist2,dist

    xmax = max(x1,x2)
    xmin = min(x1,x2)
    ymax = max(y1,y2)
    ymin = min(y1,y2)

    xmax = xmax+1000.0_dp
    ymax = ymax+1000.0_dp

    xmin = xmin-1000.0_dp
    ymin = ymin-1000.0_dp

    IF (xmax <= x0 .OR. xmin >= x0 .OR. ymax <= y0 .OR. ymin >= y0) THEN

       dist1 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1))
       dist2 = sqrt((x0-x2)*(x0-x2) + (y0-y2)*(y0-y2))
       dist = min(dist1,dist2)
    ELSE
       dist = abs((y2-y1)*x0-(x2-x1)*y0 + x2*y1 - y2*x1)/ &
            sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) )
    END IF
  END FUNCTION DistToLineBetweenTwoPoints



  !**************************************************************************

  !> Initializes a vertical damage distribution from specified 2-D damage.
  !! For floating ice only (or best for floating ice).
  !! Assumes (almost) all vertical layers are fully-damaged except
  !! a series of layers corresponding to some thickness of undamaged ice at floatation,
  !! so that both basal and surface crevasses appear (similar to effective ice thickness in
  !! Keller and Hutter, 2014). The thickness of the undamaged section is determined using
  !! a simple bisection scheme. When depth-averaged, the determined vertical damage profile
  !! should match the specified 2-D damage to some specified level of accuracy, but whether
  !! this accuracy is achieved may depend on the # of layers, the temperature distribution,
  !! and critical values of 2-D and 3-D damage.
  !! This scheme probably requires some edits for specific applications.

  SUBROUTINE InitParticleDz(Particles, No,numoflayers,btzav,invvisc,Dav,Model,Dam,dinit1)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    INTEGER :: No,numoflayers,ii,kk,last,bestiter,numiters
    REAL(KIND=dp), POINTER :: Dam(:)
    REAL(KIND=dp) :: btzav,invvisc,Dav,rhoi,rhow
    REAL(KIND=dp) :: Dmax,hw,tol,warningtol,CriticalDamage, &
         CriticalDav, zs, zb,intsize,dav2,dmax2,qqq,Yav,undamhigh,undamlow, &
         zrange,resid,newviscav,Yavnew,Yavprev,minresid,bestYav,H,critormax
    REAL(KIND=dp) :: newviscz(numoflayers),z(numoflayers),zlow(numoflayers), &
         zhigh(numoflayers)
    REAL(KIND=dp) :: minrlow,minrhigh,tlow,thigh,dmaxav
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    LOGICAL :: Visited = .FALSE., Found,dinit1

    SAVE :: Visited, SolverName, DmaxAv,hw, tol, warningtol, &
         numiters,CriticalDamage,CriticalDav,rhoi,rhow,critormax,Dmax


    IF (.NOT. Visited) THEN

       WRITE(SolverName, '(A)') 'InitParticleDz'

       rhow = Particles % rhow
       rhoi = Particles % rhoi

       IF (Particles % useisodam) THEN
          DmaxAv = Particles % initDmax
       ELSE
          DmaxAv  = Particles % DavDmaxI
       END IF

       Dmax = Particles % DmaxI

       hw = Particles % sealevel
       tol = Particles % dinittolerance
       warningtol = Particles % warningtol
       numiters = Particles % dinitnumiters

       IF (Particles % useisodam) THEN
          CriticalDamage = Particles % isodamcritdam
          CriticalDav = Particles % isodamcritdav
       ELSE

          CriticalDamage = Particles % criticaldamage
          CriticalDav = Particles % criticaldav
       END IF

       critormax = CriticalDamage  + (Dmax - CriticalDamage)/2

       Visited = .TRUE.
    END IF


    IF (Particles % nodamregion) THEN

       IF (Particles % Coordinate(No,1)>=Particles % ndxmin .AND. &
            Particles % Coordinate(No,1)<=Particles % ndxmax) THEN
          IF (Particles % Coordinate(No,2)>=Particles % ndymin .AND. &
               Particles % Coordinate(No,2)<=Particles % ndymax) THEN
             Particles % damage(No,:,:) = 0.0_dp
             Particles % dD(No,:,:) = 0.0_dp
             Particles % Dav(No,:) = 0.0_dp
             Dam = 0.0_dp
             Dav = 0.0_dp
             RETURN
          END IF
       END IF
    END IF

    IF (Particles % useriftdmax) THEN
       IF (Dav>Particles % riftdmax) THEN
          Dav = Particles % riftdmax
       END IF
    END IF


    IF (Particles % simpleinitdam) THEN
       Particles % Damage(No,:,1:3) = Dav
       Particles % Dav(No,4) = 0.0_dp
       Particles % Dav(No,1:3) = Dav
       Dam(:) = Dav
    END IF



    H = Particles % H(No)
    zs = H*(1.0_dp-rhoi/rhow)
    zb=zs-H


    intsize = H/(DBLE(numoflayers)-1.0_dp)

    DO ii = 0,(numoflayers-1)
       z(ii+1) = zb+ii*intsize
    END DO

    zlow = z-intsize/2.0_dp
    zlow(1) = z(1)
    zhigh = z+intsize/2.0_dp
    zhigh(numoflayers) = z(numoflayers)

    dav2 = 1.0_dp-Dav
    dmax2 = 1.0_dp-DmaxAv
    qqq= (dav2-1.0_dp)/(dmax2-dav2) + 1.0_dp
    Yav = H/qqq

    minrlow = 0.0_dp
    minrhigh = 0.0_dp
    tlow = 9999.0_dp
    thigh = 9999.0_dp

    last = 0

    DO kk = 1,numiters

       undamhigh = (hw+(1.0_dp-rhoi/rhow)*Yav)
       undamlow = ( hw- rhoi*Yav/rhow )
       !Dam is the damage with depth
       Dam(:) = 0.0_dp
       DO ii = 1,numoflayers
          IF (zhigh(ii) <= undamlow) THEN
             Dam(ii) = Dmax
          END IF
          IF (zlow(ii) >= ( undamhigh)) THEN
             Dam(ii) = Dmax
          END IF
          IF ((zlow(ii) <= undamlow) .AND. (zhigh(ii) >= undamlow)) THEN
             zrange = zhigh(ii)-zlow(ii)
             Dam(ii) = Dmax*((undamlow-zlow(ii))/zrange)
             ! IF (Dam(ii) < 0.0_dp) Dam(ii) = 0.0_dp
             ! IF (Dam(ii) > Dmax) Dam(ii) = Dmax
             IF (Dam(ii) > CriticalDamage) THEN
                IF (Dam(ii) > critormax) THEN
                   Dam(ii) = Dmax
                ELSE
                   Dam(ii) = CriticalDamage
                END IF
             END IF
          END IF

          IF ((zlow(ii) <= undamhigh) .AND. (zhigh(ii) >= undamhigh)) THEN
             zrange = zhigh(ii)-zlow(ii)
             Dam(ii) = Dmax*((zhigh(ii)-undamhigh)/zrange)
             ! IF (Dam(ii) < 0.0_dp) Dam(ii) = 0.0_dp
             ! IF (Dam(ii) > Dmax) Dam(ii) = Dmax
             IF (Dam(ii) > CriticalDamage)THEN
                IF (Dam(ii) > critormax) THEN
                   Dam(ii) = Dmax
                ELSE
                   Dam(ii) = CriticalDamage
                END IF
             END IF
          END IF
       END DO

       newviscz = (1.0_dp-Dam(:)) * Particles % bz(No,:) * Particles % EF(No)**(-1.0_dp/3.0_dp)
       newviscav = (((SUM(newviscz)-newviscz(1)/2.0_dp-newviscz(numoflayers)/2.0_dp)&
            *H/(DBLE(numoflayers)-1.0_dp)))/H


       IF (dinit1) THEN
          resid = newviscav/( (1.0_dp-Dav)*Particles % binit(No) * Particles % EF(No)**(-1.0_dp/3.0_dp)) -1.0_dp
       ELSE
          resid = newviscav/invvisc - 1.0_dp
       END IF


       IF (ABS(resid) > tol) THEN

          IF (resid>0.0_dp) THEN
             !your guess is too high

             IF (minrlow == 0.0_dp) THEN
                Yavnew = Yav - ABS(Yav)*0.5_dp
                !/2.0_dp
             ELSE
                Yavnew = (Yav + tlow)*0.5_dp
                !/2.0_dp
             END IF

             IF ((minrhigh == 0.0_dp) .OR. (resid<minrhigh)) THEN
                minrhigh = resid
                thigh = Yav
             END IF
          ELSE

             IF (minrhigh == 0.0_dp) THEN
                Yavnew = Yav + ABS(Yav)*0.5_dp
                !/2.0_dp
             ELSE
                Yavnew = (Yav+thigh)*0.5_dp
                !/2.0_dp
             END IF

             IF ((minrlow == 0.0_dp) .OR. (resid > minrlow)) THEN
                minrlow = resid
                tlow = Yav
             END IF
          END IF
       END IF

       Yavprev = Yav
       Yav = MIN(Yavnew,DmaxAv)
       Yav = MAX(Yavnew,0.0_dp)

       resid = ABS(resid)

       IF (kk==1) THEN
          minresid = resid
          bestYav = Yavprev
          bestiter = kk
       END IF

       IF (resid<minresid) THEN
          minresid = resid
          bestYav = Yavprev
          bestiter = kk
       END IF

       IF (resid<tol) THEN
          EXIT
       END IF
    END DO

    Yav = bestYav

    undamhigh = (hw+(1.0_dp-rhoi/rhow)*Yav)
    undamlow = ( hw- rhoi*Yav/rhow )
    !Dam is the damage with depth
    Dam(:) = 0.0_dp
    DO ii = 1,numoflayers
       IF (zhigh(ii) <= undamlow) THEN
          Dam(ii) = Dmax
       END IF
       IF (zlow(ii) >= ( undamhigh)) THEN
          Dam(ii) = Dmax
       END IF

       IF ((zlow(ii) <= undamlow) .AND. (zhigh(ii) >= undamlow)) THEN
          zrange = zhigh(ii)-zlow(ii)
          Dam(ii) = Dmax*((undamlow-zlow(ii))/zrange)
          ! IF (Dam(ii) < 0.0_dp) Dam(ii) = 0.0_dp
          ! IF (Dam(ii) > Dmax) Dam(ii) = Dmax
          IF (Dam(ii) > CriticalDamage) THEN
             IF (Dam(ii) > critormax) THEN
                Dam(ii) = Dmax
             ELSE
                Dam(ii) = CriticalDamage
             END IF
          END IF
       END IF

       IF ((zlow(ii) <= undamhigh) .AND. (zhigh(ii) >= undamhigh)) THEN
          zrange = zhigh(ii)-zlow(ii)
          Dam(ii) = Dmax*((zhigh(ii)-undamhigh)/zrange)
          ! IF (Dam(ii) < 0.0_dp) Dam(ii) = 0.0_dp
          ! IF (Dam(ii) > Dmax) Dam(ii) = Dmax
          IF (Dam(ii) > CriticalDamage) THEN
             IF (Dam(ii) > critormax) THEN
                Dam(ii) = Dmax
             ELSE
                Dam(ii) = CriticalDamage
             END IF
          END IF
       END IF
    END DO

    newviscz = (1.0_dp-Dam(:)) * Particles % bz(No,:) * Particles % EF(No)**(-1.0_dp/3.0_dp)
    newviscav = (((SUM(newviscz)-newviscz(1)/2.0_dp-newviscz(numoflayers)/2.0_dp)*H/(DBLE(numoflayers)-1.0_dp)))/H



    IF (dinit1) THEN
       resid = newviscav/( (1.0_dp-Dav)*Particles % binit(No) * Particles % EF(No)**(-1.0_dp/3.0_dp)) -1.0_dp
    ELSE
       resid = (newviscav/invvisc) - 1.0_dp
    END IF


    IF (resid>warningtol) THEN
       CALL Warn(SolverName,'Didnt even converge to warning tolerance!')
       PRINT *,'  '
       PRINT *,'Particle No: ', No
       PRINT *,'resid: ',resid
       PRINT *,'bestiter: ',bestiter
       PRINT *,'iter: ',kk
    END IF

    DO ii = 1,3
       Particles % Damage(No,:,ii) = Dam(:)
    END DO

    Dav = 1.0_dp - (newviscav/btzav)

    Particles % Dav(No,:) = 0.0_dp
    Particles % Dav(No,1:3) = Dav
    ! 1.0_dp - (newviscav/btzav)

  END SUBROUTINE InitParticleDz

  !**************************************************************************

  !> Depth-average the 3-D damage field for a particle
  !! Accounts for the influence of vertically-varying viscosity
  !! from other sources (e.g. temperature)
  SUBROUTINE VertIntDamFromVisc(Particles, No,layers,Model)

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No,layers,ii,dof
    REAL(KIND=dp), POINTER :: H, D(:),Dam(:)
    REAL(KIND=dp) :: newviscz(layers),btzav,newviscav,Dmax,critdav,damh,half=0.5_dp
    REAL(KIND=dp) :: NewDav,OldDav,invvisc, denom

    denom = 1.0_dp/(DBLE(layers)-1.0_dp)

    DO ii = 1,4
       D => Particles % Damage(No,1:layers,ii)

       newviscz(:) = Particles % Bz(No,:)
       btzav = (SUM(newviscz)-half*(newviscz(1)+newviscz(layers)) ) * denom

       newviscz(:) = D * Particles % Bz(No,:)
       newviscav = (SUM(newviscz)-half*(newviscz(1)+newviscz(layers)) ) * denom

       Particles % Dav(No,ii) = newviscav/btzav
    END DO

  END SUBROUTINE VertIntDamFromVisc

  !**************************************************************************

  !> Subroutine allocates particles before launching them.
  SUBROUTINE AllocateParticles(Particles, Model,NoParticles)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles

    INTEGER :: NoParticles, layers, templayers
    TYPE(Model_t) :: Model
    INTEGER :: PrevNoParticles, dofs, No, n, dim, TimeOrder, n1, n2,bufferparticles
    INTEGER, ALLOCATABLE :: Perm(:)
    REAL (KIND=dp) :: anisoparam

    LOGICAL, POINTER :: Static(:),UseInterpElem(:) !,isodam(:)
    REAL(KIND=dp), POINTER :: Coordinate(:,:),NextCoordinate(:,:)
    INTEGER, POINTER :: ElementIndex(:),InterpElem(:),Status(:) ,test(:),&
         test2(:),OrigNo(:),DamStatus(:)
    REAL(KIND=dp), POINTER :: dbassis(:),damstrain(:)
    REAL(KIND=dp), POINTER :: GradVel(:,:),Velocity(:,:),GridVelocity(:,:),GradZs(:,:),GradH(:,:)
    REAL(KIND=dp), POINTER :: Damage(:,:,:),dD(:,:,:),F(:,:),Dav(:,:),Length(:,:)
    REAL(KIND=dp), POINTER :: OrigLength(:,:),GMask(:),GVolume(:),Strain(:,:)
    REAL(KIND=dp), POINTER :: FAlpha(:,:),MB(:),xpic(:,:),PVolume(:),Bedrock(:),Ee(:)
    REAL(KIND=dp), POINTER :: FVolume(:),Tracer(:)
    REAL(KIND=dp), POINTER :: FP(:),Binit(:),H(:),EF(:),Bz(:,:),Mass(:)
    LOGICAL :: Visited = .FALSE.,GotIt


    SAVE :: Visited,layers,templayers,bufferparticles

    IF (.NOT. Visited ) THEN

       layers = Particles % numberofparticlelayers
       templayers = Particles % numberoftemperaturelayers

       bufferparticles = Particles % buffer

       Visited = .TRUE.
    END IF


    IF( NoParticles <= Particles % MaxNumberOfParticles ) THEN
       PRINT *,'AllocateParticles','There are already enough particles'
       RETURN
    ELSE
       WRITE(Message,'(A,I0)') 'Allocating number of particles: ',NoParticles
       CALL Info('AllocateParticles',Message)
       WRITE(Message,'(A,I0)') &
            'Also allocating an addition number of buffer particles: ',bufferparticles
       CALL Info('AllocateParticles',Message)
       NoParticles = NoParticles + bufferparticles
    END IF

    dim = Particles % dim
    dofs = dim


    PrevNoParticles = Particles % NumberOfParticles


    ! If appending particles to existing particles, we essentially have to make a copy
    ! of the old allocation, reallocate the old allocation, and then fill the reallocated
    ! memory with the copied allocation for the existing particles. This effectively
    ! means we are at least doubling the memory use temporarily, which can cause a crash.
    ! We deal with this two ways:
    ! 1. Add "buffer particles" that increase initial memory, but reduce the need
    !    to reallocate as often (reallocations only needed when allocating more particles
    !    than NoParticles + bufferparticles
    !    are exceeded.
    ! 2. Deallocate non-history variables to reduce the size of the new allocation.
    !    We retreive new values for these particles afterwards from the mesh


    ! Set pointers to the old stuff, these are needed
    ! if growing an already existing list of particles.

    IF( PrevNoParticles > 0 ) THEN

       n1 = 1
       n2 = PrevNoParticles

       Bedrock => Particles % Bedrock
       Binit => Particles % Binit
       BZ => Particles % Bz
       Coordinate => Particles % Coordinate
       Damage => Particles % Damage
       Dav => Particles % Dav
       dD => Particles % dD
       EF => Particles % EF
       F => Particles % F
       FP => Particles % FP
       GMask => Particles % GMask
       GradH => Particles % GradH
       GradVel => Particles % GradVel
       GradZs => Particles % GradZs
       GridVelocity => Particles % GridVelocity
       GVolume => Particles % GVolume
       H => Particles % H
       Length => Particles % Length
       Mass => Particles % Mass
       MB => Particles % MB
       NextCoordinate => Particles % NextCoordinate
       OrigLength => Particles % OrigLength
       PVolume => Particles % PVolume

       IF (Particles % trackstrain) THEN
          Strain => Particles % Strain
       END IF

       IF (Particles % usetracer) THEN
          Tracer => Particles % Tracer
       END IF


       Velocity => Particles % Velocity
       xpic => Particles % xpic
       OrigNo => Particles % OrigNo

       IF (Particles % outputdbassis) THEN
          dbassis => Particles % dbassis
       END IF


       ElementIndex => Particles % ElementIndex
       InterpElem => Particles % InterpElem
       Status => Particles % Status

       Static => Particles % Static
       DamStatus => Particles % DamStatus
       UseInterpElem => Particles % UseInterpElem
    END IF

    ALLOCATE( Particles % Bedrock(NoParticles) )
    ALLOCATE( Particles % Binit(NoParticles) )
    ALLOCATE( Particles % Bz(NoParticles,layers) )
    ALLOCATE( Particles % Coordinate(NoParticles,dofs) )
    ALLOCATE( Particles % Damage(NoParticles,layers,4) )
    ALLOCATE( Particles % Dav(NoParticles,4) )
    ALLOCATE( Particles % dD(NoParticles,layers,4) )
    ALLOCATE( Particles % EF(NoParticles) )
    ALLOCATE( Particles % F(NoParticles,4) )
    ALLOCATE( Particles % FP(NoParticles) )
    ALLOCATE( Particles % GMask(NoParticles) )
    ALLOCATE( Particles % GradH(NoParticles,dofs) )
    ALLOCATE( Particles % GradVel(NoParticles,4) )
    ALLOCATE( Particles % GradZs(NoParticles,dofs) )
    ALLOCATE( Particles % GridVelocity(NoParticles,dofs) )
    ALLOCATE( Particles % GVolume(NoParticles) )
    ALLOCATE( Particles % H(NoParticles) )
    ALLOCATE( Particles % Length(NoParticles,2) )
    ALLOCATE( Particles % Mass(NoParticles) )
    ALLOCATE( Particles % MB(NoParticles) )
    ALLOCATE( Particles % NextCoordinate(NoParticles,dofs) )
    ALLOCATE( Particles % OrigLength(NoParticles,2) )


    IF (Particles % outputdbassis) THEN
       ALLOCATE( Particles % dbassis(NoParticles) )
    END IF



    ALLOCATE( Particles % PVolume(NoParticles) )
    IF (Particles % trackstrain) THEN
       ALLOCATE( Particles % Strain(NoParticles,2) )
    END IF

    IF (Particles % usetracer) THEN
       ALLOCATE (Particles % Tracer(NoParticles))
    END IF


    ALLOCATE( Particles % Velocity(NoParticles,dofs) )
    ALLOCATE( Particles % xpic(NoParticles,6) )

    ALLOCATE( Particles % ElementIndex(NoParticles) )
    ALLOCATE( Particles % InterpElem(NoParticles) )
    ALLOCATE( Particles % Status(NoParticles) )

    ALLOCATE( Particles % Static(NoParticles) )
    ALLOCATE( Particles % DamStatus(NoParticles))
    ALLOCATE( Particles % UseInterpElem(NoParticles) )

    ALLOCATE( Particles % OrigNo(NoParticles) )


    !---------------------------------------------!
    IF( PrevNoParticles > 0 ) THEN
       n = 0
       ALLOCATE( Perm( PrevNoParticles ) )
       Perm = 0

       DO No=1,PrevNoParticles
          IF ( Status(No) == PARTICLE_LOST ) CYCLE
          IF ( Status(No) == PARTICLE_ALLOCATED ) CYCLE
          n = n+1
          Perm(n) = No
       END DO

       WRITE(Message,'(A,I0)') 'Number of old active particles: ',n
       CALL Info('AllocateParticles',Message,Level=4)

       IF( n < PrevNoParticles ) THEN
          WRITE(Message,'(A,I0)') 'Number of deleted particles: ',PrevNoParticles-n
          CALL Info('AllocateParticles',Message,Level=4)
       END IF

       n1 = 1
       n2 = n


       Particles % Bedrock(n1:n2) = Bedrock(Perm(n1:n2))
       Particles % Binit(n1:n2) = Binit(Perm(n1:n2))
       Particles % Bz(n1:n2,:) = BZ(Perm(n1:n2),:)
       Particles % Coordinate(n1:n2,:) = Coordinate(Perm(n1:n2),:)
       Particles % Damage(n1:n2,:,:) = Damage(Perm(n1:n2),:,:)
       Particles % Dav(n1:n2,:) = Dav(Perm(n1:n2),:)
       Particles % dD(n1:n2,:,:) = dd(Perm(n1:n2),:,:)
       Particles % EF(n1:n2) = EF(Perm(n1:n2))
       Particles % F(n1:n2,:) = F(Perm(n1:n2),:)
       Particles % FP(n1:n2) = FP(Perm(n1:n2))
       Particles % GMask(n1:n2) = GMask(Perm(n1:n2))
       Particles % GradH(n1:n2,:) = GradH(Perm(n1:n2),:)
       Particles % GradVel(n1:n2,:) = GradVel(Perm(n1:n2),:)
       Particles % GradZs(n1:n2,:) = GradZs(Perm(n1:n2),:)
       Particles % GridVelocity(n1:n2,:) = GridVelocity(Perm(n1:n2),:)
       Particles % GVolume(n1:n2) = GVolume(Perm(n1:n2))
       Particles % H(n1:n2) = H(Perm(n1:n2))
       Particles % Length(n1:n2,:) = Length(Perm(n1:n2),:)
       Particles % Mass(n1:n2) = Mass(Perm(n1:n2))
       Particles % MB(n1:n2) = MB(Perm(n1:n2))
       Particles % NextCoordinate(n1:n2,:) = NextCoordinate(Perm(n1:n2),: )
       Particles % OrigLength(n1:n2,:) = OrigLength(Perm(n1:n2),:)


       IF (Particles % outputdbassis) THEN
          Particles % dbassis(n1:n2) = dbassis(Perm(n1:n2))
       END IF


       Particles % PVolume(n1:n2) = PVolume(Perm(n1:n2))

       IF (Particles % trackstrain) THEN
          Particles % Strain(n1:n2,:) = Strain(Perm(n1:n2),:)
       END IF

       IF (Particles % usetracer) THEN
          Particles % Tracer(n1:n2) = Tracer(Perm(n1:n2))
       END IF


       Particles % Velocity(n1:n2,:) = Velocity(Perm(n1:n2),:)
       Particles % xpic(n1:n2,:) = xpic(Perm(n1:n2),:)
       Particles % ElementIndex(n1:n2) = ElementIndex(Perm(n1:n2))
       Particles % InterpElem(n1:n2) = InterpElem(Perm(n1:n2))
       Particles % Status(n1:n2) = Status(Perm(n1:n2))
       Particles % Static(n1:n2) = Static(Perm(n1:n2))
       Particles % DamStatus(n1:n2) = DamStatus(Perm(n1:n2))
       Particles % UseInterpElem(n1:n2) = UseInterpElem(Perm(n1:n2))
       Particles % OrigNo(n1:n2) = OrigNo(Perm(n1:n2))

       PrevNoParticles = n
       Particles % NumberOfParticles = n


       IF (ASSOCIATED(OrigNo)) DEALLOCATE(OrigNo)
       IF (ASSOCIATED(UseInterpElem)) DEALLOCATE(UseInterpElem)
       IF (ASSOCIATED(Static)) DEALLOCATE(Static)
       IF (ASSOCIATED(DamStatus)) DEALLOCATE(DamStatus)
       IF (ASSOCIATED(Status)) DEALLOCATE(Status)
       IF (ASSOCIATED(InterpElem )) DEALLOCATE(InterpElem )
       IF (ASSOCIATED(ElementIndex )) DEALLOCATE(ElementIndex )
       IF (ASSOCIATED(xpic )) DEALLOCATE(xpic )
       IF (ASSOCIATED(Velocity )) DEALLOCATE(Velocity )

       IF (Particles % usetracer) THEN
          IF (ASSOCIATED(Tracer)) DEALLOCATE(Tracer)
       END IF


       IF (Particles % trackstrain) THEN
          IF (ASSOCIATED(Strain )) DEALLOCATE(Strain )
       END IF

       IF (ASSOCIATED(PVolume )) DEALLOCATE(PVolume )

       IF (Particles % outputdbassis) THEN
          IF (ASSOCIATED(dbassis ) ) DEALLOCATE(dbassis)
       END IF

       IF (ASSOCIATED(OrigLength )) DEALLOCATE(OrigLength )
       IF (ASSOCIATED(NextCoordinate )) DEALLOCATE(NextCoordinate )
       IF (ASSOCIATED(MB )) DEALLOCATE(MB )
       IF (ASSOCIATED(Mass )) DEALLOCATE(Mass )
       IF (ASSOCIATED(Length )) DEALLOCATE(Length )
       IF (ASSOCIATED(H )) DEALLOCATE(H )
       IF (ASSOCIATED(GVolume )) DEALLOCATE(GVolume )
       IF (ASSOCIATED(GridVelocity )) DEALLOCATE(GridVelocity )
       IF (ASSOCIATED(GradZs )) DEALLOCATE(GradZs )
       IF (ASSOCIATED(GradVel )) DEALLOCATE(GradVel )
       IF (ASSOCIATED(GradH )) DEALLOCATE(GradH )
       IF (ASSOCIATED(GMask )) DEALLOCATE(GMask )
       IF (ASSOCIATED(FP )) DEALLOCATE(FP )
       IF (ASSOCIATED(F )) DEALLOCATE(F )
       IF (ASSOCIATED(EF )) DEALLOCATE(EF )
       IF (ASSOCIATED(dD )) DEALLOCATE(dD )
       IF (ASSOCIATED(Dav )) DEALLOCATE(Dav )
       IF (ASSOCIATED(Damage )) DEALLOCATE(Damage )
       IF (ASSOCIATED(Coordinate )) DEALLOCATE(Coordinate )
       IF (ASSOCIATED(Bz )) DEALLOCATE(Bz )
       IF (ASSOCIATED(Binit )) DEALLOCATE(Binit )
       IF (ASSOCIATED(Bedrock )) DEALLOCATE(Bedrock )
    END IF

    ! Initialize the newly allocated particles with default values
    !-------------------------------------------------------------
    n1 = PrevNoParticles+1
    n2 = NoParticles

    Particles % Bedrock(n1:n2) = 0.0_dp
    Particles % Binit(n1:n2) = 0.0_dp
    Particles % Bz(n1:n2,:) = 0.0_dp
    Particles % Coordinate(n1:n2,:) = 0.0_dp
    Particles % Damage(n1:n2,:,:) = 0.0_dp
    Particles % Dav(n1:n2,:) = 0.0_dp
    Particles % dD(n1:n2,:,:) = 0.0_dp
    Particles % EF(n1:n2) = 0.0_dp
    Particles % ElementIndex(n1:n2) = 0
    Particles % F(n1:n2,:) = 0.0_dp
    Particles % FP(n1:n2) = 0.0_dp
    Particles % GMask(n1:n2) = 0.0_dp
    Particles % GradH(n1:n2,:) = 0.0_dp
    Particles % GradVel(n1:n2,:) = 0.0_dp
    Particles % GradZs(n1:n2,:) = 0.0_dp
    Particles % GridVelocity(n1:n2,:) = 0.0_dp
    Particles % GVolume(n1:n2) = 0.0_dp
    Particles % H(n1:n2) = 0.0_dp
    Particles % InterpElem(n1:n2) = 0
    Particles % Length(n1:n2,:) = 0.0_dp
    Particles % Mass(n1:n2) = 0.0_dp
    Particles % MB(n1:n2) = 0.0_dp
    Particles % NextCoordinate(n1:n2,:) = 0.0_dp
    Particles % OrigLength(n1:n2,:) = 0.0_dp


    IF (Particles % outputdbassis) THEN
       Particles % dbassis(n1:n2) = 0.0_dp
    END IF


    Particles % PVolume(n1:n2) = 0.0_dp
    Particles % Static(n1:n2) = .FALSE.
    Particles % DamStatus(n1:n2) = 0
    Particles % Status(n1:n2) = PARTICLE_ALLOCATED

    IF (Particles % trackstrain) THEN
       Particles % Strain(n1:n2,:) = 0.0_dp
    END IF

    IF (Particles % usetracer) THEN
       Particles % Tracer(n1:n2)  = 0.0_dp
    END IF


    Particles % Velocity(n1:n2,:) = 0.0_dp
    Particles % xpic(n1:n2,:) = 0.0_dp
    Particles % OrigNo(n1:n2) = 0
    Particles % UseInterpElem(n1:n2) = .FALSE.

    !----------------------------------------------------------------!

    Particles % MaxNumberOfParticles = NoParticles

    PRINT *,'MaxNoParticles',Particles % MaxNumberOfParticles

    IF( PrevNoParticles > 0 ) THEN
       CALL Info('AllocateParticles','Deallocating particle permutation',Level=13)
       DEALLOCATE( Perm )
    END IF


  END SUBROUTINE AllocateParticles

  !**************************************************************************

  !> Subroutine deletes lost particles that have exited the computational domain
  !! TODO: once the MPM code is parallelized, this routine should also be used for
  !! particles that go to a neighboring partition, as in ParticleUtils.F90
  SUBROUTINE DeleteLostParticles(Particles)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No, n, PrevNoParticles, n1, n2
    INTEGER, ALLOCATABLE :: Perm(:)
    LOGICAL :: highestnoparticleslost

    PrevNoParticles = Particles % NumberOfParticles
    IF( PrevNoParticles == 0 ) RETURN

    ALLOCATE( Perm( PrevNoParticles ) )
    Perm = 0

    n = 0
    n1 = 0
    highestnoparticleslost = .FALSE.

    DO No=1,PrevNoParticles
       IF (Particles % Status(No) == PARTICLE_LOST .OR. &
            Particles % Status(No) == PARTICLE_ALLOCATED) THEN
          CYCLE
       END IF

       n = n + 1
       IF( n1 == 0 .AND. n /= No ) n1 = n
       Perm(n) = No
    END DO

    IF (n1 == 0 .AND. n/=PrevNoParticles) highestnoparticleslost = .TRUE.

    n2 = n
    CALL Info('DeleteLostParticles','Number of active particles: '&
         //TRIM(I2S(n2)),Level=1)

    IF(n1 == 0 ) THEN
       IF (.NOT. highestnoparticleslost) THEN
          CALL Info('DeleteLostParticles','No particles need to be deleted',Level=1)
          RETURN
       ELSE
          CALL Info('DeleteLostParticles','Some particles with highest No lost',Level=1)
       END IF
    ELSE
       CALL Info('DeleteLostParticles','First particle with changed permutation: '&
            //TRIM(I2S(n1)),Level=1)
    END IF


    Particles % Bedrock(n1:n2) = Particles % Bedrock(Perm(n1:n2))
    Particles % Binit(n1:n2) = Particles % Binit(Perm(n1:n2))
    Particles % Bz(n1:n2,:) = Particles % Bz(Perm(n1:n2),:)
    Particles % Coordinate(n1:n2,:) = Particles % Coordinate(Perm(n1:n2),:)
    Particles % Damage(n1:n2,:,:) = Particles % Damage(Perm(n1:n2),:,:)
    Particles % Dav(n1:n2,:) = Particles % Dav(Perm(n1:n2),:)
    Particles % dD(n1:n2,:,:) = Particles % dD(Perm(n1:n2),:,:)
    Particles % EF(n1:n2) = Particles % EF(Perm(n1:n2))
    Particles % ElementIndex(n1:n2) = Particles % ElementIndex(Perm(n1:n2))
    Particles % F(n1:n2,:) = Particles % F(Perm(n1:n2),:)
    Particles % FP(n1:n2) = Particles % FP(Perm(n1:n2))
    Particles % GMask(n1:n2) = Particles % GMask(Perm(n1:n2))
    Particles % GradH(n1:n2,:) = Particles % GradH(Perm(n1:n2),:)
    Particles % GradVel(n1:n2,:) = Particles % GradVel(Perm(n1:n2),:)
    Particles % GradZs(n1:n2,:) = Particles % GradZs(Perm(n1:n2),:)
    Particles % GridVelocity(n1:n2,:) = Particles % GridVelocity(Perm(n1:n2),:)
    Particles % GVolume(n1:n2) = Particles % GVolume(Perm(n1:n2))
    Particles % H(n1:n2) = Particles % H(Perm(n1:n2))
    Particles % InterpElem(n1:n2) = Particles % InterpElem(Perm(n1:n2))
    Particles % Length(n1:n2,:) = Particles % Length(Perm(n1:n2),:)
    Particles % Mass(n1:n2) = Particles % Mass(Perm(n1:n2))
    Particles % MB(n1:n2) = Particles % MB(Perm(n1:n2))
    Particles % NextCoordinate(n1:n2,:) = Particles % NextCoordinate(Perm(n1:n2),:)
    Particles % OrigLength(n1:n2,:) = Particles % OrigLength(Perm(n1:n2),:)

    IF (Particles % outputdbassis) THEN
       Particles % dbassis(n1:n2) = Particles % dbassis(Perm(n1:n2))
    END IF


    Particles % PVolume(n1:n2) = Particles % PVolume(Perm(n1:n2))
    Particles % DamStatus(n1:n2) = Particles % DamStatus(Perm(n1:n2))
    Particles % Static(n1:n2) = Particles % Static(Perm(n1:n2))
    Particles % Status(n1:n2) = Particles % Status(Perm(n1:n2))

    IF (Particles % trackstrain) THEN
       Particles % Strain(n1:n2,:) = Particles % Strain(Perm(n1:n2),:)
    END IF

    IF (Particles % usetracer) THEN
       Particles % tracer(n1:n2) = Particles % Tracer(Perm(n1:n2))
    END IF


    Particles % Velocity(n1:n2,:) = Particles % Velocity(Perm(n1:n2),:)
    Particles % xpic(n1:n2,:) = Particles % xpic(Perm(n1:n2),:)

    Particles % UseInterpElem(n1:n2) = Particles % UseInterpElem(Perm(n1:n2))

    Particles % OrigNo(n1:n2) = Particles % OrigNo(Perm(n1:n2))
    !----------------------------------------------------------------!

    Particles % NumberOfParticles = n2

    IF ( n2 < PrevNoParticles ) THEN

       Particles % Bedrock(n2+1:PrevNoParticles) = 0.0_dp
       Particles % Binit(n2+1:PrevNoParticles) = 0.0_dp
       Particles % Bz(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % Coordinate(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % Damage(n2+1:PrevNoParticles,:,:) = 0.0_dp
       Particles % Dav(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % dD(n2+1:PrevNoParticles,:,:) = 0.0_dp
       Particles % EF(n2+1:PrevNoParticles) = 0.0_dp
       Particles % ElementIndex(n2+1:PrevNoParticles) = 0
       Particles % F(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % FP(n2+1:PrevNoParticles) = 0.0_dp
       Particles % GMask(n2+1:PrevNoParticles) = 0.0_dp
       Particles % GradH(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % GradVel(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % GradZs(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % GridVelocity(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % GVolume(n2+1:PrevNoParticles) = 0.0_dp
       Particles % H(n2+1:PrevNoParticles) = 0.0_dp
       Particles % InterpElem(n2+1:PrevNoParticles) = 0
       Particles % Length(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % Mass(n2+1:PrevNoParticles) = 0.0_dp
       Particles % MB(n2+1:PrevNoParticles) = 0.0_dp
       Particles % NextCoordinate(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % OrigLength(n2+1:PrevNoParticles,:) = 0.0_dp

       IF (Particles % outputdbassis) THEN
          Particles % dbassis(n2+1:PrevNoParticles) = 0.0_dp
       END IF


       Particles % PVolume(n2+1:PrevNoParticles) = 0.0_dp
       Particles % damstatus(n2+1:PrevNoParticles) = 0
       Particles % Static(n2+1:PrevNoParticles) = .FALSE.
       Particles % Status(n2+1:PrevNoParticles) = PARTICLE_ALLOCATED

       IF (Particles % trackstrain) THEN
          Particles % Strain(n2+1:PrevNoParticles,:) = 0.0_dp
       END IF

       IF (Particles % usetracer) THEN
          Particles % Tracer(n2+1:PrevNoParticles) = 0.0_dp
       END IF

       Particles % Velocity(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % xpic(n2+1:PrevNoParticles,:) = 0.0_dp

       Particles % UseInterpElem(n2+1:PrevNoParticles) = .FALSE.

       Particles % OrigNo(n2+1:PrevNoParticles) = 0

    END IF

    DEALLOCATE( Perm )

  END SUBROUTINE DeleteLostParticles

  !**************************************************************************

  !> Particle to grid interpolations
  SUBROUTINE MPMParticlesToNodes( Particles, Model, whichtime)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t)   :: ElementNodes
    TYPE(Variable_t), POINTER :: HVar,V1Var,V2Var,PassVar,Var, WeightVar,MassVar,&
         invvar,BCVar, IV1Var,IV2Var,HweightVar,BVar,PM,xp1var,xp2var,surfvar,&
         DxxVar,DyyVar,DzzVar,DxyVar,maskvar,mfvar,mfwvar,dirvar,dir2var,opvar
    TYPE(Element_t), POINTER :: BulkElement
    INTEGER :: nn, nb, NoVar, ii, No, ni,t, whichtime
    CHARACTER(LEN=MAX_NAME_LEN) :: TargetVariableName
    INTEGER, POINTER :: NodeIndexes(:),HPerm(:),V1Perm(:),V2Perm(:),PassPerm(:),&
         Perm(:),WeightPerm(:),MassPerm(:),invPerm(:),BCPerm(:),&
         IV1Perm(:),IV2Perm(:),HWeightPerm(:),BPerm(:),PMPerm(:),xp1perm(:),&
         xp2perm(:),surfperm(:),DxxPerm(:),DyyPerm(:),DzzPerm(:),DxyPerm(:),&
         maskperm(:),mfperm(:),mfwperm(:),dirperm(:),dir2perm(:),opperm(:)
    LOGICAL :: Stat, Visited = .FALSE.,Found
    REAL(KIND=dp) :: detJ, rhow, sealevel, rhoi, g, h_im,scale,gridres,zs,Hf,Area,cm
    REAL(KIND=dp), POINTER :: HValues(:),V1Values(:),V2Values(:),PassValues(:),&
         Values(:),WeightValues(:),MassValues(:),invValues(:),BCValues(:),&
         IV1Values(:),IV2Values(:),HWeightValues(:),BValues(:),PMValues(:),&
         xp1values(:),xp2values(:),surfvalues(:),DxxValues(:),DyyValues(:),DzzValues(:),&
         DxyValues(:),mask(:),mf(:),mfw(:),dir(:),dir2(:),op(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    REAL(KIND=dp),ALLOCATABLE :: Basis(:),dBasisdx(:,:)
    REAL(KIND=dp) :: ddBasisddx(4,3,3)
    REAL(KIND=dp) :: xp1,xp2,thick,mass,tmf
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: yy2,yy1,x1,x2,mm,bb,rmf

    SAVE :: Mesh, nn, nb, rhow, sealevel, rhoi, g, &
         Visited,gridres,SolverName,&
         Basis,dBasisdx,cm,ElementNodes


    IF( .NOT. Visited ) THEN

       WRITE(SolverName, '(A)') 'MPMParticlesToNodes'
       Mesh => GetMesh()
       nn = Mesh % MaxElementNodes
       nb = Mesh % NumberOfBulkElements

       ALLOCATE(Basis(nn),dBasisdx(nn,3))
       ALLOCATE(ElementNodes % x(nn),ElementNodes % y(nn),ElementNodes % z(nn))

       rhow = Particles % rhow
       sealevel = Particles % sealevel
       rhoi = Particles % rhoi
       gridres = Particles % gridres

       g = ABS(Particles % gravity)


       cm = GetConstReal( Model % Constants, 'Viscosity Exponent', Found )
       IF (.NOT. Found) CALL Fatal(SolverName,&
            'Need to define "Viscosity Exponent = Real $1/n" in constants')

       Visited = .TRUE.
    END IF

    CALL INFO(Trim(SolverName), &
         '-----Interpolating Particles to Field ----',Level=5)


    WeightVar => VariableGet(Model % Mesh % Variables, 'TempVar' )
    WeightValues => WeightVar % Values
    WeightPerm => WeightVar % Perm
    WeightValues = 0.0_dp

    IV1Var => VariableGet(Model % Mesh % Variables, 'invvel 1' )
    IV1Values => IV1Var % Values
    IV1Perm => IV1Var % Perm

    IV2Var => VariableGet(Model % Mesh % Variables, 'invvel 2' )
    IV2Values => IV2Var % Values
    IV2Perm => IV2Var % Perm

    IF (whichtime == 1) THEN

       Var => VariableGet(Model % Mesh % Variables, 'btz' )
       Values => Var % Values
       Perm => Var % Perm
       values = 0.0_dp

       invVar => VariableGet(Model % Mesh % Variables, 'invvisc' )
       invValues => invVar % Values
       invPerm => invVar % Perm

    ELSEIF (whichtime == 2) THEN
       HVar => VariableGet(Model % Mesh % Variables, 'H' )
       HValues => HVar % Values
       HPerm => HVar % Perm
       HValues = 0.0_dp


       HWeightVar => VariableGet(Model % Mesh % Variables, 'Weight' )
       HWeightValues => HWeightVar % Values
       HWeightPerm => HWeightVar % Perm
       HWeightValues = 0.0_dp


       MassVar => VariableGet(Model % Mesh % Variables, 'Mass' )
       MassValues => MassVar % Values
       MassPerm => MassVar % Perm
       MassValues = 0.0_dp

       V1Var => VariableGet(Model % Mesh % Variables, 'SSAVelocity 1' )
       V1Values => V1Var % Values
       V1Perm => V1Var % Perm



       V2Var => VariableGet(Model % Mesh % Variables, 'SSAVelocity 2' )
       V2Values => V2Var % Values
       V2Perm => V2Var % Perm

       IF (.NOT. Particles % uplag) THEN
          V1Values = 0.0_dp
          V2Values = 0.0_dp
       END IF

       BVar => VariableGet(Model % Mesh % Variables, 'Particle B' )
       BValues => BVar % Values
       BPerm => BVar % Perm
       BValues = 0.0_dp


       maskvar => VariableGet(Model % Mesh % Variables, 'Mask' )
       maskPerm => maskvar % Perm
       mask => maskvar % Values

       IF (Particles % usedamage) THEN
          DxxVar => VariableGet(Model % Mesh % Variables, 'Mesh Damage 1' )
          DxxValues => DxxVar % Values
          DxxPerm => DxxVar % Perm
          DxxValues = 0.0_dp

          DyyVar => VariableGet(Model % Mesh % Variables, 'Mesh Damage 2' )
          DyyValues => DyyVar % Values
          DyyPerm => DyyVar % Perm
          DyyValues = 0.0_dp

          DzzVar => VariableGet(Model % Mesh % Variables, 'Mesh Damage 3' )
          DzzValues => DzzVar % Values
          DzzPerm => DzzVar % Perm
          DzzValues = 0.0_dp

          DxyVar => VariableGet(Model % Mesh % Variables, 'Mesh Damage 4' )
          DxyValues => DxyVar % Values
          DxyPerm => DxyVar % Perm
          DxyValues = 0.0_dp
       END IF

       PM => VariableGet( Mesh % Variables, 'surface')
       PMValues => PM % Values
       PMPerm => PM % Perm


    ELSEIF (whichtime == 3 .OR. whichtime == 8) THEN
       MassVar => VariableGet(Model % Mesh % Variables, 'Mass' )
       MassValues => MassVar % Values
       MassPerm => MassVar % Perm
       MassValues = 0.0_dp


       V1Var => VariableGet(Model % Mesh % Variables, 'SSAVelocity 1' )
       V1Values => V1Var % Values
       V1Perm => V1Var % Perm
       V1Values = 0.0_dp

       V2Var => VariableGet(Model % Mesh % Variables, 'SSAVelocity 2' )
       V2Values => V2Var % Values
       V2Perm => V2Var % Perm
       V2Values = 0.0_dp

    ELSEIF (whichtime == 4) THEN

       !Mass should already be on mesh
       MassVar => VariableGet(Model % Mesh % Variables, 'Mass' )
       MassValues => MassVar % Values
       MassPerm => MassVar % Perm

       V1Var => VariableGet(Model % Mesh % Variables, 'Vplus 1' )
       V1Values => V1Var % Values
       V1Perm => V1Var % Perm

       V2Var => VariableGet(Model % Mesh % Variables, 'Vplus 2' )
       V2Values => V2Var % Values
       V2Perm => V2Var % Perm

       V1Values = 0.0_dp
       V2Values = 0.0_dp

    ELSEIF (whichtime == 5) THEN

       HVar => VariableGet(Model % Mesh % Variables, 'Zs' )
       HValues => HVar % Values
       HPerm => HVar % Perm
       HValues = 0.0_dp

    ELSEIF (whichtime == 6) THEN

       Var => VariableGet(Model % Mesh % Variables, 'Damage' )
       Values => Var % Values
       Perm => Var % Perm
       values = 0.0_dp

    ELSEIF (whichtime == 7) THEN

       HVar => VariableGet(Model % Mesh % Variables, 'H' )
       HValues => HVar % Values
       HPerm => HVar % Perm
       HValues = 0.0_dp

       HWeightVar => VariableGet(Model % Mesh % Variables, 'Weight' )
       HWeightValues => HWeightVar % Values
       HWeightPerm => HWeightVar % Perm
       HWeightValues = 0.0_dp

    END IF



    DO ii = 1,nb

       IF ( ElemTrack(ii) % Status >= NOTFULL ) THEN

          BulkElement => Mesh % Elements( ii )
          IF(.NOT. ASSOCIATED( BulkElement ) ) CYCLE

          NodeIndexes => BulkElement % NodeIndexes

          CALL GetElementNodes(ElementNodes,BulkElement)

          ! NodeIndexes => BulkElement % NodeIndexes
          nn = BulkElement % TYPE % NumberofNodes

          Area = ElementArea(Model % Mesh,BulkElement,nn)


          DO t = 1,ABS(ElemParticles(ii) % NumberOfParticles)
             No = ElemParticles(ii) % p(t)

             IF (No>Particles % NumberOfParticles .OR. No<1 .OR. (No .NE. No)) CYCLE


             IF (Particles % Status(No) >= PARTICLE_LOST) CYCLE



             IF (Particles % ShapeFunctions == 'gimpm') THEN
                stat = GIMPMElementInfo( t,Particles, Model,BulkElement, ElementNodes, No, &
                     detJ, scale, .FALSE., Basis,dBasisdx)
             ELSE
                stat = sMPMElementInfo( BulkElement,Particles, Model, ElementNodes, No, &
                     Particles % gridres, Basis,dBasisdx)
                scale = 1.0_dp
                detJ = Particles % PVolume(No)
             END IF


             ! NOTE: using true for gimpmelementinfo takes particles with
             ! status PARTICLE_LEAVING (i.e. overlapping a bound, or passive element)
             ! basis functions for that particle as if it was moved and scaled to fit
             ! entirely within the volume it overlaps in the element being called.
             ! Multiplying by "scale" is not applicable here. "Scale" is for use with
             ! MeshToParticle interpolation (which also use the PARTICLE_LEAVING interpolation
             ! presented here) where scale is used to guarantee the volumes of the temporary particles
             ! add up to the full particle volume properly.

             scale = 1.0_dp

             IF (whichtime == 1) THEN
                Values(Perm(NodeIndexes)) = Values(Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % Binit(No) *detJ/(Area)

                IF (t== ABS(ElemParticles(ii) % NumberOfParticles)) THEN
                   WeightValues(WeightPerm(NodeIndexes)) = WeightValues(WeightPerm(NodeIndexes)) + 1.0_dp
                END IF

             ELSEIF (whichtime == 2) THEN

                MassValues(MassPerm(NodeIndexes)) = MassValues(MassPerm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % Mass(No) !*detJ/(Area)


                HWeightValues(HWeightPerm(NodeIndexes)) = HWeightValues(HWeightPerm(NodeIndexes)) + &
                     Basis(1:nn) * detJ !/ (Area)

                HValues(HPerm(NodeIndexes)) = HValues(HPerm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % H(No) * detJ !/ (Area)


                !we need viscosity on the mesh for FEM elements, which have PMVal=2.0
                !for all nodes. Only Interp for those for efficiency.

                IF (Particles % uplag .OR. (ANY(PMValues(PMPerm(NodeIndexes)) > 1.0_dp) &
                     .OR. (Particles % FEMifGrounded .AND. &
                     ((Particles % movegl .AND. ANY(mask(maskperm(nodeindexes)) > 0.0_dp)) &
                     .OR. ((.NOT. Particles % movegl) .AND. &
                     ANY(mask(maskperm(nodeindexes)) < 0.0_dp)))))) THEN

                   BValues(BPerm(NodeIndexes)) = BValues(BPerm(NodeIndexes)) + &
                        Basis(1:nn) * (Particles % Binit(No) *&
                        (Particles % EF(No) * 2.0_dp)**(-cm) ) * detJ


                   IF (Particles % usedamage .AND. ANY(Particles % Dav(No,:).NE.0.0_dp)) THEN

                      DxxValues(DxxPerm(NodeIndexes)) = DxxValues(DxxPerm(NodeIndexes)) + &
                           Basis(1:nn) * Particles % Dav(No,1) * detJ

                      DyyValues(DyyPerm(NodeIndexes)) = DyyValues(DyyPerm(NodeIndexes)) + &
                           Basis(1:nn) * Particles % Dav(No,2) * detJ

                      DzzValues(DzzPerm(NodeIndexes)) = DzzValues(DzzPerm(NodeIndexes)) + &
                           Basis(1:nn) * Particles % Dav(No,3) * detJ

                      DxyValues(DxyPerm(NodeIndexes)) = DxyValues(DxyPerm(NodeIndexes)) + &
                           Basis(1:nn) * Particles % Dav(No,4) * detJ
                   END IF

                END IF



                IF (.NOT. Particles % uplag) THEN

                   V1Values(V1Perm(NodeIndexes)) = V1Values(V1Perm(NodeIndexes)) + &
                        Basis(1:nn) * Particles % Velocity(No,1) * Particles % Mass(No)
                   V2Values(V2Perm(NodeIndexes)) = V2Values(V2Perm(NodeIndexes)) + &
                        Basis(1:nn) * Particles % Velocity(No,2)  * Particles % Mass(No)
                END IF

                IF (t==ABS(ElemParticles(ii) % NumberOfParticles)) THEN
                   WeightValues(WeightPerm(NodeIndexes)) = WeightValues(WeightPerm(NodeIndexes)) + Area/nn
                END IF

             ELSEIF (whichtime == 3) THEN

                V1Values(V1Perm(NodeIndexes)) = V1Values(V1Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % Velocity(No,1) * &
                     Particles % Mass(No)
                V2Values(V2Perm(NodeIndexes)) = V2Values(V2Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % Velocity(No,2) * &
                     Particles % Mass(No)

                MassValues(MassPerm(NodeIndexes)) = MassValues(MassPerm(NodeIndexes)) +&
                     Basis(1:nn)  * Particles % Mass(No)

                IF (t==ABS(ElemParticles(ii) % NumberOfParticles)) THEN
                   WeightValues(WeightPerm(NodeIndexes)) = WeightValues(WeightPerm(NodeIndexes)) + 1.0_dp
                END IF

             ELSEIF (whichtime == 4) THEN

                !you have stored particle vstar(r-1) on Particles % XPIC(No,1:2)
                V1Values(V1Perm(NodeIndexes)) = V1Values(V1Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % xpic(No,1) * Particles % Mass(No)

                V2Values(V2Perm(NodeIndexes)) = V2Values(V2Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % xpic(No,2) * Particles % Mass(No)

             ELSEIF (whichtime == 5) THEN

                Hf = rhow * (sealevel-Particles % bedrock(No))/rhoi
                IF (Particles % H(No) < Hf) THEN
                   Particles % Gmask(No) = 1.0_dp
                   zs = Particles % H(No) * (1.0_dp - rhoi/rhow)
                ELSE
                   Particles % Gmask(No) = -1.0_dp
                   zs = Particles % H(No) + Particles % bedrock(No)
                END IF

                HValues(HPerm(NodeIndexes)) = HValues(HPerm(NodeIndexes)) + &
                     Basis(1:nn) * zs  *detJ  / (Area)

                IF (t==ABS(ElemParticles(ii) % NumberOfParticles)) THEN
                   WeightValues(WeightPerm(NodeIndexes)) = WeightValues(WeightPerm(NodeIndexes)) + 1.0_dp
                END IF

             ELSEIF (whichtime == 6) THEN

                Values(Perm(NodeIndexes)) = Values(Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % Dav(No,1)  *detJ  / (Area)

                IF (t==ABS(ElemParticles(ii) % NumberOfParticles)) THEN
                   WeightValues(WeightPerm(NodeIndexes)) = WeightValues(WeightPerm(NodeIndexes)) + 1.0_dp
                END IF

             ELSEIF (whichtime == 7) THEN

                HWeightValues(HWeightPerm(NodeIndexes)) = HWeightValues(HWeightPerm(NodeIndexes)) + &
                     Basis(1:nn) * detJ

                HValues(HPerm(NodeIndexes)) = HValues(HPerm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % H(No)  *detJ

             ELSEIF (whichtime == 8) THEN

                V1Values(V1Perm(NodeIndexes)) = V1Values(V1Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % NextCoordinate(No,1) * &
                     Particles % Mass(No)
                V2Values(V2Perm(NodeIndexes)) = V2Values(V2Perm(NodeIndexes)) + &
                     Basis(1:nn) * Particles % NextCoordinate(No,2) * &
                     Particles % Mass(No)

                MassValues(MassPerm(NodeIndexes)) = MassValues(MassPerm(NodeIndexes)) +&
                     Basis(1:nn)  * Particles % Mass(No)

                IF (t==ABS(ElemParticles(ii) % NumberOfParticles)) THEN
                   WeightValues(WeightPerm(NodeIndexes)) = WeightValues(WeightPerm(NodeIndexes)) + 1.0_dp
                END IF

             END IF
          END DO
       END IF
    END DO


    !PassValues has a positive integer for nodes with contributing particles
    !this integer corresponds to how many elements with particles that node is associated with

    IF (whichtime == 1) THEN
       DO ii = 1,Mesh % NumberOfNodes

          IF (WeightValues(WeightPerm(ii))>0) THEN

             Values(Perm(ii)) = Values(Perm(ii))* (4.0_dp/WeightValues(WeightPerm(ii)))

             IF (Values(Perm(ii))< 0.0000001_dp) THEN
                Values(Perm(ii))= 0.0000001_dp
             END IF
          ELSE
             Values(Perm(ii)) = MAX(0.025_dp,invValues(invPerm(ii)))
             Values(Perm(ii)) = MIN(2.0_dp,Values(Perm(ii)))
          END IF
       END DO


    ELSEIF (whichtime == 2) THEN

       DO ii = 1,Mesh % NumberOfNodes

          IF (MassValues(MassPerm(ii))>0) THEN

             IF (.NOT. Particles % uplag) THEN
                V1Values(V1Perm(ii)) = V1Values(V1Perm(ii)) / MassValues(MassPerm(ii))
                V2Values(V2Perm(ii)) = V2Values(V2Perm(ii)) / MassValues(MassPerm(ii))
             END IF

             BValues(BPerm(ii)) = BValues(BPerm(ii)) / HWeightValues(HWeightPerm(ii))


             !If Particles % weighth or FEM element,
             !H will be weighted by the sum of area of particle contributions (Hweight)
             IF (Particles % weighth .OR. PMValues(PMPerm(ii)) > 1.0_dp) THEN

                HValues(HPerm(ii)) = HValues(HPerm(ii)) / HWeightValues(HWeightPerm(ii))

             ELSE
                !otherwise, H is weighted by number of elements contributing to node
                HValues(HPerm(ii)) = HValues(HPerm(ii)) / WeightValues(WeightPerm(ii))
             END IF

             IF (Particles % usedamage) THEN
                DxxValues(DxxPerm(ii)) = DxxValues(DxxPerm(ii))/ HWeightValues(HWeightPerm(ii))
                DyyValues(DyyPerm(ii)) = DyyValues(DyyPerm(ii))/ HWeightValues(HWeightPerm(ii))
                DzzValues(DzzPerm(ii)) = DzzValues(DzzPerm(ii))/ HWeightValues(HWeightPerm(ii))
                DxyValues(DxyPerm(ii)) = DxyValues(DxyPerm(ii))/ HWeightValues(HWeightPerm(ii))
             END IF

          ELSE

             HValues(HPerm(ii)) = 1.0_dp
             MassValues(MassPerm(ii)) = 0.0_dp
             V1Values(V1Perm(ii)) = 0.0_dp
             V2Values(V2Perm(ii)) = 0.0_dp
          END IF
       END DO


       ! CALL UpdateVelocityBoundsOnMesh( Model )
       IF (Particles % UseBCforPrevVel) THEN
          CALL BCVelocityUpdate( Model, V1Var, V1Perm, V1Values, V2Var, V2Perm, V2Values )
       END IF


       IF (Particles % UseHBC) THEN
          CALL BCHUpdate( Model, HVar, HPerm, HValues )
       END IF

       IF (Particles % SetYVelZero) V2Values = 0.0_dp


    ELSEIF (whichtime == 3 .OR. whichtime == 8) THEN

       DO ii = 1,Mesh % NumberOfNodes

          IF (MassValues(MassPerm(ii))>0) THEN

             V1Values(V1Perm(ii)) = V1Values(V1Perm(ii)) /MassValues(MassPerm(ii))
             V2Values(V2Perm(ii)) = V2Values(V2Perm(ii)) /MassValues(MassPerm(ii))

          ELSE
             V1Values(V1Perm(ii)) = 0.0_dp
             V2Values(V2Perm(ii)) = 0.0_dp
          END IF

       END DO

       ! CALL UpdateVelocityBoundsOnMesh( Model )
       IF (Particles % UseBCforPrevVel) THEN
          CALL BCVelocityUpdate( Model, V1Var, V1Perm, V1Values, V2Var, V2Perm, V2Values )
       END IF



       IF (Particles % SetYVelZero) V2Values = 0.0_dp


    ELSEIF (whichtime == 4) THEN
       !XPIC

       DO ii = 1,Mesh % NumberOfNodes

          IF (MassValues(MassPerm(ii))>0) THEN

             ! Mass should already be on mesh.

             V1Values(V1Perm(ii)) = V1Values(V1Perm(ii)) /MassValues(MassPerm(ii))
             V2Values(V2Perm(ii)) = V2Values(V2Perm(ii)) /MassValues(MassPerm(ii))

          ELSE
             V1Values(V1Perm(ii)) = 0.0_dp
             V2Values(V2Perm(ii)) = 0.0_dp
          END IF
       END DO

       ! CALL UpdateVelocityBoundsOnMesh( Model )
       CALL XPICBCVelocityUpdate( Model, V1Var, V1Perm, V1Values, V2Var, V2Perm, V2Values)

    ELSEIF (whichtime == 5) THEN

       DO ii = 1,Mesh % NumberOfNodes

          IF (HValues(HPerm(ii)) > 0) THEN
             HValues(HPerm(ii)) = HValues(HPerm(ii)) * (4.0_dp/WeightValues(WeightPerm(ii)))
          END IF

       END DO

    ELSEIF (whichtime == 6) THEN

       DO ii = 1,Mesh % NumberOfNodes

          IF (Values(Perm(ii)) > 0) THEN
             Values(Perm(ii)) = Values(Perm(ii)) * (4.0_dp/WeightValues(WeightPerm(ii)))
          END IF

       END DO

    ELSEIF (whichtime == 7) THEN

       DO ii = 1,Mesh % NumberOfNodes
          HValues(HPerm(ii)) = HValues(HPerm(ii)) / HWeightValues(HWeightPerm(ii))
       END DO

       IF (Particles % UseHBC) THEN
          CALL BCHUpdate( Model, HVar, HPerm, HValues )
       END IF

    END IF

    CALL INFO(Trim(SolverName), &
         '-----ALL DONE----------',Level=5)

  END SUBROUTINE MPMParticlesToNodes

  !**************************************************************************

  !> For interpolating scalar grid variables to particles
  SUBROUTINE MPMMeshScalarToParticle( Particles, Model, whichtime)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    TYPE(Variable_t), POINTER :: Zs,F,B,Mask,MB,PassVar,H,EF,Bed
    TYPE(Element_t), POINTER :: BulkElement
    TYPE(Nodes_t)   :: ElementNodes
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp),ALLOCATABLE :: Basis(:),dBasisdx(:,:)
    INTEGER, POINTER :: LocalPerm(:)
    REAL(KIND=dp), POINTER :: maskLocalField(:),ZsLocalField(:),FLocalField(:),&
         BLocalField(:),HLocalField(:),EFLocalField(:),BedLocalField(:),MBLocalField(:)
    INTEGER :: ii,jj,nn,dim,nb,No,t,whichtime,mm
    INTEGER, POINTER :: NodeIndexes(:),ZsPerm(:),HPerm(:),FPerm(:),&
         BPerm(:),MaskPerm(:),MBPerm(:),PassPerm(:),EFPerm(:),BedPerm(:)
    LOGICAL :: Visited=.FALSE.,stat,GotIt,movegl,constantef,constantmb,constantfric
    REAL(KIND=dp) :: scale,detj,tempdmask,gamma, CriticalDamage ,&
         efparam,mbparam,fricparam,g,rhoi,rhow,SqrtElementMetric
    REAL(KIND=dp), POINTER :: ZsVal(:),FVal(:),BVal(:),MaskVal(:),&
         MBVal(:),PassValues(:),HVal(:),EFVal(:),BedVal(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    INTEGER :: minstatus

    SAVE :: Mesh,nn,nb, LocalPerm,ZsLocalField,FLocalField,&
         BLocalField,MaskLocalField,dim,Visited,&
         HLocalField,EFLocalField,SolverName,movegl,BedLocalField,MBLocalField,&
         constantef,constantmb,constantfric,efparam,mbparam,fricparam,g,rhoi,rhow,&
         Basis,dBasisdx,ElementNodes

    IF( .NOT. Visited ) THEN

       WRITE(SolverName, '(A)') 'MPMMeshScalarToParticle'

       Mesh => GetMesh()
       nn = Mesh % MaxElementNodes
       nb = Mesh % NumberOfBulkElements


       ALLOCATE( Basis(nn), dBasisdx(nn, 3))
       ALLOCATE(ElementNodes % x(nn),ElementNodes % y(nn),ElementNodes % z(nn))

       ALLOCATE( LocalPerm(nn) )
       ALLOCATE( ZsLocalField(nn), FLocalField(nn), &
            BLocalField(nn),MaskLocalField(nn),&
            HLocalField(nn),EFLocalField(nn), &
            BedLocalField(nn), MBLocalField(nn))

       g = ABS(Particles % gravity)
       rhoi = Particles % rhoi
       rhow = Particles % rhow
       dim = 2

       movegl = Particles % movegl
       ConstantMB = Particles % constmb
       ConstantEF = Particles % constef
       ConstantFric = Particles % constfric

       IF (ConstantMB) THEN
          mbparam = GetConstReal( Model % Constants, 'mbparam', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Need to define "mbparam = Real $mbparam" in constants')
       END IF

       IF (ConstantEF) THEN
          efparam = GetConstReal( Model % Constants, 'efparam', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Need to define "efparam = Real $efparam" in constants')
       END IF

       IF (ConstantFric) THEN
          fricparam = GetConstReal( Model % Constants, 'fricparam', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Need to define "fricparam = Real $fricparam" in constants')
       END IF

       Visited = .TRUE.
    END IF


    IF (whichtime == 1) THEN

       H => VariableGet(Model % Mesh % Variables, 'dHdt' )
       HPerm => H % Perm
       HVal => H % Values

       HLocalField = 0.0_dp

    ELSEIF (whichtime == 2) THEN

       IF (movegl) THEN

          Bed => VariableGet(Model % Mesh % Variables, 'Bed' )
          BedPerm => Bed % Perm
          BedVal => Bed % Values

          Particles % Bedrock(:) = 0.0_dp
          Particles % Gmask(:) = 0.0_dp

       END IF

       mask => VariableGet(Model % Mesh % Variables, 'Mask' )
       maskPerm => mask % Perm
       maskVal => mask % Values

       MaskLocalField = 0.0_dp

       IF (.NOT. ConstantMB) THEN
          MB => VariableGet(Model % Mesh % Variables, 'MB' )
          MBPerm => MB % Perm
          MBVal => MB % Values

          MBLocalField = 0.0_dp
       ELSE
          Particles % MB(:) = mbparam
       END IF

       IF (.NOT. ConstantEF) THEN

          EF => VariableGet(Model % Mesh % Variables, 'EF' )
          EFPerm => EF % Perm
          EFVal => EF % Values

          EFLocalField = 0.0_dp
       ELSE
          Particles % EF(:) = efparam
       END IF

       IF (.NOT. ConstantFric) THEN
          F => VariableGet(Model % Mesh % Variables, 'FP' )
          FPerm => F % Perm
          FVal => F % Values
          FLocalField = 0.0_dp
       ELSE
          Particles % FP(:) = fricparam
       END IF

       DO No = 1, Particles % NumberOfParticles
          IF (Particles % InterpElem(No) .NE. Particles % ElementIndex(No)) CYCLE
          IF (Particles % UseInterpElem(No)) CYCLE

          IF (.NOT. ConstantEF) Particles % EF(No) = 0.0_dp
          IF (.NOT. ConstantFric) Particles % FP(No) = 0.0_dp
          IF (.NOT. ConstantMB) Particles % MB(:) = 0.0_dp
          Particles % Gmask(No) = 0.0_dp
       END DO


    ELSEIF (whichtime == 3) THEN
       !we put mass balance on the friction parameter since friction not needed for
       !the thickness solver, and will be reinstated on the friction parameter
       !before the SSASolver
       Particles % FP(:) = 0.0_dp

       F => VariableGet(Model % Mesh % Variables, 'FP' )
       FPerm => F % Perm
       FVal => F % Values

       MB => VariableGet(Model % Mesh % Variables, 'MB' )
       MBPerm => MB % Perm
       MBVal => MB % Values


       FLocalField = 0.0_dp

    ELSEIF (whichtime == 4) THEN

       Particles % GradZs(:,:) = 0.0_dp

       Zs => VariableGet(Model % Mesh % Variables, 'Zs' )
       ZsPerm => Zs % Perm
       ZsVal => Zs % Values

       ZsLocalField = 0.0_dp

    ELSEIF (whichtime == 5) THEN

       B => VariableGet(Model % Mesh % Variables, 'invvisc' )
       BPerm => B % Perm
       BVal => B % Values

       BLocalField = 0.0_dp

       DO No = 1, Particles % NumberOfParticles
          IF (Particles % InterpElem(No) .NE. Particles % ElementIndex(No)) CYCLE
          IF (Particles % UseInterpElem(No)) CYCLE
          IF (Particles % GMask(No) >= 0.9_dp) CYCLE
          IF (Particles % Gmask(No) >= 0.0_dp .AND. Particles % Gmask(No) < 0.9_dp) THEN
             IF (ANY(Particles % Dav(No,:).NE.0.0_dp)) CYCLE
          END IF

          Particles % Binit(No) = 0.0_dp
       END DO

    ELSEIF (whichtime == 6) THEN

       mask => VariableGet(Model % Mesh % Variables, 'Mask' )
       maskPerm => mask % Perm
       maskVal => mask % Values

       MaskLocalField = 0.0_dp

       Particles % Gmask(:) = 0.0_dp

    ELSEIF (Whichtime == 7) THEN

       MB => VariableGet(Model % Mesh % Variables, 'MB' )
       MBPerm => MB % Perm
       MBVal => MB % Values


       DO No = 1, Particles % NumberOfParticles

          IF (Particles % UseInterpElem(No)) CYCLE

          Particles % MB(No) = 0.0_dp

       END DO

    ELSEIF (whichtime == 8) THEN

       mask => VariableGet(Model % Mesh % Variables, 'Mask' )
       maskPerm => mask % Perm
       maskVal => mask % Values

       H => VariableGet(Model % Mesh % Variables, 'H' )
       HPerm => H % Perm
       HVal => H % Values


       Zs => VariableGet(Model % Mesh % Variables, 'Zs' )
       ZsPerm => Zs % Perm
       ZsVal => Zs % Values

       ZsLocalField = 0.0_dp

       Particles % Gmask(:) = 0.0_dp
       Particles % H(:) = 0.0_dp
       Particles % GradZs(:,:) = 0.0_dp

    ELSEIF (whichtime == 9) THEN

       H => VariableGet(Model % Mesh % Variables, 'Hinit' )
       HPerm => H % Perm
       HVal => H % Values

       HLocalField = 0.0_dp
    END IF


    LocalPerm = 0

    DO ii = 1,nb


       IF ( ElemTrack(ii) % Status >= FEM ) THEN

          BulkElement => Model % Mesh % Elements( ii )
          NodeIndexes => BulkElement % NodeIndexes
          CALL GetElementNodes(ElementNodes,BulkElement)

          IF (whichtime == 1 .OR. whichtime == 9) THEN
             !THICKNESS, add dhdt from mesh to particle H
             LocalPerm(1:nn) = HPerm(BulkElement % NodeIndexes)
             HLocalField(1:nn) = HVal(LocalPerm(1:nn))

          ELSE IF (whichtime == 2) THEN

             !MASK
             MaskLocalField(1:nn) = maskVal(MaskPerm(BulkElement % NodeIndexes(1:nn)))

             !BED
             IF (movegl) THEN
                BedLocalField(1:nn) = BedVal(BedPerm(BulkElement % NodeIndexes(1:nn)))
             END IF

             !FRICTION
             IF (.NOT. constantfric) THEN
                FLocalField(1:nn) = FVal(FPerm(BulkElement % NodeIndexes(1:nn)))
             END IF

             ! EF
             IF (.NOT. constantef) THEN
                EFLocalField(1:nn) = EFVal(EFPerm(BulkElement % NodeIndexes(1:nn)))
             END IF

             ! MB
             IF (.NOT. constantmb) THEN
                MBLocalField(1:nn) = MBVal(MBPerm(BulkElement % NodeIndexes(1:nn)))
             END IF

          ELSE IF (whichtime == 3 ) THEN
             !  Mass balance
             LocalPerm(1:nn) = MBPerm(BulkElement % NodeIndexes)
             MBLocalField(1:nn) = MBVal(LocalPerm(1:nn))

          ELSE IF (whichtime == 4) THEN

             ZsLocalField = 0.0_dp

             !ZSGRAD
             LocalPerm(1:nn) = ZsPerm(BulkElement % NodeIndexes)
             ZsLocalField(1:nn) = ZsVal(LocalPerm(1:nn))

          ELSE IF (whichtime == 5) THEN
             !BINIT
             BLocalField(1:nn) = BVal(BPerm(BulkElement % NodeIndexes(1:nn)))

          ELSE IF (whichtime == 6) THEN
             !MASK
             MaskLocalField(1:nn) = maskVal(MaskPerm(BulkElement % NodeIndexes(1:nn)))
          ELSE IF (whichtime == 8) THEN

             ZsLocalField = 0.0_dp

             !ZSGRAD
             LocalPerm(1:nn) = ZsPerm(BulkElement % NodeIndexes)
             ZsLocalField(1:nn) = ZsVal(LocalPerm(1:nn))

          END IF


          DO t = 1, ABS(ElemParticles(ii) % NumberOfParticles)
             No = ElemParticles(ii) % p(t)


             IF (Particles % ShapeFunctions == 'gimpm') THEN
                stat = GIMPMElementInfo( t,Particles, Model,BulkElement, ElementNodes, No, &
                     detJ, scale, .TRUE., Basis,dBasisdx)
             ELSE

                stat = sMPMElementInfo( BulkElement, Particles, Model, ElementNodes, No, &
                     Particles % gridres, Basis,dBasisdx)
                scale = 1.0_dp
             END IF



             IF (whichtime == 1) THEN

                Particles % H(No) = Particles % H(No) + &
                     SUM(Basis(1:nn) * HLocalField(1:nn)) * scale * Particles % dtime

             ELSE IF (whichtime == 2) THEN


                IF (Particles % InterpElem(No) .NE. Particles % ElementIndex(No)) CYCLE
                IF (Particles % UseInterpElem(No)) CYCLE

                IF (moveGL) THEN

                   Particles % Bedrock(No) = Particles % Bedrock(No) + &
                        SUM(Basis(1:nn) * BedLocalField(1:nn)) * scale

                   !if using floatation solver, mask is reversed on the mesh
                   Particles % GMask(No) = Particles % GMask(No) - &
                        SUM(Basis(1:nn) * maskLocalField(1:nn)) * scale
                ELSE
                   Particles % GMask(No) = Particles % GMask(No) + &
                        SUM(Basis(1:nn) * maskLocalField(1:nn)) * scale
                END IF



                IF (.NOT. constantmb) THEN
                   Particles % MB(No) = Particles % MB(No) + &
                        SUM(Basis(1:nn) * MBLocalField(1:nn)) * scale
                END IF

                IF (.NOT. constantfric) THEN
                   Particles % FP(No) = Particles % FP(No) + &
                        SUM(Basis(1:nn) * FLocalField(1:nn)) * scale
                END IF

                IF (.NOT. constantEF) THEN
                   Particles % EF(No) = Particles % EF(No) + &
                        SUM(Basis(1:nn) * EFLocalField(1:nn)) * scale

                   Particles % EF(No) = MAX(Particles % EF(No),0.0_dp)
                END IF


             ELSE IF (whichtime == 3) THEN

                !Use friction parameter to store mass balance for now
                Particles % FP(No) = Particles % FP(No) + &
                     SUM(Basis(1:nn) * MBLocalField(1:nn)) * scale

             ELSE IF (whichtime == 4) THEN
                !ZSGRAD
                DO jj = 1,dim
                   Particles % GradZs(No,jj) = Particles % GradZs(No,jj) + &
                        SUM(dBasisdx(1:nn,jj) * ZsLocalField(1:nn)) * scale
                END DO

             ELSE IF (whichtime == 5) THEN

                IF (Particles % InterpElem(No) .NE. Particles % ElementIndex(No)) CYCLE
                IF (Particles % UseInterpElem(No)) CYCLE

                IF (Particles % GMask(No) >= 0.9_dp) CYCLE

                IF (Particles % Gmask(No) >= 0.0_dp .AND. Particles % Gmask(No) < 0.9_dp) THEN
                   IF (ANY(Particles % Dav(No,:).NE.0.0_dp)) CYCLE
                END IF

                Particles % Binit(No) = Particles % Binit(No) + &
                     SUM(Basis(1:nn) * BLocalField(1:nn)) * scale

             ELSE IF (whichtime == 6) THEN

                Particles % GMask(No) = Particles % GMask(No) - &
                     SUM(Basis(1:nn) * maskLocalField(1:nn)) * scale

             ELSE IF (whichtime == 7) THEN

                IF (Particles % UseInterpElem(No)) CYCLE

                Particles % MB(No) = Particles % MB(No) + &
                     SUM(Basis(1:nn) * MBVal(MBPerm(NodeIndexes(1:nn)))) * scale

             ELSE IF (whichtime == 8) THEN

                IF (Particles % UseInterpElem(No)) CYCLE

                Particles % H(No) = Particles % H(No) + &
                     SUM(Basis(1:nn) * HVal(HPerm(NodeIndexes(1:nn)))) * scale

                DO jj = 1,dim
                   Particles % GradZs(No,jj) = Particles % GradZs(No,jj) + &
                        SUM(dBasisdx(1:nn,jj) * ZsLocalField(1:nn)) * scale
                END DO

                IF (moveGL) THEN
                   !you are using floatation solver.
                   !which outputs a mask opposite to our convention on the particles
                   Particles % GMask(No) = Particles % GMask(No) - &
                        SUM(Basis(1:nn) * &
                        MaskVal(MaskPerm(NodeIndexes(1:nn)))) *scale

                ELSE
                   Particles % GMask(No) = Particles % GMask(No) + &
                        SUM(Basis(1:nn) * &
                        MaskVal(MaskPerm(NodeIndexes(1:nn)))) *scale

                END IF

             ELSEIF (whichtime == 9) THEN

                IF (Particles % xpic(No,1) > 0.0_dp) THEN
                   Particles % H(No) = Particles % H(No) + &
                        SUM(Basis(1:nn) * HLocalField(1:nn)) * scale
                END IF
             END IF
          END DO
       END IF
    END DO

  END SUBROUTINE MPMMeshScalarToParticle

  !**************************************************************************

  !> For interpolating vector grid variables to particles
  SUBROUTINE MPMMeshVectorToParticle( Particles, Model, whichtime, count)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Variable_t), POINTER :: V,PV,Bz,Bzinit,gmask,PM
    TYPE(Element_t),POINTER :: BulkElement
    TYPE(Model_t) :: Model
    TYPE(Nodes_t)   :: ElementNodes
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp),ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    INTEGER, POINTER :: LocalPerm(:),PVLocalPerm(:),BzLocalPerm(:),gmaskPerm(:),&
         gmasklocalperm(:),PMLocalPerm(:)
    REAL(KIND=dp), POINTER :: LocalField(:,:),PVLocalField(:,:),BzLocalField(:,:),&
         templayervector(:),newviscz(:),Dam(:), gmasklocalfield(:),&
         PMLocalField(:)
    INTEGER :: ii,jj,kk,nn,dim,nb,No,t,dofs,whichtime,count,particlelayers,templayers
    INTEGER, POINTER :: NodeIndexes(:),VPerm(:),PVPerm(:),BzPerm(:),PMPerm(:)
    LOGICAL :: Visited=.FALSE.,stat, GotIt,interplayervisit = .FALSE.,noinitdam=.FALSE.,&
         constlintemp,useconsttemp
    REAL(KIND=dp) :: scale,detJ,Arr,Temp,rhoi,rhow,D,CriticalDav,&
         H,btzav,Dmax,binit,surftemp,basetemp,consttemp,SqrtElementMetric,rr
    REAL(KIND=dp), POINTER :: VVal(:),PVVal(:),BzVal(:),gmaskval(:),PMVal(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    INTEGER :: minstatus
    REAL(KIND=dp) :: Exx,Eyy,Exy1,Exy2,Exy,Ezz,Ee

    SAVE :: Mesh,nn,nb,LocalPerm,LocalField,dim,Visited,&
         PVLocalPerm,PVLocalField,BzLocalPerm,BzLocalField,templayers,&
         particlelayers,templayervector,rhoi,rhow,CriticalDav,Dmax,newviscz,&
         SolverName,interplayervisit,Dam,noinitdam,constlintemp,surftemp,&
         basetemp,useconsttemp,consttemp,Basis,dBasisdx,GMaskLocalPerm,&
         GMaskLocalField,ElementNodes,PMLocalPerm,PMLocalField

    IF( .NOT. Visited ) THEN
       Mesh => GetMesh()
       nn = Mesh % MaxElementNodes
       nb = Mesh % NumberOfBulkElements

       ALLOCATE( Basis(nn),dBasisdx(nn,3))
       ALLOCATE( ElementNodes % x(nn),ElementNodes % y(nn),ElementNodes % z(nn))

       dim = 2

       WRITE(SolverName, '(A)') 'MPMMeshVectorToParticle'

       templayers = Particles % numberoftemperaturelayers
       particlelayers = Particles % numberofparticlelayers
       rhow = Particles % rhow
       rhoi = Particles % rhoi
       IF (Particles % useriftdmax) THEN
          Dmax = Particles % riftdmax
       ELSE
          Dmax  = Particles % DmaxI
       END IF

       IF (Particles % useisodam) THEN
          CriticalDav = Particles % isodamcritdav
          Dmax = Particles % initdmax
       ELSE
          CriticalDav = Particles % CriticalDav
       END IF

       constlintemp = Particles % constlintemp

       useconsttemp = Particles % useconsttemp
       IF (useconsttemp) consttemp = Particles % consttemp

       noinitdam = GetLogical( Model % Constants, 'No Init Dam', GotIt )
       IF (.NOT. GotIt) THEN
          CALL Warn(SolverName, &
               'No Init Dam not defined in constants, assuming false!')
          noinitdam = .FALSE.
       END IF

       IF (constlintemp) THEN
          surftemp = GetConstReal( Model % Constants, 'surftemp', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Since using Constant Linear Temperature, need to define "surftemp = Real $" in constants')

          basetemp = GetConstReal( Model % Constants, 'basetemp', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Since using Constant Linear Temperature, need to define "basetemp = Real $" in constants')
       END IF

       ALLOCATE( templayervector(templayers))
       ALLOCATE( LocalPerm(nn), LocalField(nn,dim) )
       ALLOCATE( PVLocalPerm(nn), PVLocalField(nn,dim) )
       ALLOCATE( PMLocalPerm(nn), PMLocalField(nn) )
       ALLOCATE( BzLocalPerm(nn), BzLocalField(nn,templayers) )
       ALLOCATE( GMaskLocalPerm(nn), GMaskLocalField(nn))
       ALLOCATE( newviscz(particlelayers), Dam(particlelayers) )

       Visited = .TRUE.
    END IF

    LocalPerm = 0
    LocalField = 0.0_dp

    IF (whichtime < 4) THEN

       V => VariableGet(Model % Mesh % Variables, 'SSAVelocity' )
       VPerm => V % Perm
       VVal => V % Values

       PV => VariableGet(Model % Mesh % Variables, 'PrevVel' )
       PVPerm => PV % Perm
       PVVal => PV % Values
    END IF

    IF (whichtime == 1) THEN
       Particles % GradVel = 0.0_dp

       IF (Particles % VelocityDependentFriction) THEN
          Particles % GridVelocity = 0.0_dp
       END IF

    ELSEIF (whichtime == 2) THEN

       PVLocalPerm = 0
       PVLocalField = 0.0_dp


       IF (count < 1) THEN
          Particles % Velocity = 0.0_dp
       END IF

    ELSEIF (whichtime == 3) THEN

       Particles % GridVelocity = 0.0_dp


    ELSEIF (whichtime == 4) THEN
       Particles % xpic(:,1:3) = 0.0_dp

       V => VariableGet(Model % Mesh % Variables, 'Vplus' )
       VPerm => V % Perm
       VVal => V % Values

    ELSEIF (whichtime == 5) THEN

       IF ((.NOT. constlintemp) .AND. (.NOT. useconsttemp) &
            .AND. (.NOT. Particles % usegiveneta)) THEN
          Bz => VariableGet(Model % Mesh % Variables, 'temperature' )
          BzPerm => Bz % Perm
          BzVal => Bz % Values

       END IF

    END IF


    DO No = 1, Particles % NumberOfParticles
       IF (Particles % InterpElem(No) .NE. Particles % ElementIndex(No)) CYCLE
       IF (Particles % UseInterpElem(No)) CYCLE

       IF (whichtime == 5) THEN
          IF (count==1) THEN
             Particles % bz(No,:) = 0.0_dp
          ELSE
             IF (constlintemp) RETURN
             IF (useconsttemp) RETURN
             IF (Particles % usegiveneta) RETURN
             IF (Particles % Gmask(No) < 0.0_dp) CYCLE
             IF (Particles % Gmask(No) >= 0.0_dp .AND. Particles % Gmask(No) < 0.9_dp) THEN
                IF (ANY(Particles % Dav(No,:).NE.0.0_dp)) CYCLE
                Particles % bz(No,:) = 0.0_dp
             END IF
          END IF
       END IF

    END DO

    IF (whichtime < 5)  dofs = V % Dofs



    IF ((whichtime .NE. 5) .OR. ((.NOT. constlintemp) .AND. (.NOT. useconsttemp) &
         .AND. (.NOT. Particles % usegiveneta) )) THEN

       DO ii = 1,nb

          !only use full elements
          IF (ElemTrack(ii) % Status >= FEM) THEN

             BulkElement => Mesh % Elements( ii )
             NodeIndexes => BulkElement % NodeIndexes

             CALL GetElementNodes(ElementNodes,BulkElement)

             IF (whichtime < 5) THEN
                LocalPerm = 0.0_dp
                LocalField = 0.0_dp
                LocalPerm(1:nn) = VPerm(BulkElement % NodeIndexes)
                DO kk = 1,nn
                   DO jj = 1,dim
                      LocalField(kk,jj) = VVal(dofs*(LocalPerm(kk)-1)+jj)
                   END DO
                END DO

                !Use the previous velocity field to hold the velocity difference.
                !Then, ParticleVel(t+dt) = ParticleVel(t) + (NewVel - PrevVel)
                IF ( ((whichtime == 2) .AND. (count .NE. 0)) .OR. &
                     ((whichtime == 1  .AND. count<0))) THEN
                   PVLocalPerm(1:nn) = PVPerm(BulkElement % NodeIndexes)
                   DO kk = 1,nn
                      DO jj = 1,dim
                         PVLocalField(kk,jj) = VVal(dofs*(LocalPerm(kk)-1)+jj)&
                              - PVVal(dofs*(PVLocalPerm(kk)-1)+jj)
                      END DO
                   END DO
                END IF

             ELSE IF (whichtime == 5) THEN


                BzLocalPerm(1:nn) = BzPerm(BulkElement % NodeIndexes)
                DO kk = 1,nn
                   DO jj = 1,templayers
                      BzLocalField(kk,jj) = BzVal(templayers*(BzLocalPerm(kk)-1)+jj)
                   END DO
                END DO
             END IF


             DO t = 1, ABS(ElemParticles(ii) % NumberOfParticles)
                No = ElemParticles(ii) % p(t)

                IF (No == 0) CYCLE

                IF (whichtime == 5) THEN

                   IF (Particles % InterpElem(No) .NE. Particles % ElementIndex(No)) CYCLE
                   IF (Particles % UseInterpElem(No)) CYCLE

                   IF (count>1) THEN
                      !we only want to process particles that are starting to float.
                      IF (Particles % Gmask(No) < 0.0_dp .OR. Particles % Gmask(No) >= 0.9_dp) CYCLE
                      IF (ANY(Particles % Dav(No,:).NE.0.0_dp)) CYCLE
                   END IF
                END IF


                IF (Particles % ShapeFunctions == 'gimpm') THEN
                   stat = GIMPMElementInfo( t,Particles, Model,BulkElement, ElementNodes, No, &
                        detJ, scale, .TRUE., Basis,dBasisdx)
                ELSE

                   stat = sMPMElementInfo( BulkElement, Particles, Model, ElementNodes, No, &
                        Particles % gridres, Basis,dBasisdx)

                   scale = 1.0_dp
                END IF


                IF (whichtime == 1) THEN


                   !dvx/dx
                   Particles % GradVel(No,1) = Particles % GradVel(No,1) + &
                        SUM(dBasisdx(1:nn,1) * LocalField(1:nn,1)) * scale

                   !dvy/dy
                   Particles % GradVel(No,2) = Particles % GradVel(No,2) + &
                        SUM(dBasisdx(1:nn,2) * LocalField(1:nn,2)) * scale

                   !dvx/dy
                   Particles % GradVel(No,3) = Particles % GradVel(No,3) + &
                        SUM( dBasisdx(1:nn,2) * LocalField(1:nn,1) ) * scale

                   !dvy/dx
                   Particles % GradVel(No,4) = Particles % GradVel(No,4) + &
                        SUM( dBasisdx(1:nn,1) * LocalField(1:nn,2) ) * scale

                   IF (Particles % Static(No)) THEN
                      Particles % GradVel(No,:) = 0.0_dp
                   END IF


                   IF (Particles % VelocityDependentFriction .AND. Particles % Gmask(No)<0.0_dp) THEN
                      IF (count==-1) THEN
                         DO kk = 1,2
                            Particles % GridVelocity(No,kk) = Particles % Velocity(No,kk) + &
                                 SUM(Basis(1:nn) * PVLocalField(1:nn,kk)) * scale
                         END DO
                      ELSE
                         DO kk = 1,2
                            Particles % GridVelocity(No,kk) = Particles % GridVelocity(No,kk) + &
                                 SUM(Basis(1:nn) * LocalField(1:nn,kk)) * scale
                         END DO
                      END IF
                   END IF

                ELSE IF (whichtime == 2) THEN

                   !for updating particle velocity
                   IF (count < 1) THEN
                      DO kk = 1,2
                         Particles % Velocity(No,kk) = Particles % Velocity(No,kk) + &
                              SUM(Basis(1:nn) * LocalField(1:nn,kk)) * scale
                      END DO
                   ELSE

                      DO kk = 1,2
                         Particles % Velocity(No,kk) = Particles % Velocity(No,kk) + &
                              SUM(Basis(1:nn) * PVLocalField(1:nn,kk)) * scale
                      END DO
                   END IF

                   IF (Particles % Static(No)) THEN
                      Particles % Velocity(No,:) = 0.0_dp
                   END IF



                ELSE IF (whichtime == 3) THEN
                   !particle Velocity from mesh velocity, used for updating position

                   IF (count>=0) THEN
                      DO kk = 1,2
                         Particles % GridVelocity(No,kk) = Particles % GridVelocity(No,kk) + &
                              SUM(Basis(1:nn) * LocalField(1:nn,kk)) * scale
                      END DO
                   ELSE
                      DO kk = 1,2
                         Particles % GridVelocity(No,kk) = Particles % Velocity(No,kk) + &
                              SUM(Basis(1:nn) * PVLocalField(1:nn,kk)) * scale
                      END DO
                   END IF

                   IF (Particles % Static(No)) THEN
                      Particles % Velocity(No,:) = 0.0_dp
                   END IF

                ELSE IF (whichtime == 4) THEN

                   DO kk = 1,2
                      Particles % xpic(No,kk) = Particles % xpic(No,kk) + &
                           SUM(Basis(1:nn) * LocalField(1:nn,kk)) * scale
                   END DO

                ELSE IF (whichtime == 5) THEN


                   DO kk = 1,templayers
                      templayervector(kk) = SUM(Basis(1:nn) * BzLocalField(1:nn,kk)) * scale
                   END DO

                   !InterpLayers holds the interpolation map and functions to
                   !interpolate between the temperature layers and particle layers.
                   !InterpLayers % Map(particlelayer,1) and InterpLayers % Map(particlelayer,2) give the
                   !two temperature layers that a particle layers lies between.
                   !InterpLayers % InterpFun(particlelayer,1) and InterpLayers % InterpFun(particlelayer,2)
                   !give the two weights to interpolate between the particlelayer and the two closest temperature
                   !layers specified in InterpLayers % Map(particlelayer,1:2)

                   DO kk = 1,particlelayers

                      !assign the temperature at each particle layer to Particles % bz(No,particlelayer)
                      !we convert from temperature to bz at the end of this subroutine
                      Particles % bz(No,kk) = Particles % bz(No,kk) + &
                           (templayervector(InterpLayers % Map(kk,1)) * InterpLayers % InterpFun(kk,1) + &
                           templayervector(InterpLayers % Map(kk,2)) * InterpLayers % InterpFun(kk,2))
                   END DO
                END IF
             END DO
          END IF
       END DO
    END IF


    IF (whichtime == 5) THEN
       !convert temperature at each particle layer to B(z)

       DO No = 1,Particles % NumberOfParticles


          IF ((.NOT. constlintemp) .AND. (.NOT. useconsttemp) .AND. (.NOT. Particles % usegiveneta)) THEN
             IF (Particles % InterpElem(No) .NE. Particles % ElementIndex(No)) CYCLE
             IF (Particles % UseInterpElem(No)) CYCLE

             IF (count>1) THEN
                IF (Particles % Gmask(No) < 0.0_dp .OR. Particles % Gmask(No) >= 0.9_dp) CYCLE
                IF (ANY(Particles % Dav(No,:).NE.0.0_dp)) CYCLE
             END IF
          END IF

          DO ii = 1,particlelayers

             IF (constlintemp) THEN
                Temp = basetemp + (ii-1) * ( (surftemp-basetemp)/(particlelayers-1) )
             ELSEIF (useconsttemp) THEN
                Temp = consttemp
             ELSEIF (Particles % usegiveneta) THEN
                Temp = 0.0_dp
             ELSE
                Temp = Particles % Bz(No,ii) - 273.15_dp
             END IF

             IF (Temp<-10.0_dp) THEN
                Arr=3.985e-13_dp * EXP( -60.0e03_dp/(8.314_dp * (273.15_dp + Temp)))

             ELSEIF (Temp>0.0_dp) THEN
                Arr=1.916e03_dp * EXP( -139.0e03_dp/(8.314_dp *  (273.15_dp)))

             ELSE
                Arr = 1.916e03_dp * EXP( -139.0e03_dp/(8.314_dp *  (273.15_dp + Temp)))
             END IF

             Arr = (Arr)**(-1.0_dp/3.0_dp)
             Arr = Arr*1.0E-06_dp*31556926.0_dp**(-1.0_dp/3.0_dp)
             Particles % Bz(No,ii) = Arr

             IF (constlintemp) THEN
                Particles % Bz(:,ii) = Arr
                IF (ii == particlelayers) RETURN
             ELSE IF (useconsttemp) THEN
                Particles % Bz(:,:) = Arr
                RETURN
             ELSE IF (particles % usegiveneta) THEN
                Particles % Bz(:,:) = Particles % giveneta
             END IF

             !Particles % Bz is  in the format of input for viscosity in the Elmer SSA solver
          END DO


          IF (count>1) THEN
             H = Particles % H(No)

             newviscz = Particles % Bz(No,:)

             IF (particlelayers>1) THEN
                btzav = ((SUM(newviscz)-newviscz(1)/2.0_dp-newviscz(particlelayers)/2.0_dp) &
                     *H/(DBLE(particlelayers)-1.0_dp))/H
             ELSE
                btzav = newviscz(1)
             END IF


             !here, binit is invvisc, because we aren't "fully" floating.
             D = 1.0_dp- (Particles % Binit(No)/(btzav* Particles % EF(No)**(-1.0_dp/3.0_dp)))

             D = MIN(D,Dmax)
             D = MAX(D,0.0_dp)

             IF (noinitdam) THEN
                D = 0.0_dp
             END IF

             IF (Particles % nodaminflow) THEN
                Particles % damstatus(No) = -1
             END IF

             IF (D >= CriticalDav) THEN

                IF (Particles % useisodam .OR. Particles % nodaminflow) THEN
                   Particles % damstatus(No) = -1
                   D = Particles % initdmax
                ELSE
                   D=Dmax
                END IF

                Particles % Dav(No,1:3) = D
                Particles % Dav(No,4) = 0.0_dp
                Particles % Damage(No,:,1:3) = D
                Particles % Damage(No,:,4) = 0.0_dp
                Particles % Binit(No) = btzav
             END IF

             IF (D<=0.01_dp) THEN
                D=0.0_dp
                Particles % Dav(No,:) = 0.0_dp
                Particles % Damage(No,:,:) = 0.0_dp
                Particles % EF(No) = 1.0_dp
             ELSE
                Dam(:) = 0.0_dp
                binit = Particles % Binit(No)

                Particles % Binit(No) = btzav
                btzav = btzav* Particles % EF(No)**(-1.0_dp/3.0_dp)

                IF (particlelayers > 1) THEN

                   IF (.NOT. Particles % useisodam ) THEN

                      IF (Particles % nodaminflow) THEN
                         Particles % damstatus(No) = -1
                         Dam = D
                      ELSE
                         CALL InitParticleDz(Particles, No,particlelayers,btzav,binit,D,Model,Dam,.FALSE.)
                      END IF

                   ELSE
                      Particles % damstatus(No) = -1

                      IF (D>Particles % initdmax) D = Particles % initdmax
                      Dam = D

                   END IF

                ELSE
                   Dam = D
                END IF

                IF (Particles % nodamregion) THEN

                   IF (Particles % Coordinate(No,1)>=Particles % ndxmin .AND. &
                        Particles % Coordinate(No,1)<=Particles % ndxmax) THEN
                      IF (Particles % Coordinate(No,2)>=Particles % ndymin .AND. &
                           Particles % Coordinate(No,2)<=Particles % ndymax) THEN
                         Particles % damage(No,:,:) = 0.0_dp
                         Particles % dD(No,:,:) = 0.0_dp
                         Particles % Dav(No,:) = 0.0_dp
                         D = 0.0_dp
                         Dam = 0.0_dp
                      END IF
                   END IF
                END IF

                Particles % Dav(No,1:3) = D
                Particles % Dav(No,4) = 0.0_dp

                DO jj = 1,3
                   Particles % Damage(No,:,jj) = Dam(:)
                END DO
             END IF
          END IF
       END DO
    END IF

  END SUBROUTINE MPMMeshVectorToParticle

  !************************************************************************

  !> this solver finds the single temperature that yields the
  !! same viscosity as the vertically integrated viscosity corresponding to
  !! a column of temperatures

  SUBROUTINE VertIntMeshTemp( Model)

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    INTEGER :: templayers,ii,jj,kk
    REAL(KIND=dp) :: y2s,maxtemp
    REAL(KIND=dp) :: H,Arr,btzav,Temp,TT(2),TestArr(2)
    REAL(KIND=dp),ALLOCATABLE :: ArrVec(:)
    LOGICAL  :: GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    TYPE(Variable_t), POINTER :: T,Thick,VIT
    INTEGER, POINTER :: TPerm(:),ThickPerm(:),VITPerm(:)
    REAL(KIND=dp), POINTER :: TVal(:),ThickVal(:),VITVal(:)



    WRITE(SolverName, '(A)') 'VertIntMeshTemp'

    templayers = GetInteger( Model % Constants, 'number of temperature layers', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Fatal(SolverName, &
            'Need to define "number of temperature layers = Real $" in constants')
    END IF

    y2s = GetConstReal( Model % Constants, 'year in seconds', GotIt )
    IF (.NOT. GotIt) THEN
       CALL Fatal(SolverName, &
            'Need to define "year in seconds = Real $" in constants')
    END IF

    ALLOCATE( ArrVec(templayers))

    T => VariableGet(Model % Mesh % Variables, 'temperature' )
    TPerm => T % Perm
    TVal => T % Values

    Thick => VariableGet(Model % Mesh % Variables, 'H' )
    ThickPerm => Thick % Perm
    ThickVal => Thick % Values

    !vertically integrated temperature
    VIT => VariableGet(Model % Mesh % Variables, 'VIT' )
    VITPerm => VIT % Perm
    VITVal => VIT % Values


    DO ii = 1,Model % Mesh % NumberOfNodes

       H = ThickVal(ThickPerm(ii))
       H = MAX(H,1.0_dp)

       !get all viscosities in a vector for this node
       DO jj = 1,templayers
          Temp = TVal(templayers*(TPerm(ii)-1)+jj)-273.15_dp
          Temp = MIN(Temp,0.0_dp)

          IF (Temp<-10.0_dp) THEN
             Arr=3.985e-13_dp * EXP( -60.0e03_dp/(8.314_dp * (273.15_dp + Temp)))
          ELSEIF (Temp>0.0_dp) THEN
             Arr=1.916e03_dp * EXP( -139.0e03_dp/(8.314_dp *  (273.15_dp)))
          ELSE
             Arr = 1.916e03_dp * EXP( -139.0e03_dp/(8.314_dp *  (273.15_dp + Temp)))
          END IF

          Arr = (Arr)**(-1.0_dp/3.0_dp)
          Arr = Arr*1.0E-06_dp*y2s**(-1.0_dp/3.0_dp)

          ArrVec(jj) = Arr
       END DO


       btzav = (((SUM(ArrVec)-ArrVec(1)/2.0_dp-ArrVec(templayers)/2.0_dp)*H/(DBLE(templayers)-1.0_dp) ))/H


       !Find the two possible temps corresponding to btzav
       !Backcalculating for 0>T>-10
       Arr = (btzav)/(y2s**(-1.0_dp/3.0_dp)*1.0E-06_dp);
       Arr = Arr**(-3.0_dp);
       Arr = Arr/1.916E03_dp;
       Arr = log(Arr);
       TT(1) = ((-139.0E03_dp/Arr)/8.314_dp)-273.15_dp;

       !Backcalculating for T<-10
       Arr = (btzav)/(y2s**(-1.0_dp/3.0_dp)*1.0E-06_dp);
       Arr = Arr**(-3.0_dp);
       Arr = Arr/3.985E-13_dp;
       Arr = log(Arr);
       TT(2) = ((-60.0E03_dp/Arr)/8.314_dp)-273.15_dp;


       !backcalculate the Arrhenius factor using these new temps.
       !keep the temp with corresponding Arrhenius factor that matches btzav
       DO jj = 1,2

          Temp = TT(jj)

          IF (Temp<-10.0_dp) THEN
             Arr=3.985e-13_dp * EXP( -60.0e03_dp/(8.314_dp * (273.15_dp + Temp)))
          ELSEIF (Temp>0.0_dp) THEN
             Arr=1.916e03_dp * EXP( -139.0e03_dp/(8.314_dp *  (273.15_dp)))
          ELSE
             Arr = 1.916e03_dp * EXP( -139.0e03_dp/(8.314_dp *  (273.15_dp + Temp)))
          END IF

          Arr = (Arr)**(-1.0_dp/3.0_dp)
          Arr = Arr*1.0E-06_dp*y2s**(-1.0_dp/3.0_dp)
          TestArr(jj) = ABS(Arr - btzav)
       END DO

       VITVal(VITPerm(ii)) = TT(MINLOC(TestArr,1))+273.15_dp
    END DO


  END SUBROUTINE VertIntMeshTemp


  !************************************************************************

  !> updates ice surface height on the grid. Alternative to floatation.F90
  SUBROUTINE UpdateMeshZs( Particles,Model )

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t)   :: ElementNodes
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Variable_t), POINTER :: GM,Zs,bed,H
    TYPE(Element_t), POINTER :: BulkElement
    INTEGER :: nn, nb, ii, ni
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    INTEGER, POINTER :: NodeIndexes(:),ZsPerm(:),bedPerm(:),HPerm(:),&
         GMPerm(:)
    LOGICAL :: Stat, Visited = .FALSE.,GotIt,movegl
    REAL(KIND=dp) :: rhow,rhoi,Coord(3),BedAtPoint,dt,SqrtElementMetric
    REAL(KIND=dp), POINTER :: ZsVal(:),bedVal(:),HVal(:),GMVal(:)

    SAVE :: Visited, movegl


    IF (.NOT. Visited) THEN
       movegl = Particles % movegl
       Visited =.TRUE.
    END IF

    Mesh => GetMesh()

    rhow = Particles % rhow
    rhoi = Particles % rhoi

    nn = Mesh % NumberOfNodes

    Zs => VariableGet(Model % Mesh % Variables, 'Zs' )
    ZsPerm => Zs % Perm
    ZsVal => Zs % Values

    bed => VariableGet( Model % Mesh % Variables, 'Bed')
    bedPerm => bed % Perm
    bedVal => bed % Values

    H => VariableGet(Model % Mesh % Variables, 'H' )
    HPerm => H % Perm
    HVal => H % Values

    !grounded mask (nomovegl : positive if floating)
    !groudned mask (movegl : neg if float)
    GM => VariableGet( Mesh % Variables, 'mask')
    GMPerm => GM % Perm
    GMVal => GM % Values

    DO ii = 1,nn

       IF (movegl) THEN
          IF (HVal(HPerm(ii))>0.0_dp) THEN
             IF ((HVal(HPerm(ii)) * (1.0_dp-rhoi/rhow)) > &
                  (bedVal(bedPerm(ii)) + HVal(HPerm(ii)))) THEN

                !we are floating
                GMVal(GMPerm(ii)) = -1.0_dp
                ZsVal(ZsPerm(ii)) = &
                     HVal(HPerm(ii)) * (1.0_dp-rhoi/rhow)
             ELSE
                GMVal(GMPerm(ii)) = 1.0_dp

                ZsVal( ZsPerm (ii) ) = &
                     HVal( HPerm( ii ) ) + bedVal(bedPerm(ii))
             END IF
          END IF
       ELSE
          IF (HVal(HPerm(ii))>0.0_dp) THEN
             IF (GMVal(GMPerm(ii))>0.0_dp) THEN
                ZsVal(ZsPerm(ii)) = &
                     HVal(HPerm(ii)) * (1.0_dp-rhoi/rhow)
             ELSE
                ZsVal( ZsPerm (ii) ) = &
                     HVal( HPerm( ii ) ) + bedVal(bedPerm(ii))
             END IF
          END IF
       END IF
    END DO

  END SUBROUTINE UpdateMeshZs

  !**************************************************************************

  !> After updating particle positions and then mapping previous particle velocities to
  !! the grid, update grid velocities to satisfy Dirichlet boundary conditions
  SUBROUTINE UpdateVelocityBoundsOnMesh( Model )

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: OrigSurf,Vel1,Vel2,InitVel1,InitVel2,BCTrack
    INTEGER :: nn, ii
    INTEGER, POINTER :: OrigSurfPerm(:),Vel1Perm(:),Vel2Perm(:), &
         InitVel1Perm(:),InitVel2Perm(:),BCTrackPerm(:)
    LOGICAL :: Visited = .FALSE.
    REAL(KIND=dp), POINTER :: OrigSurfVal(:),Vel1Val(:),Vel2Val(:),InitVel1Val(:),&
         InitVel2Val(:),BCTrackVal(:)

    SAVE :: Mesh,nn,OrigSurf,OrigSurfPerm,OrigSurfVal,Vel1,Vel1Perm,Vel1Val, &
         Vel2,Vel2Perm,Vel2Val,InitVel1,InitVel1Perm,InitVel1Val,InitVel2,InitVel2Perm, &
         InitVel2Val,BCTrack,BCTrackPerm,BCTrackVal,Visited

    IF (.NOT. Visited) THEN
       Mesh => GetMesh()

       nn = Mesh % NumberOfNodes

       OrigSurf => VariableGet(Model % Mesh % Variables, 'OrigSurf' )
       OrigSurfPerm =>  OrigSurf % Perm
       OrigSurfVal =>  OrigSurf % Values

       Vel1 => VariableGet( Model % Mesh % Variables, 'SSAVelocity 1')
       Vel1Perm => Vel1 % Perm
       Vel1Val => Vel1 % Values

       Vel2 => VariableGet(Model % Mesh % Variables, 'SSAVelocity 2' )
       Vel2Perm => Vel2 % Perm
       Vel2Val => Vel2 % Values

       InitVel1 => VariableGet(Model % Mesh % Variables, 'InvVel 1')
       InitVel1Perm => InitVel1 % Perm
       InitVel1Val => InitVel1 % Values

       InitVel2 => VariableGet(Model % Mesh % Variables, 'InvVel 2')
       InitVel2Perm => InitVel2 % Perm
       InitVel2Val => InitVel2 % Values

       BCTrack=> VariableGet(Model % Mesh % Variables, 'BCTrack')
       BCTrackPerm => BCTrack % Perm
       BCTrackVal => BCTrack % Values

       Visited = .TRUE.
    END IF

    DO ii = 1,nn
       IF (BCTrackVal(BCTrackPerm(ii)) == 0.0_dp)  THEN
          CYCLE
       ELSE IF (BCTrackVal(BCTrackPerm(ii)) == 1.0_dp)  THEN
          IF (OrigSurfVal(OrigSurfPerm(ii)) > 0.0_dp) THEN
             Vel1Val(Vel1Perm(ii)) = InitVel1Val(InitVel1Perm(ii))
             Vel2Val(Vel2Perm(ii)) = InitVel2Val(InitVel2Perm(ii))
          END IF
       ELSE IF (BCTrackVal(BCTrackPerm(ii)) == -1.0_dp)  THEN
          Vel1Val(Vel1Perm(ii)) = 0.0_dp
          Vel2Val(Vel2Perm(ii)) = 0.0_dp
       END IF
    END DO

  END SUBROUTINE UpdateVelocityBoundsOnMesh

  !**************************************************************************

  !> After updating particle positions and then mapping previous particle velocities to
  !! the grid, update grid velocities to satisfy Dirichlet boundary conditions
  SUBROUTINE BCVelocityUpdate( Model, Vel1Var, Vel1Perm, Vel1Val, Vel2Var, Vel2Perm, Vel2Val )

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Variable_t), POINTER :: Vel1Var,Vel2Var
    INTEGER :: t,n,nd,i,j
    INTEGER, POINTER :: Vel1Perm(:),Vel2Perm(:)
    LOGICAL :: GotIt,MaskIceRises,Visited=.FALSE.
    REAL(KIND=dp), POINTER :: Vel1Val(:),Vel2Val(:)
    REAL(KIND=dp) :: rval1(Model % MaxElementNodes),rval2(Model % MaxElementNodes)
    TYPE(Variable_t), POINTER :: IRSol
    INTEGER, POINTER :: IRPerm(:)
    REAL(KIND=dp), POINTER :: IR(:)

    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: xmin,xmax,ymin,ymax,vel(2),rotmat(2,2),theta
    LOGICAL :: zeronormalvel

    SAVE :: Visited,MaskIceRises,IRSol,IRPerm,IR

    Particles => GlobalParticles
    IF (.NOT. Visited) THEN
       IRSol => VariableGet( Model % Mesh % Variables, 'Icerises',GotIt )
       IF (GotIt) THEN
          MaskIceRises = .TRUE.
       ELSE
          MaskIceRises = .FALSE.
       END IF

       IRPerm => IRSol % Perm
       IR => IRSol % Values

       Visited = .TRUE.
    END IF


    IF (MaskIceRises) THEN
       DO t = 1,Model % NumberOfBulkElements
          CurrentElement => Model % Elements(t)

          Model % CurrentElement => Model % Elements(t)

          n  = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes(1:n)

          DO i = 1,n
             IF (IR(IRPerm(NodeIndexes(i))) > 0.0_dp) THEN
                Vel1Val(Vel1Perm(NodeIndexes(i))) = 0.0_dp
                Vel2Val(Vel2Perm(NodeIndexes(i))) = 0.0_dp
             END IF
          END DO
       END DO
    END IF


    DO t = Model % NumberOfBulkElements + 1, &
         Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

       CurrentElement => Model % Elements(t)

       Model % CurrentElement => Model % Elements(t)

       n  = GetElementNOFNodes()
       NodeIndexes => CurrentElement % NodeIndexes(1:n)


       DO i=1,Model % NumberOfBCs
          IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN

             rval1(1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'PrevVel 1',n,NodeIndexes, gotIt )

             IF (GotIt) THEN
                Vel1Val(Vel1Perm(NodeIndexes(1:n))) = rval1(1:n)
             END IF

             rval2(1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'PrevVel 2',n,NodeIndexes, gotIt )

             IF (GotIt) THEN
                Vel2Val(Vel2Perm(NodeIndexes(1:n))) = rval2(1:n)
             END IF

             zeronormalvel = ListGetLogical( Model % BCs(i) % Values, &
                  'Zero Normal Vel',GotIt)
             IF (.NOT. GotIt) zeronormalvel = .FALSE.

             IF (zeronormalvel) THEN

                xmin = MINVAL(Model % Mesh % Nodes % x(NodeIndexes(1:n)))
                xmax = MAXVAL(Model % Mesh % Nodes % x(NodeIndexes(1:n)))
                ymin = MINVAL(Model % Mesh % Nodes % y(NodeIndexes(1:n)))
                ymax = MAXVAL(Model % Mesh % Nodes % y(NodeIndexes(1:n)))

                theta = ATAN( (ymax-ymin)/(xmax-xmin) )

                IF (ymax<0.0_dp) THEN
                   theta = -theta
                END IF


                rotmat(1,1) = COS(theta)
                rotmat(2,2) = rotmat(1,1)
                rotmat(2,1) = SIN(theta)
                rotmat(1,2) = -rotmat(2,1)

                DO j = 1,n
                   vel(2) = Vel1Val(Vel1Perm(NodeIndexes(j)))
                   vel(1) = Vel2Val(Vel2Perm(NodeIndexes(j)))

                   vel = MATMUL(rotmat,vel)
                   vel(1) = 0.0_dp
                   vel = MATMUL(TRANSPOSE(rotmat),vel)

                   Vel1Val(Vel1Perm(NodeIndexes(j))) = vel(2)
                   Vel2Val(Vel2Perm(NodeIndexes(j))) = vel(1)
                END DO
             END IF
          END IF
       END DO
    END DO
  END SUBROUTINE BCVelocityUpdate


  !**************************************************************************

  !> Update thickness on grid to satify boundary conditions
  SUBROUTINE BCHUpdate( Model, HVar, HPerm, HVal )

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Variable_t), POINTER :: HVar
    INTEGER :: t,n,nd,i
    INTEGER, POINTER :: HPerm(:)
    LOGICAL :: GotIt
    REAL(KIND=dp), POINTER :: HVal(:)
    REAL(KIND=dp) :: rval1(Model % MaxElementNodes)


    DO t = Model % NumberOfBulkElements + 1, &
         Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

       CurrentElement => Model % Elements(t)

       Model % CurrentElement => Model % Elements(t)

       n  = GetElementNOFNodes()
       NodeIndexes => CurrentElement % NodeIndexes(1:n)


       DO i=1,Model % NumberOfBCs
          IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN

             rval1(1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'H',n,NodeIndexes, gotIt )

             IF (GotIt) THEN
                HVal(HPerm(NodeIndexes(1:n))) = rval1(1:n)
             END IF
          END IF

       END DO
    END DO
  END SUBROUTINE BCHUpdate

  !**************************************************************************

  !> Update XPIC velocities to satisfy boundary conditions. Usually, this isn't really necessary
  !! and the impact is minimal. However, is useful, for example, in sensitive simulations to
  !! make sure that boundaries with zero normal velocities (symmetry boundaries)
  !! are strictly enforced during XPIC.
  SUBROUTINE XPICBCVelocityUpdate( Model, Vel1Var, Vel1Perm, Vel1Val, Vel2Var, Vel2Perm, Vel2Val )

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Variable_t), POINTER :: Vel1Var,Vel2Var
    INTEGER :: t,n,nd,i,j
    INTEGER, POINTER :: Vel1Perm(:),Vel2Perm(:)
    LOGICAL :: GotIt
    REAL(KIND=dp), POINTER :: Vel1Val(:),Vel2Val(:)
    REAL(KIND=dp) :: rval1(Model % MaxElementNodes),rval2(Model % MaxElementNodes)

    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: xmin,xmax,ymin,ymax,vel(2),rotmat(2,2),theta
    LOGICAL :: zeronormalvel

    Particles => GlobalParticles

    DO t = Model % NumberOfBulkElements + 1, &
         Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

       CurrentElement => Model % Elements(t)

       Model % CurrentElement => Model % Elements(t)

       n  = GetElementNOFNodes()
       NodeIndexes => CurrentElement % NodeIndexes(1:n)

       DO i=1,Model % NumberOfBCs
          IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN

             rval1(1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Vplus 1',n,NodeIndexes, gotIt )

             IF (GotIt) THEN
                Vel1Val(Vel1Perm(NodeIndexes(1:n))) = rval1(1:n)
             END IF

             rval2(1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Vplus 2',n,NodeIndexes, gotIt )

             IF (GotIt) THEN
                Vel2Val(Vel2Perm(NodeIndexes(1:n))) = rval2(1:n)
             END IF


             zeronormalvel = ListGetLogical( Model % BCs(i) % Values, &
                  'Zero Normal Vel',GotIt)
             IF (.NOT. GotIt) zeronormalvel = .FALSE.

             IF (zeronormalvel) THEN

                xmin = MINVAL(Model % Mesh % Nodes % x(NodeIndexes(1:n)))
                xmax = MAXVAL(Model % Mesh % Nodes % x(NodeIndexes(1:n)))
                ymin = MINVAL(Model % Mesh % Nodes % y(NodeIndexes(1:n)))
                ymax = MAXVAL(Model % Mesh % Nodes % y(NodeIndexes(1:n)))

                theta = ATAN( (ymax-ymin)/(xmax-xmin) )

                IF (ymax<0.0_dp) THEN
                   theta = -theta
                END IF

                rotmat(1,1) = COS(theta)
                rotmat(2,2) = rotmat(1,1)
                rotmat(2,1) = SIN(theta)
                rotmat(1,2) = -rotmat(2,1)

                DO j = 1,n
                   vel(2) = Vel1Val(Vel1Perm(NodeIndexes(j)))
                   vel(1) = Vel2Val(Vel2Perm(NodeIndexes(j)))

                   vel = MATMUL(rotmat,vel)
                   vel(1) = 0.0_dp
                   vel = MATMUL(TRANSPOSE(rotmat),vel)

                   Vel1Val(Vel1Perm(NodeIndexes(j))) = vel(2)
                   Vel2Val(Vel2Perm(NodeIndexes(j))) = vel(1)
                END DO
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE XPICBCVelocityUpdate


  !**************************************************************************

  !> Given the element & global coordinates returns the local coordinates.
  !! The idea of this routine is to transparently block the local coordinate
  !! search from the user by directly giving the basis function values related
  !! to a global coordinate. Sloppy tolerances are used since we *should*
  !! have already located the element.
  FUNCTION ParticleElementInfo( CurrentElement, GlobalCoord, &
       SqrtElementMetric, Basis, dBasisdx ) RESULT ( stat )

    IMPLICIT NONE
    TYPE(Element_t), POINTER :: CurrentElement
    REAL(KIND=dp) :: GlobalCoord(:), SqrtElementMetric, LocalDistance
    REAL(KIND=dp) :: Basis(:)
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)
    LOGICAL :: Stat, Debug
    INTEGER :: Misses(2) = 0

    SAVE Misses

    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: LocalCoord(3),u,v,w
    INTEGER :: n

    SAVE ElementNodes

    n = CurrentElement % TYPE % NumberOfNodes
    CALL GetElementNodes(ElementNodes,CurrentElement)

    Stat = PointInElement( CurrentElement, ElementNodes, &
         GlobalCoord, LocalCoord, GlobalEps = -1.0_dp, LocalEps = 1.0e3_dp, &
         LocalDistance = LocalDistance )

    IF( .NOT. Stat ) THEN
       Misses(1) = Misses(1) + 1

       IF( MODULO( SUM( Misses ), 101 ) == 100 ) PRINT *,'Misses:',Misses

       IF( .FALSE.) THEN
          IF( .NOT. Stat ) THEN
             CALL Warn('ParticleElementInfo','Should have found the node!')
          ELSE
             CALL Warn('ParticleElementInfo','Distance from element higher than expected!')
          END IF
          PRINT *,'LocalDistance:',LocalDistance,'Element:',CurrentElement % ElementIndex
          PRINT *,'Nodes X:',ElementNodes % x(1:n) - GlobalCoord(1)
          PRINT *,'Nodes Y:',ElementNodes % y(1:n) - GlobalCoord(2)
          PRINT *,'Nodes Z:',ElementNodes % z(1:n) - GlobalCoord(3)
       END IF
       RETURN
    END IF

    u = LocalCoord(1)
    v = LocalCoord(2)
    w = LocalCoord(3)

    stat = ElementInfo( CurrentElement, ElementNodes, U, V, W, SqrtElementMetric, &
         Basis, dBasisdx )
    IF(.NOT. Stat) Misses(2) = Misses(2) + 1

  END FUNCTION ParticleElementInfo

  !**************************************************************************

  !!> Particle splitting and reassignment of particle values for sMPM and GIMPM
  SUBROUTINE ParticleSplitting(Particles, Model, numoflayers )

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: maxlength, maxDPlength,mlength,davsplitthres,sl,oldcoord(2)
    INTEGER :: No, count, NoOld, jj, ii,curr,numoflayers
    LOGICAL :: Visited=.FALSE.,GotIt,savepasf
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName


    SAVE :: Visited, maxlength, maxDPlength,davsplitthres

    WRITE(SolverName, '(A)') 'ParticleSplitting'

    IF (.NOT. Visited ) THEN
       maxlength = Particles % maxlength
       maxDPlength = Particles % maxdplength
       davsplitthres = Particles % davsplitthres

       Visited = .TRUE.
    END IF

    ! if deactivate gradient thickness split:
    ! Particles % GradH(:,:) = 0.0_dp

    savepasf = Particles % alwayssplitfour

    jj = 0
    NoOld = Particles % NumberOfParticles

    DO No = 1, NoOld

       IF (MAXVAL(ABS(Particles % Dav(No,:))) < davsplitthres) THEN
          mlength = maxlength
       ELSE
          mlength = maxDPlength
       END IF

       !we marked on xpic if a particle was within a
       !mixed grounding/ungrounding element
       IF (Particles % xpic(No,1) == 1.0_dp) THEN
          mlength = MIN(mlength,Particles % maxGLlength)
       END IF

       IF     ( (Particles % Length(No,1) > mlength) .OR. &
            (Particles % Length(No,2) > mlength) ) THEN

          !just some testing of whether alwayssplitfour is really needed...
          IF (savepasf) THEN
             IF (Particles % Length(No,1) > mlength .AND. Particles % Length(No,2) > mlength) THEN
                Particles % alwayssplitfour = .TRUE.
             ELSEIF (Particles % Length(No,1) > mlength .AND. Particles % Length(No,1)>Particles % Length(No,2)) THEN
                IF (Particles % Length(No,1)/2.0_dp < Particles % Length(No,2)) Particles % alwayssplitfour = .TRUE.
             ELSEIF (Particles % Length(No,2) > mlength .AND. Particles % Length(No,2)>Particles % Length(No,1)) THEN
                IF (Particles % Length(No,2)/2.0_dp < Particles % Length(No,1)) Particles % alwayssplitfour = .TRUE.
             ELSE
                Particles % alwayssplitfour = .FALSE.
             END IF
          END IF

          IF ((Particles % Length(No,1) > mlength .AND. &
               Particles % Length(No,2) > mlength) .OR. &
               Particles % alwayssplitfour) THEN

             !split in both directions, 3 new particles
             jj = jj + 3
          ELSE
             !split in one direction, 1 new particle
             jj = jj + 1
          END IF
       END IF
    END DO

    WRITE( Message, * ) 'Number of new particles from splitting: ', jj
    CALL Info( SolverName, Message, Level=1 )


    IF (jj>0) THEN

       !Note: You should be running this subroutine after deletelostparticles
       !the reason being that allocateparticles will delete lost particles as well,
       !and this will mess up your bookkeeping here.

       !Particles % NumberOfParticles is the number of active particles, but
       !might be lower than the total number of allocated particles, which is at least
       !the number of particles you started with and is given by
       !Particles % MaxNumberOfParticles. This number can grow, but not decrease
       !(deletelostparticles just repacks the active particles from 1:n where n is
       !the new Particles % NumberOfParticles.  The "lost" particles are put in
       !n:Particles % MaxNumberOfParticles with zero values for all variables).
       !AllocateParticles will return immediately if jj+old is less than
       !MaxNumberOfParticles.  Otherwise, it will increase the allocation to jj+old.
       ! CALL AllocateParticles( Model, Particles, (jj+NoOld) )

       PRINT *, 'Allocating Particles to split'
       CALL AllocateParticles( Particles, Model, (jj+NoOld) )
       PRINT *, 'Done Allocating'

       Particles % NumberOfParticles = jj+NoOld

       Count = NoOld+1

       DO No = 1,NoOld

          IF ( MAXVAL(ABS(Particles % Dav(No,:))) < davsplitthres) THEN
             mlength = maxlength
          ELSE
             mlength = maxDPlength
          END IF

          !we marked on xpic if a particle was within a
          !mixed grounding/ungrounding element
          IF (Particles % xpic(No,1) == 1.0_dp) THEN
             mlength = MIN(mlength,Particles % maxGLlength)
          END IF


          IF   ( (Particles % Length(No,1) > mlength) .OR. &
               (Particles % Length(No,2) > mlength) ) THEN

             !just some testing of whether alwayssplitfour is really needed...
             IF (savepasf) THEN
                IF (Particles % Length(No,1) > mlength .AND. Particles % Length(No,2) > mlength) THEN
                   Particles % alwayssplitfour = .TRUE.
                ELSEIF (Particles % Length(No,1) > mlength .AND. Particles % Length(No,1)>Particles % Length(No,2)) THEN
                   IF (Particles % Length(No,1)/2.0_dp < Particles % Length(No,2)) Particles % alwayssplitfour = .TRUE.
                ELSEIF (Particles % Length(No,2) > mlength .AND. Particles % Length(No,2)>Particles % Length(No,1)) THEN
                   IF (Particles % Length(No,2)/2.0_dp < Particles % Length(No,1)) Particles % alwayssplitfour = .TRUE.
                ELSE
                   Particles % alwayssplitfour = .FALSE.
                END IF
             END IF


             IF ( (Particles % Length(No,1) > mlength .AND. &
                  Particles % Length(No,2) > mlength) .OR. Particles % alwayssplitfour) THEN

                !DUAL SPLIT
                ! PRINT *,'Double split'

                Particles % pvolume(No) = Particles % pvolume(No)/4.0_dp
                Particles % Gvolume(No) = Particles % Gvolume(No)/4.0_dp
                Particles % mass(No) = Particles % mass(No)/4.0_dp


                IF (Particles % ShapeFunctions .EQ. 'gimpm') THEN

                   Particles % OrigLength(No,:) = Particles % OrigLength(No,:)*0.5_dp

                ELSE

                   Particles % Strain(No,:) = ((Particles % Strain(No,:) + 1.0_dp)/ &
                        (2.0_dp*sqrt(0.5_dp))) - 1.0_dp

                   Particles % OrigLength(No,:) = sqrt(0.5_dp) * Particles % OrigLength(No,:)

                END IF


                DO jj = 1,3

                   curr = Count + jj - 1

                   IF (Particles % usetracer) THEN
                      Particles % Tracer(curr) = Particles % Tracer(No)
                   END IF

                   Particles % Status(curr) = Particles % Status(No)
                   Particles % Coordinate(curr,:) = Particles % Coordinate(No,:)
                   Particles % ElementIndex(curr) = Particles % ElementIndex(No)
                   Particles % InterpElem(curr) = Particles % InterpElem(No)
                   Particles % Bz(curr,:) = Particles % Bz(No,:)
                   Particles % EF(curr) = Particles % EF(No)
                   Particles % MB(curr) = Particles % MB(No)
                   Particles % Bedrock(curr) = Particles % Bedrock(No)
                   Particles % Binit(curr) = Particles % Binit(No)
                   Particles % FP(curr) = Particles % FP(No)
                   Particles % H(curr) = Particles % H(No)
                   Particles % Dav(curr,:) = Particles % Dav(No,:)
                   Particles % GradH(curr,:) = Particles % GradH(No,:)
                   Particles % GradVel(curr,:) =  Particles % GradVel(No,:)
                   Particles % GradZs(curr,:) = Particles % GradZs(No,:)
                   Particles % Velocity(curr,:) = Particles % Velocity(No,:)
                   Particles % GridVelocity(curr,:) = Particles % GridVelocity(No,:)
                   Particles % NextCoordinate(curr,:) = Particles % NextCoordinate(No,:)
                   Particles % Length(curr,1:2) =  Particles % Length(No,1:2)
                   Particles % F(curr,:) = Particles % F(No,:)
                   Particles % OrigLength(curr,:) =  Particles % OrigLength(No,:)
                   Particles % Damage(curr,:,:) = Particles % Damage(No,:,:)
                   Particles % dD(curr,:,:) = Particles % dD(No,:,:)
                   Particles % xpic(curr,:) = Particles % xpic(No,:)
                   Particles % gmask(curr) = Particles % gmask(No)
                   Particles % pvolume(curr) = Particles % pvolume(No)
                   Particles % mass(curr) = Particles % mass(No)
                   Particles % damstatus(curr) = Particles % damstatus(No)
                   Particles % Static(curr) = Particles % Static(No)
                   Particles % GVolume(curr) = Particles % GVolume(No)

                   IF (Particles % trackstrain) THEN
                      Particles % Strain(curr,:) = Particles % Strain(No,:)
                   END IF


                   Particles % OrigNo(curr) = Particles % OrigNo(No)

                   Particles % UseInterpElem(curr) = Particles % UseInterpElem(No)

                END DO

                oldcoord = Particles % Coordinate(No,:)

                !For No:
                !move x and y coord in positive direction for the splitting coord
                !(x,y)(No): xold+(Lx/4),yold+(Ly/4)
                Particles % Coordinate(No,1:2) = Particles % Coordinate(No,1:2) + &
                     Particles % Length(No,1:2)/4.0_dp

                Particles % H(No) = Particles % H(No) + &
                     Particles % GradH(No,1)*(Particles % Coordinate(No,1)-oldcoord(1)) + &
                     Particles % GradH(No,2)*(Particles % Coordinate(No,2)-oldcoord(2))

                IF (Particles % H(No) < 1.0_dp) Particles % H(No) = 1.0_dp

                !...and fix LX,LY, and vol (vol may change again
                !during particles in elem if on bound
                Particles % Length(No,:) = Particles % Length(No,:)/2.0_dp

                IF (Particles % ShapeFunctions == 'gimpm') THEN
                   Particles % GVolume(No) = Particles % Length(No,1) * &
                        Particles % Length(No,2)
                END IF


                !do the same for the rest: (x,y)(count): xold-(Lx/4),yold-(Ly/4)
                Particles % Coordinate(Count,1:2) = Particles % Coordinate(Count,1:2) - &
                     Particles % Length(Count,1:2)/4.0_dp

                Particles % H(Count) = Particles % H(Count) + &
                     Particles % GradH(Count,1)*(Particles % Coordinate(Count,1)-oldcoord(1)) + &
                     Particles % GradH(Count,2)*(Particles % Coordinate(Count,2)-oldcoord(2))

                IF (Particles % H(Count) < 1.0_dp) Particles % H(Count) = 1.0_dp

                !(x,y)(count+1): xold+(Lx/4),yold-(Ly/4)
                curr = Count + 1

                Particles % Coordinate(curr,1) = Particles % Coordinate(curr,1) + &
                     Particles % Length(curr,1)/4.0_dp
                Particles % Coordinate(curr,2) = Particles % Coordinate(curr,2) - &
                     Particles % Length(curr,2)/4.0_dp

                Particles % H(curr) = Particles % H(curr) + &
                     Particles % GradH(curr,1)*(Particles % Coordinate(curr,1)-oldcoord(1)) + &
                     Particles % GradH(curr,2)*(Particles % Coordinate(curr,2)-oldcoord(2))

                IF (Particles % H(curr) < 1.0_dp) Particles % H(curr) = 1.0_dp

                !(x,y)(count+2): xold-(Lx/4),yold+(Ly/4)
                curr = Count + 2

                Particles % Coordinate(curr,1) = Particles % Coordinate(curr,1) - &
                     Particles % Length(curr,1)/4.0_dp
                Particles % Coordinate(curr,2) = Particles % Coordinate(curr,2) + &
                     Particles % Length(curr,2)/4.0_dp


                Particles % H(curr) = Particles % H(curr) + &
                     Particles % GradH(curr,1)*(Particles % Coordinate(curr,1)-oldcoord(1)) + &
                     Particles % GradH(curr,2)*(Particles % Coordinate(curr,2)-oldcoord(2))

                IF (Particles % H(curr) < 1.0_dp) Particles % H(curr) = 1.0_dp

                DO jj = 1,3
                   curr = Count + jj - 1
                   Particles % Length(curr,:) = Particles % Length(No,:)

                   IF (Particles % ShapeFunctions == 'gimpm') THEN
                      Particles % GVolume(curr) = Particles % Length(curr,1)*Particles % Length(curr,2)
                   END IF
                END DO

                Count = Count+3

             ELSE

                !SINGLE SPLIT

                !only one coordinate to split, and only one new particle
                IF (Particles % Length(No,1) > mlength ) THEN
                   ii = 1
                   jj = 2
                ELSE
                   ii = 2
                   jj = 1
                END IF

                Particles % GVolume(No) = Particles % GVolume(No)/2.0_dp
                Particles % pvolume(No) = Particles % pvolume(No)/2.0_dp
                Particles % mass(No) = Particles % mass(No)/2.0_dp


                IF (Particles % ShapeFunctions .EQ. 'gimpm') THEN

                   Particles % OrigLength(No,ii) = Particles % OrigLength(No,ii)*0.5_dp

                ELSE

                   Particles % Strain(No,ii) = ((Particles % Strain(No,ii) + 1.0_dp)/ &
                        (2.0_dp*sqrt(0.5_dp))) - 1.0_dp

                   Particles % Strain(No,jj) = ((Particles % Strain(No,jj) + 1.0_dp)/ &
                        (sqrt(0.5_dp))) - 1.0_dp

                   Particles % OrigLength(No,:) = sqrt(0.5_dp) * Particles % OrigLength(No,:)
                END IF

                IF (Particles % usetracer) THEN
                   Particles % tracer(Count) = Particles % Tracer(No)
                END IF

                Particles % Status(Count) = Particles % Status(No)
                Particles % Coordinate(Count,:) = Particles % Coordinate(No,:)
                Particles % ElementIndex(Count) = Particles % ElementIndex(No)
                Particles % InterpElem(Count) = Particles % InterpElem(No)
                Particles % Bz(Count,:) = Particles % Bz(No,:)
                Particles % EF(Count) = Particles % EF(No)
                Particles % MB(Count) = Particles % MB(No)
                Particles % Bedrock(Count) = Particles % Bedrock(No)
                Particles % Binit(Count) = Particles % Binit(No)
                Particles % FP(Count) = Particles % FP(No)
                Particles % H(Count) = Particles % H(No)
                Particles % F(Count,:) = Particles % F(No,:)
                Particles % Dav(Count,:) = Particles % Dav(No,:)
                Particles % GradH(Count,:) = Particles % GradH(No,:)
                Particles % GradVel(Count,:) =  Particles % GradVel(No,:)
                Particles % GradZs(Count,:) = Particles % GradZs(No,:)
                Particles % Velocity(Count,:) = Particles % Velocity(No,:)
                Particles % GridVelocity(Count,:) = Particles % GridVelocity(No,:)
                Particles % NextCoordinate(Count,:) = Particles % NextCoordinate(No,:)
                Particles % Length(Count,:) =  Particles % Length(No,:)
                Particles % Damage(Count,:,:) = Particles % Damage(No,:,:)
                Particles % dD(Count,:,:) = Particles % dD(No,:,:)
                Particles % xpic(Count,:) = Particles % xpic(No,:)
                Particles % gmask(Count) = Particles % gmask(No)
                Particles % pvolume(Count) = Particles % pvolume(No)
                Particles % mass(Count) = Particles % mass(No)
                Particles % damstatus(Count) = Particles % damstatus(No)
                Particles % Static(Count) = Particles % Static(No)
                Particles % GVolume(Count) = Particles % GVolume(No)

                Particles % UseInterpElem(Count) = Particles % UseInterpElem(No)

                IF (Particles % trackstrain) THEN
                   Particles % Strain(Count,:) = Particles % Strain(No,:)
                END IF


                Particles % OrigNo(Count) = Particles % OrigNo(No)

                oldcoord = Particles % Coordinate(No,:)

                !For No:
                !move x or y coord in positive direction for the splitting coord
                !i.e. (x)(No): xold+(Lx/4)
                Particles % Coordinate(No,ii) = Particles % Coordinate(No,ii) + &
                     Particles % Length(No,ii)/4.0_dp

                !...and fix LX,LY,and vol (vol may change again during
                !particles in elem if on bound
                Particles % Length(No,ii) = Particles % Length(No,ii)/2.0_dp


                !do the same for the new particle, but in the negative direction:
                !i.e. (x)(count): xold-(Lx/4)
                Particles % Coordinate(Count,ii) = Particles % Coordinate(Count,ii) - &
                     Particles % Length(Count,ii)/4.0_dp

                Particles % Length(Count,:) = Particles % Length(No,:)
                Particles % OrigLength(Count,:) = Particles % OrigLength(No,:)


                Particles % H(No) = Particles % H(No) + &
                     Particles % GradH(No,ii)*(Particles % Coordinate(No,ii)-oldcoord(ii))

                IF (Particles % H(No) < 1.0_dp) Particles % H(No) = 1.0_dp

                Particles % H(Count) = Particles % H(Count) + &
                     Particles % GradH(Count,ii)*(Particles % Coordinate(Count,ii)-oldcoord(ii))

                IF (Particles % H(Count) < 1.0_dp) Particles % H(Count) = 1.0_dp

                IF (Particles % ShapeFunctions == 'gimpm') THEN
                   Particles % GVolume(No) = Particles % Length(No,1) * &
                        Particles % Length(No,2)

                   Particles % GVolume(Count) = Particles % Length(Count,1) * &
                        Particles % Length(Count,2)
                END IF

                Count = Count + 1
             END IF
          END IF
       END DO
    END IF

    Particles % alwayssplitfour = savepasf

  END SUBROUTINE ParticleSplitting

  !**************************************************************************

  !> Simple sMPM-style interpolation of grid vector variables and gradients
  !! to any point within an element (i.e. for grid to particle intepolation).
  !! Not as efficient as MPMMeshVectorToParticle
  SUBROUTINE GetVectorFieldInMesh(Var, CurrentElement, Basis, Velo, dBasisdx, GradVelo )

    IMPLICIT NONE
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t) :: CurrentElement
    REAL(KIND=dp) :: Basis(:), Velo(:)
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:), GradVelo(:,:)

    TYPE(Valuelist_t), POINTER :: Params
    INTEGER, POINTER :: LocalPerm(:)
    REAL(KIND=dp), POINTER :: LocalVelo(:,:)
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: VeloFieldDofs
    REAL(KIND=dp) :: SumBasis
    INTEGER :: i,j,k,n,npos,ind,dim
    LOGICAL :: GotIt, InterfaceNodes
    LOGICAL :: Visited=.FALSE.


    SAVE :: Visited, Dim, LocalVelo, LocalPerm, InterfaceNodes

    IF(.NOT. Visited ) THEN
       Mesh => GetMesh()
       Params => ListGetSolverParams()
       n = Mesh % MaxElementNodes
       ALLOCATE( LocalPerm(n), LocalVelo(n,3) )

       InterfaceNodes = GetLogical( Params,'Interface Nodes',GotIt)

       LocalPerm = 0
       LocalVelo = 0.0_dp
       Dim = Mesh % MeshDim
       Visited = .TRUE.
    END IF

    Velo = 0.0_dp
    IF( PRESENT( GradVelo ) ) GradVelo = 0.0_dp

    n = CurrentElement % TYPE % NumberOfNodes
    LocalPerm(1:n) = Var % Perm( CurrentElement % NodeIndexes )
    npos = COUNT ( LocalPerm(1:n) > 0 )

    IF( npos == 0 ) RETURN

    !-----------------------------------------------------------------
    ! compute the velocity also for case when the particle
    ! has just crossed the boundary. For example, its floating on the
    ! fluid boundary. This is a little bit fishy and could perhaps
    ! only be done conditionally....
    ! Can't really determine the gradient here
    !-----------------------------------------------------------------
    VeloFieldDofs = Var % Dofs
    IF( npos == n ) THEN
       DO i=1,n
          j = LocalPerm(i)
          DO k=1,dim
             LocalVelo(i,k) = Var % Values( VeloFieldDofs*(j-1)+k)
          END DO
       END DO
    ELSE
       IF(.NOT. InterfaceNodes ) RETURN

       SumBasis = 0.0_dp
       DO i=1,n
          j = LocalPerm(i)
          IF( j > 0 ) THEN
             SumBasis = SumBasis + Basis(i)
             DO k=1,dim
                LocalVelo(i,k) = Var % Values( VeloFieldDofs*(j-1)+k)
             END DO
          ELSE
             Basis(i) = 0.0_dp
             LocalVelo(i,1:dim) = 0.0_dp
          END IF
       END DO
    END IF

    DO i=1,dim
       Velo(i) = SUM( Basis(1:n) * LocalVelo(1:n,i) )
       IF( PRESENT( GradVelo ) ) THEN
          DO j=1,dim
             GradVelo(i,j) = SUM( dBasisdx(1:n,j) * LocalVelo(1:n,i) )
          END DO
       END IF
    END DO

    IF( npos < n ) THEN
       Velo(1:dim) = Velo(1:dim) / SumBasis
       IF( PRESENT( GradVelo ) ) THEN
          GradVelo(:,1:dim) = GradVelo(:,1:dim) / SumBasis
       END IF
    END IF

  END SUBROUTINE GetVectorFieldInMesh


  !> Simple sMPM-style interpolation of grid scalar variables and gradients
  !! to any point within an element (i.e. for grid to particle intepolation).
  !! Not as efficient as MPMMeshScalarToParticle
  SUBROUTINE GetScalarFieldInMesh(Var, CurrentElement, Basis, Pot, dBasisdx, GradPot )

    IMPLICIT NONE
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t) :: CurrentElement
    REAL(KIND=dp) :: Basis(:), Pot
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:), GradPot(:)

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: LocalPerm(:)
    REAL(KIND=dp), POINTER :: LocalField(:)
    INTEGER :: i,j,n,dim
    LOGICAL :: Visited=.FALSE.

    SAVE :: Visited, Dim, n

    IF(.NOT. Visited ) THEN
       Mesh => GetMesh()
       n = Mesh % MaxElementNodes
       Dim = Mesh % MeshDim
       Visited = .TRUE.
    END IF

    ALLOCATE( LocalPerm(n), LocalField(n) )

    LocalPerm = 0
    LocalField = 0.0_dp

    Pot = 0.0_dp
    IF( PRESENT( GradPot ) ) GradPot = 0.0_dp

    IF(.NOT. ASSOCIATED( Var ) ) RETURN

    n = CurrentElement % TYPE % NumberOfNodes
    IF( ASSOCIATED( Var % Perm ) ) THEN
       LocalPerm(1:n) = Var % Perm( CurrentElement % NodeIndexes )
       IF( .NOT. ALL ( LocalPerm(1:n) > 0 )) RETURN
       LocalField(1:n) = Var % Values( LocalPerm(1:n) )
    ELSE
       ! Some variables do not have permutation, most importantly the node coordinates
       LocalField(1:n) = Var % Values( CurrentElement % NodeIndexes )
    END IF

    Pot = SUM( Basis(1:n) * LocalField(1:n) )

    IF( PRESENT( GradPot ) ) THEN
       DO i=1,dim
          GradPot(i) = SUM( dBasisdx(1:n,i) * LocalField(1:n) )
       END DO
    END IF

    DEALLOCATE( LocalPerm, LocalField )

  END SUBROUTINE GetScalarFieldInMesh

  !**************************************************************************

  !> Finds the particle in the mesh using octree based search.
  !! This could be preferred in the initial finding of the correct elements.
  !! The major downside of the method is that there is no controlled face
  !! detection needed for wall interaction, for example.
  SUBROUTINE LocateParticleInMeshOctree( ElementIndex, GlobalCoords, &
       LocalCoords )

    USE Lists
    USE Interpolation
    USE DefUtils

    IMPLICIT NONE
    INTEGER :: ElementIndex
    REAL(KIND=dp) :: GlobalCoords(3)
    REAL(KIND=dp), OPTIONAL :: LocalCoords(3)

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Hit, Stat
    INTEGER :: i,j,k,n
    TYPE(Nodes_t), SAVE :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Quadrant_t), POINTER, SAVE :: RootQuadrant =>NULL(), LeafQuadrant
    REAL(KIND=dp) :: BoundingBox(6), eps2, eps1, uvw(3), localdumy(3)


    Mesh => GetMesh()

    ! Check that the previous hit is not hit even now
    IF( ElementIndex > 0 ) THEN

       Element => Mesh % Elements( ElementIndex )

       n = GetElementNOFNodes(Element)

       CALL GetElementNodes(ElementNodes,Element)

       IF ( PRESENT( LocalCoords) ) localdumy = localcoords
       IF ( PointInElement( Element, ElementNodes, &
            GlobalCoords, localdumy ) ) THEN !RETURN
          IF ( PRESENT( LocalCoords) )  localcoords=localdumy
          RETURN
       END IF
    END IF

    ! Find the right element using an octree search
    ! This is optimal when the particles are searched only once.
    !-----------------------------------------------------------
    IF ( .NOT.ASSOCIATED(Mesh % RootQuadrant) ) THEN
       BoundingBox(1) = MINVAL( Mesh % Nodes % x )
       BoundingBox(2) = MINVAL( Mesh % Nodes % y )
       BoundingBox(3) = MINVAL( Mesh % Nodes % z )
       BoundingBox(4) = MAXVAL( Mesh % Nodes % x )
       BoundingBox(5) = MAXVAL( Mesh % Nodes % y )
       BoundingBox(6) = MAXVAL( Mesh % Nodes % z )

       eps1 = 1.0e-3
       eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
       BoundingBox(1:3) = BoundingBox(1:3) - eps2
       BoundingBox(4:6) = BoundingBox(4:6) + eps2

       CALL BuildQuadrantTree( Mesh,BoundingBox,Mesh % RootQuadrant)

    END IF
    RootQuadrant => Mesh % RootQuadrant

    Element => NULL()
    ElementIndex = 0
    CALL FindLeafElements(GlobalCoords, Mesh % MeshDim, RootQuadrant, LeafQuadrant)
    IF ( ASSOCIATED(LeafQuadrant) ) THEN
       DO i = 1, LeafQuadrant % NElemsInQuadrant
          j = LeafQuadrant % Elements(i)
          Element => Mesh % Elements(j)

          n = GetElementNOFNodes( Element )
          CALL GetElementNodes( ElementNodes, Element)

          IF ( PointInElement( Element, ElementNodes, GlobalCoords, uvw ) ) THEN
             IF( PRESENT( LocalCoords) ) LocalCoords = uvw
             ElementIndex = j

             RETURN
          END IF
       END DO
    END IF

  END SUBROUTINE LocateParticleInMeshOctree

  !**************************************************************************

  SUBROUTINE EditParticleVolume(Particles, No,jj,N,S,E,W,xmax,xmin,ymax,ymin)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No,jj
    REAL(KIND=dp) :: N,S,E,W,xmax,xmin,ymax,ymin,s1,s2

    IF ((jj==1) .OR. (jj==2) .OR. (jj==3)) THEN
       s1 = N-ymax
    ELSE IF ((jj==4) .OR. (jj==6)) THEN
       s1 = MIN(N,ymax)-MAX(S,ymin)
    ELSE IF ((jj==7) .OR. (jj==8) .OR. (jj==9)) THEN
       s1 = ymin-S
    END IF

    IF ((jj==1) .OR. (jj==4) .OR. (jj==7)) THEN
       s2 = xmin-W
    ELSE     IF ((jj==2) .OR. (jj==8)) THEN
       s2 = MIN(E,xmax)-MAX(W,xmin)
    ELSE IF ((jj==3) .OR. (jj==6) .OR. (jj==9)) THEN
       s2 = E-xmax
    END IF

    !particle in an element, subtract volume
    !that overlaps non-associated elements
    Particles % GVolume(No) = Particles % GVolume(No) - ABS(s1*s2)

  END SUBROUTINE EditParticleVolume


  !**************************************************************************

  !> Saves particles in unstructured XML VTK format (VTU) to an external file.
  SUBROUTINE ParticleOutputVtu( Particles,Model )

    USE DefUtils
    USE MeshUtils
    USE ElementDescription
    USE AscBinOutputUtils

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    TYPE(ValueList_t),POINTER :: Params
    INTEGER, SAVE :: nTime = 0
    LOGICAL :: GotIt, Parallel, FixedMeshend,SinglePrec
    CHARACTER(MAX_NAME_LEN), SAVE :: FilePrefix
    CHARACTER(MAX_NAME_LEN) :: VtuFile, PvtuFile
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i, j, k, Partitions, Part, ExtCount, FileindexOffSet, iTime, &
         Status, MinSaveStatus, MaxSaveStatus, PrecBits, PrecSize, IntSize
    REAL(KIND=dp) :: SaveNodeFraction
    REAL(KIND=dp),POINTER :: LocalVal(:)
    LOGICAL :: BinaryOutput,AsciiOutput,Found,Visited = .FALSE.,SaveFields,SaveAll,floatingonly
    REAL(KIND=dp) :: DoubleWrk
    REAL(KIND=dp) :: SingleWrk,anisoparam,minx,maxx
    CHARACTER(MAX_NAME_LEN) :: Str,Dir
    INTEGER :: NumberOfNodes, ParallelNodes, Dim,layers,maxdofs,templayers,damdofs

    SAVE :: MinSaveStatus, MaxSaveStatus,LocalVal,Visited,layers,&
         templayers,maxdofs,floatingonly

    IF (.NOT. Visited ) THEN
       damdofs = 4
       layers = Particles % numberofparticlelayers
       templayers = Particles % numberoftemperaturelayers

       maxdofs = MAX(damdofs,layers)
       maxdofs = MAX(maxdofs,templayers)
       maxdofs = MAX(maxdofs,3)

       ALLOCATE( LocalVal( maxdofs) )
       Params => ListGetSolverParams()
       FloatingOnly = GetLogical( Params,'Floating Only',GotIt)
       IF( FloatingOnly ) THEN
          CALL Info('VtuOutputSolver','Saving Floating Particles Only!',Level=7)
       END IF

       Visited = .TRUE.
    END IF

    PRINT *,'Saving particles'

    Params => ListGetSolverParams()
    Mesh => GetMesh()

    ExtCount = ListGetInteger( Params,'Output Count',GotIt)
    IF( GotIt ) THEN
       nTime = ExtCount
    ELSE
       nTime = nTime + 1
    END IF
    FileIndexOffset = ListGetInteger( Params,'Fileindex offset',GotIt)
    iTime = nTime + FileIndexOffset

    ! IF ( nTime == 1 ) THEN
    FilePrefix = ListGetString( Params,'Filename Prefix')
    CALL Info('ParticleOutputVtu','Saving in VTK XML unstructured format to file: ' &
         //TRIM(FilePrefix)//'.vtu',Level=1)

    Dir = ListGetString( Params,'Filename Directory')

    MinSaveStatus = PARTICLE_ACTIVE
    MaxSaveStatus = PARTICLE_LOST

    BinaryOutput = GetLogical( Params,'Binary Output',GotIt)
    IF( GotIt ) THEN
       AsciiOutput = .NOT. BinaryOutput
    ELSE
       AsciiOutput = GetLogical( Params,'Ascii Output',GotIt)
       BinaryOutput = .NOT. AsciiOutput
    END IF

    SaveFields = GetLogical( Params,'Save Fields',GotIt)
    IF(.NOT. GotIt) SaveFields = .TRUE.

    SaveAll = GetLogical( Params,'Save All',GotIt)
    IF (.NOT. GotIt) SaveAll = .FALSE.

    SinglePrec = GetLogical( Params,'Single Precision',GotIt)
    IF( SinglePrec ) THEN
       CALL Info('VtuOutputSolver','Using single precision arithmetics in output!',Level=7)
    END IF


    minx = GetCReal( Params,'Min X To Save Particle',GotIt)
    IF( .NOT. GotIt ) THEN
       minx = -HUGE(1.0_dp)
    END IF

    maxx = GetCReal( Params,'Max X To Save Particle',GotIt)
    IF( .NOT. GotIt ) THEN
       maxx = HUGE(1.0_dp)
    END IF

    IF( SinglePrec ) THEN
       PrecBits = 32
       PrecSize = KIND( SingleWrk )
    ELSE
       PrecBits = 64
       PrecSize = KIND( DoubleWrk )
    END IF
    IntSize = KIND(i)

    Dim = Particles % dim

    NumberOfNodes = 0
    DO i=1,Particles % NumberOfParticles
       IF ( Particles % Status(i) > MaxSaveStatus .OR. &
            Particles % Status(i) < MinSaveStatus )  CYCLE

       IF (Particles % Coordinate(i,1) < minx) CYCLE
       IF (Particles % Coordinate(i,1) > maxx) CYCLE

       IF (FloatingOnly) THEN
          IF (Particles % Gmask(i)<0) CYCLE
       END IF


       NumberOfNodes = NumberOfNodes + 1
    END DO

    SaveNodeFraction = ListGetCReal( Params,'Particle Save Fraction',GotIt)
    IF(GotIt) THEN
       NumberOfNodes = NINT( SaveNodeFraction * NumberOfNodes )
    ELSE
       i = ListGetInteger( Params,'Particle Save Number',GotIt)
       IF( GotIt ) THEN
          NumberOfNodes = MIN(i,NumberOfNodes)
       END IF
    END IF




    IF( iTime < 10000 ) THEN
       WRITE( VtuFile,'(A,A,I4.4,".vtu")' ) TRIM(Dir),TRIM(FilePrefix),iTime
    ELSE
       WRITE( VtuFile,'(A,A,I0,".vtu")' ) TRIM(Dir),TRIM(FilePrefix),iTime
    END IF


    CALL Info('ParticleOutputVtu','Saving particles to file: '//TRIM(VtuFile),Level=8)
    CALL WriteVtuFile( VtuFile,Model )

    !---------------------- SUBROUTINE PARTICLEOUTPUTVTU -----------------------!

  CONTAINS

    SUBROUTINE WriteVtuFile( VtuFile, Model )
      IMPLICIT NONE
      TYPE(Model_t) :: Model
      CHARACTER(LEN=*), INTENT(IN) :: VtuFile
      INTEGER, PARAMETER :: VtuUnit = 58
      TYPE(Variable_t), POINTER :: Var, Solution
      CHARACTER(LEN=512) :: str
      INTEGER :: i,j,k,dofs,Rank,cumn,n,vari,sdofs,IsVector,Offset,PartDim,&
           layers,templayers,damdofs,lay
      CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
           FieldName2, BaseString, OutStr
      CHARACTER :: lf
      LOGICAL :: ScalarsExist, VectorsExist, Found, ParticleMode, ComponentVector, &
           ComplementExists, ThisOnly, Stat, Visited=.FALSE.,endsim = .TRUE.
      LOGICAL :: WriteData, WriteXML, Buffered,calcstress
      INTEGER, POINTER :: Perm(:), Perm2(:), Indexes(:)
      INTEGER, ALLOCATABLE :: ElemInd(:),ElemInd2(:)
      REAL(KIND=dp), POINTER :: Values(:),Values2(:),&
           Values3(:),VecValues(:,:),Basis(:),VecValues3(:,:,:),Val1
      REAL(KIND=dp) :: x,y,z,u,v,w,DetJ,anisoparam,val,Peff,psrten(2,2),&
           psr(2),exx,eyy,exy,ezz,SR(4),e1,e2,angle,fx,fy,vnorm,PD(3),EV(3,3),DDD(2,2),SRm(2,2)
      REAL(KIND=dp) :: EigValues(2),EigenVec(2,2)
      REAL(KIND=dp) :: pSR2(3,3),EeExp,EFExp,Ee,RHS,D(4),exxd1m1,eyyd2m1,ezzd3m1,exyd4,&
           taurhsmultthird,Tau(3,3)
      REAL(KIND=dp) :: ETau(3,3),one=1.0_dp,denom,t1d2m1,t2d1m1,t3od3m1,d4t4
      REAL(KIND=dp) :: half = 0.5_dp,onethird=1.0_dp/3.0_dp,two=2.0_dp

      REAL(KIND=dp) :: difference,tmf

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Element
      TYPE(Variable_t), POINTER :: ParticleVar

      SAVE :: layers, templayers,damdofs,Visited

      IF (.NOT. Visited ) THEN

         damdofs = 4
         layers = Particles % numberofparticlelayers
         templayers = Particles % numberoftemperaturelayers

         Visited = .TRUE.
      END IF

      ! Initialize the auxiliary module for buffered writing
      !--------------------------------------------------------------
      CALL AscBinWriteInit( AsciiOutput, SinglePrec, VtuUnit, NumberOfNodes )

      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

      n = Mesh % MaxElementDOFS
      ALLOCATE( ElemInd(n), ElemInd2(n) )

      ThisOnly = .TRUE.

      ParticleMode = .TRUE. !.NOT. ASSOCIATED( Particles % UVW )

      ! Linefeed character
      !-----------------------------------
      lf = CHAR(10)
      dim = 3

      PartDim = Particles % Dim

      WriteXML = .TRUE.
      WriteData = AsciiOutput
      Params => ListGetSolverParams()
      Buffered = .TRUE.

      ! This is a hack to ensure that the streamed saving will cover the whole file
      !----------------------------------------------------------------------------
      IF(.TRUE.) THEN
         OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'formatted', STATUS='unknown' )
         WRITE( VtuUnit,'(A)') ' '
         CLOSE( VtuUnit )
      END IF

      ! This format works both for ascii and binary output
      !-------------------------------------------------------------------------
      OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'unformatted', ACCESS = 'stream', STATUS='unknown' )

      WRITE( OutStr,'(A)') '<?xml version="1.0"?>'//lf
      CALL AscBinStrWrite( OutStr )

      IF ( LittleEndian() ) THEN
         OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
      ELSE
         OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'//lf
      END IF
      CALL AscBinStrWrite( OutStr )
      WRITE( OutStr,'(A)') '  <UnstructuredGrid>'//lf
      CALL AscBinStrWrite( OutStr )
      WRITE( OutStr,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NumberOfNodes,&
           '" NumberOfCells="',NumberOfNodes,'">'//lf
      CALL AscBinStrWrite( OutStr )


      ! Header information for nodewise data
      !---------------------------------------------------------------------
      ScalarFieldName = ListGetString( Params,'Scalar Field 1',ScalarsExist)
      VectorFieldName = ListGetString( Params,'Vector Field 1',VectorsExist)
      IF( .NOT. ( ScalarsExist .OR. VectorsExist) ) THEN
         CALL Warn('WriteVtuFile','Are there really no scalars or vectors?')
      END IF

      WRITE( OutStr,'(A)') '      <PointData>'//lf
      CALL AscBinStrWrite( OutStr )


      ! do the scalars & vectors
      !--------------------------------- -----------------------------------
100   Offset = 0

      calcstress = .FALSE.


      DO Vari = 1,999
         BaseString = 'Vector Field'
         WRITE(Txt,'(A)') TRIM(BaseString)//' '//TRIM(I2S(Vari))
         FieldName = ListGetString( Params, TRIM(Txt), Found )
         IF(.NOT. Found) EXIT

      END DO


      IF( SaveFields ) THEN

         DO IsVector = 0, 1

            DO Vari = 1, 999

               IF( IsVector == 0 ) THEN
                  BaseString = 'Scalar Field'
               ELSE
                  BaseString = 'Vector Field'
               END IF

               IF (SaveAll) THEN

                  IF (IsVector == 1) THEN
                     SELECT CASE (Vari)
                     CASE(1)
                        FieldName = 'velocity'
                     CASE(2)
                        FieldName = 'gridvelocity'
                     CASE(3)
                        FieldName = 'gradvel'
                     CASE(4)
                        FieldName = 'gradzs'
                     CASE(5)
                        FieldName = 'length'
                     CASE(6)
                        FieldName = 'dav'
                     CASE(7)
                        FieldName = 'nextcoordinate'
                     CASE(8)
                        FieldName = 'strain'
                     CASE(9)
                        FieldName = 'principal_strain_rates'
                     CASE(10)
                        FieldName = 'principal_damage'
                     CASE(11)
                        FieldName = 'pde_one'
                     CASE(12)
                        FieldName = 'pde_two'
                     CASE(13)
                        FieldName = 'pde_three'
                     CASE(14)
                        FieldName = 'psre_one'
                     CASE(15)
                        FieldName = 'psre_two'
                     CASE(16)
                        FieldName = 'psre_three'
                     CASE(17)
                        FieldName = 'f'
                     CASE(18)
                        FieldName = 'gradh'
                     CASE(19)
                        FieldName = 'principal_deviatoric_stresses'
                     CASE(20)
                        FieldName = 'deviatoric_stresses'
                     CASE(21)
                        FieldName = 'sr_in_pdam_dir'
                     CASE(22)
                        EXIT
                     END SELECT
                  ELSE
                     SELECT CASE (Vari)
                     CASE(1)
                        FieldName = 'fp'
                     CASE(2)
                        FieldName = 'h'
                     CASE(3)
                        FieldName = 'gvolume'
                     CASE(4)
                        FieldName = 'pvolume'
                     CASE(5)
                        FieldName = 'mass'
                     CASE(6)
                        FieldName = 'elementindex'
                     CASE(7)
                        FieldName = 'binit'
                     CASE(8)
                        FieldName = 'status'
                     CASE(9)
                        FieldName = 'interpelem'
                     CASE(10)
                        FieldName = 'ef'
                     CASE(11)
                        FieldName = 'gmask'
                     CASE(12)
                        FieldName = 'viscosity'
                     CASE(13)
                        FieldName = 'bedrock'
                     CASE(14)
                        FieldName = 'static'
                     CASE(15)
                        FieldName = 'mb'
                     CASE(16)
                        FieldName = 'no'
                     CASE(17)
                        FieldName = 'origno'
                     CASE(18)
                        FieldName = 'useinterpelem'
                     CASE(19)
                        FieldName = 'damstatus'
                     CASE(20)
                        FieldName = 'tracer'
                     CASE(21)
                        EXIT
                     END SELECT
                  END IF
               ELSE

                  WRITE(Txt,'(A)') TRIM(BaseString)//' '//TRIM(I2S(Vari))
                  FieldName = ListGetString( Params, TRIM(Txt), Found )
                  IF(.NOT. Found) EXIT
               END IF



               ! Get the values assuming particle mode
               !---------------------------------------------------------------------
               IF( ParticleMode ) THEN
                  IF( IsVector == 1) THEN
                     IF( FieldName == 'velocity' ) THEN
                        dofs = 2
                     ELSE IF( FieldName == 'gridvelocity' ) THEN
                        dofs = 2
                     ELSE IF( FieldName == 'gradvel') THEN
                        dofs = 4
                     ELSE IF( FieldName == 'gradzs') THEN
                        dofs = 2
                     ELSE IF (FieldName == 'gradh') THEN
                        dofs = 2
                     ELSE IF( FieldName == 'length') THEN
                        dofs = 2
                     ELSE IF (FieldName == 'origlength') THEN
                        dofs = 2
                     ELSE IF( FieldName == 'dav') THEN
                        dofs = damdofs
                     ELSE IF (FieldName == 'damage') THEN
                        !returns D(z) if damdofs == 1, otherwise, returns DI(z)
                        dofs = layers
                     ELSE IF (FieldName == 'damageii') THEN
                        !returns DII(z)
                        dofs = layers
                     ELSE IF (FieldName == 'damageiii') THEN
                        !returns DII(z)
                        dofs = layers
                     ELSE IF (FieldName == 'xdamage') THEN
                        !returns D(z) if damdofs == 1, otherwise, returns Dxx(z)
                        dofs = layers
                     ELSE IF (FieldName == 'ydamage') THEN
                        !returns D(z) if damdofs == 1, otherwise, returns Dxx(z)
                        dofs = layers
                     ELSE IF (FieldName == 'zdamage') THEN
                        !returns D(z) if damdofs == 1, otherwise, returns Dxx(z)
                        dofs = layers
                     ELSE IF (FieldName == 'xydamage') THEN
                        !returns D(z) if damdofs == 1, otherwise, returns Dxx(z)
                        dofs = layers
                     ELSE IF (FieldName == 'dd') THEN
                        !returns dD(z) if damdofs == 1, otherwise, returns dDxx(z)
                        dofs = layers
                     ELSE IF (FieldName == 'bz') THEN
                        dofs = layers
                     ELSE IF (FieldName == 'nextcoordinate') THEN
                        dofs = 2
                     ELSE IF (FieldName == 'xpic') THEN
                        dofs = 6
                     ELSE IF (FieldName == 'strain') THEN
                        dofs = 2
                     ELSE IF (FieldName == 'principal_strain_rates') THEN
                        dofs = 2
                     ELSE IF (FieldName == 'principal_damage') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'sr_in_pdam_dir') THEN
                        dofs = 2
                     ELSE IF (FieldName == 'pde_one') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'pde_two') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'pde_three' .OR. FieldName == 'eff_pds') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'pdse_two' .OR. FieldName == 'eff_pdse_two') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'psre_one') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'psre_two') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'psre_three') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'pde_two_30d') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'pde_two_20d') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'pde_two_10d') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'f') THEN
                        dofs = 4
                     ELSE IF (FieldName == 'principal_deviatoric_stresses') THEN
                        dofs = 3
                     ELSE IF (FieldName == 'deviatoric_stresses') THEN
                        dofs = 4
                     ELSE
                        WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                        CALL Warn('WriteVtuXMLFile', Txt)
                        CYCLE
                     END IF
                  ELSE
                     dofs = 1

                     IF (FieldName == 'particle time') THEN
                        dofs = 1
                     ELSEIF (FieldName == 'particle dt') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'fp' ) THEN
                        dofs = 1
                     ELSE IF( FieldName == 'h') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'gvolume') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'pvolume') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'mass') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'elementindex') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'binit') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'status') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'interpelem') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'ef') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'gmask') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'viscosity') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'bedrock') THEN
                        dofs = 1
                     ELSE IF( FieldName == 'static') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'mb') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'no') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'origno') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'useinterpelem') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'damstatus') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'isodav') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'dbassis') THEN
                        dofs = 1
                     ELSE IF (FieldName == 'tracer') THEN
                        dofs = 1
                     ELSE
                        WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                        CALL Warn('WriteVtuXMLFile', Txt)
                        CYCLE
                     END IF
                  END IF
               END IF

               ! Finally save the field values for scalars
               !---------------------------------------------------------------------
               IF( dofs == 1 ) THEN
                  sdofs = 1
               ELSE
                  sdofs = MAX(dofs,3)
               END IF

               sdofs = dofs


               IF( WriteXML ) THEN
                  WRITE( OutStr,'(A,I0,A)') '        <DataArray type="Float',PrecBits,'" Name="'//TRIM(FieldName)
                  CALL AscBinStrWrite( OutStr )

                  WRITE( OutStr,'(A,I0,A)') '" NumberOfComponents="',sdofs,'"'
                  CALL AscBinStrWrite( OutStr )

                  IF( AsciiOutput ) THEN
                     WRITE( OutStr,'(A)') ' format="ascii">'//lf
                     CALL AscBinStrWrite( OutStr )
                  ELSE
                     WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
                     CALL AscBinStrWrite( OutStr )
                  END IF
               END IF


               IF( BinaryOutput ) THEN
                  k = NumberOfNodes * PrecSize * sdofs
                  Offset = Offset + IntSize + k
               END IF

               j = 0
               IF( WriteData ) THEN
                  IF( BinaryOutput ) WRITE( VtuUnit ) k

                  DO i = 1, Particles % NumberOfParticles

                     IF ( Particles % Status(i) > MaxSaveStatus .OR. &
                          Particles % Status(i) < MinSaveStatus )  CYCLE

                     IF (Particles % Coordinate(i,1) < minx) CYCLE
                     IF (Particles % Coordinate(i,1) > maxx) CYCLE

                     IF (FloatingOnly) THEN
                        IF (Particles % Gmask(i)<0) CYCLE
                     END IF

                     j = j + 1

                     LocalVal = 0.0_dp



                     IF( ParticleMode ) THEN
                        IF( IsVector == 1) THEN

                           IF( FieldName == 'velocity' ) THEN
                              LocalVal(1:sdofs) = Particles % Velocity(i,1:sdofs)
                           ELSE IF( FieldName == 'gridvelocity' ) THEN
                              LocalVal(1:sdofs) = Particles % GridVelocity(i,1:sdofs)
                           ELSE IF( FieldName == 'gradvel') THEN
                              LocalVal(1:sdofs) = Particles % GradVel(i,1:sdofs)
                           ELSE IF( FieldName == 'gradh') THEN
                              LocalVal(1:sdofs) = Particles % GradH(i,1:sdofs)
                           ELSE IF( FieldName == 'gradzs') THEN
                              LocalVal(1:sdofs) = Particles % GradZs(i,1:sdofs)
                           ELSE IF( FieldName == 'length') THEN
                              LocalVal(1:sdofs) = Particles % Length(i,1:sdofs)
                           ELSE IF (FieldName == 'origlength') THEN
                              LocalVal(1:sdofs) = Particles % OrigLength(i,1:sdofs)

                           ELSE IF( FieldName == 'dav') THEN
                              LocalVal(1:sdofs) = Particles % Dav(i,1:sdofs)

                              DO  lay = 1,sdofs
                                 IF (LocalVal(lay) .NE. LocalVal(lay)) LocalVal(lay) = -1
                              END DO
                           ELSE IF (FieldName == 'damage') THEN

                              IF (Particles % gamma>0.0_dp) THEN
                                 DO lay=1,sdofs
                                    PD=0.0_dp
                                    CALL PrincipalDamage(Particles%Damage(i,lay,:),PD)
                                    LocalVal(lay) = PD(1)
                                 END DO

                              ELSE

                                 LocalVal(1:sdofs) = Particles % Damage(i,1:sdofs,1)
                              END IF
                           ELSE IF (FieldName == 'damageii') THEN
                              DO lay=1,sdofs
                                 PD=0.0_dp
                                 CALL PrincipalDamage(Particles%Damage(i,lay,:),PD)
                                 LocalVal(lay) = PD(2)
                              END DO
                           ELSE IF (FieldName == 'damageiii') THEN
                              DO lay=1,sdofs
                                 PD=0.0_dp
                                 CALL PrincipalDamage(Particles%Damage(i,lay,:),PD)
                                 LocalVal(lay) = PD(3)

                              END DO
                           ELSE IF (FieldName == 'xdamage') THEN

                              LocalVal(1:sdofs) = Particles % Damage(i,1:sdofs,1)
                           ELSE IF (FieldName == 'ydamage') THEN

                              LocalVal(1:sdofs) = Particles % Damage(i,1:sdofs,2)
                           ELSE IF (FieldName == 'zdamage') THEN

                              LocalVal(1:sdofs) = Particles % Damage(i,1:sdofs,3)
                           ELSE IF (FieldName == 'xydamage') THEN

                              LocalVal(1:sdofs) = Particles % Damage(i,1:sdofs,4)

                           ELSE IF (FieldName == 'dd') THEN
                              LocalVal(1:sdofs) = Particles % dD(i,1:sdofs,1)
                           ELSE IF (FieldName == 'bz') THEN
                              LocalVal(1:sdofs) = Particles % bz(i,1:sdofs)
                           ELSE IF (FieldName == 'nextcoordinate') THEN
                              LocalVal(1:sdofs) = Particles % NextCoordinate(i,1:sdofs)
                           ELSE IF (FieldName == 'xpic') THEN
                              LocalVal(1:sdofs) = Particles % xpic(i,1:sdofs)
                           ELSE IF (FieldName == 'strain') THEN
                              IF (Particles % trackstrain) THEN
                                 LocalVal(1:sdofs) = Particles % strain(i,1:sdofs)
                              ELSE
                                 LocalVal(1:sdofs) = -9999.0_dp
                              END IF

                           ELSE IF (FieldName == 'principal_strain_rates') THEN


                              exx = Particles % GradVel(i,1)
                              eyy = Particles % GradVel(i,2)
                              exy = 0.5_dp*(Particles % GradVel(i,3)  + Particles % GradVel(i,4) )

                              psrten(1,1) = exx
                              psrten(2,2) = eyy
                              psrten(1,2) = exy
                              psrten(2,1) = exy

                              CALL Eigen2DSym_TryGenFirst(psrten,EigValues,EigenVec)

                              e1 = EigValues(2)
                              e2 = EigValues(1)

                              LocalVal(1) = e1
                              LocalVal(2) = e2
                           ELSE IF (FieldName == 'principal_deviatoric_stresses') THEN


                              pSR2(1,1) = Particles % GradVel(i,1)
                              pSR2(2,1) = 0.5_dp*(Particles % Gradvel(i,3) + Particles % GradVel(i,4))
                              pSR2(1,2) = pSR2(2,1)
                              pSR2(2,2) = Particles % GradVel(i,2)
                              pSR2(3,3) = -pSR2(1,1)-pSR2(2,2)

                              EeExp =  (1.0_dp-3.0_dp)/(2.0_dp * 3.0_dp)
                              EFexp = -1.0_dp/3.0_dp


                              Ee = 0.5_dp*(pSR2(1,1)*pSR2(1,1) + pSR2(2,2)*pSR2(2,2) + &
                                   pSR2(3,3)*pSR2(3,3)) + pSR2(1,2)*pSR2(1,2)

                              IF (Ee<(Particles % criticalshearrate*Particles % criticalshearrate)) THEN
                                 Ee = Particles % criticalshearrate*Particles % criticalshearrate
                              END IF


                              RHS = Particles % Binit(i) * (Ee**EeExp) * (Particles % EF(i)**EFexp )

                              D(1) = Particles % Dav(i,1)
                              D(2) = Particles % Dav(i,2)
                              D(3) = Particles % Dav(i,3)
                              D(4) = Particles % Dav(i,4)

                              exxd1m1 = psr2(1,1)*(D(1)-1.0_dp)
                              eyyd2m1 = psr2(2,2)*(D(2)-1.0_dp)
                              ezzd3m1 = psr2(3,3)*(D(3)-1.0_dp)
                              exyd4 = psr2(1,2)*D(4)
                              taurhsmultthird = RHS/3.0_dp


                              Tau(1,1) = taurhsmultthird * (-2.0_dp*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
                              Tau(2,2) = taurhsmultthird * (exxd1m1-2.0_dp*eyyd2m1+ezzd3m1-exyd4)
                              Tau(3,3) = taurhsmultthird * (exxd1m1+eyyd2m1-2.0_dp*ezzd3m1+2.0_dp*exyd4)
                              Tau(1,2) = -RHS * 0.5_dp * (psr2(1,2)*(D(1)+D(2)-2.0_dp) + D(4)*(psr2(1,1)+psr2(2,2)))
                              Tau(2,1) = Tau(1,2)


                              CALL Eigen2DSym_TryGenFirst(Tau(1:2,1:2),EigValues,EigenVec)

                              e1 = EigValues(2)
                              e2 = EigValues(1)


                              LocalVal(1) = e1
                              LocalVal(2) = e2
                              LocalVal(3) = Tau(3,3)



                           ELSE IF (FieldName == 'deviatoric_stresses') THEN

                              pSR2(1,1) = Particles % GradVel(i,1)
                              pSR2(2,1) = 0.5_dp*(Particles % Gradvel(i,3) + Particles % GradVel(i,4))
                              pSR2(1,2) = pSR2(2,1)
                              pSR2(2,2) = Particles % GradVel(i,2)
                              pSR2(3,3) = -pSR2(1,1)-pSR2(2,2)

                              EeExp =  (1.0_dp-3.0_dp)/(2.0_dp * 3.0_dp)
                              EFexp = -1.0_dp/3.0_dp


                              Ee = 0.5_dp*(pSR2(1,1)*pSR2(1,1) + pSR2(2,2)*pSR2(2,2) + &
                                   pSR2(3,3)*pSR2(3,3)) + pSR2(1,2)*pSR2(1,2)

                              IF (Ee<(Particles % criticalshearrate*Particles % criticalshearrate)) THEN
                                 Ee = Particles % criticalshearrate*Particles % criticalshearrate
                              END IF


                              RHS = Particles % Binit(i) * (Ee**EeExp) * (Particles % EF(i)**EFexp )

                              D(1) = Particles % Dav(i,1)
                              D(2) = Particles % Dav(i,2)
                              D(3) = Particles % Dav(i,3)
                              D(4) = Particles % Dav(i,4)

                              exxd1m1 = psr2(1,1)*(D(1)-1.0_dp)
                              eyyd2m1 = psr2(2,2)*(D(2)-1.0_dp)
                              ezzd3m1 = psr2(3,3)*(D(3)-1.0_dp)
                              exyd4 = psr2(1,2)*D(4)
                              taurhsmultthird = RHS/3.0_dp


                              Tau(1,1) = taurhsmultthird * (-2.0_dp*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
                              Tau(2,2) = taurhsmultthird * (exxd1m1-2.0_dp*eyyd2m1+ezzd3m1-exyd4)
                              Tau(3,3) = taurhsmultthird * (exxd1m1+eyyd2m1-2.0_dp*ezzd3m1+2.0_dp*exyd4)
                              Tau(1,2) = -RHS * 0.5_dp * (psr2(1,2)*(D(1)+D(2)-2.0_dp) + D(4)*(psr2(1,1)+psr2(2,2)))
                              Tau(2,1) = Tau(1,2)

                              LocalVal(1) = Tau(1,1)
                              LocalVal(2) = Tau(2,2)
                              LocalVal(3) = Tau(3,3)
                              LocalVal(4) = Tau(1,2)



                           ELSE IF (FieldName == 'eff_pdse_two' .OR. FieldName == 'eff_pds') THEN


                              pSR2(1,1) = Particles % GradVel(i,1)
                              pSR2(2,1) = 0.5_dp*(Particles % Gradvel(i,3) + Particles % GradVel(i,4))
                              pSR2(1,2) = pSR2(2,1)
                              pSR2(2,2) = Particles % GradVel(i,2)
                              pSR2(3,3) = -pSR2(1,1)-pSR2(2,2)

                              EeExp =  (1.0_dp-3.0_dp)/(2.0_dp * 3.0_dp)
                              EFexp = -1.0_dp/3.0_dp


                              Ee = 0.5_dp*(pSR2(1,1)*pSR2(1,1) + pSR2(2,2)*pSR2(2,2) + &
                                   pSR2(3,3)*pSR2(3,3)) + pSR2(1,2)*pSR2(1,2)

                              IF (Ee<(Particles % criticalshearrate*Particles % criticalshearrate)) THEN
                                 Ee = Particles % criticalshearrate*Particles % criticalshearrate
                              END IF


                              RHS = Particles % Binit(i) * (Ee**EeExp) * (Particles % EF(i)**EFexp )

                              D(1) = Particles % Dav(i,1)
                              D(2) = Particles % Dav(i,2)
                              D(3) = Particles % Dav(i,3)
                              D(4) = Particles % Dav(i,4)

                              exxd1m1 = psr2(1,1)*(D(1)-1.0_dp)
                              eyyd2m1 = psr2(2,2)*(D(2)-1.0_dp)
                              ezzd3m1 = psr2(3,3)*(D(3)-1.0_dp)
                              exyd4 = psr2(1,2)*D(4)
                              taurhsmultthird = RHS/3.0_dp


                              Tau(1,1) = taurhsmultthird * (-2.0_dp*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
                              Tau(2,2) = taurhsmultthird * (exxd1m1-2.0_dp*eyyd2m1+ezzd3m1-exyd4)
                              Tau(3,3) = taurhsmultthird * (exxd1m1+eyyd2m1-2.0_dp*ezzd3m1+2.0_dp*exyd4)
                              Tau(1,2) = -RHS * 0.5_dp * (psr2(1,2)*(D(1)+D(2)-2.0_dp) + &
                                   D(4)*(psr2(1,1)+psr2(2,2)))
                              Tau(2,1) = Tau(1,2)


                              !effective principal deviatoric stress
                              denom = one/( D(4)*D(4) + D(1) + D(2) -D(1)*D(2) - one)
                              t1d2m1 = Tau(1,1)*(D(2)-one)*denom
                              t2d1m1 = Tau(2,2)*(D(1)-one)*denom
                              t3od3m1 = Tau(3,3)/(D(3)-one)
                              d4t4 = Tau(1,2)*D(4)*denom

                              Etau(1,1) = onethird*(two*t1d2m1 - t2d1m1 +t3od3m1 -d4t4)
                              Etau(2,2) = onethird*(-t1d2m1 + two*t2d1m1 +t3od3m1 -d4t4)
                              Etau(3,3) = onethird*(-t1d2m1 - t2d1m1 - two*t3od3m1 + two*d4t4)
                              Etau(1,2) = half*denom*(tau(1,2)*(D(1)+D(2)-two) - D(4)*(Tau(1,1)+Tau(2,2)))
                              Etau(2,1) = Etau(1,2)


                              SR(1) = ETau(1,1)
                              SR(2) = ETau(2,2)
                              SR(3) = ETau(3,3)
                              SR(4) = ETau(1,2)



                              IF (FieldName == 'eff_pdse_two' ) THEN

                                 CALL PrincipalEigenVec(SR,EV)
                                 LocalVal(1) = EV(1,2)
                                 LocalVal(2) = EV(2,2)
                                 LocalVal(3) = EV(3,2)

                              END IF


                              IF (FieldName == 'eff_pds') THEN
                                 CALL PrincipalDamage(SR,PD)

                                 LocalVal(1) = PD(1)
                                 LocalVal(2) = PD(2)
                                 LocalVal(3) = PD(3)

                              END IF

                           ELSE IF (FieldName == 'pdse_two') THEN


                              pSR2(1,1) = Particles % GradVel(i,1)
                              pSR2(2,1) = 0.5_dp*(Particles % Gradvel(i,3) + Particles % GradVel(i,4))
                              pSR2(1,2) = pSR2(2,1)
                              pSR2(2,2) = Particles % GradVel(i,2)
                              pSR2(3,3) = -pSR2(1,1)-pSR2(2,2)

                              EeExp =  (1.0_dp-3.0_dp)/(2.0_dp * 3.0_dp)
                              EFexp = -1.0_dp/3.0_dp


                              Ee = 0.5_dp*(pSR2(1,1)*pSR2(1,1) + pSR2(2,2)*pSR2(2,2) + &
                                   pSR2(3,3)*pSR2(3,3)) + pSR2(1,2)*pSR2(1,2)

                              IF (Ee<(Particles % criticalshearrate*Particles % criticalshearrate)) THEN
                                 Ee = Particles % criticalshearrate*Particles % criticalshearrate
                              END IF


                              RHS = Particles % Binit(i) * (Ee**EeExp) * (Particles % EF(i)**EFexp )

                              D(1) = Particles % Dav(i,1)
                              D(2) = Particles % Dav(i,2)
                              D(3) = Particles % Dav(i,3)
                              D(4) = Particles % Dav(i,4)

                              exxd1m1 = psr2(1,1)*(D(1)-1.0_dp)
                              eyyd2m1 = psr2(2,2)*(D(2)-1.0_dp)
                              ezzd3m1 = psr2(3,3)*(D(3)-1.0_dp)
                              exyd4 = psr2(1,2)*D(4)
                              taurhsmultthird = RHS/3.0_dp


                              Tau(1,1) = taurhsmultthird * (-2.0_dp*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
                              Tau(2,2) = taurhsmultthird * (exxd1m1-2.0_dp*eyyd2m1+ezzd3m1-exyd4)
                              Tau(3,3) = taurhsmultthird * (exxd1m1+eyyd2m1-2.0_dp*ezzd3m1+2.0_dp*exyd4)
                              Tau(1,2) = -RHS * 0.5_dp * (psr2(1,2)*(D(1)+D(2)-2.0_dp) + D(4)*(psr2(1,1)+psr2(2,2)))
                              Tau(2,1) = Tau(1,2)


                              SR(1) = Tau(1,1)
                              SR(2) = Tau(2,2)
                              SR(3) = Tau(3,3)
                              SR(4) = Tau(1,2)

                              CALL PrincipalEigenVec(SR,EV)

                              LocalVal(1) = EV(1,2)
                              LocalVal(2) = EV(2,2)
                              LocalVal(3) = EV(3,2)

                           ELSE IF (FieldName == 'principal_damage') THEN

                              PD=0.0_dp

                              IF (Particles % gamma>0.0_dp .or. Particles%nodzz ) THEN
                                 CALL PrincipalDamage(Particles%Dav(i,:),PD)

                                 LocalVal(1) = PD(1)
                                 LocalVal(2) = PD(2)
                                 LocalVal(3) = PD(3)

                                 IF (Particles % useriftdmax) THEN
                                    IF (Particles % DamStatus(i)==1) THEN
                                       LocalVal(1:3) = Particles%Dav(i,1)
                                    END IF
                                 END IF

                              ELSE
                                 LocalVal(1:3) = Particles%Dav(i,1)
                              END IF


                              IF (Particles % dmaxII_dom .AND. &
                                   Particles % DamStatus(i)==1) THEN

                                 LocalVal(2) = PD(1)
                                 LocalVal(1) = PD(2)

                              END IF


                              IF (LocalVal(1) .NE. LocalVal(1)) LocalVal(1) = -1
                              IF (LocalVal(2) .NE. LocalVal(2)) LocalVal(2) = -1
                              IF (LocalVal(3) .NE. LocalVal(3)) LocalVal(3) = -1


                           ELSE IF (FieldName == 'pde_one') THEN

                              EV=0.0_dp

                              CALL PrincipalEigenVec(Particles%Dav(i,:),EV)

                              LocalVal(1) = EV(1,1)
                              LocalVal(2) = EV(2,1)
                              LocalVal(3) = EV(3,1)


                              IF (Particles % dmaxII_dom .AND. &
                                   Particles % DamStatus(i)==1) THEN

                                 LocalVal(1) = EV(1,2)
                                 LocalVal(2) = EV(2,2)
                                 LocalVal(3) = EV(3,2)
                              END IF


                           ELSE IF (FieldName == 'pde_two') THEN

                              EV=0.0_dp

                              CALL PrincipalEigenVec(Particles%Dav(i,:),EV)

                              LocalVal(1) = EV(1,2)
                              LocalVal(2) = EV(2,2)
                              LocalVal(3) = EV(3,2)


                              IF (Particles % dmaxII_dom .AND. &
                                   Particles % DamStatus(i)==1) THEN

                                 LocalVal(1) = EV(1,1)
                                 LocalVal(2) = EV(2,1)
                                 LocalVal(3) = EV(3,1)
                              END IF



                           ELSE IF (FieldName == 'pde_three') THEN

                              EV=0.0_dp

                              CALL PrincipalEigenVec(Particles%Dav(i,:),EV)


                              LocalVal(1) = EV(1,3)
                              LocalVal(2) = EV(2,3)
                              LocalVal(3) = EV(3,3)

                           ELSE IF (FieldName == 'sr_in_pdam_dir') THEN


                              EV=0.0_dp

                              CALL PrincipalEigenVec(Particles%Dav(i,:),EV)


                              EV(1:3,3) = EV(1:3,1)
                              EV(1:3,1) = EV(1:3,2)
                              EV(1:3,2) = EV(1:3,3)

                              exx = Particles % GradVel(i,1)
                              eyy = Particles % GradVel(i,2)
                              exy = 0.5_dp*(Particles % GradVel(i,3)  + Particles % GradVel(i,4) )

                              psrten(1,1) = exx
                              psrten(2,2) = eyy
                              psrten(1,2) = exy
                              psrten(2,1) = exy

                              psrten = MATMUL(MATMUL(TRANSPOSE(EV(1:2,1:2)),psrten),EV(1:2,1:2))

                              LocalVal(1) = psrten(2,2)
                              LocalVal(2) = psrten(1,1)

                           ELSE IF (FieldName == 'psre_one') THEN

                              EV=0.0_dp

                              SR(1) = Particles % GradVel(i,1)
                              SR(2) = Particles % GradVel(i,2)
                              SR(3) = -SR(1) - SR(2)
                              SR(4) = 0.5_dp*(Particles % GradVel(i,3)  + Particles % GradVel(i,4) )

                              CALL PrincipalEigenVec(SR,EV)


                              LocalVal(1) = EV(1,1)
                              LocalVal(2) = EV(2,1)
                              LocalVal(3) = EV(3,1)

                           ELSE IF (FieldName == 'psre_two') THEN

                              EV=0.0_dp

                              SR(1) = Particles % GradVel(i,1)
                              SR(2) = Particles % GradVel(i,2)
                              SR(3) = -SR(1) - SR(2)
                              SR(4) = 0.5_dp*(Particles % GradVel(i,3)  + Particles % GradVel(i,4) )


                              CALL PrincipalEigenVec(SR,EV)


                              LocalVal(1) = EV(1,2)
                              LocalVal(2) = EV(2,2)
                              LocalVal(3) = EV(3,2)

                           ELSE IF (FieldName == 'psre_three') THEN

                              EV=0.0_dp

                              SR(1) = Particles % GradVel(i,1)
                              SR(2) = Particles % GradVel(i,2)
                              SR(3) = -SR(1) - SR(2)
                              SR(4) = 0.5_dp*(Particles % GradVel(i,3)  + Particles % GradVel(i,4) )


                              CALL PrincipalEigenVec(SR,EV)


                              LocalVal(1) = EV(1,3)
                              LocalVal(2) = EV(2,3)
                              LocalVal(3) = EV(3,3)


                           ELSE IF (FieldName == 'pde_two_30d') THEN

                              EV=0.0_dp

                              CALL PrincipalEigenVec(Particles%Damage(i,30,:),EV)


                              LocalVal(1) = EV(1,2)
                              LocalVal(2) = EV(2,2)
                              LocalVal(3) = EV(3,2)


                           ELSEIF (FieldName == 'pde_two_20d') THEN

                              EV=0.0_dp

                              CALL PrincipalEigenVec(Particles%Damage(i,20,:),EV)

                              !old
                              LocalVal(1) = EV(1,2)
                              LocalVal(2) = EV(2,2)
                              LocalVal(3) = EV(3,2)

                           ELSEIF (FieldName == 'pde_two_10d') THEN

                              EV=0.0_dp

                              CALL PrincipalEigenVec(Particles%Damage(i,10,:),EV)

                              !old
                              LocalVal(1) = EV(1,2)
                              LocalVal(2) = EV(2,2)
                              LocalVal(3) = EV(3,2)

                           ELSEIF (FieldName == 'f') THEN
                              LocalVal(1:4) = Particles % F(i,1:4)
                           END IF

                        ELSE
                           sdofs = 1

                           IF (FieldName == 'particle time') THEN
                              LocalVal(1) = Particles % time
                           ELSEIF (FieldName == 'particle dt') THEN
                              LocalVal(1) = Particles % dtime
                           ELSE IF( FieldName == 'fp' ) THEN
                              LocalVal(1) = MAX(Particles % FP(i),0.0_dp)
                           ELSE IF( FieldName == 'h') THEN
                              IF (Particles % H(i) .NE. Particles % H(i)) Particles % H(i) = -999.0_dp
                              LocalVal(1) = Particles % H(i)
                           ELSE IF( FieldName == 'gvolume') THEN
                              LocalVal(1) = Particles % Gvolume(i)
                           ELSE IF( FieldName == 'pvolume') THEN
                              LocalVal(1) = Particles % pVolume(i)
                           ELSE IF (FieldName == 'mass') THEN
                              LocalVal(1) = Particles % mass(i)
                           ELSE IF( FieldName == 'elementindex') THEN
                              LocalVal(1) = Particles % ElementIndex(i)
                           ELSE IF( FieldName == 'binit') THEN
                              LocalVal(1) = Particles % binit(i)
                           ELSE IF( FieldName == 'status') THEN
                              LocalVal(1) = Particles % Status(i)
                           ELSE IF( FieldName == 'interpelem') THEN
                              LocalVal(1) = Particles % InterpElem(i)
                           ELSE IF( FieldName == 'ef') THEN
                              LocalVal(1) = Particles % EF(i)
                           ELSE IF( FieldName == 'gmask') THEN
                              LocalVal(1) = Particles % Gmask(i)
                           ELSE IF( FieldName == 'viscosity') THEN
                              LocalVal(1) = (1.0_dp - Particles % Dav(i,1)) * &
                                   Particles % Binit(i) * &
                                   ((Particles % EF(i))**(-1.0_dp/3.0_dp))
                           ELSE IF (FieldName == 'bedrock') THEN
                              LocalVal(1) = Particles % Bedrock(i)
                           ELSE IF( FieldName == 'static') THEN
                              IF (Particles % Static(i)) THEN
                                 LocalVal(1) = 1.0_dp
                              ELSE
                                 LocalVal(1) = -1.0_dp
                              END IF
                           ELSE IF (FieldName == 'mb') THEN
                              LocalVal(1) = Particles % MB(i)
                           ELSE IF (FieldName == 'no') THEN
                              LocalVal(1) = i
                           ELSE IF (FieldName == 'origno') THEN
                              LocalVal(1) = Particles % origno(i)
                           ELSE IF (FieldName == 'useinterpelem') THEN
                              IF (Particles % UseInterpElem(i)) THEN
                                 LocalVal(1) = 1.0_dp
                              ELSE
                                 LocalVal(1) = -1.0_dp
                              END IF
                           ELSE IF (FieldName == 'damstatus') THEN

                              LocalVal(1) = Particles % damstatus(i)
                           ELSE IF (FieldName == 'isodav') THEN
                              LocalVal(1) = Particles % dav(i,1)
                           ELSE IF (FieldName == 'dbassis') THEN
                              LocalVal(1) = Particles % dbassis(i)

                           ELSE IF (FieldName == 'tracer') THEN
                              IF (Particles % usetracer) THEN
                                 LocalVal(1) = Particles % tracer(i)
                              ELSE
                                 LocalVal(1) = -1.0_dp
                              END IF
                           END IF
                        END IF
                     END IF

                     DO k=1,sdofs
                        IF (LocalVal(k) .NE. LocalVal(k)) THEN
                           LocalVal(k) = -9999.0_dp
                        END IF

                        CALL AscBinRealWrite( LocalVal(k) )
                     END DO

                     IF( j == NumberOfNodes ) EXIT
                  END DO

                  CALL AscBinRealWrite( 0.0_dp, .TRUE.)
               END IF

               IF( AsciiOutput ) THEN
                  WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
                  CALL AscBinStrWrite( OutStr )
               END IF
            END DO

         END DO
      END IF

      IF( WriteXML ) THEN
         WRITE( OutStr,'(A)') '      </PointData>'//lf
         CALL AscBinStrWrite( OutStr )

         WRITE( OutStr,'(A)') '      <CellData>'//lf
         CALL AscBinStrWrite( OutStr )
         WRITE( OutStr,'(A)') '      </CellData>'//lf
         CALL AscBinStrWrite( OutStr )
      END IF


      ! Coordinates of each point
      !-------------------------------------
      IF( WriteXML ) THEN
         WRITE( OutStr,'(A)') '      <Points>'//lf
         CALL AscBinStrWrite( OutStr )

         WRITE( OutStr,'(A,I0,A,I0,A)') '        <DataArray type="Float',PrecBits,'" NumberOfComponents="',dim,'"'
         CALL AscBinStrWrite( OutStr )

         IF( AsciiOutput ) THEN
            WRITE( OutStr,'(A)') ' format="ascii">'//lf
            CALL AscBinStrWrite( OutStr )
         ELSE
            WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
            CALL AscBinStrWrite( OutStr )
         END IF
      END IF

      IF( BinaryOutput ) THEN
         k = dim * NumberOfNodes * PrecSize
         Offset = Offset + IntSize + k
      END IF

      IF( WriteData ) THEN
         IF( BinaryOutput ) WRITE( VtuUnit ) k

         LocalVal = 0.0_dp
         j = 0
         DO i = 1, Particles % NumberOfParticles

            IF ( Particles % Status(i) > MaxSaveStatus .OR. &
                 Particles % Status(i) < MinSaveStatus )  CYCLE

            IF (Particles % Coordinate(i,1) < minx) CYCLE
            IF (Particles % Coordinate(i,1) > maxx) CYCLE

            IF (FloatingOnly) THEN
               IF (Particles % Gmask(i)<0) CYCLE
            END IF
            j = j + 1

            IF( ParticleMode ) THEN
               DO k=1,PartDim
                  !DO k=1,dim
                  LocalVal(k) = Particles % Coordinate(i,k)
               END DO

            END IF

            CALL AscBinRealWrite( LocalVal(1) )
            CALL AscBinRealWrite( LocalVal(2) )
            CALL AscBinRealWrite( LocalVal(3) )

            IF( j == NumberOfNodes ) EXIT
         END DO

         CALL AscBinRealWrite( 0.0_dp, .TRUE.)
      END IF

      IF( AsciiOutput ) THEN
         WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
         CALL AscBinStrWrite( OutStr )
      END IF
      IF( WriteXML ) THEN
         WRITE( OutStr,'(A)') '      </Points>'//lf
         CALL AscBinStrWrite( OutStr )
      END IF


      ! Write out the mesh
      !-------------------------------------
      IF( WriteXML ) THEN
         WRITE( OutStr,'(A)') '      <Cells>'//lf
         CALL AscBinStrWrite( OutStr )

         WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="connectivity"'
         CALL AscBinStrWrite( OutStr )

         IF( AsciiOutput ) THEN
            WRITE( OutStr,'(A)') ' format="ascii">'//lf
            CALL AscBinStrWrite( OutStr )
         ELSE
            WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
            CALL AscBinStrWrite( OutStr )
         END IF
      END IF

      IF( BinaryOutput ) THEN
         ! The offset needs to be summed over all nodes
         k = NumberOfNodes * IntSize
         Offset = Offset + k + IntSize
      END IF

      IF( WriteData ) THEN
         IF( BinaryOutput ) WRITE( VtuUnit ) k
         DO i = 1, NumberOfNodes
            CALL AscBinIntegerWrite( i - 1)
         END DO
         CALL AscBinIntegerWrite( 0, .TRUE. )
      END IF

      IF( AsciiOutput ) THEN
         WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
         CALL AscBinStrWrite( OutStr )
      END IF

      ! Offsets for element indexes
      !-------------------------------------------------------------------
      IF( WriteXML ) THEN
         WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="offsets"'
         CALL AscBinStrWrite( OutStr )

         IF( AsciiOutput ) THEN
            WRITE( OutStr,'(A)') ' format="ascii">'//lf
            CALL AscBinStrWrite( OutStr )
         ELSE
            WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
            CALL AscBinStrWrite( OutStr )
         END IF
      END IF

      IF( BinaryOutput ) THEN
         k = NumberOfNodes * IntSize
         Offset = Offset + IntSize + k
      END IF

      IF( WriteData ) THEN
         IF( BinaryOutput ) WRITE( VtuUnit ) k
         DO i = 1, NumberOfNodes
            CALL AscBinIntegerWrite( i )
         END DO
         CALL AscBinIntegerWrite( 0, .TRUE.)
      END IF

      IF( AsciiOutput ) THEN
         WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
         CALL AscBinStrWrite( OutStr )
      END IF
      IF( WriteXML ) THEN
         WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="types"'
         CALL AscBinStrWrite( OutStr )

         IF( AsciiOutput ) THEN
            WRITE( OutStr,'(A)') ' FORMAT="ascii">'//lf
            CALL AscBinStrWrite( OutStr )
         ELSE
            WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
            CALL AscBinStrWrite( OutStr )
         END IF
      END IF

      IF( BinaryOutput ) THEN
         k = NumberOfNodes * IntSize
         Offset = Offset + IntSize + k
      END IF

      IF( WriteData ) THEN
         IF( BinaryOutput ) WRITE( VtuUnit ) k
         ! elementtype is fixed to single nodes (==1)
         DO i = 1, NumberOfNodes
            CALL AscBinIntegerWrite( 1 )
         END DO
         CALL AscBinIntegerWrite( 0, .TRUE. )
      END IF

      IF( AsciiOutput ) THEN
         WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
         CALL AscBinStrWrite( OutStr )
      END IF
      IF( WriteXml ) THEN
         WRITE( OutStr,'(A)') '      </Cells>'//lf
         CALL AscBinStrWrite( OutStr )
         WRITE( OutStr,'(A)') '    </Piece>'//lf
         CALL AscBinStrWrite( OutStr )
         WRITE( OutStr,'(A)') '  </UnstructuredGrid>'//lf
         CALL AscBinStrWrite( OutStr )
      END IF


      IF( BinaryOutput ) THEN
         IF( WriteXML ) THEN
            WRITE( OutStr,'(A)') '<AppendedData encoding="raw">'//lf
            CALL AscBinStrWrite( OutStr )
            WRITE( VtuUnit ) '_'

            WriteXML = .FALSE.
            WriteData = .TRUE.
            GOTO 100
         ELSE
            WRITE( OutStr,'(A)') lf//'</AppendedData>'//lf
            CALL AscBinStrWrite( OutStr )
         END IF
      END IF


      WRITE( OutStr,'(A)') '</VTKFile>'//lf
      CALL AscBinStrWrite( OutStr )

      WRITE( OutStr,'(A)') ' '
      CALL AscBinStrWrite( OutStr )

      CLOSE( VtuUnit )


      CALL AscBinWriteFree()

      DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( ElemInd, ElemInd2 )

    END SUBROUTINE WriteVtuFile

  END SUBROUTINE ParticleOutputVtu

  !**************************************************************************

  SUBROUTINE GetMaxdDPrincipalDamageVert(Particles,No,layers,pddmax)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: D(2,2),T(3,3),pDold(3),pDnew(3)
    REAL(KIND=dp) :: Dnew(layers),newviscz(layers)
    REAL(KIND=dp) :: btzav,newviscav,oolmo
    REAL(KIND=dp) :: EigVals(2),EigVec(2,2)
    REAL(KIND=dp) :: CriticalDamage,DMax,pddmax,DavNew(4),eigdmax
    REAL(KIND=dp) :: one = 1.0_dp, half=1.0_dp/2.0_dp
    INTEGER :: infor,No,ii,layers,damdofs

    !get what your new vert int damage would be given your new dD
    oolmo = one/(DBLE(layers)-one)

    DO ii = 1,4

       Dnew = Particles % Damage(No,:,ii) + Particles % dD(No,:,ii)
       newviscz(:) = Particles % Bz(No,:)
       btzav = (SUM(newviscz)-half*(newviscz(1)+newviscz(layers))) * oolmo
       newviscz(:) = Dnew * Particles % Bz(No,:)
       newviscav = (SUM(newviscz)-half*(newviscz(1)+newviscz(layers))) * oolmo
       DavNew(ii) = newviscav/btzav
    END DO

    !old damage eigs
    D(1,1) = Particles % Dav(No,1)
    D(2,1) = Particles % Dav(No,4)
    D(1,2) = Particles % Dav(No,4)
    D(2,2) = Particles % Dav(No,2)

    CALL Eigen2DSym_TryGenFirst(D,EigVals,EigVec)

    pDold(1:2) = EigVals
    pDold(3) = Particles % Dav(No,3)

    !new damage eigs
    D(1,1) = DavNew(1)
    D(2,1) = DavNew(4)
    D(1,2) = DavNew(4)
    D(2,2) = DavNew(2)

    CALL Eigen2DSym_TryGenFirst(D,EigVals,EigVec)

    pDnew(1:2) = EigVals
    pDnew(3) = DavNew(3)

    !max difference between old and new damage
    pddmax = MAXVAL(ABS(pDnew-pDold))

  END SUBROUTINE GetMaxdDPrincipalDamageVert

  !**************************************************************************

  SUBROUTINE PrincipalDamage(Din,Dout)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: Din(4), D(2,2),Dout(3)
    REAL(KIND=dp) :: EigVals(2),EigVec(2,2)
    INTEGER :: infor,no,ii

    D(1,1) = Din(1)
    D(1,2) = Din(4)
    D(2,1) = Din(4)
    D(2,2) = Din(2)

    CALL Eigen2DSym_TryGenFirst(D,EigVals,EigVec)

    !since output of eigvalues2d is in ascending order
    Dout(1) = EigVals(2)
    Dout(2) = EigVals(1)
    Dout(3) = Din(3)

  END SUBROUTINE PrincipalDamage

  !**************************************************************************

  SUBROUTINE PrincipalEigenVec(Din,EVout)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: Din(4), D(2,2),EVout(3,3)
    REAL(KIND=dp) :: EigVals(2),EigVec(2,2)
    REAL(KIND=dp) :: one=1.0_dp,none=-1.0_dp,zero=0.0_dp

    D(1,1) = Din(1)
    D(2,1) = Din(4)
    D(1,2) = Din(4)
    D(2,2) = Din(2)

    CALL Eigen2DSym_TryGenFirst(D,EigVals,EigVec)


    IF (EigVals(1) == EigVals(2)) THEN
       !for the horizontal components:
       !if there is no preferred eigendirection,
       !return [0 0 -1]
       !This doesn't affect the vertical component, which is
       !always [0 0 1]
       EVout= zero
       EVout(3,1:2) = none
       EVout(3,3) = one
    ELSE
       Evout = zero
       EVout(1:2,2) = EigVec(1:2,1)
       EVout(1:2,1) = EigVec(1:2,2)
       EVout(3,3) = one
    END IF

  END SUBROUTINE PrincipalEigenVec

  !**************************************************************************

  SUBROUTINE FixPrincipalDamageVertInt(No,Model)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: layers
    REAL(KIND=dp) :: D(2,2)
    REAL(KIND=dp) :: CriticalDamage,DMax,CriticalDav
    REAL(KIND=dp) :: DMaxI,DMaxII,DMaxIII,RiftDmax
    REAL(KIND=dp) :: DavDMaxI,DavDMaxII,DavDMaxIII
    REAL(KIND=dp) :: TT,DD,lambda(3),en,sqrteig
    REAL(KIND=dp) :: quart=0.25_dp,half=0.5_dp,thres=0.0001_dp,zero=0.0_dp
    REAL(KIND=dp) :: EigValues(3),EigValues2d(2),WORK(68),EigenVec(2,2),EigenVec2(2,2)
    INTEGER :: infor,no,ii,jj,eigperm(3),maxloc1,midloc,minloc1,kk
    LOGICAL :: Visited = .FALSE.,rupt(3),useparticleeig
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: w,x,y,z

    SAVE :: Particles,DmaxI,DMaxII,DmaxIII,CriticalDav,&
         CriticalDamage,layers,Visited,riftdmax,DavDmaxI,DavDmaxII,&
         DavDmaxIII


    IF (.NOT. Visited) THEN

       Particles => GlobalParticles

       DmaxI = Particles % DmaxI
       DmaxII = Particles % DmaxII
       DmaxIII = Particles % DmaxIII

       DavDmaxI = Particles % DavDmaxI
       DavDmaxII = Particles % DavDmaxII
       DavDmaxIII = Particles % DavDmaxIII

       RiftDmax = Particles % RiftDmax

       CriticalDav = Particles % CriticalDav

       layers =  Particles % numberofparticlelayers
       Visited = .TRUE.
    END IF

    CriticalDamage = Particles % CriticalDamage

    D(1,1) = Particles % Dav(No,1)
    D(2,1) = Particles % Dav(No,4)
    D(1,2) = Particles % Dav(No,4)
    D(2,2) = Particles % Dav(No,2)

    TT = D(1,1)+D(2,2)
    DD = D(1,1)*D(2,2)-D(1,2)*D(1,2)

    sqrteig = quart*TT*TT-DD
    IF (sqrteig<0.0_dp) sqrteig = 0.0_dp
    sqrteig = sqrt(sqrteig)
    TT = half*TT
    lambda(1)=TT+sqrteig
    lambda(2)=TT-sqrteig
    lambda(3) = Particles % Dav(No,3)

    IF (.NOT. (ANY(lambda < 0.0_dp) .OR. ANY(lambda > CriticalDav) )) RETURN

    CALL Eigen2DSym_TryGenFirst_VecOnly(D,lambda(1:2),EigValues2d,EigenVec)

    EigValues(1:2) = EigValues2d
    EigValues(3) = Particles % Dav(No,3)

    !since eigvalues2d returns in order of smallest to biggest, max principal
    !damage is eigvalues2d(2) (assuming the horizontal D is bigger than the vertical)

    rupt(:) = .FALSE.

    DO ii = 1,3
       IF (EigValues(ii) > CriticalDav) rupt(ii) = .TRUE.
    END DO

    IF ((Particles % isorift .OR. Particles % damstatus(No) == 1)) THEN
       IF (ANY(rupt) .OR. Particles % damstatus(No) == 1) THEN
          rupt(:) = .TRUE.
       END IF
    END IF


    IF (ALL(rupt) .AND. Particles % useriftdmax) THEN
       !(DmaxI==DmaxII)) THEN

       Particles % Damage(No,:,:) = 0.0_dp
       Particles % Damage(No,:,1:2) = DMaxI
       Particles % Damage(No,:,3) = DmaxIII
       Particles % Dav(No,:) = 0.0_dp

       ! IF (Particles % useriftdmax) THEN
       Particles % Dav(No,1:3) = RiftDmax
       ! ELSE
       !    Particles % Dav(No,1:2) = DmaxI
       !    Particles % Dav(No,3) = DmaxIII
       ! END IF

       Particles % damstatus(No)=1

       RETURN
    END IF

    IF (ANY(rupt(1:2))) THEN

       DO jj = 1,layers

          D(1,1) = Particles % Damage(No,jj,1)
          D(2,2) = Particles % Damage(No,jj,2)
          D(1,2) = Particles % Damage(No,jj,4)
          D(2,1) = Particles % Damage(No,jj,4)


          IF (ANY(D>0.0_dp)) THEN
             !use the layer's eigenvectors to determine
             !which direction to rupt.

             CALL Eigen2DSym_TryGenFirst(D,EigValues2d,EigenVec2)


             IF (EigValues2d(1) == EigValues2d(2)) THEN
                UseParticleEig = .TRUE.
             ELSE

                D = 0.0_dp
                D(1,1) = EigValues2d(1)
                D(2,2) = EigValues2d(2)

                UseParticleEig = .FALSE.
             END IF

          ELSE
             UseParticleEig = .TRUE.
          END IF


          IF (UseParticleEig) THEN
             !use the particle eigenvectors to determine
             !which direction to rupt.
             EigenVec2 = EigenVec

             !eigenvalues ATA'
             !(according to the DAv EigenVec)
             D = MATMUL(TRANSPOSE(EigenVec2),D)
             D = MATMUL(D,EigenVec2)
          END IF


          IF (rupt(1)) D(1,1) = DMaxII
          IF (rupt(2)) D(2,2) = DMaxI

          !rotate back A'TA
          ! D = MATMUL(EigenVec2,D)
          ! D = MATMUL(D,TRANSPOSE(EigenVec2))

          !same as above, but faster
          w = D(1,1)*EigenVec2(1,1) + D(2,1)*EigenVec2(1,2)
          x = D(1,2)*EigenVec2(1,1) + D(2,2)*EigenVec2(1,2)
          y = D(1,1)*EigenVec2(2,1) + D(2,1)*EigenVec2(2,2)
          z = D(1,2)*EigenVec2(2,1) + D(2,2)*EigenVec2(2,2)

          D(1,1) = EigenVec2(1,1)*w + EigenVec2(1,2)*x
          D(2,2) = EigenVec2(2,1)*y + EigenVec2(2,2)*z
          D(1,2) = EigenVec2(2,1)*w + EigenVec2(2,2)*x


          Particles % Damage(No,jj,1) = D(1,1)
          Particles % Damage(No,jj,2) = D(2,2)
          Particles % Damage(No,jj,4) = D(1,2)
       END DO

       !then you need to vertically integrate again
       CALL VertIntDamFromVisc(Particles, No, layers,Model)

    END IF

    IF (rupt(3)) THEN
       Particles % Damage(No,:,3) = DmaxIII
       Particles % Dav(No,3) = DavDmaxIII
    END IF

    IF (ANY(rupt)) THEN
       Particles % DamStatus(No) = 1

       IF (Particles % useriftdmax) THEN
          !rupture all averaged components equally
          Particles % Dav(No,1:3) = RiftDmax
          Particles % Dav(No,4) = 0.0_dp
       ELSE

          IF (Particles % isorift) THEN
             !rupture all averaged components
             D(1,1) = Particles % Dav(No,1)
             D(2,2) = Particles % Dav(No,2)
             D(1,2) = Particles % Dav(No,4)
             D(2,1) = Particles % Dav(No,4)

             CALL Eigen2DSym_TryGenFirst(D,EigValues2d,EigenVec2)

             D = 0.0_dp
             D(1,1) = DavDMaxII
             D(2,2) = DavDMaxI

             w = D(1,1)*EigenVec2(1,1)
             x = D(2,2)*EigenVec2(1,2)
             y = D(1,1)*EigenVec2(2,1)
             z = D(2,2)*EigenVec2(2,2)

             Particles % Dav(No,1) = EigenVec2(1,1)*w + EigenVec2(1,2)*x
             Particles % Dav(No,2) = EigenVec2(2,1)*y + EigenVec2(2,2)*z
             Particles % Dav(No,4) = EigenVec2(2,1)*w + EigenVec2(2,2)*x
             Particles % Dav(No,3) = DavDMaxIII

          END IF
       END IF
    END IF


  END SUBROUTINE FixPrincipalDamageVertInt

  !**************************************************************************

  SUBROUTINE CheckPrincipalDamage(No)

    IMPLICIT NONE
    INTEGER :: No,ii
    REAL(KIND=dp) :: D(2,2),EigValues(2),EigenVec(2,2),DavEigenVec(2,2)
    TYPE(Particle_t), POINTER :: Particles

    Particles => GlobalParticles

    D(1,1) = Particles % Dav(No,1)
    D(2,1) = Particles % Dav(No,4)
    D(1,2) = Particles % Dav(No,4)
    D(2,2) = Particles % Dav(No,2)

    CALL Eigen2DSym_TryGenFirst(D,EigValues,DavEigenVec)

    PRINT *,''
    PRINT *,'PARTICLE NO',No
    PRINT *,'DAV EIGVALUES',EigValues
    PRINT *,'DAV EIGENVEC 1',DavEigenVec(:,2)
    PRINT *,''

    DO ii = 1,Particles % numberofparticlelayers

       IF (ANY(ABS(Particles % damage(No,ii,:)) > 0.0_dp)) THEN

          D(1,1) = Particles % damage(No,ii,1)
          D(2,1) = Particles % damage(No,ii,4)
          D(1,2) = Particles % damage(No,ii,4)
          D(2,2) = Particles % damage(No,ii,2)

          CALL Eigen2DSym_TryGenFirst(D,EigValues,EigenVec)


          PRINT *,'layer',ii
          PRINT *,'D',Particles % damage(No,ii,:)
          PRINT *,'mag',SUM(EigValues) + Particles % damage(No,ii,3)
          PRINT *,'eigvalues',EigValues
          PRINT *,'eigenvec 1',EigenVec(:,2)
          PRINT *,'diff',EigenVec(:,2)-DavEigenVec(:,2)
          PRINT *,''
       END IF
    END DO

  END SUBROUTINE CheckPrincipalDamage

  !**************************************************************************

  SUBROUTINE FixPrincipalDavAndLayers(Particles, No, layers,Model)

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No,layers,ii,jj
    REAL(KIND=dp), POINTER :: H
    REAL(KIND=dp) :: newviscz(layers,2),btzav,newviscav,Dmax,critdav,viscav3(layers)
    REAL(KIND=dp) :: D(2,2),Dav(2,2),Dav3,davvec(4)
    LOGICAL :: Visited = .FALSE.,rupt
    REAL(KIND=dp) :: lminval,lmidval,lmaxval,half=1.0_dp/2.0_dp,denom,one=1.0_dp,zero=0.0_dp
    REAL(KIND=dp) :: EigValues(2),DavEigValues(2),WORK(68),EigenVec(2,2),DavEigenVec(2,2)
    INTEGER :: infor
    REAL(KIND=dp) :: w,x,y,z

    SAVE :: critdav,visited

    IF (.NOT. Visited) THEN

       critdav = Particles % criticaldav
       Visited = .TRUE.
    END IF


    ! -----------   1.  Get principal directions for DAv ------------
    ! HERE, DAV IS NOT PARTICLE DAV. IT IS VERTICALLY-AVERAGED DAMAGE WITHOUT
    ! TEMPERATURE EFFECTS. THE POINT OF IGNORING TEMPERATURE EFFECTS HERE
    ! IS THAT WE ARE TRYING TO GET DAMAGE DIRECTIONS WEIGHTED BY
    ! DAMAGE MAGNITUDE, NOT ALSO WEIGHTED BY TEMPERATURE EFFECTS

    ! Dav(1,1) = Particles % Dav(No,1)
    ! Dav(2,1) = Particles % Dav(No,4)
    ! Dav(1,2) = Dav(2,1)
    ! Dav(2,2) = Particles % Dav(No,2)


    denom = one/(DBLE(layers)-one)

    DO ii = 1,4
       davvec(ii) = (SUM(Particles % damage(No,:,ii))-half* &
            (Particles % damage(No,1,ii)+Particles % damage(No,layers,ii))) * denom
    END DO

    Dav(1,1) = davvec(1)
    Dav(2,1) = davvec(4)
    Dav(1,2) = Dav(2,1)
    Dav(2,2) = davvec(2)

    CALL Eigen2DSym_TryGenFirst(Dav,DavEigValues,DavEigenVec)



    ! dont worry about rupturing Dav components at this point, as
    ! fixprincipaldamagevertint is called at the end, which will take care of it.

    ! ------ 2.  Get eigenvectors and eigenvalues for Damage on each layer ------
    ! Add eigenvalues to vector that vertically integrates the damage layers for Dav
    ! Rotate that layer back to cartesian using DavEigenVec


    DO ii = 1,layers

       !--- get principal values for the layer ---!

       D(1,1) = Particles % Damage(No,ii,1)
       D(2,1) = Particles % Damage(No,ii,4)
       D(1,2) = D(2,1)
       D(2,2) = Particles % Damage(No,ii,2)

       !CALL Eigen2D(D,EigValues,EigenVec)
       CALL Eigen2DSym_TryGenFirst(D,EigValues,EigenVec)



       ! --- rotate back using daveigenvec ---!
       ! D = MATMUL(DavEigenVec,D); D = MATMUL(D,TRANSPOSE(DavEigenVec))

       !same as above, but faster
       w = EigValues(1)*DavEigenVec(1,1)
       x = EigValues(2)*DavEigenVec(1,2)
       y = EigValues(1)*DavEigenVec(2,1)
       z = EigValues(2)*DavEigenVec(2,2)

       D(1,1) = DavEigenVec(1,1)*w + DavEigenVec(1,2)*x
       D(2,2) = DavEigenVec(2,1)*y + DavEigenVec(2,2)*z
       D(1,2) = DavEigenVec(2,1)*w + DavEigenVec(2,2)*x
       D(2,1) = D(1,2)

       IF (D(1,1) < zero) D(1,1) = zero
       IF (D(2,2) < zero) D(2,2) = zero

       Particles % Damage(No,ii,1) = D(1,1)
       Particles % Damage(No,ii,2) = D(2,2)
       Particles % Damage(No,ii,4) = D(1,2)

       ! --- layer viscosity for 3. --- !
       DO jj = 1,2
          newviscz(ii,jj) = EigValues(jj) * Particles % Bz(No,ii)
       END DO
    END DO

    ! ----- 3. Vert int to get Dav, fix principal values for Dav ---

    btzav = (SUM(Particles % Bz(No,:))-half* &
         (Particles % Bz(No,1)+Particles % Bz(No,layers))) * denom

    Dav = 0.0_dp

    DO jj = 1,2
       newviscav = (SUM(newviscz(:,jj))-half* &
            (newviscz(1,jj)+newviscz(layers,jj))) * denom

       Dav(jj,jj) = newviscav/btzav
       IF (Dav(jj,jj)<0.0_dp) Dav(jj,jj) = 0.0_dp
    END DO

    viscav3 = (Particles % Damage(No,:,3)) * Particles % Bz(No,:)
    Dav3 = (SUM(viscav3)-half*(viscav3(1)+viscav3(layers)) ) * denom
    Dav3 = Dav3/btzav
    IF (Dav3 < 0.0_dp) Dav3 = 0.0_dp

    IF ((ANY((Dav)>critdav) .OR. Dav3 > critdav) .AND. Particles % gamma==0.0_dp ) THEN

       Particles % Dav(No,1:2) = Particles % DMaxI
       Particles % Dav(No,3) = Particles % DmaxIII
       Particles % Dav(No,4) = 0.0_dp
       Particles % Damage(No,:,1:2) = Particles % DMaxI
       Particles % Damage(No,:,3) = Particles % DmaxIII
       Particles % Damage(No,:,4) = 0.0_dp
       RETURN
    END IF

    IF (ANY(Dav > critdav) .OR. Dav3 > critdav) THEN
       rupt = .TRUE.
    ELSE
       rupt = .FALSE.
    END IF

    ! Dav = MATMUL(DavEigenVec,Dav); Dav = MATMUL(Dav,TRANSPOSE(DavEigenVec))
    ! same as above, but faster
    w = Dav(1,1)*DavEigenVec(1,1)
    x = Dav(2,2)*DavEigenVec(1,2)
    y = Dav(1,1)*DavEigenVec(2,1)
    z = Dav(2,2)*DavEigenVec(2,2)

    Dav(1,1) = DavEigenVec(1,1)*w + DavEigenVec(1,2)*x
    Dav(2,2) = DavEigenVec(2,1)*y + DavEigenVec(2,2)*z
    Dav(1,2) = DavEigenVec(2,1)*w + DavEigenVec(2,2)*x
    Dav(2,1) = Dav(1,2)

    Particles % Dav(No,1) = Dav(1,1)
    Particles % Dav(No,2) = Dav(2,2)
    Particles % Dav(No,3) = Dav3
    Particles % Dav(No,4) = Dav(1,2)

    IF (rupt) THEN
       CALL FixPrincipalDamageVertInt(No,Model)
    END IF

  END SUBROUTINE FixPrincipalDavAndLayers

  !**************************************************************************

  SUBROUTINE GetParticleLayerStressesforEllipse(Particles,No,numoflayers,pstressdir,pstress,&
       groundbasalwaterp,ellipsesthres)

    IMPLICIT NONE

    TYPE(Particle_t),  POINTER :: Particles
    REAL(KIND=dp) :: pstressdir(numoflayers,2,2),pstress(numoflayers,2),z
    REAL(KIND=dp) :: n,Eeexp,EFExp,zsrhs,MinSRInvSquared,rhowtimesgravity,rhoitimesgravity,&
         oneovernumoflayersminus1,Identity(3,3),psr(3,3),rhs,zs,ID(3,3),&
         ESR(3,3),Tau(2,2),IDn1(3,3),ETau(2,2),peff,eigvalues(2),eigenvec(2,2)
    REAL(KIND=dp),allocatable :: zref(:)
    REAL(KIND=dp) :: one=1.0_dp,zero=0.0_dp,onethird=1.0_dp/3.0_dp,half=0.5_dp,&
         three=3.0_dp,onepfive=1.5_dp,two=2.0_dp,quart=0.25_dp, n3 = -3.0_dp
    REAL(KIND=dp) :: determ,dxx,dyy,dzz,dxy,ee,denom,ellipsesthres
    REAL(KIND=dp) :: TT,DD
    INTEGER :: numoflayers,ii,No,whichsurf,start,finish,step
    LOGICAL :: Visited = .FALSE., groundbasalwaterp
    REAL(KIND=dp) :: eyyd2m1,exxd1m1,ezzd3m1,exyd4,taurhsmultthird
    REAL(KIND=dp) :: t1d2m1,t2d1m1,t3od3m1,d4t4
    REAL(KIND=dp) :: k1,k2,maxtracesigma,Tau33,tracesigma,voutprod(3,3)
    REAL(KIND=dp) :: scalevec1(numoflayers)
    REAL(KIND=dp) :: kf,maxvec(3)
    REAL(KIND=dp) :: TraceEtau, chi, ah, td1, Bh,pcont

    SAVE :: n,Eeexp,EFExp,zsrhs,minsrinvsquared,rhowtimesgravity,rhoitimesgravity,&
         oneovernumoflayersminus1,zref,maxtracesigma,k1,k2,Visited ,ah,Bh,pcont

    IF (.NOT. Visited) THEN
       ! various variable shortcuts for damage calculations
       n = one/Particles % Viscosityexponent
       EeExp =  (one-n)/(two * n)
       EFexp = -one/n
       zsRHS = (one - Particles % rhoi/Particles % rhow)
       MinSRInvSquared = Particles % criticalshearrate*Particles % criticalshearrate
       rhowtimesgravity = Particles % rhow*ABS(Particles % gravity)
       rhoitimesgravity = Particles % rhoi*ABS(Particles % gravity)
       oneovernumoflayersminus1 = one/DBLE(numoflayers-1)

       ALLOCATE(zref(numoflayers))

       DO ii = 1,numoflayers
          zref(ii) = -one + DBLE(ii-1)*oneovernumoflayersminus1
       END DO

       k1 = Particles % k1
       k2 = Particles % k2

       IF (k2 .NE. 0.0_dp) THEN
          !to cap k_sigma at 30:
          maxtracesigma = (30.0_dp-k1)/k2
       ELSE
          maxtracesigma = HUGE(1.0_dp)
       END IF

       ah = Particles % ah
       Bh = Particles % bh
       pcont = one-ah-Bh

       Visited = .TRUE.
    END IF

    ! scalevec1(:) = 1.0_dp
    pstress = 0.0_dp
    pstressdir = 0.0_dp

    pSR(1,1) = Particles % GradVel(No,1)
    pSR(2,1) = half*(Particles % Gradvel(No,3) + Particles % GradVel(No,4))
    pSR(1,2) = pSR(2,1)
    pSR(2,2) = Particles % GradVel(No,2)
    pSR(3,3) = -pSR(1,1)-pSR(2,2)
    Ee = half*(pSR(1,1)*pSR(1,1) + pSR(2,2)*pSR(2,2) + &
         pSR(3,3)*pSR(3,3)) + pSR(1,2)*pSR(1,2)
    IF (Ee < MinSRInvSquared) THEN
       Ee = MinSRInvSquared
    END IF
    RHS = (Ee**EeExp) * (Particles % EF(No)**EFexp)
    IF (Particles % Gmask(No) < zero) THEN
       zs = Particles % H(No)+Particles % Bedrock(No)
    ELSE
       zs = Particles % H(No)*zsRHS
    END IF

    DO whichsurf = 1,2

       IF (whichsurf == 1) THEN
          !starting from bottom
          start = 1
          finish = INT(Particles % xpic(No,5))
          step = 1

          IF (finish == 0) CYCLE
       ELSE
          !starting from top
          start = numoflayers
          finish = INT(Particles % xpic(No,6))
          step = -1

          IF (finish>start) EXIT
          IF (finish==0) EXIT
       END IF

       DO ii = start,finish,step

          z = zs + Particles % H(No) * zref(ii)

          IF (whichsurf==2 .OR. z>=0.0_dp .OR. &
               (Particles % GMask(No)<zero .AND. (.NOT. groundbasalwaterp))) THEN
             !NO BASAL WATER PRESSURE
             Particles % pressure1 = rhoitimesgravity * (zs-z)
          ELSE
             !BASAL WATER PRESSURE
             Particles % pressure1 = rhoitimesgravity*(zs-z) + rhowtimesgravity*z
          END IF

          Particles % RHS = RHS*Particles % Bz(No,ii)
          Dxx = Particles % damage(No,ii,1) + Particles % dD(No,ii,1)
          Dyy = Particles % damage(No,ii,2) + Particles % dD(No,ii,2)
          Dzz = Particles % damage(No,ii,3) + Particles % dD(No,ii,3)
          Dxy = Particles % damage(No,ii,4) + Particles % dD(No,ii,4)


          !faster way, also isn't it the actual deviatoric stress we want, not the
          !effective dev stress?

          exxd1m1 = psr(1,1)*(Dxx-one)
          eyyd2m1 = psr(2,2)*(Dyy-one)
          ezzd3m1 = psr(3,3)*(Dzz-one)
          exyd4 = psr(1,2)*Dxy
          taurhsmultthird = Particles % RHS * onethird
          Tau(1,1) = taurhsmultthird * (-two*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
          Tau(2,2) = taurhsmultthird * (exxd1m1-two*eyyd2m1+ezzd3m1-exyd4)
          Tau(1,2) = -Particles % RHS * half * (psr(1,2)*(Dxx+Dyy-two) + Dxy*(psr(1,1)+psr(2,2)))
          Tau(2,1) = Tau(1,2)

          !------- effective pressure (Peff)------
          Peff = Particles % pressure1 - (Tau(1,1)+Tau(2,2))

          !------------uncomment if using actual stress
          Tau(1,1) = Tau(1,1) - Peff
          Tau(2,2) = Tau(2,2) - Peff
          CALL Eigen2DSym_TryGenFirst(Tau,EigValues,EigenVec)


          !WHERE (0.0_dp <= EigValues .AND. EigValues < 1.0e-20_dp) EigValues = 1.0e-20_dp
          !WHERE (0.0_dp >= EigValues .AND. EigValues > -1.0e-20_dp) EigValues = -1.0e-20_dp

          WHERE (EigValues >= 0.0_dp .AND. EigValues < Particles % mindam) EigValues = Particles % mindam
          WHERE (EigValues <= 0.0_dp .AND. EigValues > -Particles % mindam) EigValues = -Particles % mindam

          !------- uncomment if using actual stress ---------
          ! this prevents compressive stresses from contributing to the ellipse shape.
          !  WHERE (EigValues < Particles % mindam) EigValues = Particles % mindam

          !EigenVec(1:2,2) is the max eigenvector
          pstressdir(ii,1,:) = EigenVec(1:2,2)
          pstressdir(ii,2,:) = EigenVec(1:2,1)
          pstress(ii,1) = EigValues(2)
          pstress(ii,2) = EigValues(1)


          ! ----- using actual stress? -----
          ! denom = one/( Dxy*Dxy + Dxx + Dyy -Dxx*Dyy - one)
          ! t1d2m1 = Tau(1,1)*(Dyy-one)*denom
          ! t2d1m1 = Tau(2,2)*(Dxx-one)*denom
          ! t3od3m1 = taurhsmultthird * (exxd1m1+eyyd2m1-two*ezzd3m1+two*exyd4)/(Dzz-one)
          ! d4t4 = Tau(1,2)*Dxy*denom

          ! Etau(1,1) = onethird*(two*t1d2m1 - t2d1m1 +t3od3m1 -d4t4)
          ! Etau(2,2) = onethird*(-t1d2m1 + two*t2d1m1 +t3od3m1 -d4t4)
          ! !Etau(3,3) = onethird*(-t1d2m1 - t2d1m1 - two*t3od3m1 + two*d4t4)
          ! Etau(1,2) = half*denom*(tau(1,2)*(Dxx+Dyy-two) - Dxy*(Tau(1,1)+Tau(2,2)))
          ! Etau(2,1) = Etau(1,2)

          ! ETau(1,1) = ETau(1,1) - Peff
          ! ETau(2,2) = ETau(2,2) - Peff
          ! CALL Eigen2DSym_TryGenFirst(ETau,EigValues,EigenVec)

          ! WHERE (EigValues < Particles % mindam) EigValues = Particles % mindam

          ! !EigenVec(1:2,2) is the max eigenvector
          ! pstressdir(ii,1,:) = EigenVec(1:2,2)
          ! pstressdir(ii,2,:) = EigenVec(1:2,1)
          ! pstress(ii,1) = EigValues(2)
          ! pstress(ii,2) = EigValues(1)

       END DO
    END DO

  END SUBROUTINE GetParticleLayerStressesforEllipse

  !**************************************************************************

  SUBROUTINE nonlocalsurroundelems(Particles,Mesh,lc)
    TYPE(Particle_t),  POINTER :: Particles
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: lc,hgres,dist,t,hdist
    REAL(KIND=dp) :: mcoords(2),tcoords(3)
    INTEGER :: telems,count,ElementIndNew, ii, jj, kk
    REAL(KIND=dp), allocatable :: scoordsx(:),scoordsy(:)
    INTEGER, allocatable :: temp(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Visited=.FALSE.

    SAVE :: Visited

    IF (.NOT. Visited) THEN
       Visited = .TRUE.
    ELSE
       RETURN
    END IF

    !for each element,
    !this subroutine returns the surrounding elements
    !within the nonlocal length scale, lc,
    !and puts them in ElemParticles( currentelem ) % NLSelems()
    !and saves the total numer of these surrounding elems in
    !ElemParticles( currentelem ) % numNLSelems.

    IF (lc < Particles % gridres) THEN
       ElemParticles(:) % numNLSelems = 0
       RETURN
    END IF

    tcoords(3) = 0.0_dp
    hgres = 0.5_dp * Particles % gridres
    dist = CEILING(lc/Particles % gridres)


    !total number of nonlocal elems per elem
    t = (dist*2.0_dp + 1.0_dp)
    telems = INT(t)
    ALLOCATE(scoordsx(telems))
    ALLOCATE(scoordsy(telems))
    ALLOCATE(temp(telems*telems))


    !half dist across nonlocal square
    dist = dist* Particles % gridres
    !hdist = 0.5_dp*dist


    DO ii = 1,Mesh % NumberOfBulkElements

       Element => Mesh % Elements(ii)
       NodeIndexes => Element % NodeIndexes

       !middle coords of the current element
       mcoords(1) = MINVAL(Mesh % Nodes % x(NodeIndexes(1:4))) + hgres
       mcoords(2) = MINVAL(Mesh % Nodes % y(NodeIndexes(1:4))) + hgres

       CALL linspace(mcoords(1)-dist,mcoords(1)+dist,scoordsx)
       CALL linspace(mcoords(2)-dist,mcoords(2)+dist,scoordsy)

       temp = 0
       count = 0
       DO jj = 1,telems

          !x test coords
          tcoords(1) = scoordsx(jj)

          !dont include elements within 1 grid cell as they are already in
          !elemparticles(ii) % surroundingelems()
          ! IF (ABS(tcoords(1)-mcoords(1)) < Particles % gridres*1.1_dp) CYCLE
          DO kk = 1,telems

             !y test coords
             tcoords(2) = scoordsy(kk)
             ! IF (ABS(tcoords(1)-mcoords(1)) < Particles % gridres*1.1_dp) CYCLE
             ElementIndNew = 0
             CALL LocateParticleInMeshOctree( ElementIndNew, tcoords )

             IF (ElementIndNew == 0) THEN
                CYCLE
             ELSE
                count = count + 1
                temp(count) = ElementIndNew
             END IF
          END DO
       END DO

       ElemParticles(ii) % numNLSelems = count
       IF (count > 0) THEN
          ALLOCATE (ElemParticles(ii) % NLSelems(count) )
          ElemParticles(ii) % NLSelems(1:count) = temp(1:count)
       END IF
    END DO


    IF (ALLOCATED(scoordsx)) DEALLOCATE(scoordsx)
    IF (ALLOCATED(scoordsy)) DEALLOCATE(scoordsy)
    IF (ALLOCATED(temp)) DEALLOCATE(temp)


  END SUBROUTINE nonlocalsurroundelems

  !**************************************************************************

  SUBROUTINE nonlocalintegraldDellipseRobust(Particles, numoflayers, &
       count, lc, gaussk, gridres,vertlc, groundbasalwaterp,ellipsesthres,&
       justusegaussian,Mesh)

    IMPLICIT NONE

    TYPE(Particle_t),  POINTER :: Particles
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: count,numoflayers,extracount,curcount
    REAL(KIND=dp) :: lc,gaussk,gridres,vertlc,ellipsesthres,fsquared
    LOGICAL ::  groundbasalwaterp,justusegaussian
    REAL(KIND=dp) :: num(numoflayers,4)
    REAL(KIND=dp) :: Coord(2),CompCoord(2)
    REAL(KIND=dp) :: dnew(numoflayers,4),denomvec(numoflayers,4)
    REAL(KIND=dp) :: TTT(numoflayers),DDD(numoflayers),eig1(numoflayers)
    REAL(KIND=dp) :: pstress(numoflayers,2),diff(numoflayers,2),phivec(numoflayers)
    REAL(KIND=dp) :: cosphi(numoflayers),sinphi(numoflayers),rho(numoflayers)
    REAL(KIND=dp),ALLOCATABLE :: nonlocdD(:,:,:),elldenom(:,:,:)
    REAL(KIND=dp),ALLOCATABLE :: ruptstat(:,:)
    INTEGER, POINTER :: eligperm(:)=>NULL()
    INTEGER,TARGET :: eligref(numoflayers)
    REAL(KIND=dp) :: ruptstat2(numoflayers),match2(numoflayers)
    REAL(KIND=dp) :: match(numoflayers),lmaxval, edgescale
    INTEGER :: ElementInd,basalref,surfref,eliglayers,smoothlayers
    INTEGER :: No, CompNo, ii, jj, kk, j, ind, ind2, p, q, m, mm
    REAL(KIND=dp) :: phiparam, tot, origdam, dist, phi,denom
    REAL(KIND=dp) :: psr(3,3),gaussksquared,minrho,tottime
    REAL(KIND=dp) :: y(numoflayers),dblesmoothlayers,num1(4),lcsquared
    LOGICAL ::  Visited=.FALSE.
    REAL(KIND=dp) :: pstressdir(numoflayers,2,2),D(2,2),EigValues(2),EigenVec(2,2)
    REAL(KIND=dp) :: one=1.0_dp,zero=0.0_dp
    LOGICAL :: edge,edge2
    REAL(KIND=dp) :: phi2,Coordreflecty


    !the nonlocal regularization from duddu 2013 (gaussian) and giry 2011 (ellipse)

    fsquared = ellipsesthres*ellipsesthres
    lcsquared = lc*lc

    gaussksquared = gaussk*gaussk
    minrho = Particles % gridres * Particles % gridres

    Particles % xpic(:,5:6) = Particles % xpic(:,3:4)

    ! -------------------------------------------------------------------
    ! 1: Particles with local dD > 0 already marked (xpic(No,1)>0)
    !    Find the particles with local dD = 0 that will gain nonlocal dD
    ! -------------------------------------------------------------------

    origdam = DBLE(count)
    tot = origdam

    ! save on the elements the max number of layers from bottom and top of
    ! surrounding elems within +/- lc from all current elem nodes
    ElemTrack(:) % ddlayfrombottom1 = 0
    ElemTrack(:) % ddlayfromtop1 = numoflayers

    ! for each element, get max layer from top and bottom of all of the element's particles
    DO ind = 1, Mesh % NumberOfBulkElements
       IF (ElemParticles(ind) % NumberOfParticles < 1) CYCLE
       eligperm => ElemParticles(ind) % p(1:ElemParticles(ind) % NumberOfParticles)
       ElemTrack(ind) % ddlayfrombottom1 = MAXVAL(Particles % xpic(eligperm,5))
       ElemTrack(ind) % ddlayfromtop1 = MINVAL(Particles % xpic(eligperm,6))
    END DO

    ! now set each element max layers to max of surrounding elements within nonlocal range
    ElemTrack(:) % ddlayfrombottom2 = ElemTrack(:) % ddlayfrombottom1
    ElemTrack(:) % ddlayfromtop2 = ElemTrack(:) % ddlayfromtop1

    DO ind = 1,Mesh % NumberOfBulkElements
       IF (ElemParticles(ind) % NumberOfParticles < 1) CYCLE
       eligperm => ElemParticles(ind) % NLSelems(:)
       ElemTrack(ind) % ddlayfrombottom2 = MAXVAL(ElemTrack(eligperm) % ddlayfrombottom1)
       ElemTrack(ind) % ddlayfromtop2 = MINVAL(ElemTrack(eligperm) % ddlayfromtop1)
    END DO

    ElemTrack(:) % ddlayfrombottom1 = ElemTrack(:) % ddlayfrombottom2
    ElemTrack(:) % ddlayfromtop1 = ElemTrack(:) % ddlayfromtop2

    ! repeat this operation once more for elements that layfrombottom = 0 or layfromtop = numoflayers
    ! this is needed because particles in these elements can end up contributing to newly marked
    ! particles in below routine
    DO ind = 1,Mesh % NumberOfBulkElements
       IF (ElemParticles(ind) % NumberOfParticles < 1) CYCLE
       IF (ElemTrack(ind) % ddlayfrombottom2 == zero .OR. ElemTrack(ind) % ddlayfromtop2 == numoflayers) THEN
          eligperm => ElemParticles(ind) % NLSelems(:)
          ElemTrack(ind) % ddlayfrombottom2 = MAXVAL(ElemTrack(eligperm) % ddlayfrombottom1)
          ElemTrack(ind) % ddlayfromtop2 = MINVAL(ElemTrack(eligperm) % ddlayfromtop1)
       END IF
    END DO


    ! now go through and find the new particles that might contribute to
    ! dD particles, or the extra particles that contrib to the new particles...
    ! count = 1
    extracount = 0
    DO No = 1, Particles % NumberOfParticles

       Particles % xpic(No,5) = ElemTrack(Particles % InterpElem(No)) % ddlayfrombottom2
       Particles % xpic(No,6) = ElemTrack(Particles % InterpElem(No)) % ddlayfromtop2


       IF (Particles % Status(No) >= PARTICLE_LOST) CYCLE
       IF (Particles % xpic(No,1) .NE. zero) CYCLE

       IF (Particles % xpic(No,5) > zero .OR. Particles % xpic(No,6) < numoflayers) THEN
          Coord = Particles % Coordinate(No,1:2)
          ElementInd = Particles % InterpElem(No)

          DO ii = 1,ElemParticles(ElementInd) % numNLSelems
             ind = ElemParticles(ElementInd) % NLSelems(ii)
             DO jj = 1,ElemParticles(ind) % NumberOfParticles
                CompNo = ElemParticles(ind) % p(jj)
                IF (Particles % xpic(CompNo,1) <= zero) CYCLE
                IF (Particles % xpic(CompNo,1) <= origdam) THEN
                   ! CompCoord = GetParticleCoord( Particles, CompNo)
                   CompCoord = Particles % Coordinate(CompNo,1:2)

                   dist = (CompCoord(1)-Coord(1))*(CompCoord(1)-Coord(1)) + &
                        (CompCoord(2)-Coord(2))*(CompCoord(2)-Coord(2))
                   !mark particle
                   IF (dist < lcsquared) THEN
                      tot = tot + one
                      Particles % xpic(No,1) = tot
                      EXIT
                   END IF
                END IF
             END DO
             IF (Particles % xpic(No,1)>zero) EXIT
          END DO
          IF (Particles % xpic(No,1) == zero) THEN
             Particles % xpic(No,1) = -1
             extracount = extracount + 1
          END IF
       END IF
    END DO

    ! -------------------------------------------------------------------
    ! 2: allocate space to calculate nonlocal contributions
    ! -------------------------------------------------------------------
    IF (tot > zero) THEN

       ALLOCATE(nonlocdD(INT(tot),numoflayers,4))
       ALLOCATE(elldenom(INT(tot),numoflayers,4))
       ALLOCATE(ruptstat(INT(tot),numoflayers))


       nonlocdD = zero
       elldenom = zero
       ruptstat = zero

    ELSE
       RETURN
    END IF

    ! -------------------------------------------------------------------
    ! 3: Particle loop: nonlocal contribution
    !    from No (current) to CompNo (surrounding)
    ! -------------------------------------------------------------------


    DO No = 1, Particles % NumberOfParticles

       count = INT(Particles % xpic(No,1))
       IF (count <= 0) CYCLE

       DO ii = 1,MIN(INT(Particles % xpic(No,5)),numoflayers)
          CALL MaxPFour(Particles % damage(No,ii,1:4),lmaxval)
          IF (lmaxval < Particles % criticaldamage) THEN
             ruptstat(count,ii) = 1.0_dp
          ELSE
             ruptstat(count,ii) = 0.0_dp
          END IF
       END DO

       DO ii = numoflayers,MAX(INT(Particles % xpic(No,6)),1),-1
          CALL MaxPFour(Particles % damage(No,ii,1:4),lmaxval)
          IF (lmaxval < Particles % criticaldamage) THEN
             ruptstat(count,ii) = 1.0_dp
          ELSE
             ruptstat(count,ii) = 0.0_dp
          END IF
       END DO
    END DO


    DO No = 1, Particles % NumberOfParticles

       curcount = INT(Particles % xpic(No,1))
       IF (curcount == 0) CYCLE
       ! IF (Particles % xpic(No,1) == zero) CYCLE
       ! IF (Particles % damstatus(No) == -1) CYCLE

       !ADDED 2021---
       !For 'extracount' particles (xpic(No,1)==-1), ruptstat is determined here
       !and saved on ruptstat2
       IF (curcount < 0 ) then
          ruptstat2=zero
          DO ii = 1,MIN(INT(Particles % xpic(No,5)),numoflayers)
             CALL MaxPFour(Particles % damage(No,ii,1:4),lmaxval)
             IF (lmaxval < Particles % criticaldamage) ruptstat2(ii) = 1.0_dp
          END DO

          DO ii = numoflayers,MAX(INT(Particles % xpic(No,6)),1),-1
             CALL MaxPFour(Particles % damage(No,ii,1:4),lmaxval)
             IF (lmaxval < Particles % criticaldamage) ruptstat2(ii) = 1.0_dp
          END DO
       ENDIF
       !---

       IF (Particles % xpic(No,5) >= Particles % xpic(No,6)) THEN
          ! basal and surface crevasses overlap.
          Particles % xpic(No,5) = numoflayers
          Particles % xpic(No,6) = 0

          eligref = (/(j,j=1,numoflayers)/)
          eligperm => eligref

          eliglayers = numoflayers
       ELSE
          eligref(1) = 0
          basalref = INT(Particles % xpic(No,5))
          surfref = INT(Particles % xpic(No,6))
          IF (basalref + numoflayers-surfref + 1 >= numoflayers) THEN
             basalref = numoflayers
             surfref = 0
          END IF
          IF (basalref>0) THEN
             eligref(1:basalref) = (/(j,j=1,basalref)/)
          END IF
          IF (surfref >0) THEN
             eligref(basalref+1:basalref+1+(numoflayers-surfref)) = &
                  (/(j,j=surfref,numoflayers)/)
             eliglayers = basalref+(numoflayers-surfref+1)
          ELSE
             eliglayers = basalref
          END IF
          IF (eliglayers > numoflayers) eliglayers = numoflayers  !just to be sure
          eligperm => eligref(1:eliglayers)
          IF (eligref(1)==0) THEN
             PRINT *,'eligref(1) is 0'           !this shouldnt happen
             CYCLE
          END IF
       END IF




       ! ElementInd = GetParticleElement( Particles, No )

       ElementInd = Particles % InterpElem(No)

       ! IF (ElementInd == 0) THEN
       !    ElementInd = Particles % InterpElem(No)
       IF (ElementInd == 0) CYCLE
       ! END IF

       ! Coord = GetParticleCoord( Particles, No )
       Coord = Particles % Coordinate(No,1:2)
       ! IF (ElementInd == 0) CYCLE

       IF (.NOT. justusegaussian) THEN

          pstress = zero       !stresses for your particle
          CALL GetParticleLayerStressesForEllipse(Particles,No,numoflayers,&
               pstressdir,pstress,groundbasalwaterp,ellipsesthres) !,eliglayers,eligperm)

       END IF


       edge = .FALSE.
       IF (Particles % mismip) THEN
          IF (Coord(2) <= lc .OR. Coord(2) >= 40000.0_dp-lc) THEN
             !potential edge effects
             edge = .TRUE.
             IF (Coord(2) <= lc) THEN
                Coordreflecty = -Coord(2)
             ELSE
                Coordreflecty = 80000.0_dp-Coord(2)
             END IF

          END IF
       END IF


       ! loop through CompNo and accumulate contributions from No

       DO kk = 1,ElemParticles(ElementInd) % numNLSelems
          ind2 = ElemParticles(ElementInd) % NLSelems(kk)

          IF (ind2 == 0) CYCLE
          DO jj = 1,ElemParticles(ind2) % NumberOfParticles
             CompNo = ElemParticles(ind2) % p(jj)
             IF (CompNo > Particles % NumberOfParticles) CYCLE
             IF (Particles % Status(CompNo) > 4) CYCLE

             IF (Particles % xpic(CompNo,2) == No) CYCLE !already taken care of
             Particles % xpic(CompNo,2) = No
             IF (Particles % xpic(CompNo,1) <= zero) CYCLE
             !  IF (Particles % damstatus(No) .NE. Particles % damstatus(CompNo)) CYCLE
             ! IF (Particles % damstatus(CompNo) .EQ. -1) CYCLE

             IF (Particles % damstatus(CompNo) .NE. Particles % damstatus(No)) THEN
                IF (Particles % damstatus(CompNo) .EQ. 1 .OR. Particles % damstatus(No) .EQ. 1) CYCLE
             END IF


             ! CompCoord = GetParticleCoord( Particles, CompNo)
             CompCoord = Particles % Coordinate(CompNo,1:2)
             phi = (Coord(1)-CompCoord(1))*(Coord(1)-CompCoord(1)) + &
                  (Coord(2)-CompCoord(2))*(Coord(2)-CompCoord(2))
             IF (phi >= lcsquared) CYCLE

             edge2 = .FALSE.

             IF (Particles % mismip .AND. edge) THEN
                IF (CompCoord(2) <= lc .OR. CompCoord(2) >= 40000.0_dp-lc) THEN

                   !potential edge effects
                   phi2 = (Coord(1)-CompCoord(1))*(Coord(1)-CompCoord(1)) + &
                        (Coordreflecty-CompCoord(2))*(Coordreflecty-CompCoord(2))

                   IF (phi2 < lcsquared) THEN
                      edge2 = .TRUE.
                   END IF
                END IF
             END IF

             match = 0.0_dp

             count = INT(Particles % xpic(CompNo,1))

             !CHANGED 2021---
             !WHERE (ruptstat(curcount,:)== ruptstat(count,:)) match = 1.0_dp
             IF (curcount>0) then
                WHERE (ruptstat(curcount,:)== ruptstat(count,:)) match = 1.0_dp
             ELSE
                match2 = ruptstat(count,:)
                WHERE (match2(:)==ruptstat2(:)) match = 1.0_dp
             ENDIF
             !---

             IF (.NOT. justusegaussian) THEN

                edge2 = .FALSE.

                !---------ELLIPSE---------

                IF (phi == zero .OR. CompNo == No) THEN
                   phi = zero
                   rho(eligperm) = one !lcsquared
                ELSE
                   diff(eligperm,1) = CompCoord(1)-Coord(1)
                   diff(eligperm,2) = CompCoord(2)-Coord(2)

                   cosphi(eligperm) = SUM(diff(eligperm,:) * pstressdir(eligperm,1,1:2),DIM=2)
                   sinphi(eligperm) = SUM(diff(eligperm,:) * pstressdir(eligperm,2,1:2),DIM=2)

                   cosphi(eligperm) = cosphi(eligperm)*cosphi(eligperm)/phi
                   sinphi(eligperm) = sinphi(eligperm)*sinphi(eligperm)/phi

                   !rho squared
                   rho(eligperm) = one / &
                        (fsquared * &
                        (cosphi(eligperm)/&
                        (pstress(eligperm,1)*pstress(eligperm,1)) + &
                        sinphi(eligperm)/&
                        (pstress(eligperm,2)*pstress(eligperm,2))))

                   WHERE (rho(eligperm)>one) rho(eligperm)=one

                   rho(eligperm) = lcsquared * rho(eligperm)

                   WHERE (rho(eligperm) < minrho) &
                        rho(eligperm) = minrho
                END IF
             ELSE

                !---------GAUSSIAN---------

                rho(eligperm) = lcsquared
             END IF

             phivec(eligperm) = exp(-(gaussksquared*phi/rho(eligperm))) * match(eligperm)

             IF (edge2) THEN
                phivec(eligperm) = phivec(eligperm) + &
                     exp(-(gaussksquared*phi2/rho(eligperm))) * match(eligperm)
             END IF

             nonlocdD(count,eligperm,1) = nonlocdD(count,eligperm,1) + &
                  phivec(eligperm) * Particles % dD(No,eligperm,1)

             nonlocdD(count,eligperm,2) = nonlocdD(count,eligperm,2) + &
                  phivec(eligperm) * Particles % dD(No,eligperm,2)

             nonlocdD(count,eligperm,3) = nonlocdD(count,eligperm,3) + &
                  phivec(eligperm) * Particles % dD(No,eligperm,3)

             nonlocdD(count,eligperm,4) = nonlocdD(count,eligperm,4) + &
                  phivec(eligperm) * Particles % dD(No,eligperm,4)

             elldenom(count,eligperm,1) = elldenom(count,eligperm,1) + &
                  phivec(eligperm)


          END DO
       END DO
    END DO

    ! -------------------------------------------------------------------
    ! 4: Update Particles % dD based on nonlocal
    ! -------------------------------------------------------------------

    elldenom(:,:,2) = elldenom(:,:,1)
    elldenom(:,:,3) = elldenom(:,:,1)
    elldenom(:,:,4) = elldenom(:,:,1)


    WHERE (elldenom .NE. zero) nonlocdD = nonlocdD / elldenom

    DO No = 1, Particles % NumberOfParticles
       IF (Particles % xpic(No,1) <= zero) CYCLE
       !  IF (Particles % damstatus(No) == -1) CYCLE
       count = INT(Particles % xpic(No,1))
       Particles % dD(No,:,1:4) = nonlocdD(count,:,1:4)

    END DO


    IF (ALLOCATED(ruptstat)) DEALLOCATE(ruptstat)
    IF (ALLOCATED(elldenom)) DEALLOCATE(elldenom)
    IF (ALLOCATED(nonlocdD)) DEALLOCATE(nonlocdD)
    IF (ASSOCIATED(eligperm)) NULLIFY(eligperm)

    ! -------------------------------------------------------------------
    ! 5: Vertical smoothing
    ! -------------------------------------------------------------------

    IF (vertlc > 0) THEN
       lcsquared = vertlc*vertlc
       phiparam = -gaussksquared/lcsquared
       DO No = 1, Particles % NumberOfParticles
          IF (Particles % xpic(No,1) <= zero) CYCLE
          y = (/(m,m=0,numoflayers-1,1)/)
          y = y * (Particles % H(No)/(numoflayers-one))
          smoothlayers = INT(CEILING(vertlc/ ( y(2)-y(1) ) + 2.0_dp))
          IF (smoothlayers<=1) CYCLE
          dnew = zero
          DO ii = 1,numoflayers
             num1 = zero
             denom = zero
             DO jj = MAX(ii-smoothlayers,1),MIN(ii+smoothlayers,numoflayers)
                dist = y(ii)-y(jj)
                dist = dist*dist
                IF (dist < lcsquared) THEN
                   phi = EXP(dist*phiparam)
                   num1 = num1 + (phi * Particles % dD(No,jj,:))
                   denom = phi + denom
                END IF
             END DO
             dnew(ii,:) = num1/denom
          END DO
          Particles % dD(No,:,:) = dnew(:,:)
       END DO
       PRINT *,'done'
    END IF


  END SUBROUTINE nonlocalintegraldDellipserobust

  !**************************************************************************

  SUBROUTINE nonlocalintegraldD(Particles, numoflayers, countin, lc, gaussk, gridres,vertlc)

    IMPLICIT NONE

    TYPE(Particle_t),  POINTER :: Particles
    INTEGER,INTENT(in) :: countin
    REAL(KIND=dp) :: num(numoflayers,4), Coord(3),CompCoord(3),dnew(numoflayers,4)
    REAL(KIND=dp),allocatable :: nonlocdD(:,:,:)
    INTEGER :: No, CompNo, ii, jj, kk, ElementInd, ind, ind2, numoflayers
    INTEGER :: m, smoothlayers,count
    REAL(KIND=dp) :: lc,gaussk,phiparam, tot, origdam, denom, dist, phi, gridres
    REAL(KIND=dp) :: y(numoflayers),dblesmoothlayers,num1(4),lcsquared,vertlc

    ! the nonlocal regularization from duddu 2013.

    ! NOTE:: nonlocal scheme currently restricted to the particles within 2 surrounding grid cells
    lcsquared = lc*lc
    phiparam = -gaussk/lcsquared

    ! first, determine which particles with dD = 0 currently will gain  damage
    origdam = DBLE(countin)
    tot = DBLE(countin)
    count = 1

    DO No = 1, Particles % NumberOfParticles

       IF (Particles % Status(No) >= PARTICLE_LOST) CYCLE
       IF (Particles % xpic(No,1) == DBLE(count)) THEN
          count = count + 1

          ElementInd = GetParticleElement( Particles, No )
          Coord = GetParticleCoord( Particles, No )

          IF (ElementInd == 0) CYCLE
          DO ii = 1,9
             ind = ElemParticles(ElementInd) % SurroundingElems(ii)
             IF (ind == 0) CYCLE
             DO kk = 1,9
                ind2 = ElemParticles(ind) % SurroundingElems(kk)
                IF (ind2 == 0) CYCLE
                DO jj = 1,ElemParticles(ind2) % NumberOfParticles
                   CompNo = ElemParticles(ind2) % p(jj)

                   !
                   !
                   !
                   !IF (Particles % DmaxII_dom) THEN
                   IF (Particles % damStatus(No) .NE. Particles % damStatus(CompNo)) CYCLE
                   !END IF

                   !dont mix isodam particles and regular particles
                   !IF (Particles % damStatus(No) == -1) THEN
                   !   IF (Particles % damStatus(CompNo) .NE. -1) CYCLE
                   !ELSE
                   !   IF (Particles % damStatus(CompNo) == -1) CYCLE
                   !END IF
                   !
                   !
                   !
                   IF (Particles % xpic(CompNo,1) == 0.0_dp) THEN

                      CompCoord = GetParticleCoord( Particles, CompNo)

                      dist = (CompCoord(1)-Coord(1))*(CompCoord(1)-Coord(1)) + &
                           (CompCoord(2)-Coord(2))*(CompCoord(2)-Coord(2))
                      !mark particle
                      IF (dist <= lcsquared) THEN
                         tot = tot + 1.0_dp
                         Particles % xpic(CompNo,1) = tot
                      END IF
                   END IF
                END DO
             END DO
          END DO
       END IF
    END DO

    IF (tot > 0.0_dp) THEN
       ALLOCATE(nonlocdD(INT(tot),numoflayers,4))
       nonlocdD = 0.0_dp
    ELSE
       RETURN
    END IF

    DO No = 1, Particles % NumberOfParticles

       IF (Particles % xpic(No,1) == 0.0_dp) CYCLE

       num = 0.0_dp
       denom = 0.0_dp
       count = INT(Particles % xpic(No,1))

       ElementInd = GetParticleElement( Particles, No )
       Coord = GetParticleCoord( Particles, No )

       IF (ElementInd == 0) CYCLE

       DO ii = 1,9
          ind = ElemParticles(ElementInd) % SurroundingElems(ii)
          IF (ind == 0) CYCLE
          DO kk = 1,9
             ind2 = ElemParticles(ind) % SurroundingElems(kk)
             IF (ind2 == 0) CYCLE
             DO jj = 1,ElemParticles(ind2) % NumberOfParticles

                CompNo = ElemParticles(ind2) % p(jj)
                !
                !
                !
                !IF (Particles % DmaxII_dom) THEN
                IF (Particles % damStatus(No) .NE. Particles % damStatus(CompNo)) CYCLE
                !END IF
                !dont mix isodam particles and regular particles
                !IF (Particles % damStatus(No) == -1) THEN
                !   IF (Particles % damStatus(CompNo) .NE. -1) CYCLE
                !ELSE
                !   IF (Particles % damStatus(CompNo) == -1) CYCLE
                !END IF
                !
                !
                !
                !nonlocal between No and CompNo has already been taken care of.
                IF (Particles % xpic(CompNo,2) == No) CYCLE

                !origdam is the original amount of particles with dD coming into the routine.
                !so it is excluding the  particles that have no local dD, but are gaining dD
                !due to the nonlocal scheme.

                IF (Particles % xpic(No,1) <= origdam) THEN

                   !if xpic(compno,1) is zero, it did not have local dD and
                   !was determined to not contribute to nonlocal dD in the first set of
                   !loops here.
                   IF (Particles % xpic(CompNo,1) == 0.0_dp) CYCLE

                   CompCoord = GetParticleCoord( Particles, CompNo)
                   phi = (CompCoord(1)-Coord(1))*(CompCoord(1)-Coord(1)) + &
                        (CompCoord(2)-Coord(2))*(CompCoord(2)-Coord(2))
                ELSE
                   CompCoord = GetParticleCoord( Particles, CompNo)
                   phi = (CompCoord(1)-Coord(1))*(CompCoord(1)-Coord(1)) + &
                        (CompCoord(2)-Coord(2))*(CompCoord(2)-Coord(2))
                END IF
                !  IF (sqrt(phi) > lc) CYCLE
                IF (phi > lcsquared) CYCLE
                phi = exp(phi*phiparam)
                num = num + (phi * Particles % dD(CompNo,:,:))
                denom = phi + denom
                Particles % xpic(CompNo,2) = No
             END DO
          END DO
       END DO

       nonlocdD(count,:,:) = num/denom
    END DO

    Particles % dD = 0.0_dp
    DO No = 1, Particles % NumberOfParticles
       IF (Particles % xpic(No,1) == 0.0_dp) CYCLE
       count = INT(Particles % xpic(No,1))
       Particles % dD(No,:,:) = nonlocdD(count,:,:)
    END DO

    !vertical smoothing
    IF (vertlc > 0) THEN
       !  lc = vertlc
       lcsquared = vertlc*vertlc
       phiparam = -gaussk/lcsquared
       ! dblesmoothlayers = DBLE(smoothlayers)
       DO No = 1, Particles % NumberOfParticles
          IF (Particles % xpic(No,1) == 0.0_dp) CYCLE
          y = (/(m,m=0,numoflayers-1,1)/)
          y = y * (Particles % H(No)/(numoflayers-1.0_dp))
          !lc = y(2)-y(1)*dblesmoothlayers
          smoothlayers = INT(CEILING(vertlc/ ( y(2)-y(1) ) + 2.0_dp))
          IF (smoothlayers<=1) CYCLE
          dnew = 0.0_dp
          DO ii = 1,numoflayers
             num1 = 0.0_dp
             denom = 0.0_dp
             DO jj = MAX(ii-smoothlayers,1),MIN(ii+smoothlayers,numoflayers)
                dist = y(ii)-y(jj)
                dist = dist*dist
                IF (dist <= lcsquared) THEN
                   phi = EXP(dist*phiparam)
                   num1 = num1 + (phi * Particles % dD(No,jj,:))
                   denom = phi + denom
                END IF
             END DO
             dnew(ii,:) = num1/denom
          END DO
          Particles % dD(No,:,:) = dnew(:,:)
       END DO
    END IF
    DEALLOCATE(nonlocdD)
  END SUBROUTINE nonlocalintegraldD

  !**************************************************************************

  SUBROUTINE smoothrupth(Particles, Mesh, lc, smoothiters)

    IMPLICIT NONE

    TYPE(Particle_t),  POINTER :: Particles
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: No,ii,kk,jj,ind2,CompNo,ElementInd,smoothiters
    REAL(KIND=dp) :: Coord(2),CompCoord(2),phi,lc,gaussk=2.0_dp,lcs,gks,pm
    LOGICAL :: edge,edge2
    !,restarted=.FALSE.
    REAL(KIND=dp) :: Coordreflecty,phi2


    CALL nonlocalsurroundelems(Particles,Mesh,lc)


    Particles % xpic = 0.0_dp
    Particles % dd = 0.0_dp

    lcs = lc*lc
    gks = gaussk*gaussk
    pm = gks/lcs


    edge = .FALSE.
    edge2 = .FALSE.

    ! IF (Particles % larcmelfracmod) THEN
    ! DO No = 1, Particles % NumberOfParticles
    !    IF (Particles % damstatus(No) .EQ. -1 .AND. &
    !         Particles % Dav(No,1) == Particles % initdmax) THEN
    !       Particles % DamStatus(No) = 2
    !    END IF
    ! END DO
    ! END IF


    DO ii = 1,smoothiters
       DO No = 1, Particles % NumberOfParticles

          IF (Particles % damstatus(No)<1) THEN
             Particles % dD(No,3,1) = Particles % H(No)
             CYCLE
          END IF

          Particles % dD(No,3:4,1) = 0.0_dp

          ElementInd = Particles % InterpElem(No)
          Coord = Particles % Coordinate(No,1:2)


          DO kk = 1,ElemParticles(ElementInd) % numNLSelems
             ind2 = ElemParticles(ElementInd) % NLSelems(kk)


             IF (ind2 == 0) CYCLE


             DO jj = 1,ElemParticles(ind2) % NumberOfParticles
                CompNo = ElemParticles(ind2) % p(jj)

                IF (CompNo==0) CYCLE

                !already taken care of
                IF (Particles % xpic(CompNo,1) == No) CYCLE

                Particles % xpic(CompNo,1) = No

                CompCoord = Particles % Coordinate(CompNo,1:2)
                phi = (Coord(1)-CompCoord(1))*(Coord(1)-CompCoord(1)) + &
                     (Coord(2)-CompCoord(2))*(Coord(2)-CompCoord(2))

                IF (Particles % mismip) THEN

                   edge = .FALSE.
                   edge2 = .FALSE.

                   IF (Coord(2) <= lc .OR. Coord(2) >= 40000.0_dp-lc) THEN
                      !potential edge effects
                      edge = .TRUE.
                      IF (Coord(2) <= lc) THEN
                         Coordreflecty = -Coord(2)
                      ELSE
                         Coordreflecty = 80000.0_dp-Coord(2)
                      END IF
                   END IF

                   IF (edge) THEN

                      IF (CompCoord(2) <= lc .OR. CompCoord(2) >= 40000.0_dp-lc) THEN

                         !potential edge effects
                         phi2 = (Coord(1)-CompCoord(1))*(Coord(1)-CompCoord(1)) + &
                              (Coordreflecty-CompCoord(2))*(Coordreflecty-CompCoord(2))

                         IF (phi2 < lcs) THEN
                            edge2 = .TRUE.
                         END IF
                      END IF
                   END IF
                END IF


                IF (phi>lcs .AND. (.NOT. edge2)) THEN
                   CYCLE
                ELSE
                   IF (phi<=lcs) THEN
                      ! phi = exp(-(gaussk*gaussk*phi/(lcs)))
                      phi = exp(-pm*phi)
                      Particles % dD(No,3,1) = Particles % dD(No,3,1) + &
                           phi*Particles % H(CompNo)
                      Particles % dD(No,4,1) = Particles % dD(No,4,1) + phi
                   END IF

                   IF (edge2) THEN
                      ! phi2 = exp(-(gaussk*gaussk*phi2/(lc*lc)))
                      phi2 = exp(-pm*phi)

                      Particles % dD(No,3,1) = Particles % dD(No,3,1) + &
                           phi2*Particles % H(CompNo)

                      Particles % dD(No,4,1) = Particles % dD(No,4,1) + phi2
                   END IF
                END IF
             END DO
          END DO

          IF (Particles % dD(No,4,1)>0.0_dp) THEN

             Particles % dD(No,3,1) = Particles % dD(No,3,1)/Particles % dD(No,4,1)

          ELSE

             Particles % dD(No,3,1) = Particles % H(No)
          END IF

       END DO

       Particles % H(:) = Particles % dD(:,3,1)
       Particles % dD(:,3:4,1) = 0.0_dp
       Particles % xpic(:,1) = 0.0_dp
    END DO


    ! IF (Particles % larcmelfracmod) THEN
    !    WHERE (Particles % damstatus(:) > 1) Particles % damstatus(:) = -1
    ! END IF


  END SUBROUTINE smoothrupth

  !**************************************************************************

  SUBROUTINE GetElemParticles_GIMPM( Particles, Model )

    implicit none

    TYPE(Particle_t), POINTER :: Particles
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Model_t) :: Model
    TYPE(Element_t), POINTER :: Element, Edge,NextElement
    TYPE(Variable_t), POINTER :: PM
    INTEGER :: in(9), Elems(9),ElemVec(4),change,count,OldNoParticles
    INTEGER :: L,R,U,D,ii, nn, jj, kk,mm,pp,ne, No, Status,ElementInd,&
         ElementIndNew,ElementIndNew2,&
         istat,maxindex,tot,AddDam,whichelem,surroundelem,fulldamcount,nind,iter
    INTEGER, POINTER :: NodeIndexes(:)=>NULL(),NextNodeIndexes(:)=>NULL(),&
         PMPerm(:)=>NULL(),NodeIndexes2(:)=>NULL()
    TYPE(Nodes_t) :: ElementNodes
    INTEGER, allocatable :: elemcount(:),fulldam(:),frontedit(:),Old(:),Perm(:)
    INTEGER, POINTER :: N1, N2
    REAL(KIND=dp) :: LX,LY,N,S,E,W,xmin,xmax,ymin,ymax,yy,xx,AddVol,ElemVol,SumVol, &
         CriticalDamage, xc, yc, gridres,dmax,tracker,distance,mindist,nx,ny,&
         vol,damvol,elemdamav,damh,frac,CriticalDav,anisoparam,detj,scale
    REAL(KIND=dp) :: Coord(3),GlobalCoords(3),ExtraCoords(3), TestCoords(3)
    REAL(KIND=dp), POINTER :: x1,x2,y1,y2
    REAL(KIND=dp), POINTER :: PMVal(:)=>NULL()
    LOGICAL :: VISITED = .FALSE., GotIt, &
         firsttime,nodamfront,stat
    REAL(KIND=dp),ALLOCATABLE :: Basis(:),dBasisdx(:,:)
    INTEGER :: match,No2
    REAL(KIND=dp) :: dist
    LOGICAL :: CalvingFront
    INTEGER :: noofn
    TYPE(Element_t), POINTER :: P1,P2,BoundaryElement
    INTEGER,ALLOCATABLE :: ptable(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

    SAVE :: VISITED,Mesh,ne,elemcount,PMPerm,PMVal,PM, CriticalDamage,&
         gridres,CriticalDav,nodamfront,Basis,dBasisdx,ElementNodes,Solvername

    firsttime = .FALSE.

    IF (.NOT. VISITED) THEN

       WRITE(SolverName, '(A)') 'GetElemParticles_GIMPM'
       ALLOCATE(Basis(4),dBasisdx(4,3))
       ALLOCATE(ElementNodes % x(4),ElementNodes % y(4),ElementNodes % z(4))

       Mesh => GetMesh()
       ne = Mesh % NumberOfBulkElements

       nodamfront = Particles % nodamfront
       CriticalDamage = Particles % CriticalDamage
       CriticalDav = Particles % criticaldav
       gridres = Particles % gridres
       frac = Particles % elementfraction

       IF (Particles % mix49) THEN
          DO No = 1,Particles % numberofparticles

             IF (Particles % Coordinate(No,1) > 300000.0_dp) THEN
                Particles % Length(No,:) = gridres/3.0_dp
             ELSE
                Particles % Length(No,:) = gridres/2.0_dp
             END IF
          END DO

          Particles % OrigLength = Particles % Length
       ELSE


          IF (frac == 1.0_dp) THEN
             Particles % Length(:,:) = gridres
             Particles % OrigLength(:,:) = gridres
          ELSE IF (frac == 4.0_dp) THEN
             Particles % Length(:,:) = gridres/2.0_dp
             Particles % OrigLength(:,:) = gridres/2.0_dp
          ELSE IF (frac == 9.0_dp) THEN
             Particles % Length(:,:) = gridres/3.0_dp
             Particles % OrigLength(:,:) = gridres/3.0_dp
          ELSE IF (frac == 16.0) THEN
             Particles % Length(:,:) = gridres/4.0_dp
             Particles % OrigLength(:,:) = gridres/4.0_dp
          ELSE
             CALL Fatal(Solvername, &
                  'Particle Element Fraction can currently only be 9,4, or 1')
          END IF
       END IF

       Particles % GVolume(:) = Particles % Length(:,1) * &
            Particles % Length(:,2)
       Particles % pvolume(:) = Particles % GVolume(:)

       CALL Info(Solvername,'allocating elemcount and ElemParticles',Level=4)


       !we restrict particle size to one element in length or less
       !Therefore, a particle can contribute to at most 4 elements
       !there are 9 possible elements that a particle can contribute to

       ! -------------
       ! | 1 | 2 | 3 |
       ! |-----------|
       ! | 4 | 5 | 6 |     !where the particle is in element 5
       ! |-----------|
       ! | 7 | 8 | 9 |
       ! -------------

       IF (ALLOCATED(ElemParticles)) DEALLOCATE(ElemParticles)

       ALLOCATE(ElemParticles(ne))

       IF (ALLOCATED(ElemTrack)) DEALLOCATE(ElemTrack)

       ALLOCATE(ElemTrack(ne))

       ElemTrack(:) % InitialSurface = -1

       !Assign surrounding elements for each element
       !by finding center of those elements and using octree search

       TestCoords = 0.0_dp


       DO ii = 1,ne

          ALLOCATE(ElemParticles(ii) % SurroundingElems(9))

          Element => Mesh % Elements(ii)

          !get center of element
          nn = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
          xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
          ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
          ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))

          xc = (xmin+xmax)/2.0_dp
          yc = (ymin+ymax)/2.0_dp

          N = yc + gridres
          S = yc - gridres
          E = xc + gridres
          W = xc - gridres

          DO jj = 1,9
             SELECT CASE (jj)
             CASE (1)
                TestCoords(1) = W
                TestCoords(2) = N
             CASE (2)
                TestCoords(1) = xc
                TestCoords(2) = N
             CASE (3)
                TestCoords(1) = E
                TestCoords(2) = N
             CASE (4)
                TestCoords(1) = W
                TestCoords(2) = yc
             CASE (5)
                ElemParticles(ii) % SurroundingElems(jj) = ii
                CYCLE
             CASE (6)
                TestCoords(1) = E
                TestCoords(2) = yc
             CASE (7)
                TestCoords(1) = W
                TestCoords(2) = S
             CASE (8)
                TestCoords(1) = xc
                TestCoords(2) = S
             CASE (9)
                TestCoords(1) = E
                TestCoords(2) = S
             END SELECT

             ElementIndNew = 0

             CALL LocateParticleInMeshOctree( ElementIndNew, TestCoords )

             IF (ElementIndNew > 0) THEN
                ElemParticles(ii) % SurroundingElems(jj) = ElementIndNew
             ELSE
                ElemParticles(ii) % SurroundingElems(jj) = 0
             END IF
          END DO
       END DO

       IF (ALLOCATED(elemcount)) DEALLOCATE (elemcount)
       ALLOCATE(elemcount(ne))

       VISITED = .TRUE.

       firsttime = .TRUE.

       WHERE (Particles % ElementIndex(:) == 0) Particles % Coordinate(:,1) = -HUGE(1.0_dp)

       !Particles % Coordinate(:,1) = Particles % Coordinate(:,1) + 225.0_dp
    END IF


    elemcount = 0
    ElemTrack(:) % Volume = 0
    ElemTrack(:) % Status = EMPTY

    PRINT *,'Assigning particles to elements'
    PRINT *,'no of particles: ', Particles % NumberOfParticles


    PM => VariableGet( Mesh % Variables, 'surface')
    IF (.NOT. ASSOCIATED(PM)) Call Fatal (Solvername,'passive mask does not exist')
    PMVal => PM % Values
    PMPerm => PM % Perm

    PMVal = -1.0_dp

    fulldamcount = 0
    !Assign particle to new element
    !Set ElemTrack % Status


    DO No = 1, Particles % NumberOfParticles

       Particles % Status(No) = PARTICLE_ACTIVE
       Coord = GetParticleCoord( Particles, No)

       ElementInd = GetParticleElement( Particles, No )


       IF (ElementInd == 0) THEN
          ElementIndNew = ElementInd
          CALL LocateParticleInMeshOctree( ElementIndNew, Coord(1:3) )

          IF (ElementIndNew > 0) THEN
             !double check
             Element => Mesh % Elements(ElementIndNew)
             nn = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
             xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
             ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
             ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
             IF ((Coord(1) <= xmax) .AND. (Coord(1)>=xmin) .AND. (Coord(2)>=ymin) .AND. (Coord(2)<=ymax)) THEN
                ElementIndNew = ElementIndNew
             ELSE
                ElementIndNew = 0
             END IF
          ELSE
             ElementIndNew = 0
          END IF

          !
       ELSE


          Element => Mesh % Elements(ElementInd)

          !test if particle is still in that element
          nn = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
          xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
          ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
          ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))

          !only test if particle is in surrounding elements if not in previous elements
          IF ((Coord(1) <= xmax) .AND. (Coord(1)>=xmin) .AND. (Coord(2)>=ymin) .AND. (Coord(2)<=ymax)) THEN
             ElementIndNew = ElementInd
          ELSE

             !these should all be positive if the particle is in the element
             !example: E represents the difference between xmax and the particle x coord
             !if positive, the particle is not in surrounding elements 3,6,or 9.
             !continue until the correct element is found

             N = ymax - Coord(2)
             S = Coord(2) - ymin
             E = xmax - Coord(1)
             W = Coord(1) - xmin

             !we create test cases so that
             !we can which surrounding element (whichelem) we are in similarly to using matrix indices
             !relative to the previous element, which is spot 22

             !whichelem equals one of these spots:

             !     11 - 12 - 13
             !     21 - 22 - 23
             !     31 - 32 - 33

             whichelem = 22

             IF (N < 0.0_dp) whichelem = whichelem-10
             IF (S < 0.0_dp) whichelem = whichelem+10
             IF (W < 0.0_dp) whichelem = whichelem-1
             IF (E < 0.0_dp) whichelem = whichelem+1

             IF (N < (0.0_dp - gridres)) whichelem = 99
             IF (S < (0.0_dp - gridres)) whichelem = 99
             IF (E < (0.0_dp - gridres)) whichelem = 99
             IF (W < (0.0_dp - gridres)) whichelem = 99

             SELECT CASE (whichelem)
             CASE (11)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(1)
             CASE (12)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(2)
             CASE (13)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(3)
             CASE (21)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(4)
             CASE (22)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(5)
             CASE (23)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(6)
             CASE (31)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(7)
             CASE (32)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(8)
             CASE (33)
                ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(9)
             CASE (99)
                ElementIndNew = ElementInd
                CALL LocateParticleInMeshOctree( ElementIndNew, Coord )
             CASE DEFAULT
                ElementIndNew = ElementInd
                CALL LocateParticleInMeshOctree( ElementIndNew, Coord )
             END SELECT

             IF (ElementIndNew > 0) THEN
                !double check
                Element => Mesh % Elements(ElementIndNew)
                nn = Element % TYPE % NumberOfNodes
                NodeIndexes => Element % NodeIndexes
                xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
                xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
                ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
                ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
                IF ((Coord(1) <= xmax) .AND. (Coord(1)>=xmin) .AND. (Coord(2)>=ymin) .AND. (Coord(2)<=ymax)) THEN
                   ElementIndNew = ElementIndNew
                ELSE
                   ElementIndNew = 0
                END IF
             ELSE
                ElementIndNew = 0
             END IF

             !
          END IF
       END IF

       Particles % ElementIndex(No) = ElementIndNew
    END DO


    DO ii = 1,2

       DO No = 1, Particles % NumberOfParticles

          ElementInd = GetParticleElement( Particles, No )

          Coord = GetParticleCoord( Particles, No)

          IF (Particles % Status(No) < PARTICLE_LOST) THEN

             LX = Particles % Length(No,1)
             LY = Particles % Length(No,2)

             N = Coord(2) + (LY/2.0_dp)
             S = Coord(2) - (LY/2.0_dp)
             E = Coord(1) + (LX/2.0_dp)
             W = Coord(1) - (LX/2.0_dp)

             in(1:9) = 0

             !the corresponding elements
             Elems(1:9) = 0
             Elems(5) = ElementInd

             IF (ElementInd > 0)  THEN

                Element => Mesh % Elements(ElementInd)

                in(5) = 1

                nn = Element % TYPE % NumberOfNodes
                NodeIndexes => Element % NodeIndexes
                xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
                xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))

                ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
                ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))

                IF(N > ymax) in(2)=1
                IF(S < ymin) in(8)=1
                IF(E > xmax) in(6)=1
                IF(W < xmin) in(4)=1

                IF (in(2) ==1 .AND. in(4)==1) in(1)=1
                IF (in(2) ==1 .AND. in(6)==1) in(3)=1
                IF (in(8) ==1 .AND. in(4)==1) in(7)=1
                IF (in(8) ==1 .AND. in(6)==1) in(9)=1
             END IF


             IF (Elems(5) > 0) THEN

                DO jj = 1,9
                   IF (in(jj) == 1) THEN
                      surroundelem = ElemParticles(ElementInd) % SurroundingElems(jj)

                      IF (surroundelem > 0) THEN

                         Elems(jj) = ElemParticles(ElementInd) % SurroundingElems(jj)
                      ELSE
                         Elems(jj) = 0
                      END IF
                   END IF
                END DO

             ELSE

                !Test if particles outside the domain are overlapping and elements

                DO jj = 1,9

                   GlobalCoords(1:3) = 0.0_dp
                   ElementIndNew = ElementInd

                   SELECT CASE (jj)
                   CASE (1)
                      GlobalCoords(1) = W
                      GlobalCoords(2) = N
                   CASE (2)
                      GlobalCoords(1) = Coord(1)
                      GlobalCoords(2) = N
                   CASE (3)
                      GlobalCoords(1) = E
                      GlobalCoords(2) = N
                   CASE (4)
                      GlobalCoords(1) = W
                      GlobalCoords(2) = Coord(2)
                   CASE (5)
                      CYCLE
                   CASE (6)
                      GlobalCoords(1) = E
                      GlobalCoords(2) = Coord(2)
                   CASE (7)
                      GlobalCoords(1) = W
                      GlobalCoords(2) = S
                   CASE (8)
                      GlobalCoords(1) = Coord(1)
                      GlobalCoords(2) = S
                   CASE (9)
                      GlobalCoords(1) = E
                      GlobalCoords(2) = S
                   END SELECT

                   CALL LocateParticleInMeshOctree( ElementIndNew, GlobalCoords )

                   !double check
                   IF (ElementIndNew > 0) THEN
                      Element => Mesh % Elements(ElementIndNew)
                      nn = Element % TYPE % NumberOfNodes
                      NodeIndexes => Element % NodeIndexes
                      xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
                      xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
                      ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
                      ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))

                      IF ((GlobalCoords(1) < xmax) .AND. (GlobalCoords(1)>xmin) &
                           .AND. (GlobalCoords(2)>ymin) .AND. (GlobalCoords(2)<ymax)) THEN
                         ElementIndNew = ElementIndNew
                      ELSE
                         ElementIndNew = 0
                      END IF
                      IF (ElementIndNew < 1) ElementIndNew = 0
                   ELSE
                      ElementIndNew = 0
                   END IF


                   !assign particle to element, avoiding repeats
                   IF ((ANY(Elems(1:jj) == ElementIndNew)) .OR. (ElementIndNew < 1)) THEN
                      Elems(jj) = 0
                   ELSE
                      Elems(jj) = ElementIndNew
                   END IF
                END DO
             END IF

             IF (Particles % mismipmelt2) THEN
                IF (Particles % time <= 100.0_dp) THEN

                   IF (E>480000.0_dp) THEN
                      Elems(3) = 0
                      Elems(6) = 0
                      Elems(9) = 0
                   END IF

                   IF ( Particles % Coordinate(No,1) > 480000.0_dp) THEN
                      Elems(2) = 0
                      Elems(5) = 0
                      Elems(8) = 0
                   END IF

                   IF (W>480000.0_dp) Elems = 0
                END IF
             END IF

             IF (ii == 1) THEN
                tot = 0
                DO jj = 1,9

                   IF ( Elems(jj) > 0 .AND. Elems(jj) <= ne) THEN
                      tot = tot+1
                      elemcount(Elems(jj)) = elemcount(Elems(jj)) + 1
                   END IF

                   !particle is in an element, but part of its domain isnt
                   IF (in(jj) == 1 .AND. Elems(jj) == 0) THEN
                      Particles % Status(No) = PARTICLE_LEAVING
                   END IF

                END DO

                IF (Elems(5) <= 0) THEN
                   IF (tot<1) THEN
                      Particles % Status(No) = PARTICLE_LOST
                   ELSE
                      Particles % Status(No) = PARTICLE_LEAVING
                   END IF
                END IF

                IF (Particles % Status(No) == PARTICLE_LEAVING) THEN
                   Particles % GVolume(No) = 0.0_dp
                ELSE
                   Particles % GVolume(No) = LX*LY
                END IF


             ELSE
                IF  (Particles % Status(No) == PARTICLE_LEAVING) Particles % GVolume(No) = 0.0_dp

                DO jj = 1,9

                   IF ( Elems(jj) > 0  .AND. Elems(jj) <= ne) THEN

                      elemcount(Elems(jj)) = elemcount(Elems(jj)) + 1
                      ElemParticles(Elems(jj)) % p(elemcount(Elems(jj))) = No

                      Element => Mesh % Elements(Elems(jj))

                      vol = GetParticleVolumeInElement(Particles,No,Element,Model)

                      ElemTrack(Elems(jj)) % Volume = ElemTrack(Elems(jj)) % Volume + vol

                      IF (ElemTrack(Elems(jj)) % Volume >= (Particles % fillpercent * gridres * gridres) ) THEN
                         ElemTrack(Elems(jj)) % Status = FULL
                      ELSE
                         ElemTrack(Elems(jj)) % Status = NOTFULL
                      END IF

                      !add volume for leaving particles
                      IF  (Particles % Status(No) == PARTICLE_LEAVING) THEN
                         Particles % GVolume(No) = Particles % GVolume(No) + vol
                      END IF
                   END IF
                END DO

             END IF
          END IF
       END DO


       IF (ii==1) THEN
          DO jj = 1,ne

             IF (ALLOCATED(ElemParticles(jj)%p)) THEN
                DEALLOCATE(ElemParticles(jj)%p)
             END IF

             ALLOCATE(ElemParticles(jj) % p(MAX(elemcount(jj),1) ))

             ElemParticles(jj) % p = 0

             ElemParticles(jj) % NumberOfParticles = elemcount(jj)

          END DO
          elemcount = 0
       END IF
    END DO


    DO No = 1, Particles % NumberOfParticles
       IF (Particles % GVolume(No) < 0.1_dp) THEN
          Particles % Status(No) = PARTICLE_LOST
       END IF
    END DO



    DO ii = 1,ne

       !first some volume corrections if lost particles ended up in an element
       !this happens when locate particle in mesh octree is a little off and locates
       !a leaving particle as being slightly within an element incorrectly.
       tot = 0

       DO jj = 1,ElemParticles(ii) % NumberOfParticles

          No = ElemParticles(ii) % p(jj)

          IF (Particles % Status(No) == PARTICLE_LOST) THEN
             Element => Mesh % Elements(ii)
             vol = GetParticleVolumeInElement(Particles,No,Element,Model)
             ElemTrack(ii) % Volume = ElemTrack(ii) % Volume - vol

             ElemParticles(ii) % p(jj) = 0
             tot = tot + 1
          END IF
       END DO

       IF (tot > 0) THEN

          IF (tot >= ElemParticles(ii) % NumberOfParticles) THEN
             !no particles, empty the particle list for the element
             ElemParticles(ii) % p(:) = 0
             ElemTrack(ii) % Status = EMPTY
             ElemParticles(ii) % NumberOfParticles = 0
          ELSE

             !adjust/reorder particle list accordingly
             IF (ALLOCATED(Ptable)) DEALLOCATE(Ptable)
             ALLOCATE(Ptable(ElemParticles(ii) % NumberOfParticles))
             Ptable(1:ElemParticles(ii) % NumberOfParticles) = &
                  ElemParticles(ii) % p(1:ElemParticles(ii) % NumberOfParticles)
             ElemParticles(ii) % p(:) = 0
             kk = 0

             DO jj = 1,ElemParticles(ii) % NumberOfParticles
                No = Ptable(jj)
                IF (No > 0) THEN
                   kk = kk+1
                   ElemParticles(ii) % p(kk) = No
                END IF
             END DO
             ElemParticles(ii) % NumberOfParticles = kk
             DEALLOCATE(Ptable)
          END IF
       END IF

       !next, if the element is FULL, set PMVal to 1.

       !If NOTFULL, test if all particles overlap a FULL element.
       !If they do overlap a FULL, then subtract their
       !volume from the NOTFULL element and set them to PARTICLE_LEAVING.
       !If the element isn't next to a full element, it is set to an IGNORE status.
       !This shouldn't really happen if the front is smooth, but
       !for IGNORE, subtract the element from particle volumes, set them to PARTICLE_LEAVING, or
       !particle lost if they have no volume.

       !If they do not overlap a FULL, the element will be an FEM element (regular FEM integration
       !using values interpolated from the mesh).


       IF (ElemTrack(ii) % Status == NOTFULL) THEN

          !make sure the element is next to a FULL
          tot = 0
          DO jj = 1,9

             IF (jj == 5) CYCLE
             ElementIndNew = ElemParticles(ii) % SurroundingElems(jj)

             IF (ElementIndNew <= 0) CYCLE

             IF (ElemTrack(ElementIndNew) % Status == FULL) THEN
                tot = tot + 1
             END IF
          END DO

          !if not next to a FULL, set to IGNORE
          IF (tot == 0) THEN
             ElemTrack(ii) % Status = IGNORE
          ELSE
             !the element is next to a FULL element
             !test if there are any particles that do not overlap a FULL element
             !if all overlap, this will be a passive NOTFULL element, and all
             !particles will become PARTICLE_LEAVING and subtract the volume they have
             !in the element from their GVolume

             !If any particle does not overlap a full element, this will become an FEM element.

             DO jj = 1,ElemParticles(ii) % NumberOfParticles

                No = ElemParticles(ii) % p(jj)

                tot = 0
                DO kk = 1,9

                   ElementIndNew = ElemParticles(ii) % SurroundingElems(kk)

                   IF (ElementIndNew <= 0) CYCLE

                   IF (ElemTrack(ElementIndNew) % Status == FULL)  THEN

                      IF (ANY(ElemParticles(ElementIndNew)%p == No)) THEN
                         tot = 1
                         EXIT
                      END IF
                   END IF
                END DO

                IF (tot == 0) EXIT
             END DO


             IF (tot == 0 .OR. Particles % femforallnotfull) THEN

                IF (tot == 0) THEN
                   ElemTrack(ii) % Status = FEM

                   Element => Mesh % Elements(ii)
                   IF (.NOT. ASSOCIATED(Element) ) CYCLE
                   NodeIndexes => Element % NodeIndexes

                   PMVal(PMPerm(NodeIndexes)) = 2.0_dp

                   !the following ELSE is new.  It should work out
                   !that now, if this NOTFULL element is at the front
                   !(next to a passive element), it will gain FEM status below.
                   !If it is interior, don't treat as FEM.
                   !Before, without the else, all were set to FEM if notfull.

                   !so if there is an overlap with a full element, then all good.
                   !don't set to FEM.
                ELSE
                   Element => Mesh % Elements(ii)
                   IF (.NOT. ASSOCIATED(Element) ) CYCLE
                   NodeIndexes => Element % NodeIndexes

                   PMVal(PMPerm(NodeIndexes)) = 1.0_dp
                END IF
             END IF
          END IF
       END IF


       IF (ElemTrack(ii) % Status == FULL) THEN

          Element => Mesh % Elements(ii)
          IF (.NOT. ASSOCIATED(Element) ) CYCLE
          NodeIndexes => Element % NodeIndexes
          nn = Element % TYPE % NumberOfNodes

          DO jj = 1,nn
             PMVal(PMPerm(NodeIndexes(jj))) = MAX(1.0_dp,PMVal(PMPerm(NodeIndexes(jj))))
          END DO
       END IF
    END DO



    !now go back and look for elements with any negative passive mask nodes.
    !subtract their contribution to particle volume.
    !Also, if any element has all positive nodes, but is not marked full, it
    !gains FEM status.

    !also set interpelem to true for particles overlapping part of the non-original domain
    Particles % UseInterpElem(:) = .FALSE.

    DO ii = 1,ne

       Element => Mesh % Elements(ii)
       IF (.NOT. ASSOCIATED(Element) ) CYCLE
       NodeIndexes => Element % NodeIndexes

       !set useinterpelem to true if a particle overlaps part of non-original domain
       IF (.NOT. firsttime) THEN
          IF (ElemTrack(ii) % InitialSurface .NE. 1 ) THEN
             DO jj = 1,ElemParticles(ii) % NumberOfParticles
                No = ElemParticles(ii) % p(jj)
                Particles % UseInterpElem(No) = .TRUE.
             END DO
          END IF
       END IF

       IF (ANY(PMVal(PMPerm(NodeIndexes)) < 0.0_dp)) THEN

          DO jj = 1,ElemParticles(ii) % NumberOfParticles

             No = ElemParticles(ii) % p(jj)
             vol = GetParticleVolumeInElement(Particles,No,Element,Model)
             Particles % GVolume(No) = Particles % GVolume(No) - vol
             Particles % Status(No) = PARTICLE_LEAVING


             IF (Particles % GVolume(No) <= 0.0_dp) THEN
                Particles % Status(No) = PARTICLE_LOST
             END IF
          END DO

          !if any positive passive mask vals, then it is next to a front elem
          !if no damage at front, then remove damage from front particles (2 cells deep)

          IF (nodamfront) THEN

             DO jj = 1,9
                ElementIndNew = ElemParticles(ii) % SurroundingElems(jj)

                IF (ElementIndNew <=0) CYCLE
                DO kk = 1,9

                   ElementIndNew2 = ElemParticles(ElementIndNew) % SurroundingElems(kk)

                   IF (ElementIndNew2 <=0) CYCLE

                   DO mm = 1,ElemParticles(ElementIndNew2) % NumberOfParticles

                      No = ElemParticles(ElementIndNew2) % p(mm)
                      Particles % Damage(No,:,:) = 0.0_dp
                      Particles % Dav(No,:) = 0.0_dp
                   END DO
                END DO
             END DO

          END IF


       ELSE


          IF (firsttime) THEN
             IF (ElemParticles(ii) % NumberOfParticles > 0) THEN
                ElemTrack(ii) % InitialSurface = 1
             END IF
          END IF

          IF (ALL(PMVal(PMPerm(NodeIndexes)) >= 1.0_dp)) THEN

             !set to FEM if:
             !doesn't contain any particles
             !any surrounding element is a passive element

             IF (ElemParticles(ii) % NumberOfParticles < 1) THEN
                ElemTrack(ii) % Status = FEM
                PMVal(PMPerm(NodeIndexes)) = 2.0_dp
                CYCLE
             END IF


             DO kk = 1,9

                ElementIndNew = ElemParticles(ii) % SurroundingElems(kk)
                IF (ElementIndNew <= 0) CYCLE

                Element => Mesh % Elements(ElementIndNew)
                IF (.NOT. ASSOCIATED(Element) ) CYCLE
                NodeIndexes2 => Element % NodeIndexes

                IF (ANY(PMVal(PMPerm(NodeIndexes2))<0.0_dp)) THEN

                   ElemTrack(ii) % Status = FEM
                   PMVal(PMPerm(NodeIndexes)) = 2.0_dp
                   EXIT

                END IF
             END DO
          END IF
       END IF
    END DO

    !for a pre-defined calving front (not based on passive nodes),
    !an option to set all front elements to FEM no matter what.
    IF (Particles % alwaysfemfront) THEN

       DO ii = Model % NumberOfBulkElements + 1, &
            Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

          BoundaryElement => Model % Elements(ii)

          Model % CurrentElement => Model % Elements(ii)
          noofn  = GetElementNOFNodes()

          NodeIndexes => BoundaryElement % NodeIndexes(1:noofn)

          DO jj=1,Model % NumberOfBCs
             IF ( BoundaryElement % BoundaryInfo % Constraint == &
                  Model % BCs(jj) % Tag ) THEN

                CalvingFront = ListGetLogical( Model % BCs(jj) % Values, &
                     'Calving Front',GotIt)
                IF (.NOT. GotIt) CYCLE

                IF (CalvingFront) THEN

                   P1 => BoundaryElement % BoundaryInfo % Left
                   P2 => BoundaryElement % BoundaryInfo % Right

                   IF (ASSOCIATED(P1)) THEN
                      IF (ElemTrack(P1 % ElementIndex) % Status > IGNORE) THEN
                         ElemTrack(P1 % ElementIndex) % Status = FEM
                         NodeIndexes => P1 % NodeIndexes
                         PMVal(PMPerm(NodeIndexes)) = 2.0_dp
                      END IF
                   END IF

                   IF (ASSOCIATED(P2)) THEN
                      IF (ElemTrack(P2 % ElementIndex) % Status > IGNORE) THEN
                         ElemTrack(P2 % ElementIndex) % Status = FEM
                         NodeIndexes => P2 % NodeIndexes
                         PMVal(PMPerm(NodeIndexes)) = 2.0_dp
                      END IF
                   END IF
                END IF
             END IF
          END DO
       END DO
    END IF


    DO No=1,Particles % NumberOfParticles
       ii = Particles % ElementIndex(No)

       IF (.NOT. ii <= 0) THEN
          IF (ElemTrack(ii) % InitialSurface == 1) THEN
             Particles % InterpElem(No) = ii
          END IF
       END IF
    END DO


    !delete lost particles:

    OldNoParticles = Particles % NumberOfParticles

    !first, make a perm linking old particle No (the location in Perm)
    !to new No (value nn in perm)
    nn = 0

    ALLOCATE( Perm( Particles % NumberOfParticles))
    Perm = 0

    DO No = 1,Particles % NumberOfParticles

       IF ( ( Particles % Status(No) > PARTICLE_LEAVING ) .OR. &
            ( Particles % Status(No) == PARTICLE_ALLOCATED )) CYCLE

       nn = nn+1
       Perm(No) = nn
    END DO

    !then, delete
    CALL DeleteLostParticles( Particles )


    IF (OldNoParticles .NE. Particles % NumberOfParticles) THEN
       !reassign particles to elements using the perm and their new locs
       DO ii = 1, ne
          nn = 0
          DO jj = 1,ElemParticles(ii) % NumberOfParticles
             No = ElemParticles(ii) % p(jj)

             ElemParticles(ii) % p(jj) = Perm(No)

             IF (Perm(No) > 0) nn = nn+1
          END DO

          IF (ElemParticles(ii) % NumberofParticles .NE. nn) THEN

             IF (nn>0) THEN
                ALLOCATE(Old(nn))

                kk = 0

                DO jj = 1, ElemParticles(ii) % NumberOfParticles

                   IF (ElemParticles(ii) % p(jj) == 0) CYCLE

                   kk = kk + 1

                   Old(kk) = ElemParticles(ii) % p(jj)

                END DO

                IF (kk>0) THEN
                   IF (ALLOCATED(ElemParticles(ii)%p)) THEN
                      DEALLOCATE(ElemParticles(ii)%p)
                   END IF

                   ElemParticles(ii) % NumberOfParticles = kk

                   ALLOCATE(ElemParticles(ii) % p(kk))

                   ElemParticles(ii) % p(1:kk) = Old(1:kk)

                ELSE
                   PRINT *,'WARNING: Element has no particles'
                   ElemTrack(ii) % Status = FEM
                   Element => Mesh % Elements(ii)
                   IF (.NOT. ASSOCIATED(Element) ) CYCLE
                   NodeIndexes => Element % NodeIndexes
                   PMVal(PMPerm(NodeIndexes)) = 2.0_dp
                END IF
                DEALLOCATE(Old)
             ELSE
                PRINT *,'WARNING 2: Element has no particles'
                ElemTrack(ii) % Status = FEM
                Element => Mesh % Elements(ii)
                IF (.NOT. ASSOCIATED(Element) ) CYCLE
                NodeIndexes => Element % NodeIndexes
                PMVal(PMPerm(NodeIndexes)) = 2.0_dp
             END IF
          END IF
       END DO
    END IF

    DEALLOCATE(Perm)

    !save basis functions for each element

    IF (Particles % UseSavedBasis) THEN

       DO ii = 1,ne

          IF (ALLOCATED(ElemParticles(ii) % Basis)) DEALLOCATE(ElemParticles(ii) % Basis)
          IF (ALLOCATED(ElemParticles(ii) % dBasisdx)) DEALLOCATE(ElemParticles(ii) % dBasisdx)

          kk = ElemParticles(ii) % NumberOfParticles
          IF (kk < 1) CYCLE

          ALLOCATE(ElemParticles(ii) % Basis(kk,4))
          ALLOCATE(ElemParticles(ii) % dBasisdx(kk,4,3))

          Element => Model % Mesh % Elements(ii)
          NodeIndexes => Element % NodeIndexes

          CALL GetElementNodes(ElementNodes,Element)

          DO jj = 1,kk

             No = ElemParticles(ii) % p(jj)

             stat = GIMPMElementInfo( 0, Particles, Model,Element, ElementNodes, No, &
                  detJ, scale, .FALSE., Basis,dBasisdx)

             ElemParticles(ii) % Basis(jj,1:4) = Basis(1:4)
             ElemParticles(ii) % dBasisdx(jj,1:4,1:3) = dBasisdx(1:4,1:3)

          END DO
       END DO

    END IF


    firsttime = .FALSE.

    PRINT *,'Particles Assigned To Elements'

  END SUBROUTINE GetElemParticles_GIMPM



  !**************************************************************************

  SUBROUTINE GetElemParticles_sMPM( Particles, Model )

    implicit none

    TYPE(Particle_t), POINTER :: Particles
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Model_t) :: Model
    TYPE(Element_t), POINTER :: Element, Edge,NextElement
    TYPE(Variable_t), POINTER :: PM
    INTEGER :: in(9), Elems(9),ElemVec(4),change,count,OldNoParticles
    INTEGER :: L,R,U,D,ii, nn, jj, kk,mm,pp,ne, No, Status,ElementInd,&
         ElementIndNew,istat,maxindex,tot,AddDam,whichelem,surroundelem,fulldamcount,nind,iter
    INTEGER, POINTER :: NodeIndexes(:)=>NULL(),NextNodeIndexes(:)=>NULL(),PMPerm(:)=>NULL()
    INTEGER, allocatable :: elemcount(:),fulldam(:),frontedit(:),Old(:),Perm(:)
    INTEGER, POINTER :: N1, N2
    REAL(KIND=dp) :: LX,LY,N,S,E,W,xmin,xmax,ymin,ymax,yy,xx,AddVol,ElemVol,SumVol, &
         CriticalDamage, xc, yc, gridres,dmax,tracker,distance,mindist,nx,ny,&
         vol,damvol,elemdamav,damh,frac,CriticalDav,anisoparam
    REAL(KIND=dp) :: Coord(3),GlobalCoords(3),ExtraCoords(3), TestCoords(3)
    REAL(KIND=dp), POINTER :: x1,x2,y1,y2
    REAL(KIND=dp), POINTER :: PMVal(:)=>NULL()
    LOGICAL :: VISITED = .FALSE., GotIt, firsttime,nodamfront,inelement
    LOGICAL :: CalvingFront
    INTEGER :: noofn
    TYPE(Element_t), POINTER :: P1,P2,BoundaryElement
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

    SAVE :: VISITED,Mesh,ne,elemcount,PMPerm,PMVal,PM, CriticalDamage,&
         gridres,CriticalDav,nodamfront,Solvername

    firsttime = .FALSE.

    IF (.NOT. VISITED) THEN

       WRITE(SolverName, '(A)') 'GetElemParticles_sMPM'
       Mesh => GetMesh()
       ne = Mesh % NumberOfBulkElements

       nodamfront = Particles % nodamfront
       CriticalDamage = Particles % CriticalDamage
       CriticalDav = Particles % criticaldav
       !  Dmax = Particles % dmax
       gridres = Particles % gridres
       frac = Particles % elementfraction

       IF (.NOT. Particles % hoop) THEN
       IF (frac == 1.0_dp) THEN
          Particles % Length(:,:) = gridres
          Particles % OrigLength(:,:) = gridres
       ELSE IF (frac == 4.0_dp) THEN
          Particles % Length(:,:) = gridres/2.0_dp
          Particles % OrigLength(:,:) = gridres/2.0_dp
       ELSE IF (frac == 9.0_dp) THEN
          Particles % Length(:,:) = gridres/3.0_dp
          Particles % OrigLength(:,:) = gridres/3.0_dp
       ELSEIF (frac == 16.0_dp) THEN
          Particles % Length(:,:) = gridres/4.0_dp
          Particles % OrigLength(:,:) = gridres/4.0_dp
       ELSE
          CALL Fatal(Solvername, &
               'Particle Element Fraction can currently only be 16, 9,4, or 1')
       END IF

       Particles % GVolume(:) = Particles % Length(:,1) * &
            Particles % Length(:,2)
       Particles % pvolume(:) = Particles % GVolume(:)
    END IF


       CALL Info(Solvername,'allocating elemcount and ElemParticles',Level=4)


       !we restrict particle size to one element in length or less
       !Therefore, a particle can contribute to at most 4 elements
       !there are 9 possible elements that a particle can contribute to

       ! -------------
       ! | 1 | 2 | 3 |
       ! |-----------|
       ! | 4 | 5 | 6 |     !where the particle is in element 5
       ! |-----------|
       ! | 7 | 8 | 9 |
       ! -------------

       IF (ALLOCATED(ElemParticles)) DEALLOCATE(ElemParticles)

       ALLOCATE(ElemParticles(ne))

       IF (ALLOCATED(ElemTrack)) DEALLOCATE(ElemTrack)

       ALLOCATE(ElemTrack(ne))

       ElemTrack(:) % InitialSurface = 1


       !Assign surrounding elements for each element
       !by finding center of those elements and using octree search

       TestCoords = 0.0_dp


       DO ii = 1,ne

          ALLOCATE(ElemParticles(ii) % SurroundingElems(9))


          Element => Mesh % Elements(ii)

          !get center of element
          nn = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          xmin = MINVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
          xmax = MAXVAL(Mesh % Nodes % x(NodeIndexes(1:nn)))
          ymin = MINVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))
          ymax = MAXVAL(Mesh % Nodes % y(NodeIndexes(1:nn)))

          xc = (xmin+xmax)/2.0_dp
          yc = (ymin+ymax)/2.0_dp

          N = yc + gridres
          S = yc - gridres
          E = xc + gridres
          W = xc - gridres

          DO jj = 1,9
             SELECT CASE (jj)
             CASE (1)
                TestCoords(1) = W
                TestCoords(2) = N
             CASE (2)
                TestCoords(1) = xc
                TestCoords(2) = N
             CASE (3)
                TestCoords(1) = E
                TestCoords(2) = N
             CASE (4)
                TestCoords(1) = W
                TestCoords(2) = yc
             CASE (5)
                ElemParticles(ii) % SurroundingElems(jj) = ii
                CYCLE
             CASE (6)
                TestCoords(1) = E
                TestCoords(2) = yc
             CASE (7)
                TestCoords(1) = W
                TestCoords(2) = S
             CASE (8)
                TestCoords(1) = xc
                TestCoords(2) = S
             CASE (9)
                TestCoords(1) = E
                TestCoords(2) = S
             END SELECT

             ElementIndNew = 0

             CALL LocateParticleInMeshOctree( ElementIndNew, TestCoords )

             IF (ElementIndNew > 0) THEN

                ElemParticles(ii) % SurroundingElems(jj) = ElementIndNew

             ELSE
                ElemParticles(ii) % SurroundingElems(jj) = 0

             END IF

          END DO
       END DO

       IF (ALLOCATED(elemcount)) DEALLOCATE (elemcount)
       ALLOCATE(elemcount(ne))

       VISITED = .TRUE.

       firsttime = .TRUE.
    END IF


    elemcount = 0
    ElemTrack(:) % Volume = 0
    ElemTrack(:) % Status = EMPTY

    PRINT *,'Assigning particles to elements'
    PRINT *,'no of particles: ', Particles % NumberOfParticles


    PM => VariableGet( Mesh % Variables, 'surface')
    IF (.NOT. ASSOCIATED(PM)) Call Fatal (Solvername,'passive mask does not exist')
    PMVal => PM % Values
    PMPerm => PM % Perm

    PMVal = 1.0_dp

    fulldamcount = 0
    !Assign particle to new element
    !Set ElemTrack % Status

    DO No = 1, Particles % NumberOfParticles

       Particles % Status(No) = PARTICLE_ACTIVE
       Coord = GetParticleCoord( Particles, No)

       ElementInd = GetParticleElement( Particles, No )

       ! IF (ElementInd > 0) THEN
       !    inelement = coordsinelement(Mesh, ElementInd,Coord)
       !    IF (.NOT. inelement) THEN
       !       DO jj = 1,9
       !          IF (jj == 5) CYCLE
       !          ElementIndNew = ElemParticles(ElementInd) % SurroundingElems(jj)
       !          IF (ElementIndNew > 0) THEN
       !             inelement = coordsinelement(Mesh, ElementIndNew, Coord)
       !             IF (inelement) THEN
       !                ElementInd = ElementIndNew
       !                EXIT
       !             END IF
       !          ELSE
       !             inelement = .FALSE.
       !          END IF
       !          IF ((.NOT. inelement) .AND. (jj == 9)) THEN
       !             CALL LocateParticleInMeshOctree( ElementIndNew, Coord(1:3) )
       !             ElementInd = ElementIndNew
       !          END IF
       !       END DO
       !    END IF
       ! ELSE
       CALL LocateParticleInMeshOctree( ElementInd, Coord(1:3) )
       ! END IF

       IF (ElementInd < 0) THEN
          ElementInd = 0
       END IF

       Particles % ElementIndex(No) = ElementInd

       IF (ElementInd < 1) THEN
          Particles % Status(No) = PARTICLE_LOST
       END IF
    END DO

    ! must delete lost particles before assigning the particles to their new elements
    CALL DeleteLostParticles( Particles )

    ! determine how many particles in each element
    DO No = 1,Particles % NumberOfParticles
       ElementInd = GetParticleElement( Particles, No )

       IF (ElementInd < 1) CYCLE
       elemcount(ElementInd) = elemcount(ElementInd) + 1
    END DO

    ! allocate ElemParticles to make list of particles in each element
    DO jj = 1,ne

       IF (ALLOCATED(ElemParticles(jj)%p)) THEN
          DEALLOCATE(ElemParticles(jj)%p)
       END IF

       ALLOCATE(ElemParticles(jj) % p(MAX(elemcount(jj),1) ))

       ElemParticles(jj) % p = 0
       ElemParticles(jj) % NumberOfParticles = elemcount(jj)

       IF (elemcount(jj) > 0) THEN
          ElemTrack(jj) % Status = FULL
       END IF
    END DO

    ! assign particles to ElemParticle list, determine volume of particles in the element
    elemcount = 0

    DO No = 1, Particles % NumberOfParticles
       ElementInd = GetParticleElement( Particles, No )
       IF (ElementInd < 1) CYCLE
       elemcount(ElementInd) = elemcount(ElementInd) + 1
       ElemTrack(ElementInd) % Volume = ElemTrack(ElementInd) % Volume + Particles % PVolume(No)
       ElemParticles(ElementInd) % p(elemcount(ElementInd)) = No
    END DO

    !Set nodes of elements without particles to PMVal = -1

    DO ii = 1,ne
       IF (ElemTrack(ii) % Status > EMPTY) CYCLE
       Element => Mesh % Elements(ii)
       IF (.NOT. ASSOCIATED(Element) ) CYCLE
       NodeIndexes => Element % NodeIndexes
       PMVal(PMPerm(NodeIndexes)) = -1.0_dp

       IF (firsttime) THEN
          ElemTrack(ii) % InitialSurface = -1
       END IF
    END DO

    IF (Particles % alwaysfemfront) THEN

       DO ii = Model % NumberOfBulkElements + 1, &
            Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

          BoundaryElement => Model % Elements(ii)

          Model % CurrentElement => Model % Elements(ii)
          noofn  = GetElementNOFNodes()

          NodeIndexes => BoundaryElement % NodeIndexes(1:noofn)

          DO jj=1,Model % NumberOfBCs
             IF ( BoundaryElement % BoundaryInfo % Constraint == &
                  Model % BCs(jj) % Tag ) THEN

                CalvingFront = ListGetLogical( Model % BCs(jj) % Values, &
                     'Calving Front',GotIt)
                IF (.NOT. GotIt) CYCLE

                IF (CalvingFront) THEN

                   P1 => BoundaryElement % BoundaryInfo % Left
                   P2 => BoundaryElement % BoundaryInfo % Right

                   IF (ASSOCIATED(P1)) THEN
                      IF (ElemTrack(P1 % ElementIndex) % Status > IGNORE) THEN
                         ElemTrack(P1 % ElementIndex) % Status = FEM
                      END IF
                   END IF

                   IF (ASSOCIATED(P2)) THEN
                      IF (ElemTrack(P2 % ElementIndex) % Status > IGNORE) THEN
                         ElemTrack(P2 % ElementIndex) % Status = FEM
                      END IF
                   END IF
                END IF
             END IF
          END DO
       END DO
    END IF

    !front scheme: set to FEM if there are a mix of -1 and 1 PMVals on the element
    !or if specified using alwaysfemfront

    DO ii = 1,ne
       IF (ElemTrack(ii) % Status == EMPTY .OR. ElemTrack(ii) % Status == FEM) CYCLE

       Element => Mesh % Elements(ii)
       IF (.NOT. ASSOCIATED(Element) ) CYCLE

       NodeIndexes => Element % NodeIndexes

       Model % CurrentElement => Model % Elements(ii)
       noofn  = GetElementNOFNodes()

       IF (ABS(SUM(PMVal(PMPerm(NodeIndexes)))) .NE.  DBLE(noofn)) THEN
          ElemTrack(ii) % Status = FEM
       END IF
    END DO

    DO ii = 1,ne
       IF (ElemTrack(ii) % Status == FEM) THEN
          Element => Mesh % Elements(ii)
          IF (.NOT. ASSOCIATED(Element) ) CYCLE

          NodeIndexes => Element % NodeIndexes
          PMVal(PMPerm(NodeIndexes)) = 2.0_dp
       END IF
    END DO


    IF (firsttime) firsttime = .FALSE.

    DO No=1,Particles % NumberOfParticles
       ii = Particles % ElementIndex(No)
       IF (ElemTrack(ii) % InitialSurface == 1) THEN
          Particles % InterpElem(No) = ii
       END IF
    END DO


    PRINT *,'Particles Assigned To Elements'

  END SUBROUTINE GetElemParticles_sMPM


  !**************************************************************************

  SUBROUTINE XPIC(Particles,Model,m)

    IMPLICIT NONE

    TYPE(Particle_t), POINTER :: Particles
    TYPE(Variable_t), POINTER :: Vstar,Vk,H,VPlus,Vstar1,Vstar2,VelSol,PM
    TYPE(Model_t) :: Model
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: dt
    INTEGER, POINTER :: VstarPerm(:),VkPerm(:),HPerm(:),VPlusPerm(:),&
         Vstar1Perm(:),Vstar2Perm(:),VelSolPerm(:),PMPerm(:)
    REAL(KIND=dp), POINTER :: VstarVal(:),VkVal(:),HVal(:),VPlusVal(:),&
         VStar1Val(:),VStar2Val(:),VelSolVal(:),PMValues(:)
    LOGICAL :: Visited = .FALSE.,stat,GotIt
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: t,r,m,count,ii,jj,kk,nb,nn,No,dim,nodevec(9),aax,aay,minstatus
    REAL(KIND=dp),ALLOCATABLE :: Basis(:),dBasisdx(:,:)
    REAL(KIND=dp) :: sfvec(9), weight,detj,scale,mm,rr,bb
    INTEGER, POINTER :: VstarLocalPerm(:),VkLocalPerm(:),VelSolLocalPerm(:)
    REAL(KIND=dp), POINTER :: VstarLocalField(:,:),VsolLocalField(:,:),adtLocalField(:,:),&
         VklocalField(:,:)

    SAVE :: Visited,nb,nn,dim,VstarLocalPerm,VstarLocalField, &
         VkLocalPerm,VsolLocalField,VkLocalField, &
         Vstar,Vstarperm,VstarVal,&
         Vk,VkPerm,VkVal,H,Hperm,HVal, &
         Vplus, VplusPerm, VPlusVal,&
         Vstar1,Vstar1perm,Vstar1Val,&
         Vstar2,Vstar2perm,Vstar2Val,&
         Basis,dBasisdx,VelSolLocalPerm,adtLocalField, &
         VelSol,VelSolPerm,VelSolVal,ElementNodes,PM,PMValues,PMPerm

    IF (.NOT. Visited) THEN

       nb = Model % Mesh % NumberOfBulkElements
       nn = Model % Mesh % MaxElementNodes

       ALLOCATE(Basis(nn),dBasisdx(nn,3))
       ALLOCATE(ElementNodes % x(nn),ElementNodes % y(nn),ElementNodes % z(nn))

       dim = 2

       ALLOCATE( VstarLocalPerm(nn), VstarLocalField(nn,dim) )
       !  ALLOCATE( VkplusLocalPerm(nn), VkplusLocalField(nn,dim) )
       ALLOCATE( VkLocalPerm(nn), VSolLocalField(nn,dim), VKlocalField(nn,dim) )
       ALLOCATE( VelSolLocalPerm(nn),  adtLocalField(nn,dim) )

       Vstar => VariableGet(Model % Mesh % Variables, 'Vstar' )
       IF (.NOT. ASSOCIATED(Vstar)) CALL Fatal('xpic','Vstar does not exist ')
       VstarPerm => Vstar % Perm
       VstarVal => Vstar % Values
       VstarVal = 0.0_dp


       Vplus => VariableGet(Model % Mesh % Variables, 'Vplus' )
       IF (.NOT. ASSOCIATED(Vplus)) CALL Fatal('xpic','Vplus does not exist ')
       VplusPerm => Vplus % Perm
       VplusVal => Vplus % Values
       VplusVal = 0.0_dp

       VelSol => VariableGet(Model % Mesh % Variables, 'SSAVelocity' )
       IF (.NOT. ASSOCIATED(VelSol)) CALL Fatal('xpic','SSAVelocity does not exist ')
       VelSolPerm => VelSol % Perm
       VelSolVal => VelSol % Values

       Vk => VariableGet(Model % Mesh % Variables, 'PrevVel' )
       IF (.NOT. ASSOCIATED(Vk)) CALL Fatal('xpic','PrevVel does not exist ')
       VkPerm => Vk % Perm
       VkVal => Vk % Values

       H => VariableGet(Model % Mesh % Variables, 'H' )
       IF (.NOT. ASSOCIATED(H)) CALL Fatal('xpic','H does not exist ')
       HPerm => H % Perm
       HVal => H % Values

       Vstar1 => VariableGet(Model % Mesh % Variables, 'Vstar 1' )
       IF (.NOT. ASSOCIATED(Vstar1)) CALL Fatal('xpic','Vstar does not exist ')
       Vstar1Perm => Vstar1 % Perm
       Vstar1Val => Vstar1 % Values

       Vstar2 => VariableGet(Model % Mesh % Variables, 'Vstar 2' )
       IF (.NOT. ASSOCIATED(Vstar2)) CALL Fatal('xpic','Vstar does not exist ')
       Vstar2Perm => Vstar2 % Perm
       Vstar2Val => Vstar2 % Values

       PM => VariableGet( Model % Mesh % Variables, 'surface')
       PMValues => PM % Values
       PMPerm => PM % Perm

       Visited = .TRUE.
    END IF

    CALL Info('xpic','Assigning new velocities and next coords to particles...',Level=1)

    mm = DBLE(m)

    dt = Particles % dtime

    !INITIALIZATION
    !vplus starts (r = 1) by equaling prevvel...
    VplusVal(2*(VplusPerm(:)-1)+1) = VkVal(2*(VkPerm(:)-1)+1)
    VplusVal(2*(VplusPerm(:)-1)+2) = VkVal(2*(VkPerm(:)-1)+2)

    VstarVal = 0.0_dp

    Particles % xpic = 0.0_dp

    !LOOPS
    DO r = 2,m

       !multiply by the r and m function
       rr = DBLE(r)

       !vplus(r-1) from mesh  to particles. 3 is a dumy
       CALL MPMMeshVectorToParticle( Particles, Model, 4, 3)


       !vplus(r-1) from particles to nodes (momentum formulation)
       CALL MPMParticlesToNodes( Particles,Model,4)

       bb = (-1.0_dp)**rr

       rr = (mm-rr+1.0_dp)/rr

       VplusVal(2*(VplusPerm(:)-1)+1) = rr * VplusVal(2*(VplusPerm(:)-1)+1)
       VplusVal(2*(VplusPerm(:)-1)+2) = rr * VplusVal(2*(VplusPerm(:)-1)+2)

       VstarVal(2*(VstarPerm(:)-1)+1) = VstarVal(2*(VstarPerm(:)-1)+1) + bb * VplusVal(2*(VplusPerm(:)-1)+1)
       VstarVal(2*(VstarPerm(:)-1)+2) = VstarVal(2*(VstarPerm(:)-1)+2) + bb * VplusVal(2*(VplusPerm(:)-1)+2)
    END DO



    Particles % NextCoordinate(:,:) = 0.0_dp
    Particles % xpic = 0.0_dp
    Particles % dd = 0.0_dp


    !XPIC may have errors for the partially-filled elements that
    !commonly occur at the ice front. These errors propagate m nodes
    !upstream, and we mark these nodes using the first component of
    !VplusVal, which we set to 1. All other active nodes are set to 0.
    !Later, this marker is interpolated to each particle on dD(:,1,1),
    !which is used to control how much of the XPIC correction a particle receives,
    !with no XPIC correction if the particle is near the ice front, and a full correction
    !far from the ice front.

    VplusVal = 0.0_dp
    VplusVal(2*(VplusPerm(:)-1)+1) = PMValues(PMPerm(:))-1.0_dp

    DO jj = 1,m-1
       DO ii = 1,nb
          IF ( ElemTrack(ii) % Status < FEM ) CYCLE
          Element => Model % Mesh % Elements(ii)
          NodeIndexes => Element % NodeIndexes

          IF (ANY(VplusVal(2*(VplusPerm(NodeIndexes)-1)+1)==1.0_dp)) THEN
             VplusVal(2*(VplusPerm(NodeIndexes)-1)+2) = 1.0_dp
          END IF

       END DO
       VplusVal(2*(VplusPerm(:)-1)+1) = VplusVal(2*(VplusPerm(:)-1)+2)
    END DO



    DO ii = 1,nb

       !only use full elements
       IF ( ElemTrack(ii) % Status < FEM ) CYCLE

       Element => Model % Mesh % Elements(ii)
       NodeIndexes => Element % NodeIndexes

       CALL GetElementNodes(ElementNodes,Element)

       VstarLocalPerm(1:nn) = VstarPerm(Element % NodeIndexes)
       VkLocalPerm(1:nn) = VkPerm(Element % NodeIndexes)
       VelSolLocalPerm(1:nn) = VelSolPerm(Element % NodeIndexes)


       DO kk = 1,nn
          DO jj = 1,2
             VstarLocalField(kk,jj) = VstarVal(2*(VstarLocalPerm(kk)-1)+jj)
             VSolLocalField(kk,jj) = VelSolVal(2*(VelSolLocalPerm(kk)-1)+jj)
             !   VkLocalFied(kk,jj) = VkVal(2*(VkLocalPerm(kk)-1)+jj)

             adtLocalField(kk,jj) = VelSolVal(2*(VelSolLocalPerm(kk)-1)+jj) - &
                  VkVal(2*(VkLocalPerm(kk)-1)+jj)
          END DO
       END DO

       DO t = 1, ABS(ElemParticles(ii) % NumberOfParticles)
          No = ElemParticles(ii) % p(t)

          IF (Particles % ShapeFunctions == 'gimpm') THEN
             stat = GIMPMElementInfo( t,Particles, Model,Element, ElementNodes, No, &
                  detJ, scale, .TRUE., Basis,dBasisdx)

             Basis = Basis*scale
          ELSE

             stat = sMPMElementInfo( Element,Particles, Model, ElementNodes, No, &
                  Particles % gridres, Basis,dBasisdx)
             scale = 1.0_dp
          END IF


          ! IF (ElemTrack(ii) % Status > FEM) THEN
          ! IF (ALL(PMValues(PMPerm(Element % NodeIndexes)) .NE. 2)) THEN
          !S(a_ex^k)
          DO kk = 1,2
             Particles % xpic(No,kk) = Particles % xpic(No,kk) + &
                  mm*SUM(Basis(1:nn) * VstarLocalField(1:nn,kk)) - &
                  mm*SUM(Basis(1:nn) * (VsolLocalField(1:nn,kk) - adtLocalField(1:nn,kk))) + &
                  SUM(Basis(1:nn)) *  Particles % Velocity(No,kk)
          END DO

          !dD(No,1,1) here will have a value of 0 if xpic is fully used or 1 if FLIP if the particle
          !is near the ice front and FLIP is fully used
          Particles % dD(No,1,1) = Particles % dD(No,1,1) + SUM(Basis(1:nn) * VplusVal(2*(VplusPerm(NodeIndexes)-1)+1))

          !PMValues(PMPerm(Element % NodeIndexes)))
          !END IF

          !S(v^k+)
          DO kk = 3,4
             Particles % xpic(No,kk) = Particles % xpic(No,kk) + &
                  SUM(Basis(1:nn) * VsolLocalField(1:nn,kk-2))
          END DO

          !S(v^(k+) - v^(k))
          DO kk = 5,6
             Particles % xpic(No,kk) = Particles % xpic(No,kk) + &
                  SUM(Basis(1:nn) * adtLocalField(1:nn,kk-4))
          END DO
       END DO
    END DO

    !dD(No,1,1) here will have a value of 1 if xpic is fully used or 0 if FLIP if the particle
    !is near the ice front and FLIP is fully used (so the opposite of previously)
    Particles % dD(:,1,1) = ABS(Particles % dD(:,1,1) - 1.0_dp)

    !scale the xpic correction by Particles % dD(:,1,1)
    Particles % xpic(:,1) = Particles % xpic(:,1) * Particles % dD(:,1,1)
    Particles % xpic(:,2) = Particles % xpic(:,2) * Particles % dD(:,1,1)

    !update all particles:
    Particles % Velocity(:,1:2) = Particles % Velocity(:,1:2) + Particles % xpic(:,5:6) - Particles % xpic(:,1:2)
    Particles % NextCoordinate(:,1:2) = Particles % GridVelocity(:,1:2) - &
         0.5_dp*(Particles % xpic(:,5:6) + Particles % xpic(:,1:2))

    Particles % xpic = 0.0_dp

    CALL Info('xpic','Done.',Level=1)

  END SUBROUTINE XPIC

  !**************************************************************************

  SUBROUTINE FixPrincipalDamage(Din,CriticalDamage)

    IMPLICIT NONE

    REAL(KIND=dp) :: Din(4), D(2,2),D3
    REAL(KIND=dp) :: TT,DD,lambda(3),en
    REAL(KIND=dp) :: quart=0.25_dp,half=0.5_dp,thres=0.0001_dp,zero=0.0_dp,one=1.0_dp
    REAL(KIND=dp) :: CriticalDamage,diff1,diff2
    REAL(KIND=dp) :: EigVals(2),EigVec(2,2)
    REAL(KIND=dp) ::w,x,y,z
    INTEGER :: ii
    REAL(KIND=dp) :: layerdmax,layerdmax2,sqrteig
    TYPE(Particle_t), POINTER :: Particles
    LOGICAL :: Visited=.FALSE.

    SAVE :: Particles,Visited,layerdmax


    IF (.NOT. Visited) THEN
       Particles => GlobalParticles
       layerdmax = Particles % DmaxI
    END IF

    !check eigenvalues
    IF (Particles % currentGamma == zero) THEN
       IF (Din(1) >= CriticalDamage) THEN
          Din(1:2) = Particles % DmaxI
          IF (.NOT. Particles % nodzz) THEN
             Din(3) = Particles % DmaxIII
          END IF
       END IF

       IF (Din(1) < zero) Din = zero
       RETURN

    ELSE

       !check eigenvalues: don't compute eigenvectors unless you have to
       TT = Din(1)+Din(2)
       DD = Din(1)*Din(2)-Din(4)*Din(4)

       sqrteig = quart*TT*TT-DD
       IF (sqrteig < 0.0_dp) sqrteig = 0.0_dp
       sqrteig = sqrt(sqrteig)
       TT = half*TT

       lambda(1)=TT+sqrteig
       lambda(2)=TT-sqrteig
       lambda(3) = Din(3)

       IF (.NOT. (ANY(lambda < zero) .OR. ANY(lambda >= CriticalDamage))) THEN
          RETURN
       END IF
    END IF


    !get eigs and eigenvectors
    !this is now the same as Eigen2DSum_TryGenFirst

    D(1,1) = Din(1)
    D(2,1) = Din(4)
    D(1,2) = Din(4)
    D(2,2) = Din(2)

    IF (Din(4)==zero) THEN
       EigVals(1) = lambda(2)
       EigVals(2) = lambda(1)

       IF (Din(2)>Din(1)) THEN
          EigVec(1,1) = one
          EigVec(2,1) = zero
          EigVec(1,2) = zero
          EigVec(2,2) = one
       ELSE
          EigVec(1,1) = zero
          EigVec(2,1) = one
          EigVec(1,2) = one
          EigVec(2,2) = zero
       END IF
    ELSE

       IF (ABS(lambda(1)-lambda(2)) < thres .OR. ABS(Din(4)) < thres ) THEN
          CALL Eigen2D(D,EigVals,EigVec)
       ELSE
          EigVals(1) = lambda(2)
          EigVals(2) = lambda(1)

          !first eigenvector
          EigVec(1,2)=Din(4)
          EigVec(2,2)=lambda(1)-Din(1)

          !2nd eigenvector
          EigVec(1,1)=Din(4)
          EigVec(2,1)=lambda(2)-Din(1)

          DO ii=1,2
             en=SQRT(SUM(EigVec(:,ii)*EigVec(:,ii)))
             EigVec(:,ii)=EigVec(:,ii)/en
          END DO
       END IF
    END IF

    D3 = Din(3)

    IF (EigVals(2) >= CriticalDamage) THEN

       Particles % prupt = .TRUE.

       IF (CriticalDamage == layerdmax) THEN

          !DMAX CORRECTION
          !we are calling this from the rkm scheme:
          !diff 1 is the excess damage accumulated on
          !D1 that exceeds dmax
          !we subtract diff1 from D1.
          !and subtract diff2 from D2 and D3,which
          !is calculated in according to the assumption that
          !damage growth between principal directions is
          !dD2=dD3 = (1-gamma)*dD1.
          !note this assumption doesn't account for
          !spin contributions and misalignment between the effective principal stress,
          !but these effects for just one subestimate in the RKM
          !scheme are negligible (the layer is ruptured
          !and we set damage evolution f = 0 from now on).
          !Accounting for these effects during this correction
          !would be complex.

          diff1 = EigVals(2) - layerdmax
          diff2 = (1.0_dp - Particles % TempGamma) * diff1

          EigVals(2) = EigVals(2) - diff1
          EigVals(1) = EigVals(1) - diff2

          IF (.NOT. Particles % noDzz) THEN
             D3 = D3 - diff2
          END IF
       ELSE
          diff1 = layerdmax - EigVals(2)
          !  IF (diff1 > zero) THEN
          diff2 = (1.0_dp - Particles % TempGamma) * diff1
          EigVals(1) = EigVals(1)+diff2
          IF (.NOT. Particles % noDzz) THEN
             D3 = D3 + diff2
          ENDIF
          ! END IF
          EigVals(2) = layerdmax

          IF (Particles % forcedzz) D3 = Particles % DmaxIII
       END IF


    END IF


    IF (EigVals(1) > Particles % DmaxII) EigVals(1) = Particles % DmaxII
    IF (D3 > Particles % DmaxIII) D3 = Particles % DmaxIII


    IF (EigVals(1)< zero) EigVals(1) = zero
    IF (EigVals(2)< zero) EigVals(2) = zero
    IF (D3 < zero) D3 = zero


    ! D = zero; D(1,1) = EigVals(1); D(2,2) = EigVal(2)
    ! D = MATMUL(EigVec,D); D = MATMUL(D2,TRANSPOSE(EigVec))
    !Din(1) = D(1,1); Din(2) = D(2,2); Din(3) = D(3,3); Din(4) = D(1,2)

    !same as above, but faster:
    w = EigVals(1)*EigVec(1,1)
    x = EigVals(2)*EigVec(1,2)
    y = EigVals(1)*EigVec(2,1)
    z = EigVals(2)*EigVec(2,2)

    Din(1) = EigVec(1,1)*w + EigVec(1,2)*x
    Din(2) = EigVec(2,1)*y + EigVec(2,2)*z
    Din(3) = D3
    Din(4) = EigVec(2,1)*w + EigVec(2,2)*x

  END SUBROUTINE FixPrincipalDamage

  !**************************************************************************

  SUBROUTINE FixPrincipalDamageInc(Din,dD,r)

    IMPLICIT NONE

    REAL(KIND=dp) :: Din(4), D(2,2),D3,dD(4),r(4)
    REAL(KIND=dp) :: TT,DDD,lambda(3),en,scale
    REAL(KIND=dp) :: quart=0.25_dp,half=0.5_dp,thres=0.0001_dp,zero=0.0_dp,one=1.0_dp
    REAL(KIND=dp) :: CriticalDamage,diff1,diff2,discrim
    REAL(KIND=dp) :: EigVals(2),EigVec(2,2)
    REAL(KIND=dp) ::w,x,y,z,a,b,c
    INTEGER :: ii
    REAL(KIND=dp) :: Dmax,dmax2,dmax3,sqrteig
    TYPE(Particle_t), POINTER :: Particles
    LOGICAL :: Visited=.FALSE.,rupt

    SAVE :: Particles,Visited,Dmax,Dmax2,Dmax3


    IF (.NOT. Visited) THEN
       Particles => GlobalParticles
       dmax = Particles % DmaxI
       dmax2 = Particles % DmaxII
       dmax3 = Particles % DmaxIII
       Visited = .TRUE.
    END IF


    !check eigenvalues
    IF (Particles % currentGamma == zero) THEN
       IF (Din(1) >= Dmax) THEN
          Din(1:2) = Particles % DmaxI
          IF (.NOT. Particles % nodzz) THEN
             Din(3) = Particles % DmaxIII
          END IF
       END IF

       IF (Din(1) < zero) Din = zero

       Particles % bmd = Din(1)
       RETURN

    ELSE

       !check eigenvalues: don't compute eigenvectors unless you have to
       TT = Din(1)+Din(2)
       DDD = Din(1)*Din(2)-Din(4)*Din(4)

       sqrteig = quart*TT*TT-DDD
       IF (sqrteig < 0.0_dp) sqrteig = 0.0_dp
       sqrteig = sqrt(sqrteig)
       TT = half*TT

       lambda(1)=TT+sqrteig
       lambda(2)=TT-sqrteig
       lambda(3) = Din(3)


       Particles % bmd = MIN(lambda(1),Dmax)
       Particles % bmd = MAX(0.0_dp,Particles % bmd)

       IF (.NOT. (lambda(1) > Dmax .OR. ANY(lambda<zero)) ) THEN
          RETURN
       END IF

       IF (lambda(1) > Dmax) THEN
          rupt = .TRUE.
       ELSE
          rupt = .FALSE.
       END IF
    END IF

    IF (Particles % prupt) THEN

       !we have already ruptured, and are simply evolving via spin
       !so just correct principal values of damage so that
       !they do not exceed dmax

       D(1,1) = Din(1)
       D(2,1) = Din(4)
       D(1,2) = Din(4)
       D(2,2) = Din(2)

       IF (Din(4)==zero) THEN
          EigVals(1) = lambda(2)
          EigVals(2) = lambda(1)
          IF (Din(2)>Din(1)) THEN
             EigVec(1,1) = one
             EigVec(2,1) = zero
             EigVec(1,2) = zero
             EigVec(2,2) = one
          ELSE
             EigVec(1,1) = zero
             EigVec(2,1) = one
             EigVec(1,2) = one
             EigVec(2,2) = zero
          END IF
       ELSE

          IF (ABS(lambda(1)-lambda(2)) < thres .OR. ABS(Din(4)) < thres ) THEN
             !Use dsyev to be sure...general solution might run into
             !precision/overflow/underflow issues under these conditions
             CALL Eigen2D(D,EigVals,EigVec)
          ELSE
             ! CALL Eigen2D(D,EigValues2d,EigenVec)
             !general solution (returning with smallest eigvalues first)
             EigVals(1) = lambda(2)
             EigVals(2) = lambda(1)

             !first eigenvector
             EigVec(1,2)=Din(4)
             EigVec(2,2)=lambda(1)-Din(1)

             !2nd eigenvector
             EigVec(1,1)=Din(4)
             EigVec(2,1)=lambda(2)-Din(1)

             DO ii=1,2
                en=SQRT(SUM(EigVec(:,ii)*EigVec(:,ii)))
                EigVec(:,ii)=EigVec(:,ii)/en
             END DO
          END IF
       END IF

       IF (EigVals(1) > Dmax2) EigVals(1) = Dmax2
       IF (EigVals(2) > Dmax) EigVals(2) = Dmax
       IF (Din(3) > Dmax3) Din(3) = Dmax3
       IF (EigVals(1) < zero) EigVals(1) = zero
       IF (EigVals(2) < zero) EigVals(2) = zero
       IF (Din(3) < zero) Din(3) = zero

       w = EigVals(1)*EigVec(1,1)
       x = EigVals(2)*EigVec(1,2)
       y = EigVals(1)*EigVec(2,1)
       z = EigVals(2)*EigVec(2,2)

       Din(1) = EigVec(1,1)*w + EigVec(1,2)*x
       Din(2) = EigVec(2,1)*y + EigVec(2,2)*z
       Din(4) = EigVec(2,1)*w + EigVec(2,2)*x

    ELSE


       IF (rupt) THEN

          Din = Din - dD

          a = dD(1)*dD(2)-dD(4)*dD(4)
          b = -Dmax*(dD(1)+dD(2)) + dD(2)*Din(1) + dD(1)*Din(2) - 2.0_dp*Din(4)*dD(4)
          c = Dmax*Dmax - Dmax*(Din(1)+Din(2))+Din(1)*Din(2)-Din(4)*Din(4)

          IF (a == 0.0_dp) THEN
             IF (b .NE. 0.0_dp) THEN
                scale = -c/b
                IF (scale<0.0_dp) scale = 0.0_dp
                IF (scale>1.0_dp) scale = 1.0_dp
             ELSE
                scale = 0.0_dp
             END IF
          ELSE

             discrim = b*b-4.0_dp*a*c
             IF (discrim > 0.0_dp) THEN
                discrim = sqrt(b*b-4.0_dp*a*c)
             ELSE
                IF (discrim < -1.0e-3_dp) THEN
                   PRINT *,'discrim warning',discrim
                END IF

                discrim = 0.0_dp
             END IF

             scale = (-b-discrim)/(2.0_dp*a)

             IF (scale < zero .OR. scale > 1.0_dp) THEN
                scale = (-b+discrim)/(2.0_dp*a)
                !    PRINT *,'scale 2',scale
                IF (scale < zero .OR. scale > 1.0_dp) THEN
                   scale = 0.0_dp
                END IF
             END IF
          END IF

          IF (scale .NE. scale) scale = 0.0_dp

          Din = Din + dD*scale

          Particles % prupt = .TRUE.

          !for rest of step, let it continue rotating,
          !from spin tensor
          Din = Din + r*(1.0_dp-scale)

          TT = Din(1)+Din(2)
          DDD = Din(1)*Din(2)-Din(4)*Din(4)

          sqrteig = quart*TT*TT-DDD
          IF (sqrteig < 0.0_dp) sqrteig = 0.0_dp
          sqrteig = sqrt(sqrteig)
          TT = half*TT

          lambda(1)=TT+sqrteig
          lambda(2)=TT-sqrteig
          lambda(3) = Din(3)

       END IF

       D(1,1) = Din(1)
       D(2,1) = Din(4)
       D(1,2) = Din(4)
       D(2,2) = Din(2)

       IF (Din(4)==zero) THEN
          EigVals(1) = lambda(2)
          EigVals(2) = lambda(1)
          IF (Din(2)>Din(1)) THEN
             EigVec(1,1) = one
             EigVec(2,1) = zero
             EigVec(1,2) = zero
             EigVec(2,2) = one
          ELSE
             EigVec(1,1) = zero
             EigVec(2,1) = one
             EigVec(1,2) = one
             EigVec(2,2) = zero
          END IF
       ELSE

          IF (ABS(lambda(1)-lambda(2)) < thres .OR. ABS(Din(4)) < thres ) THEN
             !Use dsyev to be sure...general solution might run into
             !precision/overflow/underflow issues under these conditions
             CALL Eigen2D(D,EigVals,EigVec)
          ELSE
             ! CALL Eigen2D(D,EigValues2d,EigenVec)
             !general solution (returning with smallest eigvalues first)
             EigVals(1) = lambda(2)
             EigVals(2) = lambda(1)

             !first eigenvector
             EigVec(1,2)=Din(4)
             EigVec(2,2)=lambda(1)-Din(1)

             !2nd eigenvector
             EigVec(1,1)=Din(4)
             EigVec(2,1)=lambda(2)-Din(1)

             DO ii=1,2
                en=SQRT(SUM(EigVec(:,ii)*EigVec(:,ii)))
                EigVec(:,ii)=EigVec(:,ii)/en
             END DO
          END IF
       END IF

       IF (EigVals(1) > Dmax2) EigVals(1) = Dmax2
       IF (EigVals(2) > Dmax) EigVals(2) = Dmax
       IF (Din(3) > Dmax3) Din(3) = Dmax3
       IF (EigVals(1) < zero) EigVals(1) = zero
       IF (EigVals(2) < zero) EigVals(2) = zero
       IF (Din(3) < zero) Din(3) = zero

       w = EigVals(1)*EigVec(1,1)
       x = EigVals(2)*EigVec(1,2)
       y = EigVals(1)*EigVec(2,1)
       z = EigVals(2)*EigVec(2,2)

       Din(1) = EigVec(1,1)*w + EigVec(1,2)*x
       Din(2) = EigVec(2,1)*y + EigVec(2,2)*z
       Din(4) = EigVec(2,1)*w + EigVec(2,2)*x

    END IF


  END SUBROUTINE FixPrincipalDamageInc

  !**************************************************************************

  SUBROUTINE bassisinc(Particles,D,RHS,f,dd)

    IMPLICIT NONE
    REAL(KIND=dp) :: f12rhs,f(4),D(4),dd(4)
    REAL(KIND=dp) :: mat(2,2),EigVals(2),EigenVec(2,2)
    REAL(KIND=dp) :: ww,xx,yy,zz,RHS
    REAL(KIND=dp) :: zero=0.0_dp,half=0.0_dp
    TYPE(Particle_t), POINTER :: Particles

    IF (Particles % Gamma > zero) THEN
       f12rhs = D(4) * (Particles % dvdxmdudy)
       !spin contribution
       f(1) = -f12rhs
       f(2) = f12rhs
       f(3) = zero
       f(4) = -half*(Particles % dvdxmdudy)*(D(2)-D(1))
    ELSE
       f = zero
    END IF

    IF (Particles % prupt) THEN
       dd = 0.0_dp
       RETURN
    END IF

    IF (Particles % gamma > zero) THEN

       mat(1,1) = D(1)
       mat(2,2) = D(2)
       mat(1,2) = D(4)
       mat(2,1) = D(4)

       CALL Eigen2DSym_TryGenFirst(mat,EigVals,EigenVec)

       EigVals(1) = zero
       EigVals(2) = RHS*EigVals(2)

       ww = EigVals(1)*EigenVec(1,1)
       xx = EigVals(2)*EigenVec(1,2)
       yy = EigVals(1)*EigenVec(2,1)
       zz = EigVals(2)*EigenVec(2,2)

       dD(1) = EigenVec(1,1)*ww + EigenVec(1,2)*xx
       dD(2) = EigenVec(2,1)*yy + EigenVec(2,2)*zz
       dD(3) = EigVals(1)
       dD(4) = EigenVec(2,1)*ww + EigenVec(2,2)*xx

       dD = dD + f

    ELSE
       dD(1:3) = RHS*D(1:3)
       dD(4) = 0.0_dp
    END IF

  END SUBROUTINE bassisinc

  !**************************************************************************


  !> runge-kutta-merson method to solve a system of 1st order initial value ODEs
  !! (here, SSA ice creep damage from Huth and others, 2020)
  !! integrates from t=a to t=b with the initial conditions in y
  !! Step length automatically adjusts so that absolute error estimate <= tol
  SUBROUTINE runge_kutta_merson(Particles,Din,dD,a,b,tol,ifail,ddscale)

    !edited from:
    !Chivers and Sleightholme
    !Introduction to Programming with Fortran. 2015.
    !ch2602_rkm_module.f90
    !http://rhymneyconsulting.co.uk/fortran/third_edition/ch2602_rkm_module.f90

    !This version is edited for solving ice damage (Pralong, Hutter, and Funk 2006)

    IMPLICIT NONE

    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp), INTENT (out), DIMENSION (4) :: dD
    REAL(KIND=dp), INTENT (in) :: a, b, tol
    REAL(KIND=dp), INTENT (in) :: Din(4)
    INTEGER, INTENT (out) :: ifail
    REAL(KIND=dp), INTENT (out) :: ddscale
    REAL(KIND=dp), DIMENSION (4) :: s1, s2, s3, s4, s5,&
         r1,r2,r3,r4,r5,new_D_1, old_D_1,new_D_2, error,Dinit,D,inc,rinc
    REAL(KIND=dp), DIMENSION (2,2) :: D2, EigVec
    REAL(KIND=dp), DIMENSION (2) :: EigVal
    REAL(KIND=dp) :: t, hh, hh2, hh3, hh6, hh8, w,x,y,z,factor
    REAL(KIND=dp) :: max_error, smallest_step = 1.e-15_dp,rupttime
    INTEGER :: k,no_of_steps = 0,stoptimes = 1,ii
    LOGICAL :: Visited=.false.,torupt,rupted = .FALSE.
    REAL(KIND=dp) :: TT,DDD,sqrteig,en,EigVec2(2,2)
    REAL(KIND=dp) :: Dmax,critdam,firststep,div2,div3,div6,div8

    SAVE :: Visited, DMax, critdam,firststep,div2,div3,div6,div8,factor

    IF (.NOT. Visited) THEN

       Dmax = Particles % Dmax

       !if damage > criticaldamage during RKM, set to Dmax?
       !oe just constrain damage to be <= Dmax?
       IF (Particles % RKMCritDam) THEN
          critdam = Particles % criticaldamage
       ELSE
          critdam = Particles % DmaxI
       END IF

       firststep = 1.0_dp/100.0_dp

       div2 = 1.0_dp/2.0_dp
       div3 = 1.0_dp/3.0_dp
       div6 = 1.0_dp/6.0_dp
       div8 = 1.0_dp/8.0_dp

       factor = 1.e-2_dp !1.0_dp/32.0_dp !

       Visited = .TRUE.
    END IF

    no_of_steps = 0

    ifail = 0
    dD = 0

    ! check input parameters
    IF (a==b .OR. tol<=0.0) THEN
       ifail = 1
       ddscale = 1.0_dp
       RETURN
    END IF

    t = a
    hh = (b-a)*firststep
    Dinit = Din
    D = Din

    IF (Particles % prupt) THEN
       torupt = .FALSE.
       rupted = .TRUE.
    ELSE
       torupt = .TRUE.
       rupted = .FALSE.
    END IF

    DO
       hh2 = hh*div2
       hh3 = hh*div3
       hh6 = hh*div6
       hh8 = hh*div8

       !-------------S1--------------!
       ! calculate s1,s2,s3,s4,s5
       ! s1=f(t,D)
       CALL dDdtfast(D,s1,r1)

       IF (ALL(s1 == 0.0_dp)) THEN
          dD = 0.0_dp
          ddscale = 1.0_dp
          RETURN
       END IF

       inc = hh3*s1
       rinc = hh3*r1
       new_D_1 = D + inc


       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)

       !-------------S2--------------!
       ! s2 = f(t+h/3,D+h/3*s1)
       CALL dDdtfast(new_D_1,s2,r2)
       inc = hh6*s1 + hh6*s2
       rinc = hh6*r1 + hh6*r2
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)

       !-------------S3--------------!
       ! s3=f(t+h/3,D+h/6*s1+h/6*s2)
       CALL dDdtfast(new_D_1,s3,r3)
       inc = hh8*(s2+3.0_dp*s3)
       rinc = hh8*(r2+3.0_dp*r3)
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)

       !-------------S4--------------!
       ! s4=f(t+h/2,D+h/8*(s2+3*s3))
       CALL dDdtfast(new_D_1,s4,r4)
       inc = hh2*(s1-3.0_dp*s3+4.0_dp*s4)
       rinc = hh2*(r1-3.0_dp*r3+4.0_dp*r4)
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)

       !-------------S5--------------!
       ! s5=f(t+h,D+h/2*(s1-3*s3+4*s4))
       CALL dDdtfast(new_D_1,s5,r5)

       ! calculate values at t+h
       inc =  hh2*(s1-3.0_dp*s3+4.0_dp*s4)
       rinc = hh2*(r1-3.0_dp*r3+4.0_dp*r4)
       new_D_2 = D + inc
       IF (Particles % prupt) THEN
          CALL FixPrincipalDamageInc(new_D_2,inc,rinc)
       ELSE
          CALL FixPrincipalDamageInc(new_D_2,inc,rinc)
          Particles % prupt = .FALSE.
       END IF

       old_D_1 = new_D_1
       inc =  hh6*(s1+4.0_dp*s4+s5)
       rinc = hh6*(r1+4.0_dp*r4+r5)
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)

       IF (ANY(new_D_1 .NE. new_D_1)) THEN
          PRINT *,'NAN DETECTED'
          PRINT *,'No',Particles % currentno
          PRINT *,'Dav',Particles % Dav(Particles % currentno,:)
          PRINT *,'new_D_1',new_D_1
          PRINT *,'new_D_2',new_D_2
          PRINT *,'old_D_1',old_D_1
          PRINT *,'Din',Din
          PRINT *,'D',D
          PRINT *,'old inc',hh6*(s1+4.0_dp*s4+s5)
          PRINT *,'inc',inc
          PRINT *,'rinc',rinc
          PRINT *,'s1',s1
          PRINT *,'s2',s2
          PRINT *,'s3',s3
          PRINT *,'s4',s4
          PRINT *,'s5',s5
          PRINT *,'hh',hh
          PRINT *,'rupt?',Particles % prupt
          CALL Fatal('nan detected','rkm nan')
       END IF

       !------------- calculate error estimate ----------------!
       error = ABS(0.2_dp*(new_D_1-new_D_2))
       max_error = MAXVAL(error)

       IF (max_error > tol) THEN
          !halve step length and try again
          IF (ABS(hh2) < smallest_step) then

             !warn if it's particularly bad
             IF (max_error > 0.01_dp) THEN
                PRINT *,''
                CALL Warn('Runge-Kutta-Merson','RKM: VERY SMALL STEP TAKEN!')
                PRINT *,'error exceeds tolerance: ',max_error
                PRINT *,'smallest_step', smallest_step
                PRINT *,'D',D
                PRINT *,'new_D_1',D + hh6*(s1+4.0_dp*s4+s5)
                PRINT *,'new_D_2',D + hh2*(s1-3.0_dp*s3+4.0_dp*s4)
                PRINT *,'new_D_1_fixed',new_D_1
                PRINT *,'new_D_2_fixed',new_D_2
                PRINT *,'Dinit',Dinit
             END IF

             ifail = 2
             RETURN
          END IF

          hh = hh2

       ELSE

          ! accept approximation, overwrite D with D_new_1, and t with t+h
          D = new_D_1
          t = t + hh

          IF (torupt .AND. Particles % prupt .AND. (.NOT. rupted)) THEN
             rupttime = MIN(t,b)
             rupted = .TRUE.
          END IF

          !can next step be doubled?
          !slight change for zolochevsky 2007: factor = 1/32 instead of 1/100
          IF ((max_error*factor)<tol) THEN
             hh = hh*2.0_dp
          END IF

          ! if next step goes beyond interval end b, set h = b-t
          IF (t+hh>b) THEN
             hh = b - t
          END IF

          no_of_steps = no_of_steps+1
       END IF

       IF (t>=b) THEN

          dD = D - Dinit

          ddscale = 1.0_dp

          !added 10/2/19
          IF (ALL(Dinit == 0.0_dp)) THEN
             IF (ALL(ABS(dD)<Particles % mindam)) dD = 0.0_dp
          END IF

          Particles % rkmsteps = no_of_steps
          Particles % prupt = .FALSE.
          EXIT
       END IF
    END DO

  END SUBROUTINE runge_kutta_merson

  !**************************************************************************

  !> runge-kutta-merson method to solve a system of 1st order initial value ODEs
  !! integrates from t=a to t=b with the initial conditions in y
  !! Step length automatically adjusts so that absolute error estimate <= tol
  !! This version is edited for solving ice damage according to the zero-stress
  !! damage necking/mass balance modification (Bassis and Ma, 2015; Sun et al 2017;
  !! Huth et al 2020)
  SUBROUTINE runge_kutta_merson_bassis(Particles,Din,dD,a,b,tol,RHS,S0,ifail,Hin,divu,mb,bmb)

    !edited from:
    !Chivers and Sleightholme
    !Introduction to Programming with Fortran. 2015.
    !ch2602_rkm_module.f90
    !http://rhymneyconsulting.co.uk/fortran/third_edition/ch2602_rkm_module.f90

    IMPLICIT NONE

    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp), INTENT (out), DIMENSION (4) :: dD
    REAL(KIND=dp), INTENT (in) :: a, b, tol,RHS,S0,bmb
    REAL(KIND=dp), INTENT (in) :: Din(4),divu,mb
    REAL(KIND=dp), INTENT (inout) :: Hin
    INTEGER, INTENT (out) :: ifail
    REAL(KIND=dp) :: Sb,RHS1,h
    REAL(KIND=dp), DIMENSION (4) :: s1, s2, s3, s4, s5, new_D_1, new_D_2, error,Dinit,D,f
    REAL(KIND=dp) :: h1,h2,h3,h4,h5,hinc,hnew
    REAL(KIND=dp), DIMENSION (2,2) :: D2, EigVec
    REAL(KIND=dp), DIMENSION (2) :: EigVal
    REAL(KIND=dp) :: t, hh, hh2, hh3, hh6, hh8, w,x,y,z,factor ! = 1.e-2_dp
    REAL(KIND=dp) :: max_error, smallest_step = 1.e-15_dp
    INTEGER :: k,no_of_steps = 0,stoptimes = 1,ii
    LOGICAL :: Visited=.false.
    REAL(KIND=dp) :: TT,DDD,sqrteig,en,EigVec2(2,2)
    REAL(KIND=dp) :: Dmax,critdam,firststep,div2,div3,div6,div8,one=1.0_dp
    REAL(KIND=dp), DIMENSION (4) :: inc,rinc,r1,r2,r3,r4,r5

    SAVE :: Visited, DMax, critdam,firststep,div2,div3,div6,div8,factor

    IF (.NOT. Visited) THEN

       Dmax = Particles % Dmax

       !if damage > criticaldamage during RKM, set to Dmax?
       !oe just constrain damage to be <= Dmax?
       IF (Particles % RKMCritDam) THEN
          critdam = Particles % criticaldamage
       ELSE
          critdam = Particles % DmaxI
       END IF

       firststep = 1.0_dp/100.0_dp

       div2 = 1.0_dp/2.0_dp
       div3 = 1.0_dp/3.0_dp
       div6 = 1.0_dp/6.0_dp
       div8 = 1.0_dp/8.0_dp

       factor = 1.e-2_dp !1.0_dp/32.0_dp !

       Visited = .TRUE.
    END IF

    no_of_steps = 0
    ifail = 0
    dD = 0

    ! check input parameters
    IF (a==b .OR. tol<=0.0) THEN
       ifail = 1
       RETURN
    END IF

    t = a
    hh = (b-a)*firststep
    Dinit = Din
    D = Din
    H = Hin

    Particles % bmdsave = Particles % bmd

    DO
       hh2 = hh*div2
       hh3 = hh*div3
       hh6 = hh*div6
       hh8 = hh*div8

       !-------------S1--------------!
       ! calculate s1,s2,s3,s4,s5
       ! s1=f(t,D)
       Sb = S0 * H/(1.0_dp-Particles % bmdsave)
       RHS1 = RHS * (one-Sb) - bMB/H
       CALL bassisinc(Particles,D,RHS1,r1,s1)
       h1 = mb-H*divu

       inc = hh3*s1
       rinc = hh3*r1
       hinc = hh3*h1
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)
       hnew = MAX(H + hinc,one)

       !-------------S2--------------!
       ! s2 = f(t+h/3,D+h/3*s1)
       Sb = S0 * hnew/(1.0_dp-Particles % bmd)
       RHS1 = RHS * (one-Sb) - bMB/hnew
       CALL bassisinc(Particles,new_D_1,RHS1,r2,s2)
       h2 = mb-hnew*divu
       inc = hh6*s1 + hh6*s2
       rinc = hh6*r1 + hh6*r2
       hinc = hh6*h1 + hh6*h2
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)
       hnew = MAX(H + hinc,one)

       !-------------S3--------------!
       ! s3=f(t+h/3,D+h/6*s1+h/6*s2)
       Sb = S0 * hnew/(1.0_dp-Particles % bmd)
       RHS1 = RHS * (one-Sb) - bMB/hnew
       CALL bassisinc(Particles,new_D_1,RHS1,r3,s3)
       h3 = mb-hnew*divu
       inc = hh8*(s2+3.0_dp*s3)
       rinc = hh8*(r2+3.0_dp*r3)
       hinc = hh8*(h2+3.0_dp*h3)
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)
       hnew = MAX(H + hinc,one)

       !-------------S4--------------!
       ! s4=f(t+h/2,D+h/8*(s2+3*s3))
       Sb = S0 * hnew/(1.0_dp-Particles % bmd)
       RHS1 = RHS * (one-Sb) - bMB/hnew
       CALL bassisinc(Particles,new_D_1,RHS1,r4,s4)
       h4 = mb-hnew*divu
       inc = hh2*(s1-3.0_dp*s3+4.0_dp*s4)
       rinc = hh2*(r1-3.0_dp*r3+4.0_dp*r4)
       hinc = hh2*(h1-3.0_dp*h3+4.0_dp*h4)
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)
       hnew = MAX(H + hinc,one)

       !-------------S5--------------!
       ! s5=f(t+h,D+h/2*(s1-3*s3+4*s4))
       Sb = S0 * hnew/(1.0_dp-Particles % bmd)
       RHS1 = RHS * (one-Sb) - bMB/hnew
       CALL bassisinc(Particles,new_D_1,RHS1,r5,s5)
       h5 = mb-hnew*divu
       ! calculate values at t+h
       inc = hh2*(s1-3.0_dp*s3+4.0_dp*s4)
       rinc = hh2*(r1-3.0_dp*r3+4.0_dp*r4)
       new_D_2 = D + inc
       IF (Particles % prupt) THEN
          CALL FixPrincipalDamageInc(new_D_2,inc,rinc)
       ELSE
          CALL FixPrincipalDamageInc(new_D_2,inc,rinc)
          Particles % prupt = .FALSE.
       END IF

       inc = hh6*(s1+4.0_dp*s4+s5)
       rinc = hh6*(r1+4.0_dp*r4+r5)
       hinc = hh6*(h1+4.0_dp*h4+h5)
       new_D_1 = D + inc
       CALL FixPrincipalDamageInc(new_D_1,inc,rinc)
       hnew = MAX(H + hinc,one)


       !------------- calculate error estimate ----------------!
       error = ABS(0.2_dp*(new_D_1-new_D_2))
       max_error = MAXVAL(error)

       IF (max_error > tol) THEN
          !halve step length and try again
          IF (ABS(hh2) < smallest_step) then

             !warn if it's particularly bad
             IF (max_error > 0.01_dp) THEN
                PRINT *,''
                CALL Warn('Runge-Kutta-Merson','RKM: VERY SMALL STEP TAKEN!')
                PRINT *,'error exceeds tolerance: ',max_error
                PRINT *,'smallest_step', smallest_step
                PRINT *,'D',D
                PRINT *,'new_D_1',D + hh6*(s1+4.0_dp*s4+s5)
                PRINT *,'new_D_2',D + hh2*(s1-3.0_dp*s3+4.0_dp*s4)
                PRINT *,'new_D_1_fixed',new_D_1
                PRINT *,'new_D_2_fixed',new_D_2
                PRINT *,'Dinit',Dinit
             END IF

             ifail = 2
             RETURN
          END IF

          hh = hh2
       ELSE

          ! accept approximation, overwrite D with D_new_1, and t with t+h
          D = new_D_1
          Particles % bmdsave = Particles % bmd
          H = hnew
          t = t + hh

          !can next step be doubled?
          !slight change for zolochevsky 2007: factor = 1/32 instead of 1/100
          IF ((max_error*factor)<tol) THEN
             hh = hh*2.0_dp
          END IF

          ! if next step goes beyond interval end b, set h = b-t
          IF (t+hh>b) THEN
             hh = b - t
          END IF

          no_of_steps = no_of_steps+1
       END IF

       IF (t>=b) THEN

          dD = D - Dinit
          hin = h

          Particles % rkmsteps = no_of_steps
          EXIT
       END IF
    END DO

  END SUBROUTINE runge_kutta_merson_bassis

  !**************************************************************************


  !> This Runge-Kutta-Merson solver is for the case where a particle is
  !! fully failed (Particles % damstatus = 1) and
  !! Particles % DmaxII>Particles % DmaxI (Particles % dmaxII_dom == .TRUE.)
  !! Simply evolve the damage using the spin tensor, as there should only
  !! be a directional change in damage.
  !! CAUTION: may not work with recent changes?
  SUBROUTINE runge_kutta_merson_fulldam_dmaxIIdom(Particles,Din,dD,a,b,tol,ifail)

    IMPLICIT NONE

    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp), INTENT (out), DIMENSION (4) :: dD
    REAL(KIND=dp), INTENT (in) :: a, b, tol
    REAL(KIND=dp), INTENT (in) :: Din(4)
    INTEGER, INTENT (out) :: ifail
    REAL(KIND=dp), DIMENSION (4) :: s1, s2, s3, s4, s5, new_D_1, new_D_2, error,Dinit,D
    REAL(KIND=dp) :: t, hh, hh2, hh3, hh6, hh8, factor
    REAL(KIND=dp) :: max_error, smallest_step = 1.e-15_dp
    INTEGER :: k,no_of_steps = 0,stoptimes = 1
    LOGICAL :: Visited=.false.
    REAL(KIND=dp) :: critdam,firststep,div2,div3,div6,div8

    SAVE :: Visited, critdam,firststep,div2,div3,div6,div8,factor

    IF (.NOT. Visited) THEN

       !if damage > criticaldamage during RKM, set to Dmax?
       !oe just constrain damage to be <= Dmax?
       IF (Particles % RKMCritDam) THEN
          critdam = Particles % criticaldamage
       ELSE
          critdam = Particles % DmaxI
       END IF

       firststep = 1.0_dp/100.0_dp

       div2 = 1.0_dp/2.0_dp
       div3 = 1.0_dp/3.0_dp
       div6 = 1.0_dp/6.0_dp
       div8 = 1.0_dp/8.0_dp

       factor = 1.e-2_dp

       Visited = .TRUE.
    END IF

    no_of_steps = 0
    ifail = 0
    dD = 0

    ! check input parameters
    IF (a==b .OR. tol<=0.0) THEN
       ifail = 1
       RETURN
    END IF

    t = a
    hh = (b-a)*firststep
    Dinit = Din
    D = Din

    DO
       hh2 = hh*div2
       hh3 = hh*div3
       hh6 = hh*div6
       hh8 = hh*div8

       !-------------S1--------------!
       ! calculate s1,s2,s3,s4,s5
       ! s1=f(t,D)
       ! CALL dDdt(D,s1)

       s1(1) = -D(4)
       s1(2) =  D(4)
       s1(3) =  0.0_dp
       s1(4) = div2*(D(1)-D(2))
       s1 = s1*(Particles % dvdxmdudy)

       IF (ALL(s1 == 0.0_dp)) THEN
          dD = 0.0_dp
          RETURN
       END IF

       new_D_1 = D + hh3*s1
       CALL FixPrincipalDamage(new_D_1,critdam)

       !-------------S2--------------!
       ! s2 = f(t+h/3,D+h/3*s1)
       !CALL dDdt(new_D_1,s2)
       s2(1) = -new_D_1(4)
       s2(2) =  new_D_1(4)
       s2(3) =  0.0_dp
       s2(4) = div2*(new_D_1(1)-new_D_1(2))
       s2 = s2*(Particles % dvdxmdudy)

       new_D_1 = D + hh6*s1 + hh6*s2
       CALL FixPrincipalDamage(new_D_1,critdam)
       !-------------S3--------------!
       ! s3=f(t+h/3,D+h/6*s1+h/6*s2)
       !CALL dDdt(new_D_1,s3)
       s3(1) = -new_D_1(4)
       s3(2) =  new_D_1(4)
       s3(3) =  0.0_dp
       s3(4) = div2*(new_D_1(1)-new_D_1(2))
       s3 = s3*(Particles % dvdxmdudy)

       new_D_1 = D + hh8*(s2+3.0_dp*s3)
       CALL FixPrincipalDamage(new_D_1,critdam)
       !-------------S4--------------!
       ! s4=f(t+h/2,D+h/8*(s2+3*s3))
       !CALL dDdt(new_D_1,s4)
       s4(1) = -new_D_1(4)
       s4(2) =  new_D_1(4)
       s4(3) =  0.0_dp
       s4(4) = div2*(new_D_1(1)-new_D_1(2))
       s4 = s4*(Particles % dvdxmdudy)
       new_D_1 = D + hh2*(s1-3.0_dp*s3+4.0_dp*s4)
       CALL FixPrincipalDamage(new_D_1,critdam)

       !-------------S5--------------!
       ! s5=f(t+h,D+h/2*(s1-3*s3+4*s4))
       !CALL dDdt(new_D_1,s5)
       s5(1) = -new_D_1(4)
       s5(2) =  new_D_1(4)
       s5(3) =  0.0_dp
       s5(4) = div2*(new_D_1(1)-new_D_1(2))
       s5 = s5*(Particles % dvdxmdudy)

       ! calculate values at t+h
       new_D_1 = D + hh6*(s1+4.0_dp*s4+s5)
       new_D_2 = D + hh2*(s1-3.0_dp*s3+4.0_dp*s4)

       CALL FixPrincipalDamage(new_D_1,critdam)
       CALL FixPrincipalDamage(new_D_2,critdam)

       !------------- calculate error estimate ----------------!
       error = ABS(0.2_dp*(new_D_1-new_D_2))
       max_error = MAXVAL(error)

       IF (max_error > tol) THEN
          !halve step length and try again
          IF (ABS(hh2) < smallest_step) then

             !warn if it's particularly bad
             IF (max_error > 0.01_dp) THEN
                PRINT *,''
                CALL Warn('RKM_full_dam','RKM: VERY SMALL STEP TAKEN!')
                PRINT *,'error exceeds tolerance: ',max_error
                PRINT *,'smallest_step', smallest_step
             END IF

             ifail = 2
             RETURN
          END IF

          hh = hh2

       ELSE

          ! accept approximation, overwrite D with D_new_1, and t with t+h
          D = new_D_1
          t = t + hh

          !can next step be doubled?
          !slight change for zolochevsky 2007: factor = 1/32 instead of 1/100
          IF ((max_error*factor)<tol) THEN
             hh = hh*2.0_dp
          END IF

          ! if next step goes beyond interval end b, set h = b-t
          IF (t+hh>b) THEN
             hh = b - t
          END IF

          no_of_steps = no_of_steps+1
       END IF

       IF (t>=b) THEN
          dD = D - Dinit
          Particles % rkmsteps = no_of_steps
          EXIT
       END IF
    END DO

  END SUBROUTINE runge_kutta_merson_fulldam_dmaxIIdom

  !**************************************************************************

  SUBROUTINE dDdtfast(D,f,r)

    IMPLICIT NONE
    REAL(KIND=dp), INTENT (in), DIMENSION (4) :: D
    REAL(KIND=dp), OPTIONAL, DIMENSION (4) :: r
    REAL(KIND=dp), INTENT (out), DIMENSION (4) :: f

    REAL(KIND=dp) :: ah,Bf,Bh,k1,k2,rf,sthres,kf,gamma,mur88
    REAL(KIND=dp) :: Identity(3,3),Rd(3,3),maxvec(3),fd(3,3),voutprod(3,3)
    REAL(KIND=dp) :: ID(3,3),ESR(3,3),Tau(3,3),ETau(3,3),IDn1(3,3),ETau2(3,3)
    REAL(KIND=dp) :: TauP(3,3),TauN(3,3),tracesigma,rdg,f12rhs
    REAL(KIND=dp) :: EigValues(2),EigenVec(2,2),WORK(68)
    REAL(KIND=dp) :: EigValues3(3),sigmaeigval3,Q(4)
    REAL(KIND=dp) :: Dxx,Dyy,Dzz,Dxy,Peff,td1,ps1,ps2,ps3,lmod,chi,determ
    INTEGER :: infor,mm,nn,ii
    LOGICAL :: Visited=.FALSE.
    TYPE(Particle_t),  POINTER :: Particles

    REAL(KIND=dp) :: TT,DD,lambda(2),en,loceps = 1.0e-20_dp
    REAL(KIND=dp) :: quart=0.25_dp,half=0.5_dp,thres=0.0001_dp,zero=0.0_dp
    REAL(KIND=dp) :: n3=-3.0_dp,one=1.0_dp,pcont,onepfive =1.5_dp,three=3.0_dp
    REAL(KIND=dp) :: onethird=1.0_dp/3.0_dp,two=2.0_dp
    REAL(KIND=dp) :: eyyd2m1,exxd1m1,ezzd3m1,exyd4,denom,taurhsmultthird
    REAL(KIND=dp) :: lambdasqrt,d4t4,t1d2m1,t2d1m1,t3od3m1,TraceETau,maxtracesigma
    REAL(KIND=dp), POINTER :: psr(:,:)

    REAL(KIND=dp) :: e1, e2, e3, e1b,e2b,e3b
    REAL(KIND=dp) :: ez,tet,met,a12s
    REAL(KIND=dp) :: TraceEtau2,vo1(3,3),vo2(3,3),vo3(3,3),voutprod1(3,3)

    SAVE :: ah,Bh,k1,k2,rf,sthres,Identity,Visited,pcont,psr,maxtracesigma,&
         Tau,Etau !,voutprod


    Particles => GlobalParticles

    IF (.NOT. Visited) THEN

       Tau = zero
       ETau = zero
       !voutprod = zero

       ah = Particles % ah
       ! Bf = Particles % bf
       Bh = Particles % bh
       k1 = Particles % k1
       k2 = Particles % k2
       rf = Particles % rf
       sthres = Particles % sthres

       psr => Particles % psr

       pcont = one-ah-Bh

       Identity = zero
       Identity(1,1) = one
       Identity(2,2) = one
       Identity(3,3) = one

       IF (k2 .NE. 0.0_dp) THEN
          !to cap k_sigma at 30:
          maxtracesigma = (30.0_dp-k1)/k2
       ELSE
          maxtracesigma = HUGE(1.0_dp)
       END IF

       Visited = .TRUE.
    END IF


    IF (Particles % prupt .AND. Particles % noevolveruptlayers) THEN

       IF (Particles % CurrentGamma > zero) THEN
          f12rhs = D(4) * (Particles % dvdxmdudy)
          !spin contribution
          f(1) = -f12rhs
          f(2) = f12rhs
          f(3) = zero
          f(4) = -half*(Particles % dvdxmdudy)*(D(2)-D(1))
          r = f
       ELSE
          f = zero
          r = zero
       END IF

       RETURN
    END IF


    gamma = Particles % currentgamma

    !  Note:
    !  Dxx = D(1)
    !  Dyy = D(2)
    !  Dzz = D(3)
    !  Dxy = D(4)

    !------deviatoric stress-----~

    ! SLOW VERSION:
    ! !Identity - damage tensor
    ! ID = 0.0_dp
    ! ID(1,1) = D(1) !Dxx
    ! ID(2,2) = D(2) !Dyy
    ! ID(3,3) = D(3) !Dzz
    ! ID(1,2) = D(4) !Dxy
    ! ID(2,1) = D(4) !Dxy
    ! ID = Identity - ID
    ! !deviatoric effective strain rate tensor:
    ! ESR = (MATMUL(ID,Particles % pSR) + MATMUL(Particles % pSR,ID))
    ! ESR = half*(ESR - (Identity * ( ( ESR(1,1)+ ESR(2,2) + ESR(3,3) )/three) ))
    ! !deviatoric stress tensor:
    ! Tau = Particles % RHS * ESR

    !FAST VERSION:
    exxd1m1 = psr(1,1)*(D(1)-one)
    eyyd2m1 = psr(2,2)*(D(2)-one)
    ezzd3m1 = psr(3,3)*(D(3)-one)
    exyd4 = psr(1,2)*D(4)
    taurhsmultthird = Particles % RHS * onethird


    Tau(1,1) = taurhsmultthird * (-two*exxd1m1+eyyd2m1+ezzd3m1-exyd4)
    Tau(2,2) = taurhsmultthird * (exxd1m1-two*eyyd2m1+ezzd3m1-exyd4)
    Tau(3,3) = taurhsmultthird * (exxd1m1+eyyd2m1-two*ezzd3m1+two*exyd4)
    Tau(1,2) = -Particles % RHS * half * (psr(1,2)*(D(1)+D(2)-two) + D(4)*(psr(1,1)+psr(2,2)))
    Tau(2,1) = Tau(1,2)

    !------- effective pressure (Peff)------
    Peff = Particles % pressure1 - (Tau(1,1)+Tau(2,2))

    IF (Particles % usetruecauchydamage) THEN
       !Tau is now the true Cauchy stress, not deviatoric
        Tau(1,1) = Tau(1,1) - Peff
        Tau(2,2) = Tau(2,2) - Peff
        Tau(3,3) = Tau(3,3) - Peff

       !The vertical z component of actual effective cauchy stress, used later for
       !caulculating the new Peff
        sigmaeigval3 = (Tau(3,3))/(1.0_dp - D(3))
    ENDIF

    !---- deviatoric effective stress (ETau) -----!
    !IDn1 = (I-D)^-1, expanded here for speed:
    denom = one/( D(4)*D(4) + D(1) + D(2) -D(1)*D(2) - one)
    t1d2m1 = Tau(1,1)*(D(2)-one)*denom
    t2d1m1 = Tau(2,2)*(D(1)-one)*denom
    t3od3m1 = Tau(3,3)/(D(3)-one)
    d4t4 = Tau(1,2)*D(4)*denom

    Etau(1,1) = onethird*(two*t1d2m1 - t2d1m1 +t3od3m1 -d4t4)
    Etau(2,2) = onethird*(-t1d2m1 + two*t2d1m1 +t3od3m1 -d4t4)
    Etau(3,3) = onethird*(-t1d2m1 - t2d1m1 - two*t3od3m1 + two*d4t4)
    Etau(1,2) = half*denom*(tau(1,2)*(D(1)+D(2)-two) - D(4)*(Tau(1,1)+Tau(2,2)))
    Etau(2,1) = Etau(1,2)

    IF (Particles % usetruecauchydamage) THEN
       Peff = Etau(3,3)-sigmaeigval3
    ENDIF


    !---- eigenvalues -----!
    TT = ETau(1,1)+ETau(2,2)
    DD = ETau(1,1)*ETau(2,2)-ETau(1,2)*ETau(1,2)
    lambdasqrt = quart*TT*TT-DD

    IF (lambdasqrt<0.0_dp) lambdasqrt = 0.0_dp
    lambdasqrt = sqrt(lambdasqrt)

    lambda(1)=half*TT + lambdasqrt
    lambda(2)=half*TT - lambdasqrt

    IF (ANY(lambda .NE. lambda)) THEN
       CALL Eigen2D(ETau(1:2,1:2),EigValues,EigenVec)
       lambda(2)=EigValues(1)
       lambda(1)=EigValues(2)
    END IF

    ! e3 = ETau(3,3)

    !-----max effective principal stress------
       td1 = lambda(1)-Peff

    !no change in damage
    IF (td1 < zero ) THEN !.AND. ALL(D==zero)) THEN
       f = zero
       ! RETURN

       IF (Particles % CurrentGamma > zero) THEN
          f12rhs = D(4) * (Particles % dvdxmdudy)
          !spin contribution
          f(1) = f(1) - f12rhs
          f(2) = f(2) + f12rhs
          f(4) = f(4) - half*(Particles % dvdxmdudy)*(D(2)-D(1))
          r = f
       END IF

       RETURN
    END IF

    sigmaeigval3 = ETau(3,3)-Peff

    ! IF (gamma == zero) THEN
    !    IF (td1 < zero) THEN
    !       f=zero
    !       RETURN
    !    END IF
    ! END IF

    !------- tracesigma (for kf) ------
    tracesigma = Tau(1,1)+Tau(2,2)+Tau(3,3)-three*Peff
    IF (tracesigma<zero) tracesigma = zero
    IF (tracesigma>maxtracesigma) tracesigma = maxtracesigma


    !-----kf------
    ! tracesigma = ETau(1,1)+ETau(2,2)+ETau(3,3)-three*Peff
    ! IF (tracesigma<zero) tracesigma = zero
    ! IF (tracesigma>maxtracesigma) tracesigma = maxtracesigma
    !e.g.  kf = 3.0_dp + 6.0_dp * tracesigma
    IF (tracesigma > Particles % stresshigh) Particles % stresshigh = tracesigma
    IF (tracesigma < Particles % stresslow) Particles % stresslow = tracesigma
    kf = k1 + k2*tracesigma


    !------modified Murakami creep damage------
    ! corrects symmetry issue
    !Ganczarski and Skrzypek, 2001
    IF (gamma>0.0_dp .AND. Particles % modifiedmurakami) THEN

       EigValues3(1) = td1
       EigValues3(2) = lambda(2) - Peff
       EigValues3(3) = sigmaeigval3

       ! 2d
       ! IF (ALL(EigValues3(1:2) > 0.0_dp)) THEN
       !    lmod = 1.0_dp - MINVAL(EigValues3(1:2))/MAXVAL(EigValues3(1:2))
       !    gamma = gamma*lmod
       ! END IF

       !3d
       ps3 = MINVAL(EigValues3)

       IF (ps3 > zero) THEN
          ps1 = MAXVAL(EigValues3)
          ps2 = SUM(EigValues3) - ps1 - ps3
          lmod = (ps2+ps3)/ps1
          lmod = one - lmod*(one-(ps3/(two*ps2)))
          gamma = gamma*lmod
       END IF
    END IF


    !----------- The damage rate -------!

    !here, f is the damage rate we return.
    !see Murakami et al (1980,1988,etc) and Pralong and Funk(2005,2006)

    !f = fd - [D,W]

    !where fd = f1 * Rd   and    [D,W] = DW-WD

    !f1 is defined below
    !Rd is a linear combination of isotropic and orthotropic contributions
    !D is current damage
    !W is the spin tensor


    !--- Rd
    IF (lambda(1) .NE. lambda(2)) THEN

       CALL Eigen2DSym_TryGenFirst_VecOnly(ETau(1:2,1:2),lambda(1:2),EigValues,EigenVec)

       maxvec(1:2) = EigenVec(1:2,2)
       maxvec(3) = zero

       voutprod(:,1) = maxvec(1) * maxvec(:)
       voutprod(:,2) = maxvec(2) * maxvec(:)
       voutprod(:,3) = zero !maxvec(3) * maxvec(:) !zero

       !  IF (ALL(voutprod .NE. 1)) THEN
       ! PRINT *,'voutprod',voutprod
       !  END IF
    ELSE

       CALL Eigen2DSym_TryGenFirst(ETau(1:2,1:2),EigValues,EigenVec)

       IF (EigValues(1) .EQ. EigValues(2)) THEN

          CALL Eigen2DSym_TryGenFirst(Tau(1:2,1:2),EigValues,EigenVec)

          IF (EigValues(1) .EQ. EigValues(2)) THEN
             Particles % equaleigcount = Particles % equaleigcount + 1

             CALL Eigen2DSym_TryGenFirst(ETau(1:2,1:2),EigValues,EigenVec)

             gamma = gamma * 0.99_dp
          END IF
       END IF

       maxvec(1:2) = EigenVec(1:2,2)
       maxvec(3) = zero

       voutprod(:,1) = maxvec(1) * maxvec(:)
       voutprod(:,2) = maxvec(2) * maxvec(:)
       voutprod(:,3) = zero !maxvec(3) * maxvec(:)
    END IF


    rdg = one-gamma

    Rd = gamma * voutprod
    Rd(1,1) = Rd(1,1) + rdg
    Rd(2,2) = Rd(2,2) + rdg
    Rd(3,3) = Rd(3,3) + (one-gamma)

    IF (Particles % forcedzz) Rd(3,3) = 1.0_dp
    !Particles % CurrentGamma) !rdg
    !rd(3,3) will always be zero. And using currentgamma for z right now
    !to test the effect of only applying  modified murikami to horizontal dirs
    !  Rd(3,3) = one-gamma

    !trace(etau*etau)
    TraceEtau = ETau(1,1)*ETau(1,1) + two*ETau(1,2)*ETau(1,2) + &
         ETau(2,2)*ETau(2,2)+ETau(3,3)*ETau(3,3)


    chi = ah*td1 + Bh * &
         sqrt( onepfive * TraceEtau ) &
         + pcont * (n3*Peff)

    !don't allow damage to accumulate if:
    IF (chi<Particles % sthresmod) THEN
       chi = Particles % sthresmod
       !ELSEIF (td1 <= zero) THEN
       !   chi = Particles % sthresmod
       ! ELSEIF (tracesigma <= zero) THEN
       !    chi = Particles % sthresmod
    END IF


    !------fd = f1 * Rd
    !f1 is of the form f1 = Bf * (X)^r * (m)^kf
    !(Murakami, 1988)

    !get m
    ! voutprod = MATMUL(IDn1,voutprod)
    ! mur88 = voutprod(1,1)+ voutprod(2,2)+voutprod(3,3)
    IF (Particles % prupt) THEN

       Q = D
       CALL removemaxd(Q)
       denom = one/( Q(4)*Q(4) + Q(1) + Q(2) -Q(1)*Q(2) - one)

       mur88 = (voutprod(1,1)*(Q(2)-one) + voutprod(2,2)*(Q(1)-one) &
            - two*Q(4)*voutprod(1,2) ) * denom


    ELSE
       mur88 = (voutprod(1,1)*(D(2)-one) + voutprod(2,2)*(D(1)-one) &
            - two*D(4)*voutprod(1,2) ) * denom

       !note that "- voutprod(3,3)/(D(3)-one)" is not needed above
       !because voutprod(3,3) is always zero in the current implementation
    END IF



    fd = (Particles % bf*((chi-Particles % sthresmod)**rf)*(mur88**kf))*Rd

    IF (Particles % noDzz) fd(3,3) = zero

    !   f = fd - [D,W] reduces to:

    !1. damage accumulation (no spin)
    f(1) = fd(1,1)                !dDxx/dt
    f(2) = fd(2,2)                !dDyy/dt
    f(3) = fd(3,3)                !dDzz/dt
    f(4) = fd(1,2)                !dDxy/dt and dDyx/dt

    !2. spin contribution
    f12rhs = D(4) * (Particles % dvdxmdudy)
    r(1) = -f12rhs
    r(2) = f12rhs
    r(3) = zero
    r(4) = -half*(Particles % dvdxmdudy)*(D(2)-D(1))

    IF (Particles % nospin) r = 0.0_dp

    f = f+r

    Particles % TempGamma = gamma

  END SUBROUTINE dDdtfast

  !**************************************************************************

  !> Update particle variables before SSA solution
  SUBROUTINE SSAPrepMeshToParticles( Particles, Model)

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: nn,nb, dim,dofs
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: rhoi,rhow, mbparam, efparam, fricparam,slope1,slope2
    REAL(KIND=dp),ALLOCATABLE :: Basis(:),dBasisdx(:,:)
    LOGICAL :: gotit,Visited=.FALSE.,UpdateB,stat
    LOGICAL :: movegl, ConstantMB, ConstantEF, ConstantFP
    TYPE(Variable_t), POINTER :: mask, MB, EF, F, B, Zs, V,Bed, HVar
    INTEGER, POINTER :: maskPerm(:),MBPerm(:),EFPerm(:),FPerm(:),&
         BPerm(:),ZsPerm(:),VPerm(:),BedPerm(:),LocalPerm(:),HPerm(:)
    REAL(KIND=dp), POINTER :: maskVal(:),MBVal(:),EFVal(:),FVal(:),&
         BVal(:),ZsVal(:), VVal(:), BedVal(:),HVal(:)
    REAL(KIND=dp), POINTER :: ZsLocalField(:),FPLocalField(:),BLocalField(:),&
         maskLocalField(:),EFLocalField(:),VLocalField(:,:),&
         BedLocalField(:),MBLocalField(:),HLocalField(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    INTEGER :: No,ii,jj,kk,surelem,t
    TYPE(Element_t), POINTER :: BulkElement
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: detJ,scale, Hf,min_first_binit,min_first_ef,max_first_ef

    SAVE :: Mesh,nn,nb,dim,Basis,dBasisdx,ZsLocalField,FPLocalField,BLocalField,&
         maskLocalField,EFLocalField,VLocalField,&
         BedLocalField,MBLocalField, rhoi, rhow, movegl, ConstantMB, ConstantEF,&
         ConstantFP, mbparam,efparam,fricparam, mask, maskperm, maskval,&
         mb,mbperm,mbval,ef,efperm,efval,F,FPerm,FVal,B,BPerm,BVal,&
         Zs,ZsPerm,ZsVal,V,Vperm,VVal,Bed,BedPerm,BedVal,Visited,&
         UpdateB,SolverName,dofs,LocalPerm,min_first_binit,&
         HVar,HVal,HPerm,HlocalField,ElementNodes,&
         min_first_ef,max_first_ef


    IF( .NOT. Visited ) THEN

       WRITE(SolverName, '(A)') 'SSAPrepMeshToParticles'

       Mesh => GetMesh()
       nn = Mesh % MaxElementNodes
       nb = Mesh % NumberOfBulkElements
       dim = 2

       ALLOCATE( Basis(nn), dBasisdx(nn, 3))
       ALLOCATE( ElementNodes % x(nn),ElementNodes % y(nn), ElementNodes % z(nn))


       ALLOCATE( ZsLocalField(nn), FPLocalField(nn), &
            BLocalField(nn),MaskLocalField(nn),&
            EFLocalField(nn), VLocalField(nn,dim), &
            BedLocalField(nn), MBLocalField(nn),LocalPerm(nn),&
            HLocalField(nn))


       rhoi = Particles % rhoi
       rhow = Particles % rhow


       movegl = Particles % movegl
       ConstantMB = Particles % constmb
       ConstantEF = Particles % constef
       ConstantFP = Particles % constfric

       HVar => VariableGet(Model % Mesh % Variables, 'H' )
       HPerm => HVar % Perm
       HVal => HVar % Values

       IF (ConstantMB) THEN
          mbparam = GetConstReal( Model % Constants, 'mbparam', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Need to define "mbparam = Real $mbparam" in constants')
       END IF

       IF (ConstantEF) THEN
          efparam = GetConstReal( Model % Constants, 'efparam', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Need to define "efparam = Real $efparam" in constants')
       END IF

       IF (ConstantFP) THEN
          fricparam = GetConstReal( Model % Constants, 'fricparam', GotIt )
          IF (.NOT. GotIt) CALL Fatal(SolverName, &
               'Need to define "fricparam = Real $fricparam" in constants')
       END IF

       mask => VariableGet(Model % Mesh % Variables, 'Mask' )
       maskPerm => mask % Perm
       maskVal => mask % Values

       IF (.NOT. ConstantMB) THEN
          MB => VariableGet(Model % Mesh % Variables, 'MB' )
          MBPerm => MB % Perm
          MBVal => MB % Values
       END IF

       IF (.NOT. ConstantEF) THEN

          EF => VariableGet(Model % Mesh % Variables, 'EF' )
          EFPerm => EF % Perm
          EFVal => EF % Values
       END IF

       IF (.NOT. ConstantFP) THEN
          F => VariableGet(Model % Mesh % Variables, 'FP' )
          FPerm => F % Perm
          FVal => F % Values
       END IF

       IF ((.NOT.  Particles % constlintemp) .AND. &
            (.NOT. Particles % useconsttemp) .AND. &
            (.NOT. Particles % usegiveneta)) THEN

          B => VariableGet(Model % Mesh % Variables, 'invvisc' )
          BPerm => B % Perm
          BVal => B % Values

          UpdateB = .TRUE.
       ELSE
          UpdateB = .FALSE.
       END IF


       Zs => VariableGet(Model % Mesh % Variables, 'Zs' )
       ZsPerm => Zs % Perm
       ZsVal => Zs % Values

       V => VariableGet(Model % Mesh % Variables, 'SSAVelocity' )
       VPerm => V % Perm
       VVal => V % Values
       dofs = V % Dofs

       Bed => VariableGet(Model % Mesh % Variables, 'Bed' )
       BedPerm => Bed % Perm
       BedVal => Bed % Values

       min_first_binit = MINVAL(Particles % binit(1:Particles % numberofparticles))

       min_first_ef =  MINVAL(Particles % EF(1:Particles % numberofparticles))
       max_first_ef =  MAXVAL(Particles % EF(1:Particles % numberofparticles))

       PRINT *,'min_first_binit',min_first_binit
       PRINT *,'min_first_ef',min_first_ef
       PRINT *,'max_first_ef',max_first_ef


       Visited = .TRUE.
    END IF


    CALL Info(SolverName,&
         'Starting SSA Prep: Mesh to Particles',Level=3)

    !-----------INITIALIZATION---------------!

    Particles % xpic = 0.0_dp

    !for interpelem particles, save the old values you need on xpic
    !new binit values for (grounded and .not. interpelem)
    !are saved on xpic4 and transferred over later
    !since we dont change binit for floating particles
    !or interpelem at this time

    IF (Particles % binitlowerboundfromfirsttimestep) THEN
       Particles % xpic(:,5) = Particles % Binit(:)
    END IF



    DO No = 1, Particles % NumberOfParticles
       IF (.NOT. Particles % UseInterpElem(No) ) CYCLE
       !IF (Particles % UseInterpElem(No) .OR. Particles % Status(No)==PARTICLE_LEAVING) THEN
       Particles % xpic(No,2) = Particles % EF(No)
       Particles % xpic(No,3) = Particles % FP(No)
       Particles % xpic(No,5) = Particles % Binit(No)
       ! END IF
    END DO

    Particles % Gmask = 0.0_dp
    Particles % GradH = 0.0_dp
    Particles % GradZs = 0.0_dp
    Particles % GradVel = 0.0_dp

    !   IF (.NOT. ConstantMB) Particles % MB = 0.0_dp
    IF (.NOT. ConstantEF) Particles % EF = 0.0_dp
    IF (.NOT. ConstantFP) Particles % FP = 0.0_dp
    !IF (movegl .OR. Particles % SEP .OR. Particles % CoulombFriction)
    Particles % Bedrock = 0.0_dp
    IF (Particles % VelocityDependentFriction) Particles % GridVelocity = 0.0_dp


    !----------------------------------------------------------------------------!
    !-------------------------- SET ELEMTRACK % GSTATUS   -----------------------!
    !----------------------------------------------------------------------------!
    !we do this for efficiency...
    !   Mark elements using ElemTrack % Gstatus as floating (1), ungrounding (0), or grounded (-1)
    !   ElemTrack % Gstatus = 0 (Ungrounding) is assigned for  elements with
    !both grounded and ungrounded nodes, as well as the elements next to them.

    !set all to floating to start
    ElemTrack(:) % Gstatus = 1

    !reverse mask vals if using floatation solver
    IF (movegl) THEN
       maskval = -maskval
    END IF

    DO ii = 1,nb

       IF ( ElemTrack(ii) % Status < FEM ) CYCLE
       ! IF ( ElemTrack(ii) % Status < NOTFULL ) CYCLE

       BulkElement => Model % Mesh % Elements( ii )
       NodeIndexes => BulkElement % NodeIndexes

       IF (ANY(maskval(maskperm(nodeindexes)) <= 0.0_dp)) THEN
          !element is grounded or ungrounding

          IF (ANY(maskval(maskperm(nodeindexes)) >= 0.0_dp)) THEN
             !element is  ungrounding
             ElemTrack(ii) % GStatus = 0
             DO jj = 1,9
                surelem = ElemParticles(ii) % SurroundingElems(jj)
                IF (surelem < 1) CYCLE
                !ungrounding.
                ElemTrack(surelem) % Gstatus = 0
                DO t = 1, ElemParticles(surelem) % NumberOfParticles
                   No = ElemParticles(surelem) % p(t)
                   !mark particles as ungrounding
                   Particles % xpic(No,1) = 1.0_dp
                END DO
             END DO
          ELSE
             !element is fully grounded, unless already indicated as being
             !next to an ungrounded element, which means it may contain ungrounding particles
             IF (ElemTrack(ii) % Gstatus == 0) THEN
                !ungrounding.
                CYCLE
             ELSE
                !grounded.
                ElemTrack(ii) % GStatus = -1
             END IF
          END IF
       END IF
    END DO

    !----------------------------------------------------------------------------!
    !------------------------- MESH VARS TO PARTICLES   -------------------------!
    !----------------------------------------------------------------------------!

    DO ii = 1,nb
       IF ( ElemTrack(ii) % Status < FEM ) CYCLE
       ! IF ( ElemTrack(ii) % Status < NOTFULL ) CYCLE

       BulkElement => Model % Mesh % Elements(ii)
       NodeIndexes => BulkElement % NodeIndexes(1:nn)
       CALL GetElementNodes(ElementNodes,BulkElement)


       !---VARS FOR PARTICLES OF ANY GROUNDED STATUS---
       !mass balance
       ! IF (.NOT. ConstantMB) MBLocalField(1:nn) = MBVal(MBPerm(NodeIndexes))

       !gradvel/gridvelocity
       LocalPerm = VPerm(NodeIndexes)
       DO kk = 1,nn
          DO jj = 1,dim
             VLocalField(kk,jj) = VVal(dofs*(LocalPerm(kk)-1)+jj)
          END DO
       END DO


       !gradzs
       ZsLocalField = ZsVal(ZSPerm(NodeIndexes))

       !gradH
       HLocalField = HVal(HPerm(NodeIndexes))
       !-------------------------------------------


       IF (ElemTrack(ii) % GStatus == 1) THEN

          !--------------------------FULLY FLOATING ELEMENT-----------------------

          IF (.NOT. ConstantEF) EFLocalField = EFVal(EFPerm(NodeIndexes))
          IF (.NOT. ConstantFP) FPLocalField = FVal(FPerm(NodeIndexes))

          DO t = 1, ABS(ElemParticles(ii) % NumberOfParticles)
             No = ElemParticles(ii) % p(t)

             IF (Particles % ShapeFunctions == 'gimpm') THEN
                stat = GIMPMElementInfo( t,Particles, Model,BulkElement, ElementNodes, No, &
                     detJ, scale, .TRUE., Basis,dBasisdx)
             ELSE

                stat = sMPMElementInfo( Bulkelement, Particles, Model, ElementNodes, No, &
                     Particles % gridres, Basis,dBasisdx)
                scale = 1.0_dp
             END IF

             ! IF (.NOT. ConstantMB) THEN
             !    Particles % MB(No) = Particles % MB(No) + &
             !         SUM(Basis(1:nn) * MBLocalField(1:nn)) * scale
             ! END IF

             DO jj = 1,dim
                Particles % GradZs(No,jj) = Particles % GradZs(No,jj) + &
                     SUM(dBasisdx(1:nn,jj) * ZsLocalField(1:nn)) * scale

                Particles % GradH(No,jj) = Particles % GradH(No,jj) + &
                     SUM(dBasisdx(1:nn,jj) * HLocalField(1:nn)) * scale
             END DO

             !dvx/dx
             Particles % GradVel(No,1) = Particles % GradVel(No,1) + &
                  SUM(dBasisdx(1:nn,1) * VLocalField(1:nn,1)) * scale

             !dvy/dy
             Particles % GradVel(No,2) = Particles % GradVel(No,2) + &
                  SUM(dBasisdx(1:nn,2) * VLocalField(1:nn,2)) * scale

             !dvx/dy
             Particles % GradVel(No,3) = Particles % GradVel(No,3) + &
                  SUM( dBasisdx(1:nn,2) * VLocalField(1:nn,1) ) * scale

             !dvy/dx
             Particles % GradVel(No,4) = Particles % GradVel(No,4) + &
                  SUM( dBasisdx(1:nn,1) * VLocalField(1:nn,2) ) * scale


             IF (Particles % VelocityDependentFriction) THEN
                DO kk = 1,2
                   Particles % GridVelocity(No,kk) = Particles % GridVelocity(No,kk) + &
                        SUM(Basis(1:nn) * VLocalField(1:nn,kk)) * scale
                END DO
             END IF

             ! IF (movegl .OR. Particles % CoulombFriction .OR. Particles % SEP) THEN
             Particles % Bedrock(No) = Particles % Bedrock(No) + &
                  SUM(Basis(1:nn) * BedVal(BedPerm(NodeIndexes))) * scale
             ! END IF

             IF (Particles % xpic(No,1) .NE. 1.0_dp) THEN
                Particles % Gmask(No) = 1.0_dp
             ELSE
                MaskLocalField = maskVal(MaskPerm(NodeIndexes))
                Particles % Gmask(No) = Particles % Gmask(No) + &
                     SUM(Basis(1:nn) * maskLocalField(1:nn)) * scale
             END IF


             IF (Particles % UseInterpElem(No)) CYCLE
             !previous EF is retained if interpelem

             IF (.NOT. ConstantEF) THEN
                Particles % EF(No) = Particles % EF(No) + &
                     SUM(Basis(1:nn) * EFLocalField(1:nn)) * scale
                !Particles % EF(No) = MAX(Particles % EF(No),0.0_dp)
             END IF

             IF (.NOT. ConstantFP) THEN
                Particles % FP(No) = Particles % FP(No) + &
                     SUM(Basis(1:nn) * FPLocalField(1:nn)) * scale
             END IF

          END DO !particle loop


       ELSE IF (ElemTrack(ii) % GStatus == -1) THEN



          !------------------------FULLY GROUNDED PARTICLES-------------------

          IF (Particles % FEMifGrounded) THEN
             ElemTrack(ii) % Status = FEM
          END IF


          IF (.NOT. ConstantEF) THEN
             EFLocalField = EFVal(EFPerm(NodeIndexes))
          END IF

          IF (.NOT. ConstantFP) FPLocalField = FVal(FPerm(NodeIndexes))
          !IF (Particles % CoulombFriction)
          BedLocalField = BedVal(BedPerm(NodeIndexes))
          IF (UpdateB) BLocalField = BVal(BPerm(NodeIndexes))

          DO t = 1, ElemParticles(ii) % NumberOfParticles
             No = ElemParticles(ii) % p(t)

             IF (Particles % ShapeFunctions == 'gimpm') THEN
                stat = GIMPMElementInfo( t,Particles, Model,BulkElement, ElementNodes, No, &
                     detJ, scale, .TRUE., Basis,dBasisdx)
             ELSE
                stat = sMPMElementInfo( Bulkelement, Particles, Model, ElementNodes, No, &
                     Particles % gridres, Basis,dBasisdx)
                scale = 1.0_dp
             END IF

             ! IF (.NOT. ConstantMB) THEN
             !    Particles % MB(No) = Particles % MB(No) + &
             !         SUM(Basis(1:nn) * MBLocalField(1:nn)) * scale
             ! END IF

             DO jj = 1,dim
                Particles % GradZs(No,jj) = Particles % GradZs(No,jj) + &
                     SUM(dBasisdx(1:nn,jj) * ZsLocalField(1:nn)) * scale

                Particles % GradH(No,jj) = Particles % GradH(No,jj) + &
                     SUM(dBasisdx(1:nn,jj) * HLocalField(1:nn)) * scale
             END DO

             !dvx/dx
             Particles % GradVel(No,1) = Particles % GradVel(No,1) + &
                  SUM(dBasisdx(1:nn,1) * VLocalField(1:nn,1)) * scale

             !dvy/dy
             Particles % GradVel(No,2) = Particles % GradVel(No,2) + &
                  SUM(dBasisdx(1:nn,2) * VLocalField(1:nn,2)) * scale


             !dvx/dy
             Particles % GradVel(No,3) = Particles % GradVel(No,3) + &
                  SUM( dBasisdx(1:nn,2) * VLocalField(1:nn,1) ) * scale

             !dvy/dx
             Particles % GradVel(No,4) = Particles % GradVel(No,4) + &
                  SUM( dBasisdx(1:nn,1) * VLocalField(1:nn,2) ) * scale


             IF (Particles % VelocityDependentFriction) THEN
                DO kk = 1,2
                   Particles % GridVelocity(No,kk) = Particles % GridVelocity(No,kk) + &
                        SUM(Basis(1:nn) * VLocalField(1:nn,kk)) * scale
                END DO
             END IF


             !  IF (Particles % CoulombFriction .OR. Particles % SEP) THEN
             Particles % Bedrock(No) = Particles % Bedrock(No) + &
                  SUM(Basis(1:nn) * BedLocalField(1:nn)) * scale
             !  END IF

             IF (Particles % xpic(No,1) .NE. 1.0_dp) THEN
                Particles % Gmask(No) = -1.0_dp
             ELSE
                MaskLocalField(1:nn) = maskVal(MaskPerm(NodeIndexes))
                Particles % Gmask(No) = Particles % Gmask(No) + &
                     SUM(Basis(1:nn) * maskLocalField(1:nn)) * scale
             END IF


             !previous EF,FP, and eta are retained if interpelem
             IF (Particles % UseInterpElem(No)) CYCLE


             ! IF (Particles % xpic(No,1) .NE. 1.0_dp) THEN
             IF (.NOT. ConstantEF) THEN
                Particles % EF(No) = Particles % EF(No) + &
                     SUM(Basis(1:nn) * EFLocalField(1:nn)) * scale
             END IF
             !ELSE
             !IF (.NOT. ConstantEF) THEN
             !    Particles % EF(No) = Particles % EF(No) + &
             !         SUM(Basis(1:nn) * EFLocalField(1:nn)) * scale
             !  END IF
             !END IF

             IF (.NOT. ConstantFP) THEN
                Particles % FP(No) = Particles % FP(No) + &
                     SUM(Basis(1:nn) * FPLocalField(1:nn)) * scale
             END IF

             IF (UpdateB) THEN
                Particles % xpic(No,4) = Particles % xpic(No,4) + &
                     SUM(Basis(1:nn) * BLocalField(1:nn)) * scale
             END IF

          END DO !particle loop

       ELSE
          !-------UNGROUNDING PARTICLES--------

          IF (.NOT. ConstantEF) EFLocalField = EFVal(EFPerm(NodeIndexes))
          IF (.NOT. ConstantFP) FPLocalField = FVal(FPerm(NodeIndexes))

          ! IF (movegl .OR. Particles % CoulombFriction .OR. Particles % SEP) THEN
          BedLocalField = BedVal(BedPerm(NodeIndexes))
          !  END IF

          IF (UpdateB) BLocalField = BVal(BPerm(NodeIndexes))

          MaskLocalField = maskVal(MaskPerm(NodeIndexes))

          DO t = 1, ElemParticles(ii) % NumberOfParticles
             No = ElemParticles(ii) % p(t)

             IF (Particles % ShapeFunctions == 'gimpm') THEN
                stat = GIMPMElementInfo( t,Particles, Model,BulkElement, ElementNodes, No, &
                     detJ, scale, .TRUE., Basis,dBasisdx)
             ELSE
                stat = sMPMElementInfo( Bulkelement, Particles, Model, ElementNodes, No, &
                     Particles % gridres, Basis,dBasisdx)
                scale = 1.0_dp
             END IF

             ! IF (.NOT. ConstantMB) THEN
             !    Particles % MB(No) = Particles % MB(No) + &
             !         SUM(Basis(1:nn) * MBLocalField(1:nn)) * scale
             ! END IF

             DO jj = 1,dim
                Particles % GradZs(No,jj) = Particles % GradZs(No,jj) + &
                     SUM(dBasisdx(1:nn,jj) * ZsLocalField(1:nn)) * scale

                Particles % GradH(No,jj) = Particles % GradH(No,jj) + &
                     SUM(dBasisdx(1:nn,jj) * HLocalField(1:nn)) * scale
             END DO

             !dvx/dx
             Particles % GradVel(No,1) = Particles % GradVel(No,1) + &
                  SUM(dBasisdx(1:nn,1) * VLocalField(1:nn,1)) * scale

             !dvy/dy
             Particles % GradVel(No,2) = Particles % GradVel(No,2) + &
                  SUM(dBasisdx(1:nn,2) * VLocalField(1:nn,2)) * scale


             !dvx/dy
             Particles % GradVel(No,3) = Particles % GradVel(No,3) + &
                  SUM( dBasisdx(1:nn,2) * VLocalField(1:nn,1) ) * scale

             !dvy/dx
             Particles % GradVel(No,4) = Particles % GradVel(No,4) + &
                  SUM( dBasisdx(1:nn,1) * VLocalField(1:nn,2) ) * scale


             IF (Particles % VelocityDependentFriction) THEN
                DO kk = 1,2
                   Particles % GridVelocity(No,kk) = Particles % GridVelocity(No,kk) + &
                        SUM(Basis(1:nn) * VLocalField(1:nn,kk)) * scale
                END DO
             END IF

             ! IF (Particles % CoulombFriction .OR. Particles % SEP) THEN
             Particles % Bedrock(No) = Particles % Bedrock(No) + &
                  SUM(Basis(1:nn) * BedLocalField(1:nn)) * scale
             ! END IF

             Particles % Gmask(No) = Particles % Gmask(No) + &
                  SUM(Basis(1:nn) * maskLocalField(1:nn)) * scale

             !previous EF,FP,eta,and mask are retained if interpelem
             IF (Particles % UseInterpElem(No)) CYCLE

             IF (.NOT. ConstantEF) THEN
                Particles % EF(No) = Particles % EF(No) + &
                     SUM(Basis(1:nn) * EFLocalField(1:nn)) * scale
             END IF

             IF (.NOT. ConstantFP) THEN
                Particles % FP(No) = Particles % FP(No) + &
                     SUM(Basis(1:nn) * FPLocalField(1:nn)) * scale
             END IF

             IF (UpdateB) THEN
                !put on xpic4 for now.  Will be moved to binit if
                !gmask < 0.9
                Particles % xpic(No,4) = Particles % xpic(No,4) + &
                     SUM(Basis(1:nn) * BLocalField(1:nn)) * scale
             END IF

          END DO !particle loop
       END IF !which element grounding status
    END DO !element loop

    !SEP, Binit, InterpElemParticles corrections
    IF (Particles % SEP) THEN

       DO No = 1, Particles % NumberOfParticles

          IF (Particles % Gmask(No) > 0.99_dp) THEN
             Particles % Gmask(No) = 1.0_dp
          ELSE IF (Particles % Gmask(No) < -0.99_dp) THEN
             Particles % Gmask(No) = -1.0_dp
          ELSE
             Hf = rhow * (Particles % sealevel-Particles % bedrock(No))/rhoi
             IF (Particles % H(No) .LT. Hf) THEN
                Particles % Gmask(No) = 1.0_dp
             ELSE
                Particles % Gmask(No) = -1.0_dp
             END IF
          END IF

          IF (Particles % UseInterpElem(No)) THEN
             !IF (Particles % UseInterpElem(No) .OR. Particles % Status(No)==PARTICLE_LEAVING) THEN
             Particles % EF(No) = Particles % XPIC(No,2)
             Particles % FP(No) = MAX(Particles % XPIC(No,3),0.0_dp)
             Particles % Binit(No) = Particles % xpic(No,5)
          ELSE
             IF (.NOT. UpdateB) CYCLE
             IF (Particles % GMask(No) >= 0.9_dp) CYCLE
             IF (ANY(Particles % Dav(No,:).NE.0.0_dp)) CYCLE

             Particles % Binit(No) = Particles % xpic(No,4)
          END IF
       END DO

    ELSE

       IF ((.NOT. ConstantEF) .OR. (.NOT. ConstantFP) .OR.  (UpdateB)) THEN
          DO No = 1, Particles % NumberOfParticles

             IF (Particles % UseInterpElem(No)) THEN
                Particles % EF(No) = Particles % XPIC(No,2)
                Particles % FP(No) = MAX(Particles % XPIC(No,3),0.0_dp)
                Particles % Binit(No) = Particles % xpic(No,5)
             ELSE

                IF (.NOT. UpdateB) CYCLE
                IF (Particles % GMask(No) >= 0.9_dp) CYCLE
                IF (ANY(Particles % Dav(No,:).NE.0.0_dp)) CYCLE

                Particles % Binit(No) = Particles % xpic(No,4)
             END IF
          END DO
       END IF
    END IF

    IF (Particles % binitlowerboundfromfirsttimestep) THEN
       WHERE (Particles % binit(1:Particles % NumberOfParticles) .NE. &
            Particles % binit(1:Particles % NumberOfParticles)) &
            Particles % binit(1:Particles % NumberOfParticles) = Particles % XPIC(1:Particles % NumberofParticles,5)
       WHERE (Particles % binit(1:Particles % NumberOfParticles) .NE. &
            Particles % binit(1:Particles % NumberOfParticles)) &
            Particles % binit(1:Particles % NumberOfParticles) = min_first_binit
       WHERE (Particles % binit(1:Particles % NumberOfParticles) < min_first_binit) &
            Particles % binit(1:Particles % NumberOfParticles) = min_first_binit
    END IF

    IF (Particles % efboundsfromfirsttimestep) THEN
       WHERE (Particles % EF(1:Particles % NumberOfParticles) < min_first_ef) &
            Particles % EF(1:Particles % NumberOfParticles) = min_first_ef
       WHERE (Particles % EF(1:Particles % NumberOfParticles) > max_first_ef) &
            Particles % EF(1:Particles % NumberOfParticles) = max_first_ef
    END IF


    Particles % xpic = 0.0_dp

    !if needed, reset grounding mask from floatation solver
    IF (movegl) maskval = -maskval


    IF (Particles % VelocityDependentFriction) THEN
       IF (Particles % flipvelfric) THEN
          Particles % GridVelocity(1:Particles % NumberOfParticles,:) = &
               Particles % Velocity(1:Particles % NumberOfParticles,:)
       END IF
    END IF



    CALL Info(SolverName,&
         'Done with SSA Prep: Mesh to Particles',Level=1)

    !---------------- DONE ----------------!
  END SUBROUTINE SSAPrepMeshToParticles

  !**************************************************************************

  !> Borstad, 2016 constitutive damage
  !! Not as thoroughly tested as creep or zero-stress damage....
  SUBROUTINE BorstadDamage(Particles, No,srthres,GradVel,dDBorstad,n,k)

    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: EffSR,srthres,D,Dinit,GradVel(4)
    REAL(KIND=dp) :: dDBorstad, Exx,Eyy,Ezz,Exy,n,k

    Dinit = Particles % Dav(No,1)

    IF (SRThres < 0.0_dp) THEN
       D = 0.0_dp
    ELSE

       Exx=GradVel(1)
       Eyy=GradVel(2)
       Ezz=-Exx-Eyy

       IF (Ezz > 0.0_dp) THEN
          dDBorstad = 0.0_dp
          RETURN
       END IF

       Exy=0.5_dp*(GradVel(3) + GradVel(4))
       EffSR = sqrt(0.5_dp*(Exx*Exx + Eyy*Eyy + Ezz*Ezz) + Exy*Exy)
       D = 1.0_dp-((EffSR/SRThres)**(-1.0_dp/n)) &
            *EXP((-1.0_dp*(EffSR-SRThres))/(SRThres*(K-1.0_dp)))

       IF (D>Particles % Dmax) D=Particles % DMax
       IF (D<0.0_dp) D=0.0_dp

       IF (isnan(D)) THEN
          D = Particles % Dav(No,1)
          PRINT *,'damagenan no: ',no
          PRINT *,'EffSR:        ',EFFSR
          PRINT *,'SRThres :     ',SRThres
          PRINT *,'K:            ',K
          PRINT *,' '
       END IF

    END IF

    IF (D<Particles % Dav(No,1)) D = Particles % Dav(No,1)
    Particles % Dav(No,1) = D
    dDBorstad = (D - Dinit)

  END SUBROUTINE BorstadDamage

  !**************************************************************************

  SUBROUTINE removemaxd(D)

    implicit none
    REAL(KIND=dp) :: D(4),Q(2,2),EigValues(2),EigenVec(2,2)
    REAL(KIND=dp) :: zero=0.0_dp,w,y
    ! REAL(KIND=dp) :: x,z


    Q(1,1) = D(1)
    Q(2,2) = D(2)
    Q(1,2) = D(4)
    Q(2,1) = D(4)

    CALL Eigen2DSym_TryGenFirst(Q,EigValues,EigenVec)

    ! EigValues(2) = zero

    w = EigValues(1)*EigenVec(1,1)
    ! x = EigValues(2)*EigenVec(1,2)
    y = EigValues(1)*EigenVec(2,1)
    !  z = EigValues(2)*EigenVec(2,2)

    D(1) = EigenVec(1,1)*w !+ EigenVec(1,2)*x
    D(2) = EigenVec(2,1)*y !+ EigenVec(2,2)*z
    D(4) = EigenVec(2,1)*w !+ EigenVec(2,2)*x
    ! D(4) = D(1,2)


  END SUBROUTINE removemaxd

  !**************************************************************************

  ! Return eigenvalues (lambda) of 2x2 matrice a
  ! Return eigenvectors (columns of e) of 2x2 matrice a
  ! edited from original to return in ascending order
  SUBROUTINE Eigen2Db(a,lambdaout,eout)
    implicit none
    REAL(KIND=dp),intent(in) :: a(2,2)
    REAL(KIND=dp),intent(out) :: lambdaout(2)
    REAL(KIND=dp),intent(out) :: eout(2,2)

    REAL(KIND=dp) :: T,D,en,lambda(2),e(2,2),lvec(2,2),test(2,2),cs1,sn1
    REAL(KIND=dp),parameter :: Zero=1.0e-8  !Zero=100*AEPS
    INTEGER :: i

    REAL(KIND=dp) :: EigValues(2),EigenVec(2,2),WORK(68)
    INTEGER :: infor
    TYPE(Particle_t),  POINTER :: Particles

    Particles => GlobalParticles

    T=a(1,1)+a(2,2)
    D=a(1,1)*a(2,2)-a(1,2)*a(2,1)

    lambda(1)=0.5_dp*T+sqrt(0.25_dp*T*T-D)
    lambda(2)=0.5_dp*T-sqrt(0.25_dp*T*T-D)

    e=0.0_dp
    lvec=0.0_dp
    ! only 1 double eigenvalue (ONLY append when interscting isotropic
    ! metric??) - TO CHECK
    IF (abs((lambda(1)-lambda(2))/T).LT.Zero) THEN
       !IF (lambda(1)-lambda(2) == 0.0_dp) THEN
       e(1,1)=1.0_dp
       e(2,2)=1.0_dp
    ELSE IF (abs(a(1,2)/T).GT.Zero) THEN
       !first eigenvector
       e(1,1)=a(1,2)
       e(2,1)=lambda(1)-a(1,1)
       !2nd eigenvector
       e(1,2)=a(1,2)
       e(2,2)=lambda(2)-a(1,1)
    ELSE
       !  PRINT *,'alt'
       lambda(1)=a(1,1)
       lambda(2)=a(2,2)
       e(1,1)=1.0_dp
       e(2,2)=1.0_dp
    END IF

    DO i=1,2
       en=SQRT(SUM(e(:,i)*e(:,i)))
       e(:,i)=e(:,i)/en
    END DO

    lambdaout(1)=lambda(2)
    lambdaout(2)=lambda(1)

    eout(:,1) = e(:,2)
    eout(:,2) = e(:,1)

  END SUBROUTINE Eigen2Db

  !**************************************************************************

  ! Return eigenvalues (lambda) of 2x2 matrice a
  ! Return eigenvectors (columns of e) of 2x2 matrice a
  ! edited from original to return in ascending order
  SUBROUTINE Eigen2D(a,lambdaout,eout)
    implicit none
    REAL(KIND=dp),intent(in) :: a(2,2)
    REAL(KIND=dp),intent(out) :: lambdaout(2)
    REAL(KIND=dp),intent(out) :: eout(2,2)

    REAL(KIND=dp) :: T,D,en,lambda(2),e(2,2),lvec(2,2),test(2,2),cs1,sn1
    REAL(KIND=dp),parameter :: Zero=1.0e-6 !8  !Zero=100*AEPS
    INTEGER :: i

    REAL(KIND=dp) :: EigValues(2),EigenVec(2,2),WORK(68)
    INTEGER :: infor

    EigenVec = a
    CALL DSYEV('V', 'U', 2, EigenVec, 2, EigValues, Work, 68, infor)
    IF (infor.ne.0) CALL FATAL('Compute EigenValues', 'Failed to compute EigenValues')

    eout = EigenVec
    lambdaout = EigValues
  END SUBROUTINE Eigen2D

  !**************************************************************************

  SUBROUTINE Eigen2DSym_TryGenFirstValOnly(a,EigValues)

    IMPLICIT NONE

    REAL(KIND=dp),intent(in) :: a(2,2)
    REAL(KIND=dp),intent(out) :: EigValues(2)
    REAL(KIND=dp) :: T,D,en,sqrteig
    REAL(KIND=dp),parameter :: zero=0.0_dp, thres = 0.0001_dp, one = 1.0_dp
    REAL(KIND=dp),parameter :: half = 0.5_dp, quart = 0.25_dp
    INTEGER :: ii

    T=a(1,1)+a(2,2)
    D=a(1,1)*a(2,2)-a(1,2)*a(2,1)

    !eigenvalues:
    !2 is the biggest
    sqrteig = quart*T*T-D
    IF (sqrteig<0.0_dp) sqrteig = 0.0_dp
    sqrteig = sqrt(sqrteig)
    T = half*T
    EigValues(1)=T-sqrteig
    EigValues(2)=T+sqrteig


  END SUBROUTINE Eigen2DSym_TryGenFirstValOnly

  !**************************************************************************

  SUBROUTINE MaxPFour(a,lmaxval)

    IMPLICIT NONE

    REAL(KIND=dp),intent(in) :: a(4)
    REAL(KIND=dp),intent(out) :: lmaxval
    REAL(KIND=dp) :: T,D,en,sqrteig
    REAL(KIND=dp),parameter :: zero=0.0_dp, thres = 0.0001_dp, one = 1.0_dp
    REAL(KIND=dp),parameter :: half = 0.5_dp, quart = 0.25_dp
    INTEGER :: ii

    T=a(1)+a(2)
    D=a(1)*a(2)-a(4)*a(4)

    !eigenvalues:
    !2 is the biggest
    sqrteig = quart*T*T-D
    IF (sqrteig<0.0_dp) sqrteig = 0.0_dp
    sqrteig = sqrt(sqrteig)
    ! T = half*T
    !  EigValues(1)=T-sqrteig
    lmaxval=half*T+sqrteig


  END SUBROUTINE MaxPFour

  !**************************************************************************

  SUBROUTINE Eigen2DSym_TryGenFirst(a,EigValues,EigenVec)

    IMPLICIT NONE

    REAL(KIND=dp),intent(in) :: a(2,2)
    REAL(KIND=dp),intent(out) :: EigValues(2)
    REAL(KIND=dp),intent(out) :: EigenVec(2,2)
    REAL(KIND=dp) :: T,D,en,sqrteig
    REAL(KIND=dp),parameter :: zero=0.0_dp, thres = 0.0001_dp, one = 1.0_dp
    REAL(KIND=dp),parameter :: half = 0.5_dp, quart = 0.25_dp
    INTEGER :: ii

    T=a(1,1)+a(2,2)
    D=a(1,1)*a(2,2)-a(1,2)*a(2,1)

    !eigenvalues:
    !2 is the biggest
    sqrteig = quart*T*T-D
    IF (sqrteig<0.0_dp) sqrteig = 0.0_dp

    sqrteig = sqrt(sqrteig)
    T = half*T
    EigValues(1)=T-sqrteig
    EigValues(2)=T+sqrteig

    IF (ANY(EigValues .NE. EigValues)) THEN
       CALL Eigen2D(a,EigValues,EigenVec)
       RETURN
    END IF


    IF (a(1,2)==zero) THEN
       IF (a(2,2)>a(1,1)) THEN
          EigenVec(1,1) = one
          EigenVec(2,1) = zero
          EigenVec(1,2) = zero
          EigenVec(2,2) = one
       ELSE
          EigenVec(1,1) = zero
          EigenVec(2,1) = one
          EigenVec(1,2) = one
          EigenVec(2,2) = zero
       END IF
    ELSE
       IF ( ABS(EigValues(2)-EigValues(1)) < thres .OR. ABS(a(1,2)) < thres ) THEN

          !Use dsyev to be sure...general solution might run into
          !precision/overflow/underflow issues under these conditions
          CALL Eigen2D(a,EigValues,EigenVec)
       ELSE
          !first eigenvector
          EigenVec(1,2)=a(1,2)
          EigenVec(2,2)=EigValues(2)-a(1,1)

          !2nd eigenvector
          EigenVec(1,1)=a(1,2)
          EigenVec(2,1)=EigValues(1)-a(1,1)

          DO ii=1,2
             en=SQRT(SUM(EigenVec(:,ii)*EigenVec(:,ii)))
             EigenVec(:,ii)=EigenVec(:,ii)/en
          END DO
       END IF
    END IF


  END SUBROUTINE Eigen2DSym_TryGenFirst

  !**************************************************************************

  SUBROUTINE Eigen2DSym_TryGenFirst_VecOnly(a,lambda,EigValues,EigenVec)

    IMPLICIT NONE

    REAL(KIND=dp),intent(in) :: a(2,2),lambda(2)
    REAL(KIND=dp),intent(out) :: EigValues(2),EigenVec(2,2)
    REAL(KIND=dp),parameter :: Zero=0.0_dp, thres = 0.0001_dp, one = 1.0_dp
    REAL(KIND=dp),parameter :: half = 0.5_dp, quart = 0.25_dp
    REAL(KIND=dp) :: en
    INTEGER :: ii

    ! T=a(1,1)+a(2,2)
    ! D=a(1,1)*a(2,2)-a(1,2)*a(2,1)

    ! !eigenvalues:
    ! !2 is the biggest
    ! sqrteig = sqrt(quart*T*T-D)
    ! T = half*T
    ! EigValues(1)=T-sqrteig
    ! EigValues(2)=T+sqrteig

    EigValues(2) = lambda(1)
    EigValues(1) = lambda(2)


    IF (a(1,2)==zero) THEN
       IF (a(2,2)>a(1,1)) THEN
          EigenVec(1,1) = one
          EigenVec(2,1) = zero
          EigenVec(1,2) = zero
          EigenVec(2,2) = one
       ELSE
          EigenVec(1,1) = zero
          EigenVec(2,1) = one
          EigenVec(1,2) = one
          EigenVec(2,2) = zero
       END IF
    ELSE
       IF ( ABS(EigValues(2)-EigValues(1)) < thres .OR. ABS(a(1,2)) < thres ) THEN
          !Use dsyev to be sure...general solution might run into
          !precision/overflow/underflow issues under these conditions
          CALL Eigen2D(a,EigValues,EigenVec)
       ELSE
          !first eigenvector
          EigenVec(1,2)=a(1,2)
          EigenVec(2,2)=EigValues(2)-a(1,1)

          !2nd eigenvector
          EigenVec(1,1)=a(1,2)
          EigenVec(2,1)=EigValues(1)-a(1,1)

          DO ii=1,2
             en=SQRT(SUM(EigenVec(:,ii)*EigenVec(:,ii)))
             EigenVec(:,ii)=EigenVec(:,ii)/en
          END DO
       END IF
    END IF


  END SUBROUTINE Eigen2DSym_TryGenFirst_VecOnly


  !**************************************************************************
  !<<<<<<<<<<<<<<<<<<<<<< MODULE MPMUtils <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  !**************************************************************************

END MODULE MPMUtils

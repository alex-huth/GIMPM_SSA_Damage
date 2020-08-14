FUNCTION hoopsurf(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(2),VarOut,x,y,r,buffer,yc
  LOGICAL:: Visited,GotIt
  
  r = 70000.0_dp 
  x = ABS(VarIn(1))
  y = ABS(VarIn(2))

  IF (x>r .OR. y>r) THEN
     VarOut = -1.0_dp
  ELSE
     yc = sqrt(y*y+x*x)

     IF (yc<r) THEN
        VarOut = 1.0_dp
     ELSE
        VarOut = -1.0_dp
     END IF     
  END IF

End FUNCTION hoopsurf


!x,y,v1,v2
FUNCTION upstream_xvel_adj(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(3),VarOut,x,y,hyp,theta,mag

  x = VarIn(1)
  y = ABS(VarIn(2))
  mag = VarIn(3)

  IF (x == 0.0_dp) THEN
     VarOut = 0.0_dp
  ELSE
     theta = y/x
     VarOut = mag/sqrt(1.0_dp+theta*theta)
  END IF

End FUNCTION upstream_xvel_adj


!x,y,v1,v2
FUNCTION upstream_yvel_adj(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(3),VarOut,x,y,hyp,theta,mag

  x = VarIn(1)
  y = VarIn(2)
  mag = VarIn(3)
  
  IF (y<0.0_dp) THEN
     y = -y
     mag = -mag
  END IF
  
  IF (x == 0.0_dp) THEN
     VarOut = mag     
  ELSE
     theta = y/x
     VarOut = mag * theta/sqrt(1.0_dp+theta*theta)
  END IF

End FUNCTION upstream_yvel_adj


FUNCTION constantvarssmpm(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(2),VarOut,x,y,r1,r2,buffer,dist
  LOGICAL:: Visited,GotIt

  x = VarIn(1)
  y = VarIn(2)

  dist = sqrt(x*x+y*y)

  IF (dist > 69999.0_dp .AND. dist < 70001.0_dp) THEN
     VarOut = 1.0_dp
  ELSE
     VarOut = -1.0_dp
  END IF
  
End FUNCTION constantvarssmpm

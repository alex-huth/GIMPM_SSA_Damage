FUNCTION initsurf(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(1),VarOut

  IF (VarIn(1) < 0.0_dp) THEN
     VarOut = 1.0_dp
  ELSE
     VarOut = -1.0_dp
  END IF
 
End FUNCTION initsurf

FUNCTION initH(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(1),VarOut,H0,V0,C,Q0,termA,TermB,x,gridres,H,dH
  LOGICAL :: found  

  gridres = GetConstReal( Model % Constants,'Grid Resolution',Found )
  IF (.NOT. Found) CALL Fatal('USF_1dtest:', &
       'initH: Need to define "gridres = Real $" in constants')

  H0 = 600.0_dp
  v0 = 300.0_dp

  C = (((910.0_dp*9.81_dp/(4.0_dp*1.9E8_dp))*&
       (1.0_dp-910.0_dp/1028.0_dp))**3.0_dp)*31556926.0_dp  
  Q0 = H0*v0

  termA = 4.0_dp*C/Q0
  termB = 1.0_dp/(H0*H0*H0*H0)

  IF (VarIn(1) <= 0.0_dp) THEN
     VarOut = 600.0_dp
  ELSE
     VarOut = (termA * VarIn(1) + termB)**(-0.25_dp)
  END IF

End FUNCTION initH


FUNCTION initVel(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(1),VarOut,H0,V0,C,Q0,termA,TermB,x,gridres,vel,dvel
  LOGICAL :: found

  gridres = GetConstReal( Model % Constants,'Grid Resolution',Found )
  IF (.NOT. Found) CALL Fatal('USF_1dtest:', &
       'initvel: Need to define "gridres = Real $" in constants')

  H0 = 600.0_dp
  v0 = 300.0_dp

  C = (((910.0_dp*9.81_dp/(4.0_dp*1.9E8_dp))*&
       (1.0_dp-910.0_dp/1028.0_dp))**3.0_dp)*31556926.0_dp

  Q0 = H0*v0

  termA = 4.0_dp*C/Q0
  termB = 1.0_dp/(H0*H0*H0*H0)

  IF (VarIn(1) <= 0.0_dp) THEN
     VarOut = 300.0_dp
  ELSE
     VarOut = Q0/((termA * VarIn(1) + termB)**(-0.25_dp))
  END IF

End FUNCTION initVel

FUNCTION getpassive(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(1),VarOut

  VarOut=(VarIn(1))*(-1.0_dp)
End FUNCTION getpassive



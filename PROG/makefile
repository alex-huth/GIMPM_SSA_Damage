all: MPM

MPM: MPM_Utils.F90  MPM_Procedures.F90  MPM_SSA.F90 
	elmerf90 $^ -o $@ 

clean:
	rm -f MPM
	rm -f *.o
	rm -f *.mod

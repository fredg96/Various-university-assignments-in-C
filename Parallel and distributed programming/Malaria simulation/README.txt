To compile program: 
	1. make sure that gcc and openmpi is installed.
	2. use the makefile by writing "make" in prompt.

To run program:
	1. program has two inputs N, nr of MC experiments has to be divisible by number of PEs, and name of output file.
	   example: mpirun -n 16 ./malaria 2000 "output.txt"

To remove program:
	1. write "make clean" in prompt.
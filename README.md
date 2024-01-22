# LibForQ-v1
Fortran code for generating random numbers, random probability vectors, unitaries, and quantum states

LibForQ: A Fortran Library for Quantum Information Science
In this first version, generators for random numbers, probability vectors, unitary matrices, state vectors, 
and density matrices are provided. Some matrix functions are also included (see the files matfun.f90 and qnesses.f90).

One can use this library simply by copying all the files to the main program folder and compiling them with:
  gfortran -lblas -llapack *.f90
To run the executable file a.out generated just use:
  ./a.out

Another, more convenient, manner of using the code is by creating a library. For that purpose you may follow 
the simple steps below:
1) Extract code.zip
2) Go to the folder libforq
3) Create a library with the commands:
  gfortran -O2 -c *.f90
and
  ar cr libforq.a *.o *.mod 
To compile the main program and the tests using this library, copy libforq.a to your program's folder and use the command: 
  gfortran -lblas -llapack libforq.a tests.f90 main.f90
You can also copy it to your computer libraries folder, e.g. with:
  sudo cp libforq.a /usr/local/lib
and use it, anywhere in your computer, via
  gfortran -lblas -llapack -lforq tests.f90 main.f90 

REMARK. Regarding the tests provided in the file tests.f90, you can uncomment, in the main program, the test to be 
performed and see the details in the corresponding subroutine. Some basic Gnuplot commands to see the results are also 
included there. You can also, of course, erase the program main.f90 (and eventually also tests.f90) and use your own 
main program.

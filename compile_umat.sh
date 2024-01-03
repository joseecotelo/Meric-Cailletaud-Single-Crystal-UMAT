gfortran -c -O2 -fcheck=all umat_subroutines.f MericCailletaud_module.f math_subroutines.f
gfortran MericCailletaudTest.f -fcheck=all umat_subroutines.o  math_subroutines.o MericCailletaud_module.o
# ./a.out

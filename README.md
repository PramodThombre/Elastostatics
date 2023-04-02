# Elastostatics
This code was developed as a part of the assignments for the Introduction to Finite Element Methods course.
It uses the open source Eigen library in cpp for handling matrices, vectors. Eigen library has to be installed for running this code.

Finite element (FE) code in C++ to solve:

Part1 is a 1D elastostatics problem with the following data:
E = 1011 Pa, A = 10−4 m2,  ̃f = 106 Nm−1,  ̄f = 107 Nm−2, L = 0.1 m, g1 = 0, g2 = 0.001 m, h = 106 N
Discretiztion is done using 
1) Linear lagrange shape functions with 3 elements.
2) Linear lagrange shape functions with 100 elements.

Part2 is 2D transient heat conduction problem.

Part3 is 3D Elasto statics problem.

# isoparametric-FEM

C++ finite element method implementation of a simple beam structure under applied load using isoparametric bilinear quadrilateral elements of the Lagrange type. The constitutive equations are established from Hamilton's principle and the Galerkin method of weighted residuals is used. Numerical integration is performed using Gauss quadrature. The loads, stiffness, boundary conditions, mesh definitions and beam geometry can be applied in main.cpp. The code outputs the maximum and minimum displacements, as well as .txt files with all nodal coordinates. To visualise the nodal displacements, the C++ implementation can be compiled and run from the code.py file (Python wrapper file) which additionally plots the undeformed/deformed finite element mesh (requires ideally Python >3.7 and matplotlib installation for plotting).

Final Project for AMCS-329 at KAUST: Finite Elements
====================================================

***Implementation of a Finite Elements solver in 2D from scratch***

**Author:** *Me*

The goal is to solve Poisson's equation \lambda(u)=f(x) over a simple domain.
For the moment I will assume Dirichlet boundary conditions and the mesh will
be given a priori by some triangulation process (probably Delaunay).
The method will be standard (continuous) FE without using an specialized
library (e.g. FENICS or DEAL II). 
I will use a Lagrange basis (the degrees of freedom will be the nodal values)
and P1 elements. Maybe I'll switch to higher order in the future.
The error (maybe also rate of convergence?) will be measured with some grid 
function norm against an exact solution).

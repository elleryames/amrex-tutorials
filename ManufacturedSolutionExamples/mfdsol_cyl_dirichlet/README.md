README

This example verifies the code by method of manufactured solutions by solving the Poisson problem

    - phi_rr - phi_zz - phi_r/r = f

where

    f = 2*r*r(1 - r) + (9*r - 4)*(1 + z)*(1 - z)

on the domain (r,z) in [0,1]x[-1,1], 
and with homogeneous Dirichlet boundary conditions on all boundaries. 

The exact solution is

    phi = r*r*(1-r)*(1-z)*(1+z).


We use the ABecLaplace solver, rather than the Poisson solver for practice. 


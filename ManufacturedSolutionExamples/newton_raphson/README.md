README

This example verifies the code by method of manufactured solutions by solving the Poisson problem

    - phi_xx - phi_yy  + phi * phi = f

where

    f = x*(1.-x) * y*(1.-y)

on the domain [0,1]x[0,1], and with homogeneous boundary conditions.

The nonlinearity is treated with an outer Newton-Raphson iteration. 

We use the ABecLaplace solver. 


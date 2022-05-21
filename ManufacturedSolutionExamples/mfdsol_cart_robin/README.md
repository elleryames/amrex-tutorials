README

This example verifies the code by method of manufactured solutions by solving the Poisson problem

    - phi_xx - phi_yy = f

where

    f = (1 + pi*pi) * sin(x) * sin(pi * y)

on the domain [0,1]x[0,1], and with the boundary conditions: 

    1. Homogenous Dirichlet at x = 0, y = 0, and y = 1.
    2. Robin condition at x = 1 given by 

            robin_a * phi + robin_b * dphidx = g

            robin_a = 1
            robin_b = 1
            g = robin_a * sin(x)*sin(pi * y) 
                    + robin_b * cos(x)*sin(pi * y)

The exact solution is 

    phi = sin(x) * sin(pi * y).

We use the ABecLaplace solver, rather than the Poisson solver for practice. 


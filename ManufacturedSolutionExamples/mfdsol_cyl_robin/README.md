README

This example verifies the code by method of manufactured solutions by solving the Poisson problem

    - phi_rr - phi_zz - phi_r/r = f

where

    f = (6. * s2 - 4. *(r*r + z*z)) / s2 / s2 * exp(-(r*r+z*z)/s2)
    s2 = 0.25

on the domain (r,z) in [0,1]x[-1,1], and with the boundary conditions: 

    1. Homogenous Neumann at r = 0.
    2. Robin condition at the outer boundary given by 

            robin_a * phi + robin_b * dphidn = g

            robin_a = 1
            robin_b = 1
            g = robin_b * phi + robin_b * (r * dphidr + z * dphidz)

        Note that the derivative term in the Robin is radial 
            <r, z> . grad(phi).

The exact solution is 

    phi = exp(-(r*r + z*z)/s2).

We use the ABecLaplace solver, rather than the Poisson solver for practice. 


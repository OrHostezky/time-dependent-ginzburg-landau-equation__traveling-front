# Time-Dependent Ginzburg-Landau Equation - Traveling Front

In this work, basic simulations and analysis are carried out of a certain traveling-front case of the 1-dimensional Time-Dependent Ginzburg-Landau Equation (TDGLE). This equation is a PDE that involves diffusion and an additional reaction(s), namely
```math
\frac{\partial u}{\partial t} = D\frac{\partial^2 u}{\partial x^2}+f(u)    \quad ,
```
where $u(x,t)$ can be regarded as a particle density, $D$ is the diffusion coefficient and $f(u)$ is the reaction function, which describes additionl, density-dependent effects in the system.
Specifically, the TDGLE involves a *bistable* reaction function, which by itself has two stable fixed points, introducing two stable phases to which the solution tends at long times. 
Here, a generic cubic reaction function is considered, namely $f(u)=-k(u-u_1)(u-u_2)(u-u_3)$ for $k>0$ and $u_1 < u_2 < u_3$ (a detailed chemical intuition can be drawn for $f(u)$, see [*included paper*](Time_Dependent_Ginzburg_Landau_Equation___Traveling_Front.pdf)).

Due to the fact that there are two stable phases in the infinite problem, a competition between the two can be generated by feeding the equation with initial conditions that have certain properties, which then creates a traveling-front dynamical scheme, such that the solution is of the form $u(x,t)=u(x-ct)$ for some velocity $c$. Subsequently, one of the phases eventually takes over as $t\to\inf$. The asymptotic form and velocity of the traveling front are compared with the theoretical predictions, and the invariance of the asymptotic front to initial condition is demonstrated. In addition, both the dependencies of the wave velocity and convergence time on the reaction function are characterized. 

The dynamics are numerically solved for using *MATLAB*'s built-in PDE solver `pdepe` solver, which discretizes space and integrates the ODEs resulting from this spatial discretization to obtain the approximate solutions at the desired times.

## Code


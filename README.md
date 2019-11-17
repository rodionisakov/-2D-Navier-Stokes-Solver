# -2D-Navier-Stokes-Solver

To run: open all .m Matlab files and run Main_code.m

About:
This is an incompressible 2D Navier-Stokes solver in MATLAB with simple numerical methods. The code discretises the equations on a
collocated uniform grid distribution. Time marching is formulated with Runge-Kutta third order scheme. All derivatives are approximated by
cell-centred schemes. The cylinder is modelled by immersed-boundary volume of fluid method. 

The aim is to solve for fluid flow in 2D around and in the wake of a square cylinder with dimension of D = 20 cm
placed in an approaching stream withvelocity V = 5 cm/s and the kinematic viscosity of the fluid is v = 10^-4 m^2/s. 

The simulation box (shown in Figure 1.PNG) is a rectangle with dimensions Lx = 5 m
and Ly = 1 m, and the centre of the cylinder is located at xc = 70 cm and yc = 50 cm.

Assumed that the initial velocity components are u = 5 cm/s and v = 0. 
The baseline grid (used to discretize the computational domain) is
partitioned in small square cells of size ∆x = ∆y = 1 cm. The time step is
automatically adjusted by the code based on the CFL condition. The inlet
boundary condition is dictated by the constant approaching flow, while the
outlet obeys Orlansky convective velocity condition. The top and the bottom
of the domain is assumed to be free-slip.

Results demonstrate(Results descriptiion.pdf):
A non-dimensional form of the governing equations and a discussion on
the non-dimensional parameter that appear.

Distribution of u and v velocity components as well as the
streamlines at following time instants using default MATLAB plot
settings: t = 0, 5, 10, 20, 40, 100 s .

Strouhal number related to the vortex shedding.

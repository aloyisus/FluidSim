FluidSim
========

Fluid simulator based on Robert Bridson's Siggraph notes, which can be found at:

http://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf


This project goes about as far as Chapter Five, which is to say it simulates a single-phase, incompressible fluid with no free surface, being entirely confined to a box.

We do this by numerically solving the incompressible Navier-Stokes equations for fluid flow. The solution is broken up into an Advection step followed by a Projection step.

The advection step is the semi-Lagrangian method introduced to computer graphics by Jos Stam, and the projection step solves a system of N (where N is the number of cells in the grid) linear equations using the pre-conditioned conjugate gradients (PCG) method.

I have also added support for curved solid boundaries. So far this is either a cuboid or a sphere.

The result at each step is written out to an openVDB container, which means we can use the openVDB renderer, or we can import the cache into Houdini to view.

To solver can be run from the commandline, type 'fluidsim --help' for available options.

To view some sample renders, go to my vimeo page at https://vimeo.com/johnakelly/videos/

![screenshot](https://raw.githubusercontent.com/aloyisus/FluidSim/master/rotating_paddle-high.gif)

Note: this project makes use of the open-source Eigen library for linear algebra and matrix operations.

FluidSim
========

Fluid simulator based on Robert Bridson's Siggraph notes, which can be found at:

http://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf

This project goes about as far as Chapter Five, which is to say it simulates a single-phase, incompressible fluid with no free surface, being entirely confined to a box.

This is achieved by numerically solving the incompressible Navier-Stokes equations for fluid flow. The solution is split into an Advection step followed by a Projection step:
 - The advection step is the semi-Lagrangian method (Jos Stam)
 - the projection step solves a system of linear equations using the pre-conditioned conjugate gradients (PCG) method.

The result at each step is written out to an openVDB container, which can be imported into e.g. Houdini to view and render.

The solver can be run from the commandline (type 'fluidsim --help' for available options), which allows some limited control of the volume size and time step. Test solids can be added by modifying main.cpp directly using the setSolid method on the FluidGrid class.

To view some sample renders, go to my vimeo page at https://vimeo.com/johnakelly/videos/

![screenshot](https://raw.githubusercontent.com/aloyisus/FluidSim/master/rotating_paddle-high.gif)

UPDATE 21/8/23: I made this nine years ago as a toy project to help me learn some concepts in CFD. My original implementation used custom code for the linear algebra operations and the PCG solver, but recently I've been learning OpenVDB so I started using that library's grid representation and PCG solver instead. Further work would likely delegate the buffering and advection to OpenVDB.
Be warned, I suspect there are some serious bugs, particularly relating to the handling of solid boundary conditions. Testing would help to fix these, but this would require a good understanding of the expected state of the grid at each stage of the solve, and I'm still trying to achieve that :) 
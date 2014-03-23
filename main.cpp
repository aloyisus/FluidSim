//Copyright (c) 2014 John Kelly
//
//
//
//


#include "grid.h"




int main() {

    int w = 100;
    int h = 80;
    int d = 100;
    
    double density = 0.1;
    double timestep = 0.02;
    
    
    FluidGrid* solver = new FluidGrid(w,h,d,timestep,density);
    //solver->setSolid(new Cuboid(solver, w/2,80,d/2,64,8,32,0.0,1.0));

    solver->setSolid(new Sphere(solver, w/2,h/2,d/2,16));
    
    while (solver->getCurrtime() < 5.0) {
        
        solver->Update();
    }
    
    return 0;
}
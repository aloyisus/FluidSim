//Copyright (c) 2014 John Kelly
//
//
//
//


#include "grid.h"




int main() {

    int w = 80;
    int h = 128;
    int d = 80;
    
    double density = 0.1;
    double timestep = 0.02;
    
    FluidGrid* solver = new FluidGrid(w,h,d,timestep,density);
    SolidCuboid cuboid(40,80,40,64,8,32,0.125,0.1);
    
    solver->solid = &cuboid;
    
    while (solver->getSimtime() < 10.0) {
        
        solver->Update();
    }
    
    return 0;
}
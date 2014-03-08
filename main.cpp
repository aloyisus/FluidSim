//Copyright (c) 2014 John Kelly
//
//
//
//


#include "grid.h"




int main() {

    int w = 128;
    int h = 128;
    int d = 128;
    
    double density = 0.1;
    double timestep = 0.02;
    
    FluidGrid* solver = new FluidGrid(w,h,d,timestep,density);

    while (solver->getSimtime() < 10.0) {
        
        solver->Update();

    }
    
    return 0;
}
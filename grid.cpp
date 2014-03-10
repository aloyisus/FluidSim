 //
//  grid.cpp
//  Fluid
//
//  Created by John Kelly on 26/01/2014.
//  Copyright (c) 2014 John Kelly. All rights reserved.
//

#include "grid.h"


#define round5(n) floor(n * 100000 + 0.5)/100000

#define filepath "/Users/JohnnyK/Downloads/Fluids14/"

#define filepath "/Users/JohnnyK/Downloads/Fluids17/"


//Constructor for FluidGrid 3d
FluidGrid::FluidGrid(int wh, int ht, int dh, double tstep, double rh){
    
    cellwidth = 1.0/std::min(std::min(wh, ht),dh);

    nwidth = wh;
    nheight = ht;
    ndepth = dh;
    
    rho = rh;
    
    CFLnumber = 8;
    framenumber = 0;
    currtime = 0;
    framedeltaT = tstep;
    nextframetime = tstep;
    bool writetocache = false;
    
    pressure = new vectord(wh*ht*dh);
    u = new FluidQuantity(0,0.5,0.5,wh,ht,dh,cellwidth);
    v = new FluidQuantity(0.5,0,0.5,wh,ht,dh,cellwidth);
    w = new FluidQuantity(0.5,0.5,0,wh,ht,dh,cellwidth);
    density = new FluidQuantity(0.5,0.5,0.5,wh,ht,dh,cellwidth);
    r = new vectord(wh*ht*dh);
    A = new symmbandd(wh*ht*dh,7,1,wh,wh*ht);
    Ei = new diagonald(wh*ht*dh);
    
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();
    // Create an empty floating-point grid.
    grid = openvdb::FloatGrid::create();
    
}


//Constructor for FluidGrid 2d
FluidGrid::FluidGrid(int wh, int ht, double tstep, double rh){
    
    cellwidth = 1.0/std::min(wh, ht);
    
    nwidth = wh;
    nheight = ht;
    
    rho = rh;
    
    CFLnumber = 8;
    framenumber = 0;
    currtime = 0;
    framedeltaT = tstep;
    nextframetime = tstep;
    bool writetocache = false;
    
    pressure = new vectord(wh*ht);
    u = new FluidQuantity(0,0.5,wh,ht,cellwidth);
    v = new FluidQuantity(0.5,0,wh,ht,cellwidth);
    density = new FluidQuantity(0.5,0.5,wh,ht,cellwidth);
    r = new vectord(wh*ht);
    A = new symmbandd(wh*ht,5,1,wh,0);
    Ei = new diagonald(wh*ht);
    
    image = new unsigned char[wh*ht*4];
    
}


//destructor
FluidGrid::~FluidGrid(){
    
    delete u;
    delete v;
    delete w;
    delete density;
    delete pressure;
    delete r;
    delete Ei;
    delete A;
       
}



//This builds the pressure coefficients for a 2d grid
void FluidGrid::setupA2d(){
    
    double scale =  deltaT/(rho*cellwidth*cellwidth);
    
    int jnwidth;
    
    A->fill(0);
    
    for(int j=0; j<nheight; j++)
        for(int i=0; i<nwidth; i++){
            
            jnwidth = j*nwidth;
            
            if (i > 0)
                A->at(i+jnwidth,0) += scale;
            
            if (j > 0)
                A->at(i+jnwidth,0) += scale;
            
            if (i < nwidth-1){
                A->at(i+jnwidth,0) += scale;
                A->at(i+jnwidth,1) = -scale;
            }
            
            if (j < nheight-1){
                A->at(i+jnwidth,0) += scale;
                A->at(i+jnwidth,2) = -scale;
            }
            
        };

}


//This builds the pressure coefficients for a 3d grid
void FluidGrid::setupA3d(){
    
    double scale =  deltaT/(rho*cellwidth*cellwidth);
    
    int jnwidth,knwidthnheight;
    
    A->fill(0);
    
    for(int k=0; k<ndepth; k++)
        for(int j=0; j<nheight; j++)
            for(int i=0; i<nwidth; i++){
                
                jnwidth = j*nwidth;
                knwidthnheight = k*nwidth*nheight;
                
                if (i > 0)
                    A->at(i+jnwidth+knwidthnheight,0) += scale;
                
                if (j > 0)
                    A->at(i+jnwidth+knwidthnheight,0) += scale;
                
                if (k > 0)
                    A->at(i+jnwidth+knwidthnheight,0) += scale;
                
                if (i < nwidth-1){
                    A->at(i+jnwidth+knwidthnheight,0) += scale;
                    A->at(i+jnwidth+knwidthnheight,1) = -scale;
                }
                
                if (j < nheight-1){
                    A->at(i+jnwidth+knwidthnheight,0) += scale;
                    A->at(i+jnwidth+knwidthnheight,2) = -scale;
                }
                
                if (k < ndepth-1){
                    A->at(i+jnwidth+knwidthnheight,0) += scale;
                    A->at(i+jnwidth+knwidthnheight,3) = -scale;
                }
            
            
            };
    
}


void FluidGrid::Advect2d(){
    
    density->Advect(deltaT,*u,*v);
    u->Advect(deltaT,*u,*v);
    v->Advect(deltaT,*u,*v);
    
    u->swap();
    v->swap();
    density->swap();
    
}


void FluidGrid::Advect3d(){
    
    density->Advect(deltaT,*u,*v,*w);
    u->Advect(deltaT,*u,*v,*w);
    v->Advect(deltaT,*u,*v,*w);
    w->Advect(deltaT,*u,*v,*w);
    
    u->swap();
    v->swap();
    w->swap();
    density->swap();
    
}


//Construct RHS of pressure equations for 2d grid
void FluidGrid::setupRHS2d() {
    double scale = 1.0/cellwidth;
    
    for (int j = 0; j < nheight; j++) {
        for (int i = 0; i < nwidth; i++) {

            setr(i,j) = -scale*(u->getQuantity(i+1,j) - u->getQuantity(i,j) +
                                v->getQuantity(i,j+1) - v->getQuantity(i,j));
        }
    }
    
    
    //Solid wall boundary conditions
    for (int j = 0; j < nheight; j++){
        setr(0,j) = -scale*u->getQuantity(0,j);
        setr(nwidth-1,j) = +scale*u->getQuantity(nwidth,j);
        setr(0,j) += -scale*u->getQuantity(0,j);
        setr(nwidth-1,j) += +scale*u->getQuantity(nwidth,j);
    }
    for (int i = 0; i < nwidth; i++){
        setr(i,0) = -scale*v->getQuantity(i,0);
        setr(i,nheight-1) = +scale*v->getQuantity(i,nheight);
        setr(i,0) += -scale*v->getQuantity(i,0);
        setr(i,nheight-1) += +scale*v->getQuantity(i,nheight);
    }
    
    
}

//Construct RHS of pressure equations for 3d grid
void FluidGrid::setupRHS3d() {
    double scale = 1.0/cellwidth;
    
    //see Bridson p.45
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++)
            for (int i = 0; i < nwidth; i++) {
                
                setr(i,j,k) = -scale*(u->getQuantity(i+1,j,k) - u->getQuantity(i,j,k) +
                                      v->getQuantity(i,j+1,k) - v->getQuantity(i,j,k) +
                                      w->getQuantity(i,j,k+1) - w->getQuantity(i,j,k));
            }

    

    //Solid wall boundary conditions see Bridson p.49
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++){
            setr(0,j,k) = -scale*u->getQuantity(0,j,k);
            setr(nwidth-1,j,k) = +scale*u->getQuantity(nwidth,j,k);
            setr(0,j,k) += -scale*u->getQuantity(0,j,k);
            setr(nwidth-1,j,k) += +scale*u->getQuantity(nwidth,j,k);
        }
    for (int k = 0; k < ndepth; k++)
        for (int i = 0; i < nwidth; i++){
            setr(i,0,k) = -scale*v->getQuantity(i,0,k);
            setr(i,nheight-1,k) = +scale*v->getQuantity(i,nheight,k);
            setr(i,0,k) += -scale*v->getQuantity(i,0,k);
            setr(i,nheight-1,k) += +scale*v->getQuantity(i,nheight,k);
        }
    for (int j = 0; j < nheight; j++)
        for (int i = 0; i < nwidth; i++){
            setr(i,j,0) = -scale*w->getQuantity(i,j,0);
            setr(i,j,ndepth-1) = +scale*w->getQuantity(i,j,ndepth);
            setr(i,j,0) += -scale*w->getQuantity(i,j,0);
            setr(i,j,ndepth-1) += +scale*w->getQuantity(i,j,ndepth);
        }
   
    
    
    
}


void FluidGrid::updateVelocities2d() {
    double scale = deltaT/(rho*cellwidth);
    
    //subtract pressure gradients to update velocities
    for (int j = 0; j < nheight; j++)
        for (int i = 0; i < nwidth; i++) {
            
            double cellpressure = scale*getPressure(i,j);
            u->setQuantity(i,j) -= cellpressure;
            u->setQuantity(i+1,j) += cellpressure;
            v->setQuantity(i,j) -= cellpressure;
            v->setQuantity(i,j+1) += cellpressure;
            
        }
  
    
    //Solid wall boundary conditions
    for (int j = 0; j < nheight; j++){
        u->setQuantity(0,j) = 0.0;
        u->setQuantity(nwidth,j) = 0.0;
    }
    for (int i = 0; i < nwidth; i++){
        v->setQuantity(i,0) =  0.0;
        v->setQuantity(i,nheight) = 0.0;
    }
  
    
}


void FluidGrid::updateVelocities3d() {
    double scale = deltaT/(rho*cellwidth);
    
    //subtract pressure gradients to update velocities
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++)
            for (int i = 0; i < nwidth; i++) {
                
                double cellpressure = scale*getPressure(i,j,k);
                u->setQuantity(i,j,k) -= cellpressure;
                u->setQuantity(i+1,j,k) += cellpressure;
                v->setQuantity(i,j,k) -= cellpressure;
                v->setQuantity(i,j+1,k) += cellpressure;
                w->setQuantity(i,j,k) -= cellpressure;
                w->setQuantity(i,j,k+1) += cellpressure;
                
            }
    
    
    //Solid wall boundary conditions
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++){
            u->setQuantity(0,j,k) = 0.0;
            u->setQuantity(nwidth,j,k) = 0.0;
        }
    for (int k = 0; k < ndepth; k++)
        for (int i = 0; i < nwidth; i++){
            v->setQuantity(i,0,k) =  0.0;
            v->setQuantity(i,nheight,k) = 0.0;
        }
    for (int j = 0; j < nheight; j++)
        for (int i = 0; i < nwidth; i++){
            w->setQuantity(i,j,0) =  0.0;
            w->setQuantity(i,j,ndepth) = 0.0;
        }
    
    
}


double FluidGrid::maxDivergence2d() const{
    
    double maxdiv = 0;
    double div;
    
    for (int j=0; j<nheight; j++)
        for (int i=0; i<nwidth; i++){
            
            div = 0;
            div =  fabs(u->getQuantity(i+1,j) - u->getQuantity(i,j) + v->getQuantity(i,j+1) - v->getQuantity(i,j))/cellwidth;
            
            maxdiv = fmax(div,maxdiv);
            
        }

    return maxdiv;
    
}


double FluidGrid::maxDivergence3d() const{
    
    double maxdiv = 0;
    double div;
    
    for (int k=0; k<ndepth; k++)
        for (int j=0; j<nheight; j++)
            for (int i=0; i<nwidth; i++){
                
                div = 0;
                div =  fabs(  u->getQuantity(i+1,j,k) - u->getQuantity(i,j,k)
                            + v->getQuantity(i,j+1,k) - v->getQuantity(i,j,k)
                            + w->getQuantity(i,j,k+1) - w->getQuantity(i,j,k))/cellwidth;
                
                maxdiv = fmax(div,maxdiv);
                
            }
    
    return maxdiv;
    
}


//Add density and velocity inside a rectangle
void FluidGrid::addDensity(double x, double y, double w, double h, double dval, double uval, double vval) {
    density->addEmitter(x, y, x + w, y + h, dval);
    u->addEmitter(x, y, x + w, y + h, uval);
    v->addEmitter(x, y, x + w, y + h, vval);
    
}


//Add density and velocity inside a cuboid
void FluidGrid::addDensity(double x, double y, double z, double wh, double h, double l, double dval, double uval, double vval, double wval) {
    density->addEmitter(x, y, z, x + wh, y + h, z + l, dval);
    u->addEmitter(x, y, z, x + wh, y + h, z + l, uval);
    v->addEmitter(x, y, z, x + wh, y + h, z + l, vval);
    w->addEmitter(x, y, z, x + wh, y + h, z + l, wval);
    
}


void FluidGrid::setDeltaT2d(){
   
    double temp = std::max(u->max(),v->max());
    
    if (temp==0)
        deltaT = framedeltaT;
    else
        deltaT = CFLnumber*cellwidth / temp;
    
    if (currtime + 1.01*deltaT >= nextframetime){
        
        deltaT = nextframetime - currtime;
        writetocache = true;
        nextframetime += framedeltaT;
        
    }
    
    printf("deltaT=%f, currtime=%f, writetocache=%d, framenumber= %d, A->max()=%f\n" ,deltaT,currtime,writetocache,framenumber);
    
}


void FluidGrid::setDeltaT3d(){
    
    double maxvel = std::max(w->max(),std::max(u->max(),v->max()));
    
    if (maxvel==0)
        deltaT = framedeltaT;
    else
        deltaT = CFLnumber*cellwidth / maxvel;
    
    if (currtime + 1.01*deltaT >= nextframetime){
        
        deltaT = nextframetime - currtime;
        writetocache = true;
        nextframetime += framedeltaT;
        
    }
    
    printf("deltaT=%f, currtime=%f, writetocache=%d, framenumber= %d\n" ,deltaT,currtime,writetocache,framenumber);
    
}


void FluidGrid::Project(){
    
    setupRHS3d();
    
    //CGSolver(*A, *r, *pressure, TOL, MAXITER);
    PCGSolver(*A, *Ei, *r, *pressure, TOL, MAXITER);
    //GaussSeidelSolver(*A, *r, *pressure, TOL, MAXITER);
    
    updateVelocities3d();
    
}






void FluidGrid::Update(){

    float temp;
    
    setDeltaT3d();
    
    Advect3d();
    
    //AddForces();
    
    addDensity(0.45, 0.2, 0.45, 0.1, 0.05, 0.1, 0.75, 0.0, 2.0, 0.0);
    addDensity(0.45, 0.3, 0.45, 0.1, 0.05, 0.1, 1.0, 0.0, 2.0, 0.0);
    
    temp = maxDivergence3d();
    
    printf("max. divergence before project %f\n", temp);
   
    setupA3d();
    MIC0precon(*A, *Ei);
    
    Project();
    
    currtime += deltaT;
    
    temp = maxDivergence3d();
    
    printf("max. divergence after project %f\n", temp);
    
    //3d write to cache
    if (writetocache){
        WriteToCache();
        writetocache = false;
    }
   
}


void FluidGrid::WriteToCache(){
    
    // Define a coordinate with large signed indices.
    openvdb::Coord ijk;
    
    int& i = ijk[0];
    int& j = ijk[1];
    int& k = ijk[2];
    
    // Get an accessor for coordinate-based access to voxels.
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    
    for (k=0; k<ndepth; k++)
        for (j=0; j<nheight; j++)
            for (i=0; i<nwidth; i++)
                // Set the distance for voxel (i,j,k).
                accessor.setValue(ijk, density->getQuantity(i,j,k));
    
    char outfile[128];
    sprintf(outfile, "%sFrame%05d.vdb",filepath, framenumber++);
    
    // Create a VDB file object.
    openvdb::io::File file(outfile);
    
    // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    // Write out the contents of the container.
    file.write(grids);
    file.close();
   
}


void FluidQuantity::Advect(double timestep, const FluidQuantity& u, const FluidQuantity& v){
    
    double x,y;
    double xmid,ymid;
    
    double uvel,vvel;
    
    for (int j=0; j<ysamples; j++)
        for (int i=0; i<xsamples; i++){
            
            x = i+offsetx;
            y = j+offsety;
            
            uvel = u.InterpolateLinear(x,y) / cellwidth;
            vvel = v.InterpolateLinear(x,y) / cellwidth;
            
            
            //Forward Euler
            /*
            x -= timestep*uvel;
            y -= timestep*vvel;
            */
            
            
            
            //Runge-Kutta 2
            xmid = x - 0.5*timestep*uvel;
            ymid = y - 0.5*timestep*vvel;
            
            uvel = u.InterpolateLinear(xmid,ymid) / cellwidth;
            vvel = v.InterpolateLinear(xmid,ymid) / cellwidth;
            
            x -= timestep*uvel;
            y -= timestep*vvel;
            
            
            
            //setBuffer(i,j) = InterpolateLinear(x,y);
            setBuffer(i,j) = InterpolateCM(x,y);
            
            
        }

    
}


void FluidQuantity::Advect(double timestep, const FluidQuantity& u, const FluidQuantity& v, const FluidQuantity& w){
    
    double x,y,z;
    double xmid,ymid,zmid;
    
    double uvel,vvel,wvel;
    
    for (int k=0; k<zsamples; k++)
        for (int j=0; j<ysamples; j++)
            for (int i=0; i<xsamples; i++){
                
                x = i+offsetx;
                y = j+offsety;
                z = k+offsetz;
                
                uvel = u.InterpolateLinear(x,y,z) / cellwidth;
                vvel = v.InterpolateLinear(x,y,z) / cellwidth;
                wvel = w.InterpolateLinear(x,y,z) / cellwidth;
                
                
                //Forward Euler
                /*
                 x -= timestep*uvel;
                 y -= timestep*vvel;
                 z -= timestep*wvel;
                */
                
                
                
                //Runge-Kutta 2
                
                xmid = x - 0.5*timestep*uvel;
                ymid = y - 0.5*timestep*vvel;
                zmid = z - 0.5*timestep*wvel;
                
                uvel = u.InterpolateLinear(xmid,ymid,zmid) / cellwidth;
                vvel = v.InterpolateLinear(xmid,ymid,zmid) / cellwidth;
                wvel = w.InterpolateLinear(xmid,ymid,zmid) / cellwidth;
                
                x -= timestep*uvel;
                y -= timestep*vvel;
                z -= timestep*wvel;
                
                
                
                //setBuffer(i,j,k) = InterpolateLinear(x,y,z);
                setBuffer(i,j,k) = InterpolateCM(x,y,z);
                
                
            }
    
    
}


//Constructor for FluidQuantity
FluidQuantity::FluidQuantity(double ofx, double ofy, double ofz,  int w, int h, int d, double cw){
        
    offsetx = ofx;
    offsety = ofy;
    offsetz = ofz;
    
    if (offsetx == 0){
        xsamples = w+1;
    }
    else{
        xsamples = w;
    }
    if (offsety == 0){
        ysamples = h+1;
    }
    else{
        ysamples = h;
    }
    if (offsetz == 0){
        zsamples = d+1;
    }
    else{
        zsamples = d;
    }
    
    cellwidth = cw;
    
    quantity = new vectord(xsamples*ysamples*zsamples,0);
    buffer = new vectord(xsamples*ysamples*zsamples,0);
    
}


//Constructor for FluidQuantity
FluidQuantity::FluidQuantity(double ofx, double ofy, int w, int h, double cw){
    
    offsetx = ofx;
    offsety = ofy;
    
    if (offsetx == 0){
        xsamples = w+1;
    }
    else{
        xsamples = w;
    }
    if (offsety == 0){
        ysamples = h+1;
    }
    else{
        ysamples = h;
    }
    
    cellwidth = cw;
    
    quantity = new vectord(xsamples*ysamples,0);
    buffer = new vectord(xsamples*ysamples,0);
    
}


//Destuctor for FluidQuantity
FluidQuantity::~FluidQuantity(){
    
    delete quantity;
    delete buffer;
    
    
}


double FluidQuantity::InterpolateLinear(double x, double y) const{
    
    //transform into the correct dimensionless cellspace for the given offset
    double cx = x-offsetx;
    double cy = y-offsety;
    
    //Clamp x,y to grid
    cx = fmax(fmin(cx,xsamples-1.001),0.0);
    cy = fmax(fmin(cy,ysamples-1.001),0.0);
    
    //Get grid co-ordinates of x,y: i,j such that i,j is closer to the origin than x,y
    int i = int(cx);
    int j = int(cy);
    
    //get relative cell co-ordinates
    cx -= i;
    cy -= j;
    
    double p00 = getQuantity(i,j);
    double p10 = getQuantity(i+1,j);
    double p01 = getQuantity(i,j+1);
    double p11 = getQuantity(i+1,j+1);
    
    //linearly interpolate value of quantity at x,y
    double edge0 = (1-cx)*p00 + cx*p10;
    double edge1 = (1-cx)*p01 + cx*p11;
    
    return (1-cy)*edge0 + cy*edge1;
    
    
}


double FluidQuantity::InterpolateLinear(double x, double y, double z) const{
    
    //transform into the correct dimensionless cellspace for the given offset
    double cx = x-offsetx;
    double cy = y-offsety;
    double cz = z-offsetz;

    //Clamp x,y to grid
    cx = fmax(fmin(cx,xsamples-1.001),0.0);
    cy = fmax(fmin(cy,ysamples-1.001),0.0);
    cz = fmax(fmin(cz,zsamples-1.001),0.0);
    
    //Get grid co-ordinates of x,y,z: i,j,k such that i,j,k is closer to the origin than x,y,z
    int i = int(cx);
    int j = int(cy);
    int k = int(cz);
    
    //get relative cell co-ordinates
    cx -= i;
    cy -= j;
    cz -= k;
   
    double p000 = getQuantity(i,j,k);
    double p100 = getQuantity(i+1,j,k);
    double p010 = getQuantity(i,j+1,k);
    double p110 = getQuantity(i+1,j+1,k);
    double p001 = getQuantity(i,j,k+1);
    double p101 = getQuantity(i+1,j,k+1);
    double p011 = getQuantity(i,j+1,k+1);
    double p111 = getQuantity(i+1,j+1,k+1);
    
    //linearly interpolate value of quantity at x,y,z
    double edge0 = (1-cx)*p000 + cx*p100;
    double edge1 = (1-cx)*p010 + cx*p110;
    double face0 = (1-cy)*edge0 + cy*edge1;
    
    double edge2 = (1-cx)*p001 + cx*p101;
    double edge3 = (1-cx)*p011 + cx*p111;
    double face1 = (1-cy)*edge2 + cy*edge3;
    
    return (1-cz)*face0 + cz*face1;
    
}


double FluidQuantity::CatmullRom(double x, double f0, double f1, double f2, double f3){
    
    double x2 = x*x;
    double x3 = x*x2;
    
    return f0*(    -0.5*x +     x2 - 0.5*x3)
    + f1*(1          - 2.5*x2 + 1.5*x3)
    + f2*(     0.5*x +   2*x2 - 1.5*x3)
    + f3*(            -0.5*x2 + 0.5*x3);
    
    
    
}


double FluidQuantity::InterpolateCM(double x, double y, double z){
    
    
    //transform into the correct dimensionless cellspace for the given offset
    double cz = z-offsetz;
    
    //Clamp z to grid
    cz = fmax(fmin(cz,zsamples-1.001),0.0);
    
    //Get grid co-ordinates of z: k such that k is closer to the origin than z
    int k = int(cz);
    
    //get relative cell z co-ordinates
    cz -= k;
    
    //Need to interpolate on four 2d grids nearest x,y,z
    //making sure to clamp so they lie on the grid in terms of depth
    int mink = std::max(k-1,0);
    int maxk = std::min(k+2,zsamples-1);
    
    //Use InterpolateCMSlice function to interpolate on each of four 2d grids
    double grid1 = InterpolateCMSlice(x, y, mink);
    double grid2 = InterpolateCMSlice(x, y, k);
    double grid3 = InterpolateCMSlice(x, y, k+1);
    double grid4 = InterpolateCMSlice(x, y, maxk);
    
    //Finally, interpolate from 4 points above:
    return CatmullRom(cz, grid1, grid2, grid3, grid4);
    
}


double FluidQuantity::InterpolateCMSlice(double x, double y, int k){
    
    //transform into the correct dimensionless cellspace for the given offset
    double cx = x-offsetx;
    double cy = y-offsety;
    
    //Clamp x,y to grid
    cx = fmax(fmin(cx,xsamples-1.001),0.0);
    cy = fmax(fmin(cy,ysamples-1.001),0.0);
    
    //Get grid co-ordinates of x,y: i,j such that i,j is closer to the origin than x,y
    int i = int(cx);
    int j = int(cy);
    
    //get relative cell co-ordinates
    cx -= i;
    cy -= j;
    
    
    /* i,j are the co-ordinates of p11 in the 2d grid. We want to find sixteen points below,
       while carefully handling grid boundaries so we are never outside the grid
     
     p00   p10   p20   p30
     
     
     p01   p11   p21   p31
     
     
     p02   p12   p22   p32
     
     
     p03   p13   p23   p33
     
     */

    //now clamp the i,j of the twelve points around the central square so they lie on the grid
    int minj = std::max(j-1,0);
    int maxj = std::min(j+2,ysamples-1);
    int mini = std::max(i-1,0);
    int maxi = std::min(i+2,xsamples-1);


    //Use CatmullRom interpolation to get four more points, one for each row
    double row1 = CatmullRom(cx, getQuantity(mini,minj,k), getQuantity(i,minj,k), getQuantity(i+1,minj,k), getQuantity(maxi,minj,k));
    double row2 = CatmullRom(cx, getQuantity(mini,j,k), getQuantity(i,j,k), getQuantity(i+1,j,k), getQuantity(maxi,j,k));
    double row3 = CatmullRom(cx, getQuantity(mini,j+1,k), getQuantity(i,j+1,k), getQuantity(i+1,j+1,k), getQuantity(maxi,j+1,k));
    double row4 = CatmullRom(cx, getQuantity(mini,maxj,k), getQuantity(i,maxj,k), getQuantity(i+1,maxj,k), getQuantity(maxi,maxj,k));
    

    //Finally, interpolate from 4 points above:
    return CatmullRom(cy, row1, row2, row3, row4);
    
}


//Set FluidQuantity inside rectangle
void FluidQuantity::addEmitter(double x0, double y0, double x1, double y1, double v) {
    int ix0 = (int)(x0/cellwidth - offsetx);
    int iy0 = (int)(y0/cellwidth - offsety);
    int ix1 = (int)(x1/cellwidth - offsetx);
    int iy1 = (int)(y1/cellwidth - offsety);
    
    for (int y = std::max(iy0, 0); y < std::min(iy1, ysamples); y++)
        for (int x = std::max(ix0, 0); x < std::min(ix1, xsamples); x++)
            if (fabs(getQuantity(x,y)) < fabs(v))
                setQuantity(x,y) = v;
};


//Set FluidQuantity inside cuboid
void FluidQuantity::addEmitter(double x0, double y0, double z0, double x1, double y1, double z1, double v) {
    int ix0 = (int)(x0/cellwidth - offsetx);
    int iy0 = (int)(y0/cellwidth - offsety);
    int iz0 = (int)(z0/cellwidth - offsetz);
    int ix1 = (int)(x1/cellwidth - offsetx);
    int iy1 = (int)(y1/cellwidth - offsety);
    int iz1 = (int)(z1/cellwidth - offsetz);
    
    for (int z = std::max(iz0, 0); z < std::min(iz1, zsamples); z++)
        for (int y = std::max(iy0, 0); y < std::min(iy1, ysamples); y++)
            for (int x = std::max(ix0, 0); x < std::min(ix1, xsamples); x++)
                if (fabs(getQuantity(x,y,z)) < fabs(v))
                    setQuantity(x,y,z) = v;
};






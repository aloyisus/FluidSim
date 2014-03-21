 //
//  grid.cpp
//  Fluid
//
//  Created by John Kelly on 26/01/2014.
//  Copyright (c) 2014 John Kelly. All rights reserved.
//

#include "grid.h"


#define round5(n) floor(n * 100000 + 0.5)/100000


#define filepath "/Users/JohnnyK/Downloads/Fluids44/"



//Constructor for SolidCuboid class
//Note that position and dimensions are given in terms of cells, not space
SolidCuboid::SolidCuboid(double centrex, double centrey, double centrez, double width, double height, double depth, double alph, double omg)  : centrepos(3,0) {
    
    centrepos.at(0) = centrex;
    centrepos.at(1) = centrey;
    centrepos.at(2) = centrez;

    w = width;
    h = height;
    d = depth;
    
    alpha = 2*pi*alph;
    omega = 2*pi*omg;
    
}

//Function to determine if point i,j,k is inside this solid cuboid
bool SolidCuboid::PointIsInside(double i, double j, double k, double time){
    
    vectord cellpos(3,0);
    cellpos.at(0) = i;
    cellpos.at(1) = j;
    cellpos.at(2) = k;
    
    bool temp;
    matrixd invtransform(3,3,0);
    
    //rotation of -theta about z-axis
    invtransform.at(0,0) = cos(alpha + omega*time);
    invtransform.at(1,0) = -sin(alpha + omega*time);
    invtransform.at(0,1) = sin(alpha + omega*time);
    invtransform.at(1,1) = cos(alpha + omega*time);
    invtransform.at(2,2) = 1.0;
    
    cellpos = invtransform*(cellpos - centrepos);
    
    temp = (cellpos.at(0) >= - w/2) && (cellpos.at(0) <= + w/2) &&
           (cellpos.at(1) >= - h/2) && (cellpos.at(1) <= + h/2) &&
           (cellpos.at(2) >= - d/2) && (cellpos.at(2) <= + d/2);
    
    if (temp)
        printf("");
    
    return temp;
    
}

vectord SolidCuboid::getVelocity(double x, double y, double z){
    
    vectord temp = vectord(3,0);
    
    //vector from axis of rotation (z-axis) to x,y,z
    double rx = x - centrepos.at(0);
    double ry = y - centrepos.at(1);
    
    double distance = sqrt(rx*rx + ry*ry);
    
    //get the unit normal to vector (rx,ry)
    double normx = -ry / distance;
    double normy = rx / distance;
    
    //get magnitude of velocity at x,y, scale the norm vector by this
    normx *= distance*omega;
    normy *= distance*omega;
    
    temp.at(0) = normx;
    temp.at(1) = normy;
    temp.at(2) = 0;
    
    return temp;
    
}


vectord SolidCuboid::getNormal(double x0, double y0, double z0, double time){
    
    vectord temp = vectord(3,0);
    matrixd transform(3,3,0);
    vectord r = vectord(3,0);
    
    r.at(0) = x0;
    r.at(1) = y0;
    r.at(2) = z0;
    
    //rotation of -theta about z-axis
    transform.at(0,0) = cos(alpha + omega*time);
    transform.at(1,0) = -sin(alpha + omega*time);
    transform.at(0,1) = sin(alpha + omega*time);
    transform.at(1,1) = cos(alpha + omega*time);
    transform.at(2,2) = 1.0;
    
    r = transform*(r - centrepos);

    double x = r.at(0);
    double y = r.at(1);
    double z = r.at(2);
    
    
    if ((y >= fabs(z)*h/d) && (y >= fabs(x)*h/w)){
        
        temp.at(0) = 0;
        temp.at(1) = 1;
        temp.at(2) = 0;
    
    }
    else if ((y <= fabs(z)*h/d) && (y <= fabs(x)*h/w)){
        
        temp.at(0) = 0;
        temp.at(1) = -1;
        temp.at(2) = 0;
        
    }
    else if ((x >= fabs(z)*w/d) && (x >= fabs(y)*w/h)){
        
        temp.at(0) = 1;
        temp.at(1) = 0;
        temp.at(2) = 0;
        
    }
    else if ((x <= fabs(z)*w/d) && (x <= fabs(y)*w/h)){
        
        temp.at(0) = -1;
        temp.at(1) = 0;
        temp.at(2) = 0;
        
    }
    else if ((z >= fabs(x)*d/w) && (z >= fabs(y)*d/h)){
        
        temp.at(0) = 0;
        temp.at(1) = 0;
        temp.at(2) = 1;
        
    }
    else if ((z <= fabs(x)*d/w) && (z <= fabs(y)*d/h)){
        
        temp.at(0) = 0;
        temp.at(1) = 0;
        temp.at(2) = -1;
        
    }

    
    //rotation of +theta about z-axis
    transform.at(0,0) = cos(alpha + omega*time);
    transform.at(1,0) = sin(alpha + omega*time);
    transform.at(0,1) = -sin(alpha + omega*time);
    transform.at(1,1) = cos(alpha + omega*time);
    transform.at(2,2) = 1.0;
    
    return transform*temp;
    
    
}


//Constructor for FluidGrid 3d
FluidGrid::FluidGrid(int wh, int ht, int dh, double tstep, double rh){
    
    cellwidth = 1.0/std::min(std::min(wh, ht),dh);

    nwidth = wh;
    nheight = ht;
    ndepth = dh;
    
    rho = rh;
    
    CFLnumber = 5;
    framenumber = 0;
    currtime = 0;
    framedeltaT = tstep;
    nextframetime = tstep;
    bool writetocache = false;
    
    pressure = new vectord(wh*ht*dh);
    r = new vectord(wh*ht*dh);
   
    u = new FluidQuantity(0,0.5,0.5,wh,ht,dh,cellwidth);
    v = new FluidQuantity(0.5,0,0.5,wh,ht,dh,cellwidth);
    w = new FluidQuantity(0.5,0.5,0,wh,ht,dh,cellwidth);
    density = new FluidQuantity(0.5,0.5,0.5,wh,ht,dh,cellwidth);
 
    A = new symmbandd(wh*ht*dh,7,1,wh,wh*ht);
    Ei = new diagonald(wh*ht*dh);
    
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();
    // Create an empty floating-point grid.
    grid = openvdb::FloatGrid::create();
    
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



//This builds the system of equations Ap = b
void FluidGrid::BuildLinearSystem(){
    
    double Ascale =  deltaT/(rho*cellwidth*cellwidth);
    double rscale = 1.0/cellwidth;
    
    int jnwidth,knwidthnheight;
    
    A->fill(0);
    r->fill(0);
    
    //Setup up the matrix of pressure coefficients, weighted by volumes
    for(int k=0; k<ndepth; k++)
        for(int j=0; j<nheight; j++)
            for(int i=0; i<nwidth; i++){
                
                jnwidth = j*nwidth;
                knwidthnheight = k*nwidth*nheight;

                A->at(i+jnwidth+knwidthnheight,0) += Ascale*u->getVolume(i,j,k);

                A->at(i+jnwidth+knwidthnheight,0) += Ascale*v->getVolume(i,j,k);

                A->at(i+jnwidth+knwidthnheight,0) += Ascale*w->getVolume(i,j,k);
            
                A->at(i+jnwidth+knwidthnheight,0) += Ascale*u->getVolume(i+1,j,k);
                A->at(i+jnwidth+knwidthnheight,1) = -Ascale*u->getVolume(i+1,j,k);

                A->at(i+jnwidth+knwidthnheight,0) += Ascale*v->getVolume(i,j+1,k);
                A->at(i+jnwidth+knwidthnheight,2) = -Ascale*v->getVolume(i,j+1,k);
            
                A->at(i+jnwidth+knwidthnheight,0) += Ascale*w->getVolume(i,j,k+1);
                A->at(i+jnwidth+knwidthnheight,3) = -Ascale*w->getVolume(i,j,k+1);

                
                
                setr(i,j,k) = -rscale*(u->getVolume(i+1,j,k)*u->getQuantity(i+1,j,k) - u->getVolume(i,j,k)*u->getQuantity(i,j,k) +
                                      v->getVolume(i,j+1,k)*v->getQuantity(i,j+1,k) - v->getVolume(i,j,k)*v->getQuantity(i,j,k) +
                                      w->getVolume(i,j,k+1)*w->getQuantity(i,j,k+1) - w->getVolume(i,j,k)*w->getQuantity(i,j,k));

                setr(i,j,k) += rscale*((u->getVolume(i+1,j,k)-density->getVolume(i,j,k))*
                                      cellwidth*solid->getVelocity(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz()).at(0) -
                                      (u->getVolume(i,j,k)-density->getVolume(i,j,k))*
                                      cellwidth*solid->getVelocity(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz()).at(0) +
                                      (v->getVolume(i,j+1,k)-density->getVolume(i,j,k))*
                                      cellwidth*solid->getVelocity(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz()).at(1) -
                                      (v->getVolume(i,j,k)-density->getVolume(i,j,k))*
                                      cellwidth*solid->getVelocity(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz()).at(1) +
                                       (w->getVolume(i,j,k+1)-density->getVolume(i,j,k))*
                                      cellwidth*solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz()).at(2) -
                                      (w->getVolume(i,j,k)-density->getVolume(i,j,k))*
                                      cellwidth*solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz()).at(2));
            
            };
    
    
    //DEBUG
    /*
    double sumr = 0;
    int count = 0;
    
    //go through all connected cells (they are fluid cells) and sum r
    for(int k=0; k<ndepth; k++)
        for(int j=0; j<nheight; j++)
            for(int i=0; i<nwidth; i++)
                if (density->getVolume(i,j,k) > 0){
                    sumr += getr(i,j,k);
                    count++;
                }
    
    //take the mean of r for the connected cells
    sumr /= count;
    
    //subtract the mean of r from the connected cells
    for(int k=0; k<ndepth; k++)
        for(int j=0; j<nheight; j++)
            for(int i=0; i<nwidth; i++)
                if (density->getVolume(i,j,k) > 0){
                    setr(i,j,k) -= sumr;
                }
     */
     
    
}


void FluidGrid::Advect3d(){
    
    density->Advect(deltaT,*u,*v,*w, solid, currtime);
    u->Advect(deltaT,*u,*v,*w, solid, currtime);
    v->Advect(deltaT,*u,*v,*w, solid, currtime);
    w->Advect(deltaT,*u,*v,*w, solid, currtime);
    
    u->swap();
    v->swap();
    w->swap();
    density->swap();
    
}


//Fill up the volume samples of the grid
void FluidGrid::CalculateVolumes(){
    
    double* supervol = new double[nwidth*nheight*ndepth*8];
    
    //Perform 2x2x2 supersampling over entire grid, so each ijk cell will be sampled 8 times
    for (int k=0; k<2*ndepth; k++)
        for (int j=0; j<2*nheight; j++)
            for (int i=0; i<2*nwidth; i++){
                
                supervol[i+j*2*nwidth+k*4*nwidth*nheight] = solid->PointIsInside(0.25+i*0.5, 0.25+j*0.5, 0.25+k*0.5, getSimtime());
               
            }
   
    double posx,posy,posz;
    
    //Note that u samples at (0,j,k) and (nwidth,j,k) are unchanged as zero
    for (int k=0; k<u->getzSamples(); k++)
        for (int j=0; j<u->getySamples(); j++)
            for (int i=1; i<u->getxSamples()-1; i++){
                
                //get position of sample in 'cell space'
                posx = i+u->getOffsetx();
                posy = j+u->getOffsety();
                posz = k+u->getOffsetz();
                
                int ix0 = floor((posx-0.25)*2);
                int ix1 = ix0 + 1;
                
                int iy0 = floor((posy-0.25)*2);
                int iy1 = iy0 + 1;
                
                int iz0 = floor((posz-0.25)*2);
                int iz1 = iz0 + 1;
                
                u->setVolume(i,j,k) = 8.0;
                

                u->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz0*4*nwidth*nheight];
                u->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz0*4*nwidth*nheight];
                u->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz1*4*nwidth*nheight];
                u->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz1*4*nwidth*nheight];
                u->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz0*4*nwidth*nheight];
                u->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz0*4*nwidth*nheight];
                u->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz1*4*nwidth*nheight];
                u->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz1*4*nwidth*nheight];

                
                u->setVolume(i,j,k) = u->setVolume(i,j,k)/8.0;
                
                
            }
    //Note that v samples at (i,0,k) and (i,nheight,k) are unchanged as zero
    for (int k=0; k<v->getzSamples(); k++)
        for (int j=1; j<v->getySamples()-1; j++)
            for (int i=0; i<v->getxSamples(); i++){
                
                //get position of sample in 'cell space'
                posx = i+v->getOffsetx();
                posy = j+v->getOffsety();
                posz = k+v->getOffsetz();
                
                int ix0 = floor((posx-0.25)*2);
                int ix1 = ix0 + 1;
                
                int iy0 = floor((posy-0.25)*2);
                int iy1 = iy0 + 1;
                
                int iz0 = floor((posz-0.25)*2);
                int iz1 = iz0 + 1;
                
                v->setVolume(i,j,k) = 8.0;
                

                v->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz0*4*nwidth*nheight];
                v->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz0*4*nwidth*nheight];
                v->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz1*4*nwidth*nheight];
                v->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz1*4*nwidth*nheight];
                v->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz0*4*nwidth*nheight];
                v->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz0*4*nwidth*nheight];
                v->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz1*4*nwidth*nheight];
                v->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz1*4*nwidth*nheight];

                
                v->setVolume(i,j,k) = v->setVolume(i,j,k)/8.0;

                
            }
                
    //Note that w samples at (i,j,0) and (i,j,ndepth) are unchanged as zero
    for (int k=1; k<w->getzSamples()-1; k++)
        for (int j=0; j<w->getySamples(); j++)
            for (int i=0; i<w->getxSamples(); i++){
                
                //get position of sample in 'cell space'
                posx = i+w->getOffsetx();
                posy = j+w->getOffsety();
                posz = k+w->getOffsetz();
            
                int ix0 = floor((posx-0.25)*2);
                int ix1 = ix0 + 1;
                
                int iy0 = floor((posy-0.25)*2);
                int iy1 = iy0 + 1;
                
                int iz0 = floor((posz-0.25)*2);
                int iz1 = iz0 + 1;
                
                w->setVolume(i,j,k) = 8.0;
                

                w->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz0*4*nwidth*nheight];
                w->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz0*4*nwidth*nheight];
                w->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz0*4*nwidth*nheight];
                w->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz0*4*nwidth*nheight];
                w->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz1*4*nwidth*nheight];
                w->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz1*4*nwidth*nheight];
                w->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz1*4*nwidth*nheight];
                w->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz1*4*nwidth*nheight];

                
                w->setVolume(i,j,k) = w->setVolume(i,j,k)/8.0;

    
            }
    
    
    for (int k=0; k<density->getzSamples(); k++)
        for (int j=0; j<density->getySamples(); j++)
            for (int i=0; i<density->getxSamples(); i++){
                
                //get position of sample in 'cell space'
                posx = i+density->getOffsetx();
                posy = j+density->getOffsety();
                posz = k+density->getOffsetz();
                
                int ix0 = floor((posx-0.25)*2);
                int ix1 = ix0 + 1;
                
                int iy0 = floor((posy-0.25)*2);
                int iy1 = iy0 + 1;
                
                int iz0 = floor((posz-0.25)*2);
                int iz1 = iz0 + 1;


                density->setVolume(i,j,k) = 8.0;
                
                density->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz0*4*nwidth*nheight];
                density->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz0*4*nwidth*nheight];
                density->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz0*4*nwidth*nheight];
                density->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz0*4*nwidth*nheight];
                density->setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz1*4*nwidth*nheight];
                density->setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz1*4*nwidth*nheight];
                density->setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz1*4*nwidth*nheight];
                density->setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz1*4*nwidth*nheight];
                
                density->setVolume(i,j,k) = density->setVolume(i,j,k)/8.0;
                
                
            }
   
    delete supervol;
    
}



void FluidGrid::updateVelocities3d() {
    double scale = deltaT/(rho*cellwidth);
    
    //We don't touch the samples that are aligned with the edges of the grid here
    
    for (int k = 0; k < u->getzSamples(); k++)
        for (int j = 0; j < u->getySamples(); j++)
            for (int i = 1; i < u->getxSamples()-1; i++) {
                
                if (u->getVolume(i,j,k) > 0)
                    u->setQuantity(i,j,k) -= scale*(getPressure(i,j,k) - getPressure(i-1,j,k));
                else
                    u->setQuantity(i,j,k) = 0;
               
            }
    
    for (int k = 0; k < v->getzSamples(); k++)
        for (int j = 1; j < v->getySamples()-1; j++)
            for (int i = 0; i < v->getxSamples(); i++) {
                
                if (v->getVolume(i,j,k) > 0)
                    v->setQuantity(i,j,k) -= scale*(getPressure(i,j,k) - getPressure(i,j-1,k));
                else
                    v->setQuantity(i,j,k) = 0;
                
            }
    
    for (int k = 1; k < w->getzSamples()-1; k++)
        for (int j = 0; j < w->getySamples(); j++)
            for (int i = 0; i < w->getxSamples(); i++) {
                
                if (w->getVolume(i,j,k) > 0)
                    w->setQuantity(i,j,k) -= scale*(getPressure(i,j,k) - getPressure(i,j,k-1));
                else
                    w->setQuantity(i,j,k) = 0;
                
            }
    
    
    
    
    
    double maxvel = std::max(w->max(),std::max(u->max(),v->max()));
   
    printf("max velocity after pressure update = %f\n",maxvel);
    
    density->ExtrapolateVelocity();
    u->ExtrapolateVelocity();
    v->ExtrapolateVelocity();
    w->ExtrapolateVelocity();
    
    maxvel = std::max(w->max(),std::max(u->max(),v->max()));
    
    printf("max velocity after velocity extrapolation = %f\n",maxvel);
    
    SetSolidBoundaries();
    
    maxvel = std::max(w->max(),std::max(u->max(),v->max()));
    
    printf("max velocity after SetSolidBoundaries = %f\n",maxvel);
   
    
}


vectord FluidGrid::getVelocity(double x, double y, double z){
    
    
    vectord temp = vectord(3,0);
    
    temp.at(0) = u->InterpolateLinear(x,y,z);
    temp.at(1) = v->InterpolateLinear(x,y,z);
    temp.at(2) = w->InterpolateLinear(x,y,z);
    
    return temp;
    
    
}


void FluidGrid::SetSolidBoundaries(){

    u->CopytoBuffer();
    v->CopytoBuffer();
    w->CopytoBuffer();
    
    vectord vfluid(3,0);
    vectord vsolid(3,0);
    vectord normal(3,0);
    
    //For cells which have zero volume, set the velocity samples bordering to have zero relative velocity in the
    //direction of the normal to the solid boundary
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++)
            for (int i = 0; i < nwidth; i++) {
                
                if (density->getVolume(i,j,k) == 0)
                {
                    

                    vsolid = cellwidth*solid->getVelocity(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    vfluid = getVelocity(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    normal = solid->getNormal(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz(),currtime);
                    
                    u->setBuffer(i,j,k) =  (vfluid - ((vfluid - vsolid)*normal)*normal).at(0);
                    
                    vsolid = cellwidth*solid->getVelocity(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    vfluid = getVelocity(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    normal = solid->getNormal(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz(),currtime);
                    
                    u->setBuffer(i+1,j,k) = (vfluid - ((vfluid - vsolid)*normal)*normal).at(0);
                    
                    vsolid = cellwidth*solid->getVelocity(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz());
                    vfluid = getVelocity(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz());
                    normal = solid->getNormal(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz(),currtime);
                    
                    v->setBuffer(i,j,k) = (vfluid - ((vfluid - vsolid)*normal)*normal).at(1);
                    
                    vsolid = cellwidth*solid->getVelocity(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz());
                    vfluid = getVelocity(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz());
                    normal = solid->getNormal(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz(),currtime);
                    
                    v->setBuffer(i,j+1,k) = (vfluid - ((vfluid - vsolid)*normal)*normal).at(1);
                    
                    vsolid = cellwidth*solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz());
                    vfluid = getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz());
                    normal = solid->getNormal(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz(),currtime);
                    
                    w->setBuffer(i,j,k) = (vfluid - ((vfluid - vsolid)*normal)*normal).at(2);
                    
                    vsolid = cellwidth*solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz());
                    vfluid = getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz());
                    normal = solid->getNormal(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz(),currtime);
                    
                    w->setBuffer(i,j,k+1) = (vfluid - ((vfluid - vsolid)*normal)*normal).at(2);
                   
                    
                }
                
            }
    
    
    u->swap();
    v->swap();
    w->swap();
    
    
    //Boundary conditions for the walls of the box
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


//Add density and velocity inside a cuboid
void FluidGrid::addDensity(double x, double y, double z, double wh, double h, double l, double dval, double uval, double vval, double wval) {
    
    double x0 = int(x*nwidth);
    double x1 = int((x+wh)*nwidth);
    double y0 = int(y*nheight);
    double y1 = int((y+h)*nheight);
    double z0 = int(z*ndepth);
    double z1 = int((z+l)*ndepth);
    
    density->addEmitter(x0, y0, z0, x1, y1, z1, dval);
    u->addEmitter(x0, y0, z0, x1, y1, z1, uval);
    v->addEmitter(x0, y0, z0, x1, y1, z1, vval);
    w->addEmitter(x0, y0, z0, x1, y1, z1, wval);
    
}


void FluidGrid::setDeltaT3d(){
    
    double maxvel = std::max(w->max(),std::max(u->max(),v->max()));
    
    
    printf("max velocity before project = %f\n",maxvel);
    
    
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
    
    
    //CGSolver(*A, *r, *pressure, TOL, MAXITER);
    PCGSolver(*A, *Ei, *r, *pressure, TOL, MAXITER);
    //GaussSeidelSolver(*A, *r, *pressure, TOL, MAXITER);
    
    updateVelocities3d();
    
}


void FluidGrid::Update(){

    float temp;
    
    setDeltaT3d();

    //u->ExtrapolateVelocity();
    //v->ExtrapolateVelocity();
    //w->ExtrapolateVelocity();
    
    Advect3d();
    
    //AddForces();
    
    if (writetocache)
        addDensity(0.45, 0.2, 0.45, 0.1, 0.05, 0.2, 1.0, 0.0, 2.5, 0.0);
    
    //temp = maxDivergence3d();
    
    //printf("max. divergence before project %f\n", temp);
   
    CalculateVolumes();
    
    BuildLinearSystem();
    MIC0precon(*A, *Ei);
   
    Project();
    
    currtime += deltaT;
    
    //temp = maxDivergence3d();
    
    //printf("max. divergence after project %f\n", temp);
    
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
            for (i=0; i<nwidth; i++){
                
                //accessor.setValue(ijk, density->getQuantity(i,j,k));
                accessor.setValue(ijk, density->getQuantity(i,j,k) + 1.0 - density->getVolume(i,j,k));
    
            }
                
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


//DEBUG FOR TESTING solid
void FluidGrid::SolidDensityUpdate(SolidCuboid& solid){
    
    double dens;
    
    // Define a coordinate with large signed indices.
    openvdb::Coord ijk;
    
    int& i = ijk[0];
    int& j = ijk[1];
    int& k = ijk[2];
    
    // Get an accessor for coordinate-based access to voxels.
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    
    for (k=0; k<ndepth; k++)
        for (j=0; j<nheight; j++)
            for (i=0; i<nwidth; i++){

                //dens = solid.PointIsInside(i+0.5,j+0.5,k+0.5,currtime);
                dens = 1.0 - density->getVolume(i,j,k);
            
                accessor.setValue(ijk, dens);
            
            }
                
    char outfile[128];
    sprintf(outfile, "%sCubeTest%05d.vdb",filepath, framenumber++);
    
    // Create a VDB file object.
    openvdb::io::File file(outfile);
    
    // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    // Write out the contents of the container.
    file.write(grids);
    file.close();
    
    currtime += framedeltaT;
    
    
    
}



//if x,y is inside a solid, backproject to nearest surface point and return that
void SolidCuboid::ClosestSurface(double &x, double &y, double &z, double time) {

    matrixd transform(3,3,0);
    vectord r = vectord(3,0);
    
    r.at(0) = x;
    r.at(1) = y;
    r.at(2) = z;
    
    //rotation of -theta about z-axis
    transform.at(0,0) = cos(alpha + omega*time);
    transform.at(1,0) = -sin(alpha + omega*time);
    transform.at(0,1) = sin(alpha + omega*time);
    transform.at(1,1) = cos(alpha + omega*time);
    transform.at(2,2) = 1.0;
    
    r = transform*(r - centrepos);

    
    double dx = w - 2*fabs(x);
    double dy = h - 2*fabs(y);
    double dz = d - 2*fabs(z);
    
    if ((dx <= dy) && (dx <= dz))
        x = (fabs(x)/x)*w/2;
    else if ((dy <= dx) && (dy <= dz))
        y = (fabs(y)/y)*h/2;
    else
        z = (fabs(z)/z)*d/2;
    
   
    //rotation of +theta about z-axis
    transform.at(0,0) = cos(alpha + omega*time);
    transform.at(1,0) = sin(alpha + omega*time);
    transform.at(0,1) = -sin(alpha + omega*time);
    transform.at(1,1) = cos(alpha + omega*time);
    transform.at(2,2) = 1.0;

    r = centrepos + transform*r;
    
    x = r.at(0);
    y = r.at(1);
    z = r.at(2);
    
    
    
}


void FluidQuantity::Advect(double timestep, const FluidQuantity& u, const FluidQuantity& v, const FluidQuantity& w, SolidCuboid* gridsolid, double time ){
    
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
                
                
                //ARE WE inside the SOLID?
                if (gridsolid->PointIsInside(x, y, z, time))
                    gridsolid->ClosestSurface(x, y, z, time);
                
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
    volume = new vectord(xsamples*ysamples*zsamples,0);
    
}


//Destuctor for FluidQuantity
FluidQuantity::~FluidQuantity(){
    
    delete quantity;
    delete buffer;
    
    
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


//Set FluidQuantity inside a given volume of space
void FluidQuantity::addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double v) {

    
    for (int z = std::max(z0, 0); z < std::min(z1, zsamples); z++)
        for (int y = std::max(y0, 0); y < std::min(y1, ysamples); y++)
            for (int x = std::max(x0, 0); x < std::min(x1, xsamples); x++)
                if (fabs(getQuantity(x,y,z)) < fabs(v))
                    setQuantity(x,y,z) = v;
}




void FluidQuantity::ExtrapolateVelocity(){
    
    bool* valid = new bool[xsamples*ysamples*zsamples];
    bool* old_valid = new bool[xsamples*ysamples*zsamples];
    
    //Initialize the list of valid cells IF WEIGHT IS NON-ZERO THEN CELL IS VALID
    for(int k = 0; k < zsamples; k++)
        for(int j = 0; j < ysamples; j++)
            for(int i = 0; i < xsamples; i++)
                valid[i + j*xsamples + k*xsamples*ysamples] = getVolume(i,j,k) > 0;
    
    for(int layers = 0; layers < 5; ++layers) {
        
        for(int i =0; i<xsamples*ysamples*zsamples; i++)
            old_valid[i] = valid[i];
        
        
        //DON'T ITERATE OVER THE EDGES OF THE GRID
        for(int k = 1; k < zsamples-1; k++)
            for(int j = 1; j < ysamples-1; j++)
                for(int i = 1; i < xsamples-1; i++){
            
                    double sum = 0;
                    int count = 0;
                    
                    //IF CELL IS *NOT* VALID
                    if(!old_valid[i + j*xsamples + k*xsamples*ysamples]) {
                        
                        if(old_valid[i+1 + j*xsamples + k*xsamples*ysamples]) {
                            sum += getQuantity(i+1,j,k);
                            ++count;
                        }
                        if(old_valid[i-1 + j*xsamples + k*xsamples*ysamples]) {
                            sum += getQuantity(i-1,j,k);
                            ++count;
                        }
                        if(old_valid[i + (j+1)*xsamples + k*xsamples*ysamples]) {
                            sum += getQuantity(i,j+1,k);
                            ++count;
                        }
                        if(old_valid[i + (j-1)*xsamples + k*xsamples*ysamples]) {
                            sum += getQuantity(i,j-1,k);
                            ++count;
                        }
                        if(old_valid[i + j*xsamples + (k+1)*xsamples*ysamples]) {
                            sum += getQuantity(i,j,k+1);
                            ++count;
                        }
                        if(old_valid[i + j*xsamples + (k-1)*xsamples*ysamples]) {
                            sum += getQuantity(i,j,k-1);
                            ++count;
                        }
                        
                        //If any of neighbour cells were valid,
                        //assign the cell their average value and tag it as valid
                        if(count > 0) {
                            setQuantity(i,j,k) = sum /(double)count;
                            valid[i + j*xsamples + k*xsamples*ysamples] = true;
                        }
                        
                    }
                    
        }
        
    }
    
    delete valid;
    delete old_valid;
   
};





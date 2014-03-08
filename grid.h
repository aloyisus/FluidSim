//
//  grid.h
//  Fluid
//
//  Created by John Kelly on 26/01/2014.
//  Copyright (c) 2014 John Kelly. All rights reserved.
//

#ifndef __Fluid__grid__
#define __Fluid__grid__

#include <iostream>
#include "matrix.h"
#include "openvdb/openvdb.h"


class FluidQuantity {


    vectord* quantity;
    vectord* buffer;
    
    int xsamples,ysamples,zsamples;
    float offsetx, offsety, offsetz;
    float cellwidth;
    
    double CatmullRom(double x, double f0, double f1, double f2, double f3);
    
    double InterpolateLinear(double x, double y, double z) const;
    double InterpolateLinear(double x, double y) const;
    
    double InterpolateCMSlice(double x, double y, int k);
    
    double InterpolateCM(double x, double y, double z);
    
    double InterpolateCM(double x, double y){
        
        return InterpolateCMSlice(x,y,0);
        
    }
    
   
    
public:

    //Constructor for FluidQuantity
    FluidQuantity(double ofx, double ofy, double ofz, int w, int h, int d, double cw);
    FluidQuantity(double ofx, double ofy, int w, int h, double cw);
        
    //Destuctor for FluidQuantity
    ~FluidQuantity();

    void swap(){ std::swap(quantity, buffer); }
    
    void Advect(double timestep, const FluidQuantity& u, const FluidQuantity& v);
    void Advect(double timestep, const FluidQuantity& u, const FluidQuantity& v, const FluidQuantity& w);
    
    double getQuantity(int i, int j) const{ return quantity->at(i + j*xsamples); }
    double &setQuantity(int i, int j){ return quantity->at(i + j*xsamples); };
    double getQuantity(int i, int j, int k) const{ return quantity->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setQuantity(int i, int j, int k){ return quantity->at(i + j*xsamples + k*xsamples*ysamples); };
    
    
    double getBuffer(int i, int j) const{ return buffer->at(i + j*xsamples); }
    double &setBuffer(int i, int j){ return buffer->at(i + j*xsamples); };
    double getBuffer(int i, int j, int k) const{ return buffer->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setBuffer(int i, int j, int k){ return buffer->at(i + j*xsamples + k*xsamples*ysamples); };
    
    double sum() const { return quantity->sum(); }
    double max() const { return quantity->max(); }
   
    void addEmitter(double x0, double y0, double x1, double y1, double v);
    void addEmitter(double x0, double y0, double z0, double x1, double y1, double z1, double v);
   
};




class FluidGrid {
    
    int nwidth,nheight,ndepth;
    
    double cellwidth;
    
    double rho;
    
    int CFLnumber;

    double currtime;
    double nextframetime;
    double framedeltaT;
    double deltaT;
    int framenumber;
    bool writetocache;
    
    //2d grid is written to bitmap
    unsigned char *image;
    
    //3d grid is written to openvdb container
    openvdb::FloatGrid::Ptr grid;

    //Matrices for pressure coefficients A and MIC(0) preconditioner Ei
    symmbandd* A;
    diagonald* Ei;
   
    //vectors for pressure and rhs of matrix equation
    vectord* pressure;
    vectord *r; /* Right hand side of pressure solve */

    //Velocity components
    FluidQuantity* u;
    FluidQuantity* v;
    FluidQuantity* w;
    
    //Quantity defined on grid eg. density
    FluidQuantity* density;
    
    void setupA2d();
    void setupA3d();
    
    void setupRHS2d();
    void setupRHS3d();
    
    void Advect2d();
    void Advect3d();
    
    void updateVelocities2d();
    void updateVelocities3d();
    
    void Project();
    
    void AddForces();
    
    void setDeltaT2d();
    void setDeltaT3d();
    
    void WriteToCache();
    
    double maxDivergence2d() const;
    double maxDivergence3d() const;
    
    void addDensity(double x, double y, double w, double h, double d, double u, double v);
    void addDensity(double x, double y, double z, double wh, double h, double l, double d, double u, double v, double w);
    
    
    //Getters and Setters for pressure
    double getPressure(int i, int j) const { return pressure->at(j*nwidth + i); }
    double& setPressure( int i, int j ){ return pressure->at(j*nwidth + i); }
    double getPressure(int i, int j, int k ) const { return pressure->at(k*nwidth*nheight + j*nwidth + i); }
    double& setPressure( int i, int j, int k ){ return pressure->at(k*nwidth*nheight + j*nwidth + i); }
    
    
    //Getters and Setters for r
    double getr(int i, int j) const { return r->at(j*nwidth + i); }
    double& setr( int i, int j ){ return r->at(j*nwidth + i); }
    double getr(int i, int j, int k ) const { return r->at(k*nwidth*nheight + j*nwidth + i); }
    double& setr( int i, int j, int k ){ return r->at(k*nwidth*nheight + j*nwidth + i); }
    
    int getframenumber() const { return framenumber; }
    double getDeltaT() const { return deltaT; }
    
public:
    
    //Constructor for FluidGrid
    FluidGrid(int w, int h, int d, double tstep, double rh);
    
    //Constructor for FluidGrid
    FluidGrid(int w, int h, double tstep, double rh);
    
    //Destuctor for FluidGrid
    ~FluidGrid();

    double getSimtime() const { return currtime; }
    
    void Update();
    
    
};



#endif /* defined(__Fluid__grid__) */

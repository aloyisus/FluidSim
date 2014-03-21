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

#define pi 3.14159265359

class SolidCuboid {
    
    vectord centrepos;
    double w,h,d;
    double alpha;
    double omega;
    
public:
    
    SolidCuboid(double centre0x, double centrey, double centrez, double width, double height, double depth, double alpha, double omega);
   
    bool PointIsInside(double i, double j, double k, double time);
    
    vectord getVelocity(double x, double y, double z);
   
    vectord getNormal(double x, double y, double z, double time);
    
    void ClosestSurface(double &x, double &y, double &z, double time);
    
};


class FluidQuantity {


    vectord* quantity;
    vectord* buffer;
    vectord* volume;
    
    
    int xsamples,ysamples,zsamples;
    double offsetx, offsety, offsetz;
    double cellwidth;
    
    double CatmullRom(double x, double f0, double f1, double f2, double f3);
    

    


    
public:

    double InterpolateLinear(double x, double y, double z) const;
    
    double InterpolateCMSlice(double x, double y, int k);
    
    double InterpolateCM(double x, double y, double z);
    
    double InterpolateCM(double x, double y){
        
        return InterpolateCMSlice(x,y,0);
        
    }
    
    void CopytoBuffer(){ *buffer = *quantity; }
    
   
    int getxSamples(){ return xsamples; };
    int getySamples(){ return ysamples; };
    int getzSamples(){ return zsamples; };
    
    double getOffsetx(){ return offsetx; };
    double getOffsety(){ return offsety; };
    double getOffsetz(){ return offsetz; };
    
    //Constructor for FluidQuantity
    FluidQuantity(double ofx, double ofy, double ofz, int w, int h, int d, double cw);
    
    //Destuctor for FluidQuantity
    ~FluidQuantity();

    void swap(){ std::swap(quantity, buffer); }
    
    void Advect(double timestep, const FluidQuantity& u, const FluidQuantity& v, const FluidQuantity& w, SolidCuboid* gridsolid, double time);
    
    double getQuantity(int i, int j, int k) const{ return quantity->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setQuantity(int i, int j, int k){ return quantity->at(i + j*xsamples + k*xsamples*ysamples); };
    
    double getBuffer(int i, int j, int k) const{ return buffer->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setBuffer(int i, int j, int k){ return buffer->at(i + j*xsamples + k*xsamples*ysamples); };
    
    double getVolume(int i, int j, int k) const{ return volume->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setVolume(int i, int j, int k){ return volume->at(i + j*xsamples + k*xsamples*ysamples); };
    
    double sum() const { return quantity->sum(); }
    double max() const { return quantity->max(); }
   
    void addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double v);

    void ExtrapolateVelocity();
    
   
};




class FluidGrid {
    
    int nwidth,nheight,ndepth;
    
    double cellwidth;
    
    double rho;
    
    double CFLnumber;

    double currtime;
    double nextframetime;
    double framedeltaT;
    double deltaT;
    int framenumber;
    bool writetocache;
    
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
    
    
    //vectors for volumes of u cells, v cells, w cells and d cells
    vectord* u_volume;
    vectord* v_volume;
    vectord* w_volume;
    vectord* d_volume;
    
    
    //Quantity defined on grid eg. density
    FluidQuantity* density;
    
    void BuildLinearSystem();
    
    void Advect3d();
    
    void updateVelocities3d();
    
    void Project();
    
    void AddForces();
    
    void setDeltaT3d();
    
    void WriteToCache();
    
    double maxDivergence3d() const;
    
    void addDensity(double x, double y, double z, double wh, double h, double l, double d, double u, double v, double w);
    
    
    //Getters and Setters for pressure
    double getPressure(int i, int j, int k ) const { return pressure->at(k*nwidth*nheight + j*nwidth + i); }
    double& setPressure( int i, int j, int k ){ return pressure->at(k*nwidth*nheight + j*nwidth + i); }
    
    
    //Getters and Setters for r
    double getr(int i, int j, int k ) const { return r->at(k*nwidth*nheight + j*nwidth + i); }
    double& setr( int i, int j, int k ){ return r->at(k*nwidth*nheight + j*nwidth + i); }
    
    int getframenumber() const { return framenumber; }
    double getDeltaT() const { return deltaT; }
    

    
public:
    
    //Constructor for FluidGrid
    FluidGrid(int w, int h, int d, double tstep, double rh);
    
    //Destuctor for FluidGrid
    ~FluidGrid();

    double getSimtime() const { return currtime; }
    
    double getWidth() const {return nwidth;};
    double getHeight() const {return nheight;};
    double getDepth() const {return ndepth;};
    
    
    
    void Update();
    
    //DEBUG FOR TESTING solid
    void SolidDensityUpdate(SolidCuboid& solid);
    
    SolidCuboid* solid;
    
    void CalculateVolumes();
    
    vectord getVelocity(double x, double y, double z);
    
    void SetSolidBoundaries();
    
    
};



#endif /* defined(__Fluid__grid__) */

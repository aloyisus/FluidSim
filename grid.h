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


class FluidGrid;


//Abstract solid class
class Solid {
    
protected:
    
    vectord centrepos;
    FluidGrid* parentgrid;
    Solid() : centrepos(3,0) { };
    
public:
    
    virtual bool PointIsInside(double x, double y, double z) = 0;
    
    virtual vectord getVelocity(double x, double y, double z) = 0;
    
    virtual vectord getNormal(double x, double y, double z) = 0;
    
    virtual void ClosestSurface(double &x, double &y, double &z) = 0;
    
    
    
};




class Cuboid : public Solid {
    

    double w,h,d;
    double alpha;
    double omega;
    
public:
    
    Cuboid(FluidGrid* parent, double centrex, double centrey, double centrez, double width, double height, double depth, double alpha, double omega);
    
    bool PointIsInside(double x, double y, double z);
    vectord getVelocity(double x, double y, double z);
    vectord getNormal(double x, double y, double z);
    void ClosestSurface(double &x, double &y, double &z);
    
};




class Sphere : public Solid {

    double radius;
    
public:
    
    Sphere(FluidGrid* parent, double centrex, double centrey, double centrez, double radius);
    
    bool PointIsInside(double x, double y, double z);
    vectord getVelocity(double x, double y, double z);
    vectord getNormal(double x, double y, double z);
    void ClosestSurface(double &x, double &y, double &z);
    
};




class FluidQuantity {


    vectord* quantity;
    vectord* buffer;
    vectord* volume;
    
    FluidGrid* parentgrid;
    
    int xsamples,ysamples,zsamples;
    double offsetx, offsety, offsetz;

    
    double CatmullRom(double x, double f0, double f1, double f2, double f3);
    double InterpolateCMSlice(double x, double y, int k);
    
    double InterpolateCM(double x, double y, double z);
    
    double InterpolateCM(double x, double y){
        
        return InterpolateCMSlice(x,y,0);
        
    }
    
public:
    
    //Constructor for FluidQuantity
    FluidQuantity(FluidGrid* parent, double ofx, double ofy, double ofz);
    
    //Destuctor for FluidQuantity
    ~FluidQuantity();

    
    //Accessors
    int getxSamples(){ return xsamples; };
    int getySamples(){ return ysamples; };
    int getzSamples(){ return zsamples; };
    
    double getOffsetx(){ return offsetx; };
    double getOffsety(){ return offsety; };
    double getOffsetz(){ return offsetz; };
    
    
    double getQuantity(int i, int j, int k) const{ return quantity->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setQuantity(int i, int j, int k){ return quantity->at(i + j*xsamples + k*xsamples*ysamples); };
    
    double getBuffer(int i, int j, int k) const{ return buffer->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setBuffer(int i, int j, int k){ return buffer->at(i + j*xsamples + k*xsamples*ysamples); };
    
    double getVolume(int i, int j, int k) const{ return volume->at(i + j*xsamples + k*xsamples*ysamples); };
    double &setVolume(int i, int j, int k){ return volume->at(i + j*xsamples + k*xsamples*ysamples); };
    
    
    
    double sum() const { return quantity->sum(); }
    double max() const { return quantity->max(); }
   
    void CopyQuantitytoBuffer(){ *buffer = *quantity; }
    
    void swap(){ std::swap(quantity, buffer); }
    
    void Advect();
    
    void addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double v);

    void ExtrapolateVelocity();
    
    double InterpolateLinear(double x, double y, double z) const;
    
    void CalculateVolumes(double* supervol);
    
   
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
    
    Solid* solid;
    
    
    //Accessors
    double getPressure(int i, int j, int k ) const { return pressure->at(k*nwidth*nheight + j*nwidth + i); }
    double& setPressure( int i, int j, int k ){ return pressure->at(k*nwidth*nheight + j*nwidth + i); }
    double getr(int i, int j, int k ) const { return r->at(k*nwidth*nheight + j*nwidth + i); }
    double& setr( int i, int j, int k ){ return r->at(k*nwidth*nheight + j*nwidth + i); }
    int getframenumber() const { return framenumber; }
    

    void BuildLinearSystem();
    void Advect();
    void updateVelocities();
    void Project();
    void AddForces();
    void setDeltaT();
    void WriteToCache();
    double maxDivergence() const;
    void addDensity(double x, double y, double z, double wh, double h, double l, double d, double u, double v, double w);
    
    
    
public:
    
    //Constructor for FluidGrid
    FluidGrid(int w, int h, int d, double tstep, double rh);
    
    //Destuctor for FluidGrid
    ~FluidGrid();

    
    //Accessors
    double getDeltaT() const { return deltaT; }
    double getCellwidth() const { return cellwidth; }
    double getCurrtime() const {return currtime;}
    Solid* getSolid() const {return solid;}
    void setSolid(Solid* obj) {solid = obj;}
    int getWidth() const {return nwidth;}
    int getHeight() const {return nheight;}
    int getDepth() const {return ndepth;}
    FluidQuantity* getU() const {return u;}
    FluidQuantity* getV() const {return v;}
    FluidQuantity* getW() const {return w;}
    
    
    void Update();
    void CalculateVolumes();
    vectord getVelocity(double x, double y, double z);
    void SetSolidBoundaries();
    
    
};



#endif /* defined(__Fluid__grid__) */

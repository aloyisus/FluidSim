#ifndef __Fluid__grid__
#define __Fluid__grid__

#include "openvdb/openvdb.h"
#include <openvdb/tools/Interpolation.h>
#include <openvdb/math/ConjGradient.h>

#include <iostream>
#include <chrono>

typedef openvdb::math::pcg::SparseStencilMatrix<double, 7> SymmBandMatrix;
typedef openvdb::math::pcg::IncompleteCholeskyPreconditioner<SymmBandMatrix> CholeskyPrecondMatrix;

#define AMBIENT_TEMPERATURE 273.0
#define GRAVITY -9.8
#define TOL 1e-5
#define MAXITER 100000



class FluidGrid;


//Abstract solid class
class Solid {
    
protected:
    
    openvdb::Vec3d centrepos;
    FluidGrid* parentgrid;
    // Solid() : centrepos(3,0) { };
    
public:
    
    virtual bool PointIsInside(double x, double y, double z) = 0;
    
    virtual openvdb::Vec3d getVelocity(openvdb::Vec3d xyz) = 0;
    
    virtual openvdb::Vec3d getNormal(openvdb::Vec3d xyz) = 0;
    
    virtual void ClosestSurface(double &x, double &y, double &z) = 0;

};


class Cuboid : public Solid {
    

    double w,h,d;
    double alpha;
    double omega;
    
public:
    
    Cuboid(FluidGrid* parent, double centrex, double centrey, double centrez, double width, double height, double depth, double alpha, double omega);
    
    bool PointIsInside(double x, double y, double z);
    openvdb::Vec3d getVelocity(openvdb::Vec3d xyz);
    openvdb::Vec3d getNormal(openvdb::Vec3d xyz);
    void ClosestSurface(double &x, double &y, double &z);
    
};


class Sphere : public Solid {

    double radius;
    
public:
    
    Sphere(FluidGrid* parent, double centrex, double centrey, double centrez, double radius);
    
    bool PointIsInside(double x, double y, double z);
    openvdb::Vec3d getVelocity(openvdb::Vec3d xyz);
    openvdb::Vec3d getNormal(openvdb::Vec3d xyz);
    void ClosestSurface(double &x, double &y, double &z);
    
};


class FluidQuantity {

    openvdb::FloatGrid::Ptr quantity;
    openvdb::FloatGrid::Ptr buffer;
    openvdb::FloatGrid::Ptr volume;
    
    openvdb::FloatGrid::Accessor* q_access;
    openvdb::FloatGrid::Accessor* b_access;
    openvdb::FloatGrid::Accessor* v_access;

    FluidGrid* parentgrid;

    int xsamples, ysamples, zsamples, numsamples;
    double offsetx, offsety, offsetz;

    double CatmullRom(double x, double f0, double f1, double f2, double f3);
    double InterpolateCMSlice(double x, double y, int k);
    
public:
    
    FluidQuantity(FluidGrid* parent, double ofx, double ofy, double ofz, double value=0);
    
    //Destuctor for FluidQuantity
    ~FluidQuantity();

    //Accessors
    int getxSamples(){ return xsamples; };
    int getySamples(){ return ysamples; };
    int getzSamples(){ return zsamples; };
    
    double getOffsetx(){ return offsetx; };
    double getOffsety(){ return offsety; };
    double getOffsetz(){ return offsetz; };
    
    openvdb::math::Transform::Ptr getTransform(){return quantity->transformPtr();};

    inline double getQuantity(int i, int j, int k) const {
        return q_access->getValue(openvdb::Coord(i, j, k));
    }

    inline void setQuantity(int i, int j, int k, double value) {
        q_access->setValue(openvdb::Coord(i, j, k), value);
    }
    
    inline double getBuffer(int i, int j, int k) const {
        return b_access->getValue(openvdb::Coord(i, j, k));
    }

    inline void setBuffer(int i, int j, int k, double value) {
        b_access->setValue(openvdb::Coord(i, j, k), value);
    }
    
    inline double getVolume(int i, int j, int k) const {
        return v_access->getValue(openvdb::Coord(i, j, k));
    }

    inline void setVolume(int i, int j, int k, double value) {
        v_access->setValue(openvdb::Coord(i, j, k), value);
    }
    
    double sum() const;
    double max() const;
   
    inline void CopyQuantitytoBuffer() {
        buffer = quantity->deepCopy(); // This will copy the entire grid NOTE WE WERE USING JUST copy() before which is shallow
        delete b_access;
        b_access = new openvdb::FloatGrid::Accessor(buffer->getAccessor());
    }
    
    inline void swap() {
        std::swap(quantity, buffer);
        std::swap(q_access, b_access);
    }
    
    openvdb::Vec3d indexToWorld(int i, int j, int k);

    void Advect();
  
    void AddForces();
    
    void addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double v);

    void ExtrapolateVelocity();
    
    double InterpolateLinear(double x, double y, double z) const;
    double InterpolateCubic(double x, double y, double z) const;    
    double InterpolateCM(double x, double y, double z);
    
    void CalculateVolumes(double* supervol);
};


class FluidGrid {
    
    int nwidth, nheight, ndepth, ncells;
    
    double cellwidth;
    
    double density;
    
    double CFLnumber;

    double currtime;
    double nextframetime;
    double framedeltaT;
    double deltaT;
    int framenumber;
    bool writetocache;

    openvdb::BBoxd bbox;

    std::string outputpath;
    
    //3d grid is written to openvdb container
    openvdb::FloatGrid::Ptr grid;

    //Matrices for pressure coefficients A
    SymmBandMatrix *A;
   
    //vectors for pressure and rhs of matrix equation
    SymmBandMatrix::VectorType *pressure;
    SymmBandMatrix::VectorType *r;

    //Velocity components
    FluidQuantity* u;
    FluidQuantity* v;
    FluidQuantity* w;
    
    //Quantities defined on grid eg. smoke, temperature
    FluidQuantity* smoke;
    FluidQuantity* temperature;
    
    Solid* solid;

    //Accessors
    double getPressure(int i, int j, int k ) const { return (*pressure)[k*nwidth*nheight + j*nwidth + i]; }
    double& setPressure( int i, int j, int k ){ return (*pressure)[k*nwidth*nheight + j*nwidth + i]; }
    double getr(int i, int j, int k ) const { return (*r)[k*nwidth*nheight + j*nwidth + i]; }
    double& setr( int i, int j, int k ){ return (*r)[k*nwidth*nheight + j*nwidth + i]; }
    int getframenumber() const { return framenumber; }
    
    void addToA(int cellindex, int offset, double value);
    void BuildLinearSystem();
    void Advect();

    void updateVelocities();
    void Project();
    void AddForces();
    void setDeltaT();
    void WriteToCache();
    void addSmoke(double x, double y, double z, double wh, double h, double l, double d, double t);

public:

    //Constructor for FluidGrid
    FluidGrid(int w, int h, int d, double tstep, double rh, std::string filepath);
    
    //Destuctor for FluidGrid
    ~FluidGrid();
    
    //Accessors
    double getDeltaT() const { return deltaT; }
    double getCellwidth() const { return cellwidth; }
    openvdb::BBoxd getBBox() const {return bbox;}
    double getCurrtime() const {return currtime;}
    Solid* getSolid() const {return solid;}
    void setSolid(Solid* obj) {solid = obj;}
    int getWidth() const {return nwidth;}
    int getHeight() const {return nheight;}
    int getDepth() const {return ndepth;}
    int getSize() const {return ncells;}
    FluidQuantity* getU() const {return u;}
    FluidQuantity* getV() const {return v;}
    FluidQuantity* getW() const {return w;}
    
    
    void Update();
    void CalculateVolumes();
    openvdb::Vec3d getVelocity(openvdb::Vec3d xyz);
    void SetSolidBoundaries();
    void SetWallBoundaries();
    
    double maxDivergence() const;
    
};


#endif /* defined(__Fluid__grid__) */

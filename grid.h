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
    
    openvdb::math::Transform::Ptr getTransform() const {return quantity->transformPtr();}

    double getQuantity(int i, int j, int k) const {
        return q_access->getValue(openvdb::Coord(i, j, k));
    }

    void setQuantity(int i, int j, int k, double value) {
        q_access->setValue(openvdb::Coord(i, j, k), value);
    }
    
    double getBuffer(int i, int j, int k) const {
        return b_access->getValue(openvdb::Coord(i, j, k));
    }

    void setBuffer(int i, int j, int k, double value) {
        b_access->setValue(openvdb::Coord(i, j, k), value);
    }
    
    double getVolume(int i, int j, int k) const {
        return v_access->getValue(openvdb::Coord(i, j, k));
    }

    void setVolume(int i, int j, int k, double value) {
        v_access->setValue(openvdb::Coord(i, j, k), value);
    }
    
    double sum() const;
    double max() const;
    double sumVolumes() const;
   
    void CopyQuantitytoBuffer() {
        buffer = quantity->deepCopy(); // This will copy the entire grid NOTE WE WERE USING JUST copy() before which is shallow
        delete b_access;
        b_access = new openvdb::FloatGrid::Accessor(buffer->getAccessor());
    }
    
    void swap() {
        std::swap(quantity, buffer);
        std::swap(q_access, b_access);
    }
    
    openvdb::Vec3d indexToWorld(int i, int j, int k);

    void Advect();
  
    void AddForces();
    
    void addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double v);
   
    double InterpolateLinear(openvdb::Vec3d xyz) const;
    double InterpolateQuadratic(openvdb::Vec3d xyz) const;    
    double InterpolateCubic(openvdb::Vec3d xyz) const;

    void CalculateVolumes(double*);
    
};


class FluidGrid {
    
    public:

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
    void CalculateVolumes();

// public:

    //Constructor for FluidGrid
    FluidGrid(int w, int h, int d, double tstep, double rh, std::string filepath);
    
    //Destuctor for FluidGrid
    ~FluidGrid();
    
    //Accessors
    double getDeltaT() const { return deltaT; }
    double getCellwidth() const { return cellwidth; }
    openvdb::BBoxd getBBox() const {return bbox;}
    double getCurrtime() const {return currtime;}
    int getWidth() const {return nwidth;}
    int getHeight() const {return nheight;}
    int getDepth() const {return ndepth;}
    int getSize() const {return ncells;}
    FluidQuantity* getU() const {return u;}
    FluidQuantity* getV() const {return v;}
    FluidQuantity* getW() const {return w;}
    
    
    void Update();
    void SetWallBoundaries();
    
    double maxDivergence() const;

    void averagevoxelsaround(openvdb::Vec3d ijk);
    
};

bool isMatrixSymmetric(const SymmBandMatrix&, SymmBandMatrix::ValueType);
bool rowsSumToZero(const SymmBandMatrix&, SymmBandMatrix::ValueType);

template <typename SparseMatrixType>
void clearSparseMatrix(SparseMatrixType* matrix);

#endif /* defined(__Fluid__grid__) */

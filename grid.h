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

template <class T>
class FluidQuantity {

    typename T::Accessor* q_access;
    typename T::Accessor* b_access;

    FluidGrid* parentgrid;

    int xsamples, ysamples, zsamples;
    double offsetx, offsety, offsetz;

    bool staggered;
    
public:

    using ValueType = typename T::ValueType;

    typename T::Ptr quantity;
    typename T::Ptr buffer;
    
    FluidQuantity(FluidGrid* parent, bool staggered, ValueType value);
    
    ~FluidQuantity() {
        // No explicit cleanup required for OpenVDB grids
        delete q_access;
        delete b_access;
    }

    //Accessors
    int getxSamples(){ return xsamples; };
    int getySamples(){ return ysamples; };
    int getzSamples(){ return zsamples; };
    
    double getOffsetx(){ return offsetx; };
    double getOffsety(){ return offsety; };
    double getOffsetz(){ return offsetz; };
    
    openvdb::math::Transform::Ptr getTransform() const {return quantity->transformPtr();}

    ValueType getQuantity(int i, int j, int k) const {
        return q_access->getValue(openvdb::Coord(i, j, k));
    }

    void setQuantity(int i, int j, int k, ValueType value) {
        q_access->setValue(openvdb::Coord(i, j, k), value);
    }

    void setBuffer(int i, int j, int k, ValueType value) {
        b_access->setValue(openvdb::Coord(i, j, k), value);
    }
    
    ValueType max() const;
       
    void swap() {
        std::swap(quantity, buffer);
        std::swap(q_access, b_access);
    }
    
    void Advect();
  
    void AddForces();
    
    void addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, ValueType v);
   
    ValueType InterpolateLinear(openvdb::Vec3d xyz) const;
    ValueType StaggeredInterpolateLinear(openvdb::Vec3d xyz) const;
    ValueType InterpolateQuadratic(openvdb::Vec3d xyz) const;    
    ValueType InterpolateCubic(openvdb::Vec3d xyz) const;
    
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
    FluidQuantity<openvdb::Vec3dGrid>* velocity;

    //Quantities defined on grid eg. smoke, temperature
    FluidQuantity<openvdb::FloatGrid>* smoke;
    FluidQuantity<openvdb::FloatGrid>* temperature;

    //Accessors
    double getPressure(int i, int j, int k ) const { return (*pressure)[k*nwidth*nheight + j*nwidth + i]; }
    int getframenumber() const { return framenumber; }

    void addToA(int i, int offset, double value){
        switch (offset){
            case 0:
                A->setValue(i, i, A->getValue(i, i) + value);
                break;
            case 1:
                if (i + 1 < ncells){
                    A->setValue(i, i+1, A->getValue(i, i+1) + value);
                    A->setValue(i+1, i, A->getValue(i+1, i) + value);
                }
                break;
            case 2:
                if (i + nwidth < ncells){
                    A->setValue(i, i + nwidth, A->getValue(i, i + nwidth) + value);
                    A->setValue(i + nwidth, i, A->getValue(i + nwidth, i) + value);
                }
                break;
            case 3:
                if (i + nwidth*nheight < ncells){
                    A->setValue(i, i + nwidth*nheight, A->getValue(i, i + nwidth*nheight) + value);
                    A->setValue(i + nwidth*nheight, i, A->getValue(i + nwidth*nheight, i) + value);
                }
                break;
        }
    }


    void setA(int i, int offset, double value){
        switch (offset){
            case 0:
                A->setValue(i, i, value);
                break;
            case 1:
                if (i + 1 < ncells){
                    A->setValue(i, i+1, value);
                    A->setValue(i+1, i, value);
                }
                break;
            case 2:
                if (i + nwidth < ncells){
                    A->setValue(i, i + nwidth, value);
                    A->setValue(i + nwidth, i, value);
                }
                break;
            case 3:
                if (i + nwidth*nheight < ncells){
                    A->setValue(i, i + nwidth*nheight, value);
                    A->setValue(i + nwidth*nheight, i, value);
                }
                break;
        }
    }


    void BuildLinearSystem();
    void Advect();
    void updateVelocities();
    void Project();
    void AddForces();
    void setDeltaT();
    void WriteToCache();
    void addSmoke(double x, double y, double z, double wh, double h, double l, double d, double t);

// public:

    //Constructor for FluidGrid
    FluidGrid(int w, int h, int d, double tstep, double rh, std::string filepath);
    
    ~FluidGrid(){
        delete velocity;
        delete smoke;
        delete temperature;
        delete pressure;
        delete r;
        delete A;
    }
    
    //Accessors
    double getDeltaT() const { return deltaT; }
    double getCellwidth() const { return cellwidth; }
    openvdb::BBoxd getBBox() const {return bbox;}
    double getCurrtime() const {return currtime;}
    int getWidth() const {return nwidth;}
    int getHeight() const {return nheight;}
    int getDepth() const {return ndepth;}
    FluidQuantity<openvdb::Vec3dGrid>* getVelocity() const {return velocity;}

    void Update();
    void SetWallBoundaries();
    
    double maxDivergence() const{

        double maxdiv = 0;
        double div;

        for (int k=0; k<ndepth; k++)
            for (int j=0; j<nheight; j++)
                for (int i=0; i<nwidth; i++){
                    
                    div =  fabs(  velocity->getQuantity(i+1,j,k)[0] - velocity->getQuantity(i,j,k)[0]
                                + velocity->getQuantity(i,j+1,k)[1] - velocity->getQuantity(i,j,k)[1]
                                + velocity->getQuantity(i,j,k+1)[2] - velocity->getQuantity(i,j,k)[2])/cellwidth;
                    
                    maxdiv = fmax(div,maxdiv);
                    
                }
        return maxdiv;
    }
    
    bool sanityChecks(double eps) {

        int slice_size  = nwidth * nheight;
        int total_cells = nwidth * nheight * ndepth;
        for(int k=0; k<ndepth; k++)
            for(int j=0; j<nheight; j++)
                for(int i=0; i<nwidth; i++){

                    auto ix = i + (j * nwidth) + (k * slice_size);
                    auto row = A->getConstRow(ix);

                    int nn = 0;
                    int count = 0;
                    SymmBandMatrix::ValueType sum = 0;
                    // Iterate ONLY over the non-zero elements in this row
                    for (auto it = row.cbegin(); it; ++it) {
                        sum += *it;
                        count++;

                        int jx = it.column();
                        SymmBandMatrix::ValueType val_ij = *it;
                        if ((jx==ix) && !(val_ij > eps)){
                            std::cout << "non-positive diagonal element " << ix << std::endl;
                            return false;
                        }
                        
                        // Fetch the mirrored value (A_ji)
                        if (jx != ix){
                            SymmBandMatrix::ValueType val_ji = A->getValue(jx, ix);
                            
                            // Floating point comparison
                            if (std::abs(val_ij - val_ji) > eps) {
                                std::cout << "asymmetry detected at row " << ix << " column " << jx << " " << val_ij << " " << val_ji << std::endl;
                                return false; 
                            }
                        }
                    }

                    if (i > 0) nn += 1;
                    if (i < nwidth-1) nn += 1;
                    if (j > 0) nn += 1;
                    if (j < nheight-1) nn += 1;                
                    if (k > 0) nn += 1;
                    if (k < ndepth-1) nn += 1;

                    if (count != nn + 1){
                        std::cout << "count in row ix = " << ix << "should be " << nn+1 << " but is " << count << std::endl;
                        return false;
                    }

                    if (!count){
                        std::cout << "no non-zero elements in row " << i << std::endl; 
                        return false;
                    }

                    if (std::abs(sum) > eps){
                        std::cout << count << " elements summing to non-zero in row " << i << " total is " << std::abs(sum) << std::endl;

                        for (auto it = row.cbegin(); it; ++it) {
                            std::cout << "  col=" << it.column()
                                    << " value=" << *it
                                    << '\n';
                        }

                        return false;
                    }

        }
        return true;
    }

};


template <typename SparseMatrixType>
void clearSparseMatrix(SparseMatrixType* matrix);

#endif /* defined(__Fluid__grid__) */

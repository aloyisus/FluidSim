#ifndef __Fluid__grid__
#define __Fluid__grid__

#include <iostream>
#include "openvdb/openvdb.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>


#define AMBIENT_TEMPERATURE 273.0
#define GRAVITY -9.8;
#define TOL double(0.001)
#define MAXITER 10000


Eigen::SparseMatrix<double>* fillSymmBandMatrix(int size, int offset1, int offset2, int offset3, double value);

class FluidGrid;


//Abstract solid class
class Solid {
    
protected:
    
    Eigen::Vector3d centrepos;
    FluidGrid* parentgrid;
    // Solid() : centrepos(3,0) { };
    
public:
    
    virtual bool PointIsInside(double x, double y, double z) = 0;
    
    virtual Eigen::Vector3d getVelocity(double x, double y, double z) = 0;
    
    virtual Eigen::Vector3d getNormal(double x, double y, double z) = 0;
    
    virtual void ClosestSurface(double &x, double &y, double &z) = 0;
    
    
    
};




class Cuboid : public Solid {
    

    double w,h,d;
    double alpha;
    double omega;
    
public:
    
    Cuboid(FluidGrid* parent, double centrex, double centrey, double centrez, double width, double height, double depth, double alpha, double omega);
    
    bool PointIsInside(double x, double y, double z);
    Eigen::Vector3d getVelocity(double x, double y, double z);
    Eigen::Vector3d getNormal(double x, double y, double z);
    void ClosestSurface(double &x, double &y, double &z);
    
};




class Sphere : public Solid {

    double radius;
    
public:
    
    Sphere(FluidGrid* parent, double centrex, double centrey, double centrez, double radius);
    
    bool PointIsInside(double x, double y, double z);
    Eigen::Vector3d getVelocity(double x, double y, double z);
    Eigen::Vector3d getNormal(double x, double y, double z);
    void ClosestSurface(double &x, double &y, double &z);
    
};




class FluidQuantity {

    Eigen::VectorXd* quantity;
    Eigen::VectorXd* buffer;
    Eigen::VectorXd* volume;
    
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
    
    
    double getQuantity(int i, int j, int k) const{ return (*quantity)(i + j*xsamples + k*xsamples*ysamples); };
    double &setQuantity(int i, int j, int k){ return (*quantity)(i + j*xsamples + k*xsamples*ysamples); };
    
    double getBuffer(int i, int j, int k) const{ return (*buffer)(i + j*xsamples + k*xsamples*ysamples); };
    double &setBuffer(int i, int j, int k){ return (*buffer)(i + j*xsamples + k*xsamples*ysamples); };
    
    double getVolume(int i, int j, int k) const{ return (*volume)(i + j*xsamples + k*xsamples*ysamples); };
    double &setVolume(int i, int j, int k){ return (*volume)(i + j*xsamples + k*xsamples*ysamples); };
    
    double sum() const { return quantity->sum(); }
    double max() const { return quantity->maxCoeff(); }
   
    void CopyQuantitytoBuffer(){ *buffer = *quantity; }
    
    void swap(){ std::swap(quantity, buffer); }
    
    void Advect();
    
    void AddForces();
    
    void addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double v);

    void ExtrapolateVelocity();
    
    double InterpolateLinear(double x, double y, double z) const;
    
    void CalculateVolumes(double* supervol);
};




class FluidGrid {
    
    int nwidth,nheight,ndepth;
    
    double cellwidth;
    
    double density;
    
    double CFLnumber;

    double currtime;
    double nextframetime;
    double framedeltaT;
    double deltaT;
    int framenumber;
    bool writetocache;

    std::string outputpath;
    
    //3d grid is written to openvdb container
    openvdb::FloatGrid::Ptr grid;

    //Matrices for pressure coefficients A and MIC(0) preconditioner Ei
    Eigen::SparseMatrix<double>* A;
    Eigen::MatrixXd*  Ei;
   
    //vectors for pressure and rhs of matrix equation
    Eigen::VectorXd* pressure;
    Eigen::VectorXd *r; /* Right hand side of pressure solve */

    //Velocity components
    FluidQuantity* u;
    FluidQuantity* v;
    FluidQuantity* w;
    
    //Quantities defined on grid eg. smoke, temperature
    FluidQuantity* smoke;
    FluidQuantity* temperature;
    
    Solid* solid;

    //Accessors
    double getPressure(int i, int j, int k ) const { return (*pressure)(k*nwidth*nheight + j*nwidth + i); }
    double& setPressure( int i, int j, int k ){ return (*pressure)(k*nwidth*nheight + j*nwidth + i); }
    double getr(int i, int j, int k ) const { return (*r)(k*nwidth*nheight + j*nwidth + i); }
    double& setr( int i, int j, int k ){ return (*r)(k*nwidth*nheight + j*nwidth + i); }
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
    Eigen::Vector3d getVelocity(double x, double y, double z);
    void SetSolidBoundaries();
    void SetWallBoundaries();
    
    double maxDivergence() const;
    
};



#endif /* defined(__Fluid__grid__) */

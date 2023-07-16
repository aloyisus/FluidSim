#include "grid.h"

#include <Eigen/IterativeLinearSolvers>

#define round5(n) floor(n * 100000 + 0.5)/100000


/**
 * Create a sparse symmetric banded matrix.
 *
 * This function creates a sparse symmetric banded matrix of size `size` with the specified offsets and values.
 * The elements within the specified offsets will have the value `value`, and the elements outside the offsets will be zero.
 *
 * @param size The size of the square matrix.
 * @param offset1 The first offset for the banded structure.
 * @param offset2 The second offset for the banded structure.
 * @param offset3 The third offset for the banded structure.
 * @param value The value to assign to the elements within the offsets.
 */
 void fillSymmBandMatrix(Eigen::SparseMatrix<double>* mat, int size, int offset1, int offset2, int offset3, double value){
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < size; i++){
            triplets.push_back(Eigen::Triplet<double>(i, i, value));
            if (i + offset1 < size){
                triplets.push_back(Eigen::Triplet<double>(i, i + offset1, value));
                triplets.push_back(Eigen::Triplet<double>(i + offset1, i, value));
            }
            if (i + offset2 < size){
                triplets.push_back(Eigen::Triplet<double>(i, i + offset2, value));
                triplets.push_back(Eigen::Triplet<double>(i + offset2, i, value));
            }
            if (i + offset3 < size){
                triplets.push_back(Eigen::Triplet<double>(i, i + offset3, value));
                triplets.push_back(Eigen::Triplet<double>(i + offset3, i, value));
            }
    }
    mat->setFromTriplets(triplets.begin(), triplets.end());
}

/**
 * @brief Constructs a new Cuboid instance.
 *
 * This constructor initializes a new Cuboid instance with the specified position, dimensions,
 * and other parameters. Note that the position and dimensions are given in terms of cells,
 * not spatial coordinates.
 *
 * @param parent Pointer to the FluidGrid that this cuboid is a part of.
 * @param centrex X-coordinate of the cuboid's center in cell units.
 * @param centrey Y-coordinate of the cuboid's center in cell units.
 * @param centrez Z-coordinate of the cuboid's center in cell units.
 * @param width Width of the cuboid in cell units.
 * @param height Height of the cuboid in cell units.
 * @param depth Depth of the cuboid in cell units.
 * @param alph Initial rotation of cuboid in radians.
 * @param omg Angular velocity of cuboid in radians/s.
 */
Cuboid::Cuboid(FluidGrid* parent, double centrex, double centrey, double centrez, double width, double height, double depth, double alph, double omg) {
    
    parentgrid = parent;
    
    centrepos(0) = centrex;
    centrepos(1) = centrey;
    centrepos(2) = centrez;

    w = width;
    h = height;
    d = depth;
    
    alpha = alph;
    omega = omg;
    
}

/**
 * Determines whether a given point is located within the cuboid.
 *
 * This method takes the coordinates of a point (x, y, z) and checks if it lies 
 * within the spatial domain of the cuboid. The cuboid's rotation is taken into 
 * account based on the current time of the parent fluid grid. A rotation matrix 
 * is formed from `alpha` and `omega` parameters, which represent 
 * a rotation angle and an angular velocity, respectively. The rotation of the 
 * cuboid oscillates over time with the frequency of `omega` and the phase of 
 * `alpha`, both in radians.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return True if the point is inside the cuboid, false otherwise.
 */
bool Cuboid::PointIsInside(double x, double y, double z){
    
    Eigen::Vector3d cellpos(x,y,z);
    
    double t = parentgrid->getCurrtime();
    
    bool temp;
    Eigen::Matrix3d invtransform;
    
    //rotation of -theta about z-axis
    invtransform.setZero();
    invtransform(0,0) = cos(alpha + omega*t);
    invtransform(1,0) = -sin(alpha + omega*t);
    invtransform(0,1) = sin(alpha + omega*t);
    invtransform(1,1) = cos(alpha + omega*t);
    invtransform(2,2) = 1.0;
    
    cellpos = invtransform*(cellpos - centrepos);
    
    temp = (cellpos(0) >= - w/2) && (cellpos(0) <= + w/2) &&
           (cellpos(1) >= - h/2) && (cellpos(1) <= + h/2) &&
           (cellpos(2) >= - d/2) && (cellpos(2) <= + d/2);
    
    return temp;
    
}

/**
 * Computes the velocity at a given point in the cuboid.
 *
 * This function calculates the velocity vector for a point located at coordinates 
 * (x, y, z) within the cuboid. The cuboid is assumed to be rotating about the z-axis 
 * with angular velocity `omega`. The velocity of a point in a rotating body is 
 * proportional to its radial distance from the axis of rotation and the direction 
 * of velocity is perpendicular to the radial vector (following the right-hand rule). 
 * The magnitude of the velocity is then scaled by the width of the cells in the 
 * parent grid.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return The velocity vector at the given point (x, y, z).
 */
Eigen::Vector3d Cuboid::getVelocity(double x, double y, double z){

    //vector from axis of rotation (z-axis) to x,y,z
    double rx = x - centrepos(0);
    double ry = y - centrepos(1);

    double distance = sqrt(rx*rx + ry*ry);

    //get the unit normal to vector (rx,ry)
    double normx = -ry / distance;
    double normy = rx / distance;

    //get magnitude of velocity at x,y, scale the norm vector by this
    normx *= distance*omega;
    normy *= distance*omega;

    Eigen::Vector3d temp(normx, normy, 0);
    temp *=  parentgrid->getCellwidth();

    return temp;
}

/**
 * Computes the normal vector at a given point in the cuboid.
 *
 * This function calculates the normal vector for a point located at coordinates 
 * (x, y, z) within the cuboid. The cuboid is assumed to be rotating about the z-axis 
 * with angular velocity `omega`. The rotation is characterized by an angle `alpha` 
 * which varies with time (t = `omega` * current_time). 
 * 
 * The function first transforms the point coordinates into the rotating frame of 
 * the cuboid and determines which face of the cuboid the point is closest to. 
 * It assigns a unit normal vector corresponding to that face. The normal vector 
 * is then transformed back to the original frame by inverse rotation.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return The normal vector at the given point (x, y, z).
 */
Eigen::Vector3d Cuboid::getNormal(double x, double y, double z){

    Eigen::Vector3d temp;
    Eigen::Matrix3d transform;
    Eigen::Vector3d r(x,y,z);

    double t = parentgrid->getCurrtime();

    //rotation of -theta about z-axis
    transform.setZero();
    transform(0,0) = cos(alpha + omega*t);
    transform(1,0) = -sin(alpha + omega*t);
    transform(0,1) = sin(alpha + omega*t);
    transform(1,1) = cos(alpha + omega*t);
    transform(2,2) = 1.0;

    r = transform*(r - centrepos);

    x = r(0);
    y = r(1);
    z = r(2);

    if ((y >= fabs(z)*h/d) && (y >= fabs(x)*h/w)){
        
        temp(0) = 0;
        temp(1) = 1;
        temp(2) = 0;
    }
    else if ((y <= fabs(z)*h/d) && (y <= fabs(x)*h/w)){
        
        temp(0) = 0;
        temp(1) = -1;
        temp(2) = 0;
    }
    else if ((x >= fabs(z)*w/d) && (x >= fabs(y)*w/h)){
        
        temp(0) = 1;
        temp(1) = 0;
        temp(2) = 0;
    }
    else if ((x <= fabs(z)*w/d) && (x <= fabs(y)*w/h)){
        
        temp(0) = -1;
        temp(1) = 0;
        temp(2) = 0;
    }
    else if ((z >= fabs(x)*d/w) && (z >= fabs(y)*d/h)){
        
        temp(0) = 0;
        temp(1) = 0;
        temp(2) = 1;
    }
    else if ((z <= fabs(x)*d/w) && (z <= fabs(y)*d/h)){
        
        temp(0) = 0;
        temp(1) = 0;
        temp(2) = -1;
    }

    //rotation of +theta about z-axis
    transform(0,0) = cos(alpha + omega*t);
    transform(1,0) = sin(alpha + omega*t);
    transform(0,1) = -sin(alpha + omega*t);
    transform(1,1) = cos(alpha + omega*t);
    transform(2,2) = 1.0;
    
    temp = transform * temp;
    
    return temp;
}

/**
 * Modifies the input coordinates (x, y, z) to represent the closest surface point 
 * on the cuboid if the original point is inside the cuboid.
 *
 * This function checks if a given point (x, y, z) is located inside a rotating cuboid. 
 * If so, it transforms the point into the rotating frame of the cuboid, calculates 
 * the distances to the cuboid's surfaces in each direction, and adjusts the point 
 * coordinates to match the closest surface point.
 * 
 * The function then transforms the updated point back into the original frame. 
 * The cuboid's rotation is characterized by an angle `alpha` which varies with time 
 * (t = `omega` * current_time).
 *
 * The original point coordinates are modified in-place.
 *
 * @param x Reference to the x-coordinate of the point. 
 *          Modified to match the x-coordinate of the closest surface point if the 
 *          original point is inside the cuboid.
 * @param y Reference to the y-coordinate of the point. 
 *          Modified to match the y-coordinate of the closest surface point if the 
 *          original point is inside the cuboid.
 * @param z Reference to the z-coordinate of the point. 
 *          Modified to match the z-coordinate of the closest surface point if the 
 *          original point is inside the cuboid.
 */
void Cuboid::ClosestSurface(double &x, double &y, double &z) {

    Eigen::Matrix3d transform;
    Eigen::Vector3d r(x,y,z);

    double t = parentgrid->getCurrtime();

    //rotation of -theta about z-axis
    transform.setZero();
    transform(0,0) = cos(alpha + omega*t);
    transform(1,0) = -sin(alpha + omega*t);
    transform(0,1) = sin(alpha + omega*t);
    transform(1,1) = cos(alpha + omega*t);
    transform(2,2) = 1.0;

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
    transform(0,0) = cos(alpha + omega*t);
    transform(1,0) = sin(alpha + omega*t);
    transform(0,1) = -sin(alpha + omega*t);
    transform(1,1) = cos(alpha + omega*t);
    transform(2,2) = 1.0;

    r = centrepos + transform*r;

    x = r(0);
    y = r(1);
    z = r(2);
}

/**
 * Constructor for the `Sphere` class.
 *
 * This constructor initializes a `Sphere` object, which represents a sphere within a fluid grid.
 * The position of the sphere and its radius are given in terms of cell counts, not physical space.
 *
 * @param parent Pointer to the parent `FluidGrid` object in which this sphere exists. 
 *        The `Sphere` object uses this reference to access properties of the fluid grid.
 * @param centrex The x-coordinate of the sphere's center, in cell units.
 * @param centrey The y-coordinate of the sphere's center, in cell units.
 * @param centrez The z-coordinate of the sphere's center, in cell units.
 * @param rd The radius of the sphere, in cell units.
 */
Sphere::Sphere(FluidGrid* parent, double centrex, double centrey, double centrez, double rd) {

    parentgrid = parent;

    centrepos(0) = centrex;
    centrepos(1) = centrey;
    centrepos(2) = centrez;

    radius = rd;
}

/**
 * @brief Determines if a point is inside the sphere.
 *
 * This function checks if the given point (x, y, z) is inside the sphere. It does this by
 * calculating the squared distance from the center of the sphere to the point. If this 
 * distance is less than or equal to the square of the radius of the sphere, the point is 
 * inside the sphere.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return true if the point is inside the sphere, false otherwise.
 */
bool Sphere::PointIsInside(double x, double y, double z){

    Eigen::Vector3d cellpos(x,y,z);

    double sqrdistance = (centrepos - cellpos).dot(centrepos - cellpos);

    return (sqrdistance <= radius*radius);
}

/**
 * @brief Retrieves the velocity of the sphere at the specified point.
 *
 * As the sphere has no velocity component normal to its surface, this function always returns zero
 * regardless of the input position (x, y, z). This is represented by an Eigen::Vector3d containing three zeroes.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return An Eigen::Vector3d representing the velocity of the sphere at the point, which is always (0,0,0).
 */
Eigen::Vector3d Sphere::getVelocity(double x, double y, double z){
    Eigen::Vector3d temp(0,0,0);
    return temp;
}

/**
 * @brief Retrieves the velocity of the sphere at the specified point.
 *
 * As the sphere has no velocity component normal to its surface, this function always returns zero
 * regardless of the input position (x, y, z). This is represented by an Eigen::Vector3d containing three zeroes.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return An Eigen::Vector3d representing the velocity of the sphere at the point, which is always (0,0,0).
 */
Eigen::Vector3d Sphere::getNormal(double x, double y, double z){
    Eigen::Vector3d rpos(x,y,z);
    rpos = (rpos - centrepos);
    double mod = sqrt(rpos.dot(rpos));

    if (mod == 0)
        return rpos;

    rpos = (1/mod)*rpos;
    return rpos;
}

/**
 * @brief Backprojects a point from the sphere's interior to its surface.
 *
 * If the provided point (x, y, z) is inside the sphere, this function adjusts its coordinates 
 * to correspond to the nearest point on the sphere's surface. It uses the surface normal for the 
 * backprojection. 
 *
 * Note that the arguments are passed by reference and will be modified by this function.
 *
 * @param x The x-coordinate of the point, replaced by the x-coordinate of the nearest surface point if inside the sphere.
 * @param y The y-coordinate of the point, replaced by the y-coordinate of the nearest surface point if inside the sphere.
 * @param z The z-coordinate of the point, replaced by the z-coordinate of the nearest surface point if inside the sphere.
 */
void Sphere::ClosestSurface(double &x, double &y, double &z) {

    Eigen::Vector3d r;
    Eigen::Vector3d norm = getNormal(x,y,z);

    r = centrepos + radius*norm;

    x = r(0);
    y = r(1);
    z = r(2);
}

/**
 * @brief Constructor for FluidGrid class.
 *
 * Creates a FluidGrid object with specified dimensions, timestep, density, and file path for output. 
 * Grid dimensions are in number of cells, and the shortest dimension is normalized to 1 unit of spatial measure.
 * Each quantity of interest (velocity in three directions, smoke, and temperature) is represented by a FluidQuantity object.
 * OpenVDB library is initialized to create an empty floating-point grid.
 *
 * @param wh Width of the grid in cells.
 * @param ht Height of the grid in cells.
 * @param dh Depth of the grid in cells.
 * @param tstep Time interval between frames.
 * @param rh Density of the fluid.
 * @param filepath Path to the output file for storing simulation results.
 */
FluidGrid::FluidGrid(int wh, int ht, int dh, double tstep, double rh, std::string filepath){
    //the shortest dimension measures 1 spatial unit, so cellwidth is 1 / no. of cells in shortest dimension
    cellwidth = 1.0/std::min(std::min(wh, ht),dh);

    nwidth = wh;
    nheight = ht;
    ndepth = dh;

    density = rh;

    CFLnumber = 5;
    framenumber = 0;
    currtime = 0;
    framedeltaT = tstep;
    nextframetime = tstep;
    bool writetocache = false;

    outputpath = filepath;

    pressure = new Eigen::VectorXd(wh*ht*dh);
    r = new Eigen::VectorXd(wh*ht*dh);
    A = new Eigen::SparseMatrix<double>(wh*ht*dh, wh*ht*dh);
    u = new FluidQuantity(this, 0,0.5,0.5);
    v = new FluidQuantity(this, 0.5,0,0.5);
    w = new FluidQuantity(this, 0.5,0.5,0);
    smoke = new FluidQuantity(this, 0.5,0.5,0.5);
    temperature = new FluidQuantity(this, 0.5,0.5,0.5,AMBIENT_TEMPERATURE);

    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();
}

//destructor
FluidGrid::~FluidGrid(){
    delete u;
    delete v;
    delete w;
    delete smoke;
    delete temperature;
    delete pressure;
    delete r;
    delete Ei;
    delete A;  
}

/**
 * Add a value to the specified offsets in the sparse matrix A.
 *
 * This function adds the specified value to the elements of the matrix A based on the given offset and index.
 * The offset determines the position of the elements to modify relative to the current index.
 *
 * @param i The index of the row/column to modify.
 * @param offset The offset indicating the position of the elements to modify relative to the current index.
 * @param value The value to add to the specified elements.
 */
void FluidGrid::addToA(int i, int offset, double value){
    int size = A->rows();
    switch (offset){
        case 0:
            A->coeffRef(i, i) += value;
            break;
        case 1:
            if (i + 1 < size){
                A->coeffRef(i, i+1) += value;
                A->coeffRef(i+1, i) += value;
            }
            break;
        case 2:
            if (i + nwidth < size){
                A->coeffRef(i, i + nwidth) += value;
                A->coeffRef(i + nwidth, i) += value;
            }
            break;
        case 3:
            if (i + nwidth*nheight < size){
                A->coeffRef(i, i + nwidth*nheight) += value;
                A->coeffRef(i + nwidth*nheight, i) += value;
            }
            break;
    }
}

/**
 * @brief Constructs the system of linear equations, Ap = b, to solve for pressure.
 *
 * This function builds the matrix A of pressure coefficients and the vector b of
 * additional terms that account for external influences such as the divergence of 
 * the velocity field and solid boundaries. The pressure is calculated based on the
 * volumes weighted by the fluid quantities within each grid cell. The function 
 * constructs these matrices by iterating through the 3D fluid grid.
 */
void FluidGrid::BuildLinearSystem(){
    double Ascale =  deltaT/(density*cellwidth*cellwidth);
    double rscale = 1.0/cellwidth;
    int row, column;
    r->setZero();
    fillSymmBandMatrix(A, nwidth*nheight*ndepth, 1, nwidth, nwidth*nheight, 0); 

    //Setup up the matrix of pressure coefficients, weighted by volumes
    for(int k=0; k<ndepth; k++)
        for(int j=0; j<nheight; j++)
            for(int i=0; i<nwidth; i++){

                row = i+j*nwidth+k*nwidth*nheight;
                column = row;

                addToA(row, 0, Ascale*u->getVolume(i,j,k));

                addToA(row, 0, Ascale*v->getVolume(i,j,k));

                addToA(row, 0, Ascale*w->getVolume(i,j,k));
            
                addToA(row, 0, Ascale*u->getVolume(i+1,j,k));
                addToA(row, 1, -Ascale*u->getVolume(i+1,j,k));

                addToA(row, 0, Ascale*v->getVolume(i,j+1,k));
                addToA(row, 2, -Ascale*v->getVolume(i,j+1,k));
            
                addToA(row, 0, Ascale*w->getVolume(i,j,k+1));
                addToA(row, 3, -Ascale*w->getVolume(i,j,k+1));

                (*r)(row) = -rscale*(u->getVolume(i+1,j,k)*u->getQuantity(i+1,j,k) - u->getVolume(i,j,k)*u->getQuantity(i,j,k) +
                                      v->getVolume(i,j+1,k)*v->getQuantity(i,j+1,k) - v->getVolume(i,j,k)*v->getQuantity(i,j,k) +
                                      w->getVolume(i,j,k+1)*w->getQuantity(i,j,k+1) - w->getVolume(i,j,k)*w->getQuantity(i,j,k));

                if (solid != nullptr){
                    (*r)(row) += rscale*((u->getVolume(i+1,j,k)-smoke->getVolume(i,j,k))*
                                        solid->getVelocity(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz())(0) -
                                        (u->getVolume(i,j,k)-smoke->getVolume(i,j,k))*
                                        solid->getVelocity(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz())(0) +
                                        (v->getVolume(i,j+1,k)-smoke->getVolume(i,j,k))*
                                        solid->getVelocity(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz())(1) -
                                        (v->getVolume(i,j,k)-smoke->getVolume(i,j,k))*
                                        solid->getVelocity(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz())(1) +
                                        (w->getVolume(i,j,k+1)-smoke->getVolume(i,j,k))*
                                        solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz())(2) -
                                        (w->getVolume(i,j,k)-smoke->getVolume(i,j,k))*
                                        solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz())(2));
                }             
            };

}

/**
 * @brief Carries out advection for all fluid quantities within the grid.
 *
 * The function advances the physical properties of the fluid through the grid by 
 * a single step in time. The advection process is performed for each of the 
 * fluid's properties: temperature, smoke density, and the three components of 
 * velocity (u, v, w). Once advection has been completed, the old and new states 
 * of each fluid property are swapped, ready for the next time step.
 */
void FluidGrid::Advect(){
    
    temperature->Advect();
    smoke->Advect();
    u->Advect();
    v->Advect();
    w->Advect();
    
    u->swap();
    v->swap();
    w->swap();
    smoke->swap();
    temperature->swap();
    
}

/**
 * @brief Calculates the volume fractions of fluid quantity in each cell of its parent grid.
 *
 * This function iterates over the samples of the fluid quantity, calculates their positions in 'cell space', 
 * and uses these positions to access corresponding supervoxel fractions from the input `supervol` array. These fractions 
 * are subtracted from the initial volume (set as 8.0) to yield a volume fraction that represents the proportion of the cell 
 * not occupied by solid matter. This calculated volume fraction is stored back into the fluid quantity's volume field.
 * A calculated volume of 1.0 indicates that the cell is not occupied by any solid.
 * 
 * This function assumes that the `supervol` array is filled with supervoxel fractions that represent whether each supervoxel 
 * in the grid is inside a solid or not. This array is typically provided by `FluidGrid::CalculateVolumes()`.
 *
 * Note: The function skips the samples at the boundaries of the grid (i.e., u volume samples at (0,j,k) and (nwidth,j,k)).
 *
 * @param supervol Pointer to an array of supervoxel fractions, typically obtained from `FluidGrid::CalculateVolumes()`. 
 * This array should have dimensions twice as large as the parent grid in each dimension (hence, 2x2x2 supersampling).
 */
void FluidQuantity::CalculateVolumes(double* supervol){

    int nwidth = parentgrid->getWidth();
    int nheight = parentgrid->getHeight();

    //Note that u volume samples at (0,j,k) and (nwidth,j,k) are unchanged as zero
    //so we don't need to set them. To achieve this, we iterate from int(1-offset) to samples - int(1-offset)
    //because if offset == 0, int(1-offset) will be 1, otherwise it will be zero, so this ensures we skip over
    //volume samples for grid-aligned, offset == 0 samples
    for (int k=int(1-offsetz); k<getzSamples()-int(1-offsetz); k++)
        for (int j=int(1-offsety); j<getySamples()-int(1-offsety); j++)
            for (int i=int(1-offsetx); i<getxSamples()-int(1-offsetx); i++){

                //get position of sample in 'cell space'
                double posx = i+getOffsetx();
                double posy = j+getOffsety();
                double posz = k+getOffsetz();

                int ix0 = floor((posx-0.25)*2);
                int ix1 = ix0 + 1;

                int iy0 = floor((posy-0.25)*2);
                int iy1 = iy0 + 1;

                int iz0 = floor((posz-0.25)*2);
                int iz1 = iz0 + 1;

                setVolume(i,j,k) = 8.0;

                setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz0*4*nwidth*nheight];
                setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz0*4*nwidth*nheight];
                setVolume(i,j,k) -= supervol[ix0+iy0*2*nwidth+iz1*4*nwidth*nheight];
                setVolume(i,j,k) -= supervol[ix0+iy1*2*nwidth+iz1*4*nwidth*nheight];
                setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz0*4*nwidth*nheight];
                setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz0*4*nwidth*nheight];
                setVolume(i,j,k) -= supervol[ix1+iy0*2*nwidth+iz1*4*nwidth*nheight];
                setVolume(i,j,k) -= supervol[ix1+iy1*2*nwidth+iz1*4*nwidth*nheight];

                setVolume(i,j,k) = setVolume(i,j,k)/8.0;

            }
}

/**
 * @brief This function calculates the volume fractions of the fluid grid cells by 2x2x2 supersampling.
 * 
 * First, it creates an array of supervoxel volume fractions by checking whether each supervoxel (a 
 * sub-division of a cell in the fluid grid) is inside a solid. This is done by iterating over a grid 
 * twice as dense as the fluid grid in all dimensions (hence, 2x2x2 supersampling). It then passes this 
 * array to the `CalculateVolumes()` function of each of the u, v, w, and smoke fluid quantities.
 *
 * @param none.
 * @return void.
 */
void FluidGrid::CalculateVolumes(){
    double* supervol = new double[nwidth*nheight*ndepth*8];
    bool pointisinside;

    //Perform 2x2x2 supersampling over entire grid, so each ijk cell will be sampled 8 times
    for (int k=0; k<2*ndepth; k++)
        for (int j=0; j<2*nheight; j++)
            for (int i=0; i<2*nwidth; i++){
                pointisinside = false;
                if (solid != nullptr)
                    pointisinside = solid->PointIsInside(0.25+i*0.5, 0.25+j*0.5, 0.25+k*0.5);
                supervol[i+j*2*nwidth+k*4*nwidth*nheight] = pointisinside;
            }

    u->CalculateVolumes(supervol);
    v->CalculateVolumes(supervol);
    w->CalculateVolumes(supervol);
    smoke->CalculateVolumes(supervol);

    delete[] supervol;
}

/**
 * @brief Updates the velocity fields of the fluid grid based on the current pressure field.
 *
 * This function iterates over all the velocity components (u, v, w) stored in the grid. For each component, 
 * the function calculates the velocity update based on the pressure gradient between neighboring cells. 
 * This calculation is done according to the fluid dynamics equation: 
 * Δv = - (Δt / ρΔx) * ΔP, where Δt is the time step, ρ is the fluid density, Δx is the cell width, and ΔP is the pressure difference.
 * If the volume of the cell is non-zero, the function updates the velocity component by subtracting the calculated update. 
 * If the volume of the cell is zero, the function sets the velocity component to zero.
 * After the velocity components have been updated, the function extrapolates the velocity fields 
 * of temperature, smoke, u, v, and w and then enforces the boundary conditions.
 */
void FluidGrid::updateVelocities() {
    
    double scale = deltaT/(density*cellwidth);
    
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

    temperature->ExtrapolateVelocity();
    smoke->ExtrapolateVelocity();
    u->ExtrapolateVelocity();
    v->ExtrapolateVelocity();
    w->ExtrapolateVelocity();

    if (solid != nullptr)
        SetSolidBoundaries();
    SetWallBoundaries();
}

/**
 * @brief Retrieves the fluid velocity at a specific point in the grid.
 *
 * This function interpolates the fluid velocity from the surrounding grid points 
 * at the given position (x, y, z) using a linear method. The three components 
 * of the velocity vector (u, v, w) are each interpolated separately and then 
 * combined into the final velocity vector.
 *
 * @param x The x-coordinate of the point in the grid.
 * @param y The y-coordinate of the point in the grid.
 * @param z The z-coordinate of the point in the grid.
 * @return The interpolated fluid velocity vector at the point (x, y, z).
 */
Eigen::Vector3d FluidGrid::getVelocity(double x, double y, double z){

    Eigen::Vector3d temp;

    temp(0) = u->InterpolateLinear(x,y,z);
    temp(1) = v->InterpolateLinear(x,y,z);
    temp(2) = w->InterpolateLinear(x,y,z);

    return temp;
}

/**
 * @brief Sets the solid boundaries within the fluid grid.
 *
 * This function is responsible for properly setting the solid boundaries within the 
 * fluid grid. The boundary conditions are set in such a way that the velocity 
 * of the fluid is zero relative to the solid boundary. The function iterates over 
 * all cells in the grid, and for cells with zero volume (assumed to be solid cells),
 * it sets the velocities of the fluid in the bordering cells so that there is no
 * relative velocity in the direction of the solid boundary. This is achieved by
 * subtracting from the fluid velocity the component of the fluid velocity that 
 * is along the normal to the solid boundary.
 *
 * The function also sets the boundary conditions for the walls of the box, where the 
 * fluid velocity is set to zero.
 *
 * This method uses a buffer to store intermediate velocity values, and at the end, 
 * the buffer values are swapped with the actual fluid velocities.
 *
 * This method does not return any value.
 */
void FluidGrid::SetSolidBoundaries(){

    u->CopyQuantitytoBuffer();
    v->CopyQuantitytoBuffer();
    w->CopyQuantitytoBuffer();

    Eigen::Vector3d vfluid;
    Eigen::Vector3d vsolid;
    Eigen::Vector3d normal;
    
    //For cells which have zero volume, set the velocity samples bordering to have zero relative velocity in the
    //direction of the normal to the solid boundary
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++)
            for (int i = 0; i < nwidth; i++) {
                if (smoke->getVolume(i,j,k) == 0)
                {
                    vsolid = solid->getVelocity(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    vfluid = getVelocity(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    normal = solid->getNormal(i+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    
                    u->setBuffer(i,j,k) =  (vfluid - (vfluid - vsolid).dot(normal)*normal)(0);
                    
                    vsolid = solid->getVelocity(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    vfluid = getVelocity(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    normal = solid->getNormal(i+1+u->getOffsetx(),j+u->getOffsety(),k+u->getOffsetz());
                    
                    u->setBuffer(i+1,j,k) = (vfluid - (vfluid - vsolid).dot(normal)*normal)(0);
                    
                    vsolid = solid->getVelocity(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz());
                    vfluid = getVelocity(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz());
                    normal = solid->getNormal(i+v->getOffsetx(),j+v->getOffsety(),k+v->getOffsetz());
                    
                    v->setBuffer(i,j,k) = (vfluid - (vfluid - vsolid).dot(normal)*normal)(1);
                    
                    vsolid = solid->getVelocity(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz());
                    vfluid = getVelocity(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz());
                    normal = solid->getNormal(i+v->getOffsetx(),j+1+v->getOffsety(),k+v->getOffsetz());
                    
                    v->setBuffer(i,j+1,k) = (vfluid - (vfluid - vsolid).dot(normal)*normal)(1);
                    
                    vsolid = solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz());
                    vfluid = getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz());
                    normal = solid->getNormal(i+w->getOffsetx(),j+w->getOffsety(),k+w->getOffsetz());
                    
                    w->setBuffer(i,j,k) = (vfluid - (vfluid - vsolid).dot(normal)*normal)(2);
                    
                    vsolid = solid->getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz());
                    vfluid = getVelocity(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz());
                    normal = solid->getNormal(i+w->getOffsetx(),j+w->getOffsety(),k+1+w->getOffsetz());
                    
                    w->setBuffer(i,j,k+1) = (vfluid - (vfluid - vsolid).dot(normal)*normal)(2);
                }
            }
    u->swap();
    v->swap();
    w->swap();
}

void FluidGrid::SetWallBoundaries(){
    //Boundary conditions for the walls of the box
    //Set the normal component to the wall to zero
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



/**
 * @brief Calculates and returns the maximum divergence in the fluid grid.
 *
 * This function computes the divergence at each point in the grid, which is
 * a measure of the rate at which fluid is exiting or entering a small volume
 * around the point. The divergence is calculated based on the fluid velocities
 * in the x, y, and z directions (`u`, `v`, and `w` respectively).
 * The maximum absolute divergence value across all points in the grid is then 
 * returned. This function is useful for checking the accuracy of the simulation, 
 * as the divergence should be zero in incompressible flows.
 *
 * @return double The maximum absolute divergence value in the grid.
 */
double FluidGrid::maxDivergence() const{

    double maxdiv = 0;
    double div;

    for (int k=0; k<ndepth; k++)
        for (int j=0; j<nheight; j++)
            for (int i=0; i<nwidth; i++){
                
                div =  fabs(  u->getQuantity(i+1,j,k) - u->getQuantity(i,j,k)
                            + v->getQuantity(i,j+1,k) - v->getQuantity(i,j,k)
                            + w->getQuantity(i,j,k+1) - w->getQuantity(i,j,k))/cellwidth;
                
                maxdiv = fmax(div,maxdiv);
                
            }

    return maxdiv;
}

/**
 * @brief Adds smoke and temperature within a specified cuboid region of the fluid grid.
 *
 * This function is designed to add smoke and temperature inside a cuboid region within the fluid grid. 
 * The size and position of the cuboid are specified by the user in the form of the input parameters.
 * The added temperature and smoke are modelled as emitters in the specified region. The temperature 
 * and smoke density values for these emitters are also user-defined.
 *
 * @param x X-coordinate of the cuboid's starting point.
 * @param y Y-coordinate of the cuboid's starting point.
 * @param z Z-coordinate of the cuboid's starting point.
 * @param wh Width of the cuboid along the X-axis.
 * @param h Height of the cuboid along the Y-axis.
 * @param l Length of the cuboid along the Z-axis.
 * @param dval Smoke density value for the emitter.
 * @param tval Temperature value for the emitter.
 */
void FluidGrid::addSmoke(double x, double y, double z, double wh, double h, double l, double dval, double tval) {

    double x0 = int(x*nwidth);
    double x1 = int((x+wh)*nwidth);
    double y0 = int(y*nheight);
    double y1 = int((y+h)*nheight);
    double z0 = int(z*ndepth);
    double z1 = int((z+l)*ndepth);

    temperature->addEmitter(x0, y0, z0, x1, y1, z1, tval);;
    smoke->addEmitter(x0, y0, z0, x1, y1, z1, dval);

}

/**
 * @brief Updates the time step (`deltaT`) for the simulation based on the fluid's maximum velocity and
 * the Courant-Friedrichs-Lewy (CFL) condition.
 *
 * This method adjusts the simulation's time step based on the current maximum velocity in the fluid grid to
 * ensure numerical stability according to the CFL condition. 
 * If the maximum velocity is zero, the time step is set equal to the frame time step. Otherwise, it is set as
 * the CFL number multiplied by the cell width, divided by the maximum velocity. 
 * Additionally, if the next simulation time (current time + 1.01 * `deltaT`) is greater or equal to the next
 * frame time, the time step is adjusted to match the remaining time to the next frame, and a flag to write
 * the current simulation state to the cache is activated.
 *
 * @note This method will output diagnostic information to the console, including the maximum velocity before projection.
 */
void FluidGrid::setDeltaT(){

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

/**
 * @brief Solves the pressure Poisson equation and updates the velocity field of the fluid grid.
 *
 * This function first solves the pressure Poisson equation Ap = r, where A is the pressure coefficient matrix, 
 * p is the pressure vector, and r is the divergence vector. The solver used is the Preconditioned Conjugate Gradient (PCG) 
 * method, but other solvers such as the Conjugate Gradient or Gauss-Seidel could be used as well.
 * The solution of the equation provides the pressure field that ensures the fluid is incompressible.
 * After the pressure field is calculated, the function calls the updateVelocities() function to update the velocity fields 
 * (u, v, w) based on the newly calculated pressure field. This ensures the fluid motion follows the Navier-Stokes equations.
 */
void FluidGrid::Project(){
    //CGSolver(*A, *r, *pressure, TOL, MAXITER);
    //PCGSolver(*A, *Ei, *r, *pressure, TOL, MAXITER);
    //GaussSeidelSolver(*A, *r, *pressure, TOL, MAXITER);

    // Set up and configure the PCG solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double> > pcgSolver;
    pcgSolver.setTolerance(TOL);  // Set the convergence tolerance
    pcgSolver.setMaxIterations(MAXITER);  // Set the maximum number of iterations

    // Solve the system of equations
    pcgSolver.compute(*A);
    *pressure = pcgSolver.solve(*r);

    updateVelocities();
}

/**
 * @brief Adds force contributions to the vertical velocity of the fluid in the grid.
 * 
 * This method applies the buoyancy force to the vertical velocity field of the fluid grid.
 * The buoyancy force is determined by the distribution of smoke density and temperature 
 * in the grid. Both contribute to the buoyancy force but in opposite directions:
 * 
 * - Smoke increases the buoyancy force (alpha * smoke), leading to a rise in fluid.
 * - Temperature difference from the ambient temperature reduces the buoyancy (beta * (temperature - ambient temperature)).
 * 
 * The total buoyancy force at each point in the vertical velocity field is then scaled by the gravitational acceleration 
 * and the time step size before being added to the vertical velocity.
 * 
 * @note The method only modifies the vertical velocity field `v`.
 */
void FluidGrid::AddForces(){
    float alpha = 1;
    float beta = 1/AMBIENT_TEMPERATURE;

    for (int k = 0; k < v->getzSamples(); k++)
        for (int j = 0; j < v->getySamples(); j++)
            for (int i = 0; i < v->getxSamples(); i++){

                //get position of sample in 'cell space'
                double x = i + v->getOffsetx();
                double y = j + v->getOffsety();
                double z = k + v->getOffsetz();

                v->setQuantity(i,j,k) += getDeltaT() * (alpha*smoke->InterpolateLinear(x,y,z) - beta*(temperature->InterpolateLinear(x,y,z)-AMBIENT_TEMPERATURE)) * GRAVITY;
                //v->setQuantity(i,j,k) += getDeltaT() * GRAVITY;
            }
}

/**
 * @brief Advances the fluid simulation by one time step.
 * 
 * This method performs several operations needed to update the fluid grid for the next time step. Here's the breakdown:
 * 
 * 1. `setDeltaT()`: Computes the time step size for this iteration based on fluid velocities and the 
 *    Courant–Friedrichs–Lewy (CFL) condition.
 * 
 * 2. `Advect()`: Moves the fluid quantities (like velocity and smoke density) around the grid by the velocity field,
 *    a step known as advection.
 * 
 * 3. `AddForces()`: Adds force contributions, such as buoyancy and gravity, to the fluid's velocity.
 * 
 * 4. `addSmoke()`: Emits smoke at a specific location in the grid.
 * 
 * 5. `CalculateVolumes()`: Updates the solid volumes within each cell of the grid.
 * 
 * 6. `BuildLinearSystem()`: Sets up a linear system of equations representing the fluid flow constraints (e.g. incompressibility).
 * 
 * 7. `MIC0precon()`: Applies Modified Incomplete Cholesky preconditioner to the linear system for faster convergence of the solver.
 * 
 * 8. `Project()`: Solves the linear system to ensure the fluid velocity field remains divergence-free (i.e., incompressible).
 * 
 * 9. Updates the current time by adding the time step size.
 * 
 * 10. If it's time to save a frame of the simulation, `WriteToCache()` is called to store the current state of the simulation.
 */
void FluidGrid::Update(){

    setDeltaT();
    Advect();
    AddForces();
    addSmoke(0.25, 0.1, 0.45, 0.1, 0.05, 0.1, 0.5, 450);
    CalculateVolumes();
    BuildLinearSystem();
    Project();

    currtime += deltaT;

    if (writetocache){
        WriteToCache();
        writetocache = false;
    }
}

/**
 * @brief Writes the current state of the smoke simulation to an OpenVDB file.
 * 
 * This method uses the OpenVDB library to store the current state of the smoke
 * simulation as a volumetric grid in a file. 
 * Each cell in the fluid grid is stored as a voxel in the OpenVDB grid.
 * 
 * For each voxel, the value is determined by the smoke density and the temperature.
 * The resultant voxel values are then stored in the file.
 * 
 * The generated file is saved in a specified output path.
 * The frame number is then incremented for the next save operation.
 */
void FluidGrid::WriteToCache(){
    // Create density grid and get its accessor
    openvdb::FloatGrid::Ptr densityGrid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor densityAccessor = densityGrid->getAccessor();

    // Create temperature grid and get its accessor
    openvdb::FloatGrid::Ptr temperatureGrid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor temperatureAccessor = temperatureGrid->getAccessor();

    // Define a coordinate with large signed indices.
    openvdb::Coord ijk;

    int& i = ijk[0];
    int& j = ijk[1];
    int& k = ijk[2];

    for (k=0; k<ndepth; k++)
        for (j=0; j<nheight; j++)
            for (i=0; i<nwidth; i++){
                // densityAccessor.setValue(ijk, smoke->getQuantity(i,j,k) + 1.0 - smoke->getVolume(i,j,k));
                densityAccessor.setValue(ijk, smoke->getQuantity(i,j,k));
                temperatureAccessor.setValue(ijk, temperature->getQuantity(i,j,k));
            }

    char outfile[128];
    sprintf(outfile, "%sFrame%05d.vdb", outputpath.c_str(), framenumber++);

    // Create a VDB file object.
    openvdb::io::File file(outfile);

    densityGrid->setName("density");
    temperatureGrid->setName("temperature");

    // Add the grid pointers to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(densityGrid);
    grids.push_back(temperatureGrid);

    // Write out the contents of the container.
    file.write(grids);
    file.close();
}

/**
 * @brief Performs advection of the fluid quantity in the grid.
 *
 * This method advects a fluid quantity across the fluid grid based on 
 * the grid's current velocity field. The advection is computed using a 
 * second-order Runge-Kutta method. 
 *
 * The method works by updating the fluid quantity at each point in the 
 * grid based on the local velocity and the time step. The method also 
 * checks whether the advection has moved the fluid quantity inside a 
 * solid object in the grid, and if so, it moves the fluid quantity to 
 * the closest point on the surface of the solid object.
 *
 * The new fluid quantity values are stored in a buffer to ensure that 
 * the advection of all points is based on the same initial state, and 
 * are later swapped with the current fluid quantities.
 *
 * @note The function uses Catmull-Rom interpolation for the final 
 * advection step.
 */
void FluidQuantity::Advect(){

    double x,y,z;
    double xmid,ymid,zmid;

    double uvel,vvel,wvel;

    double timestep = parentgrid->getDeltaT();
    double cellwidth = parentgrid->getCellwidth();
    Solid* gridsolid = parentgrid->getSolid();

    FluidQuantity* u = parentgrid->getU();
    FluidQuantity* v = parentgrid->getV();
    FluidQuantity* w = parentgrid->getW();

    for (int k=0; k<zsamples; k++)
        for (int j=0; j<ysamples; j++)
            for (int i=0; i<xsamples; i++){

                x = i+offsetx;
                y = j+offsety;
                z = k+offsetz;

                uvel = u->InterpolateLinear(x,y,z) / cellwidth;
                vvel = v->InterpolateLinear(x,y,z) / cellwidth;
                wvel = w->InterpolateLinear(x,y,z) / cellwidth;

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

                uvel = u->InterpolateLinear(xmid,ymid,zmid) / cellwidth;
                vvel = v->InterpolateLinear(xmid,ymid,zmid) / cellwidth;
                wvel = w->InterpolateLinear(xmid,ymid,zmid) / cellwidth;

                x -= timestep*uvel;
                y -= timestep*vvel;
                z -= timestep*wvel;

                if (gridsolid != nullptr){
                //Are we inside the solid?
                    if (gridsolid->PointIsInside(x, y, z))
                        gridsolid->ClosestSurface(x, y, z);
                }
                
                //setBuffer(i,j,k) = InterpolateLinear(x,y,z);
                setBuffer(i,j,k) = InterpolateCM(x,y,z);

            }  
}

/**
 * @brief Constructs a new FluidQuantity instance.
 *
 * @param parent Pointer to the parent FluidGrid object.
 * @param ofx Offset along the x-dimension.
 * @param ofy Offset along the y-dimension.
 * @param ofz Offset along the z-dimension.
 * @param value Initialise with this value.
 *
 * This constructor initializes a new FluidQuantity object by setting
 * its parent grid and offset values. It determines the number of 
 * samples along each dimension (x, y, z) based on the offset values,
 * and creates and initializes new vectors for quantity, buffer, and volume.
 * 
 * If the offset in a given dimension is zero, the number of samples in
 * that dimension is incremented by 1. This ensures that for non-offset 
 * (grid-aligned) samples, the quantity includes the boundary samples, 
 * thus enabling the handling of boundary conditions.
 * 
 * The quantity and buffer vectors are initialized to value, the volume to zero.
 * Each has a size equal to the product of the number of samples along 
 * each dimension.
 */
FluidQuantity::FluidQuantity(FluidGrid* parent, double ofx, double ofy, double ofz, double value){

    offsetx = ofx;
    offsety = ofy;
    offsetz = ofz;

    parentgrid = parent;

    if (offsetx == 0){
        xsamples = parent->getWidth()+1;
    }
    else{
        xsamples = parent->getWidth();
    }
    if (offsety == 0){
        ysamples = parent->getHeight()+1;
    }
    else{
        ysamples = parent->getHeight();
    }
    if (offsetz == 0){
        zsamples = parent->getDepth()+1;
    }
    else{
        zsamples = parent->getDepth();
    }

    quantity = new Eigen::VectorXd(xsamples*ysamples*zsamples);
    buffer = new Eigen::VectorXd(xsamples*ysamples*zsamples);
    volume = new Eigen::VectorXd(xsamples*ysamples*zsamples);

    quantity->fill(value);
    buffer->fill(value);
    volume->fill(0);
}

/* Destuctor for FluidQuantity */
FluidQuantity::~FluidQuantity(){
    delete quantity;
    delete buffer;
    delete volume;
}

/**
 * Interpolates the quantity value at a given position using trilinear interpolation.
 *
 * This function uses the method of trilinear interpolation to estimate the 
 * quantity value at the specified position (x, y, z). The interpolation takes into 
 * account the relative cell coordinates, and interpolates across the three dimensions.
 * The result provides a smoother representation of the quantity field compared to using 
 * the values at the grid points directly.
 *
 * @param x The x-coordinate of the position where the quantity value is to be interpolated.
 * @param y The y-coordinate of the position where the quantity value is to be interpolated.
 * @param z The z-coordinate of the position where the quantity value is to be interpolated.
 * @return The interpolated value of the quantity at the given position.
 */
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

/**
 * Evaluates the Catmull-Rom spline at a given position.
 *
 * This function uses the Catmull-Rom spline, which is a type of interpolating 
 * spline defined by four points, to estimate the value at the specified position.
 * The Catmull-Rom spline is a cubic polynomial that ensures a smooth transition
 * between the points.
 *
 * @param x The position where the Catmull-Rom spline is to be evaluated.
 * @param f0 The value of the function at the first of the four points.
 * @param f1 The value of the function at the second of the four points.
 * @param f2 The value of the function at the third of the four points.
 * @param f3 The value of the function at the fourth of the four points.
 * @return The evaluated value of the Catmull-Rom spline at the given position.
 */
double FluidQuantity::CatmullRom(double x, double f0, double f1, double f2, double f3){
    double x2 = x*x;
    double x3 = x*x2;

    return f0*(    -0.5*x +     x2 - 0.5*x3)
    + f1*(1          - 2.5*x2 + 1.5*x3)
    + f2*(     0.5*x +   2*x2 - 1.5*x3)
    + f3*(            -0.5*x2 + 0.5*x3);
}

/**
 * Performs 3D interpolation of a fluid quantity using Catmull-Rom splines.
 *
 * This function first transforms the z-coordinate into the correct dimensionless cell space for the given offset. 
 * It then clamps the z-coordinate to ensure it's within the grid. 
 * The method subsequently calls the `InterpolateCMSlice` method on four slices of the grid (determined by the z-coordinate)
 * to get interpolated values in 2D. Finally, it performs another Catmull-Rom interpolation along the z-direction using these 
 * interpolated 2D values to obtain the final interpolated 3D value.
 *
 * @param x The x-coordinate of the position to interpolate.
 * @param y The y-coordinate of the position to interpolate.
 * @param z The z-coordinate of the position to interpolate.
 * @return The interpolated value at the given position (x, y, z).
 */
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

/**
 * Uses the Catmull-Rom spline to perform a 2D interpolation on a specified slice of the grid.
 *
 * This method first transforms the input coordinates into the dimensionless cell space 
 * according to the offsets. The input coordinates are then clamped to the grid. 
 * The Catmull-Rom interpolation is then performed in two stages: 
 * first, it's performed along the x-direction for four rows, 
 * and then the results of these interpolations are used for another 
 * Catmull-Rom interpolation along the y-direction.
 *
 * @param x The x-coordinate of the position to interpolate.
 * @param y The y-coordinate of the position to interpolate.
 * @param k The z-index of the slice to interpolate on.
 * @return The interpolated value at the given position in the specified slice.
 */
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

/**
 * Sets the FluidQuantity value within a specified rectangular volume in the grid.
 *
 * This function assigns a new value to each point inside a given volume specified by 
 * the rectangular region from (x0, y0, z0) to (x1, y1, z1). The value is only assigned 
 * if the absolute value of the current fluid quantity at the grid point is less than the 
 * absolute value of 'v'. The function also ensures that the rectangular region lies within 
 * the valid boundaries of the fluid grid.
 *
 * @param x0 The starting x-coordinate of the rectangular volume.
 * @param y0 The starting y-coordinate of the rectangular volume.
 * @param z0 The starting z-coordinate of the rectangular volume.
 * @param x1 The ending x-coordinate of the rectangular volume.
 * @param y1 The ending y-coordinate of the rectangular volume.
 * @param z1 The ending z-coordinate of the rectangular volume.
 * @param value  The fluid quantity value to assign within the given volume.
 */
void FluidQuantity::addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double value) {
    for (int z = std::max(z0, 0); z < std::min(z1, zsamples); z++)
        for (int y = std::max(y0, 0); y < std::min(y1, ysamples); y++)
            for (int x = std::max(x0, 0); x < std::min(x1, xsamples); x++)
                if (fabs(getQuantity(x,y,z)) < fabs(value))
                    setQuantity(x,y,z) = value;
}

/**
 * @brief Extrapolates the velocity of a fluid to cells neighbouring fluid cells.
 * 
 * This method extrapolates the velocity of a fluid to cells that were not originally 
 * marked as fluid cells, but neighbour cells containing fluid.
 * 
 * @details
 * The function iterates through the 3D grid of the fluid and determines for each cell whether 
 * it is a valid cell. A cell is valid if it's not entirely occupied by solid objects. 
 * 
 * The extrapolation is performed layer by layer. At each layer, the function checks the neighbouring 
 * cells of each invalid cell and if any of them are valid, it assigns the invalid cell the 
 * average quantity of its valid neighbours.
 * 
 * The process is repeated for a predefined number of layers, spreading the fluid's velocity 
 * outward, into regions that were originally non-fluid.
 * 
 * Note: The function assumes that there is a solid boundary at the edges of the grid.
 */
void FluidQuantity::ExtrapolateVelocity(){
    bool* valid = new bool[xsamples*ysamples*zsamples];
    bool* old_valid = new bool[xsamples*ysamples*zsamples];

    //Initialize the list of valid cells if weight is non-zero then cell is valid
    //non-zero volume indicates the presence of quantity in that cell
    for(int k = 0; k < zsamples; k++)
        for(int j = 0; j < ysamples; j++)
            for(int i = 0; i < xsamples; i++)
                valid[i + j*xsamples + k*xsamples*ysamples] = getVolume(i,j,k) > 0;

    for(int layers = 0; layers < 5; ++layers) {
        for(int i =0; i<xsamples*ysamples*zsamples; i++)
            old_valid[i] = valid[i];

        //Don't iterate over the edges of the grid
        for(int k = 1; k < zsamples-1; k++)
            for(int j = 1; j < ysamples-1; j++)
                for(int i = 1; i < xsamples-1; i++){
                    double sum = 0;
                    int count = 0;
                    //if the cell is not valid
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
    delete[] valid;
    delete[] old_valid;
};

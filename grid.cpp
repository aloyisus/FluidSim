#include "grid.h"


/**
 * @brief Constructs a new Cuboid instance.
 *
 * This constructor initializes a new Cuboid instance with the specified position, dimensions,
 * and other parameters. Note that the position and dimensions are given in terms of
 * world-spatial coordinates.
 *
 * @param parent Pointer to the FluidGrid that this cuboid is a part of.
 * @param centrex X-coordinate of the cuboid's center.
 * @param centrey Y-coordinate of the cuboid's center.
 * @param centrez Z-coordinate of the cuboid's center.
 * @param width Width of the cuboid.
 * @param height Height of the cuboid.
 * @param depth Depth of the cuboid.
 * @param alph Initial rotation of cuboid in radians.
 * @param omg Angular velocity of cuboid in radians/s.
 */
Cuboid::Cuboid(FluidGrid* parent, double centrex, double centrey, double centrez, double width, double height, double depth, double alph, double omg) {
    
    parentgrid = parent;
    auto bbox = parent->getBBox();

    openvdb::Vec3d size = bbox.max() - bbox.min();

    centrepos(0) = centrex * size[0];
    centrepos(1) = centrey * size[1];
    centrepos(2) = centrez * size[2];

    w = width * size[0];
    h = height * size[1];
    d = depth * size[2];
    
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
    
    openvdb::Vec3d cellpos(x,y,z);
    
    double t = parentgrid->getCurrtime();
    
    bool temp;
    openvdb::Mat3d invtransform;
    
    //rotation of -theta about z-axis
    invtransform(2,0) = 0;
    invtransform(0,2) = 0;
    invtransform(2,1) = 0;
    invtransform(1,2) = 0;
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
openvdb::Vec3d Cuboid::getVelocity(openvdb::Vec3d xyz){

    //vector from axis of rotation (z-axis) to x,y,z
    double rx = xyz[0] - centrepos(0);
    double ry = xyz[1] - centrepos(1);

    double distance = sqrt(rx*rx + ry*ry);

    //get the unit normal to vector (rx,ry)
    double normx = -ry / distance;
    double normy = rx / distance;

    //get magnitude of velocity at x,y, scale the norm vector by this
    normx *= distance*omega;
    normy *= distance*omega;

    openvdb::Vec3d temp(normx, normy, 0);

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
openvdb::Vec3d Cuboid::getNormal(openvdb::Vec3d xyz){

    openvdb::Vec3d temp;
    openvdb::Mat3d transform;

    double t = parentgrid->getCurrtime();

    //rotation of -theta about z-axis
    transform(2,0) = 0;
    transform(0,2) = 0;
    transform(2,1) = 0;
    transform(1,2) = 0;
    transform(0,0) = cos(alpha + omega*t);
    transform(1,0) = -sin(alpha + omega*t);
    transform(0,1) = sin(alpha + omega*t);
    transform(1,1) = cos(alpha + omega*t);
    transform(2,2) = 1.0;

    xyz = transform*(xyz - centrepos);

    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];

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
    transform(2,0) = 0;
    transform(0,2) = 0;
    transform(2,1) = 0;
    transform(1,2) = 0;
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

    openvdb::Mat3d transform;
    openvdb::Vec3d r(x,y,z);

    double t = parentgrid->getCurrtime();

    //rotation of -theta about z-axis
    transform(2,0) = 0;
    transform(0,2) = 0;
    transform(2,1) = 0;
    transform(1,2) = 0;
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
    transform(2,0) = 0;
    transform(0,2) = 0;
    transform(2,1) = 0;
    transform(1,2) = 0;
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

    openvdb::Vec3d cellpos(x,y,z);

    double sqrdistance = (centrepos - cellpos).dot(centrepos - cellpos);

    return (sqrdistance <= radius*radius);
}

/**
 * @brief Retrieves the velocity of the sphere at the specified point.
 *
 * As the sphere has no velocity component normal to its surface, this function always returns zero
 * regardless of the input position (x, y, z). This is represented by an openvdb::Vec3d containing three zeroes.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return An openvdb::Vec3d representing the velocity of the sphere at the point, which is always (0,0,0).
 */
openvdb::Vec3d Sphere::getVelocity(openvdb::Vec3d xyz){
    openvdb::Vec3d temp(0,0,0);
    return temp;
}

/**
 * @brief Retrieves the velocity of the sphere at the specified point.
 *
 * As the sphere has no velocity component normal to its surface, this function always returns zero
 * regardless of the input position (x, y, z). This is represented by an openvdb::Vec3d containing three zeroes.
 * 
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @return An openvdb::Vec3d representing the velocity of the sphere at the point, which is always (0,0,0).
 */
openvdb::Vec3d Sphere::getNormal(openvdb::Vec3d xyz){

    xyz = (xyz - centrepos);
    double mod = sqrt(xyz.dot(xyz));

    if (mod == 0)
        return xyz;

    xyz = (1/mod)*xyz;
    return xyz;
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

    openvdb::Vec3d r;
    openvdb::Vec3d norm = getNormal(openvdb::Vec3d(x,y,z));

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
    // make cellwidth be 1 / no. of cells in shortest dimension
    // this results in a volume of unit length when measured along its shortest side
    cellwidth = 1.0/std::min(std::min(wh, ht),dh);

    nwidth = wh;
    nheight = ht;
    ndepth = dh;
    ncells = wh*ht*dh;

    openvdb::Vec3d min(0, 0, 0);
    openvdb::Vec3d max(nwidth * cellwidth, nheight * cellwidth, ndepth * cellwidth);
    bbox = openvdb::BBoxd(min, max);

    density = rh;

    CFLnumber = 5;
    framenumber = 0;
    currtime = 0;
    framedeltaT = tstep;
    nextframetime = tstep;
    bool writetocache = false;

    outputpath = filepath;

    r = new SymmBandMatrix::VectorType(wh*ht*dh);
    pressure = new SymmBandMatrix::VectorType(wh*ht*dh);
    A = new SymmBandMatrix(wh*ht*dh);

    u = new FluidQuantity(this, 0,0.5,0.5);
    v = new FluidQuantity(this, 0.5,0,0.5);
    w = new FluidQuantity(this, 0.5,0.5,0);
    smoke = new FluidQuantity(this, 0.5,0.5,0.5);
    temperature = new FluidQuantity(this, 0.5,0.5,0.5,AMBIENT_TEMPERATURE);
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

// MODIFY THIS ONE.................................JK
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
    r->fill(0);
    A->scale(0);

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
                
                // HACK! Set zero diagonal elements to be very small instead
                // this is an attempt to pass a validation test in openvdb's PCG solver
                if (!A->getValue(row, row))
                    A->setValue(row, row, 1e-5);

                (*r)[row] = -rscale*(u->getVolume(i+1,j,k)*u->getQuantity(i+1,j,k) - u->getVolume(i,j,k)*u->getQuantity(i,j,k) +
                                      v->getVolume(i,j+1,k)*v->getQuantity(i,j+1,k) - v->getVolume(i,j,k)*v->getQuantity(i,j,k) +
                                      w->getVolume(i,j,k+1)*w->getQuantity(i,j,k+1) - w->getVolume(i,j,k)*w->getQuantity(i,j,k));

                if (solid != nullptr){
                    // This is derived from figure 4.3 on page 49 of Bridson's book
                    (*r)[row] += rscale*((u->getVolume(i+1,j,k) - smoke->getVolume(i,j,k))*
                                          solid->getVelocity(u->indexToWorld(i+1,j,k))(0) -
                                         (u->getVolume(i,j,k) - smoke->getVolume(i,j,k))*
                                          solid->getVelocity(u->indexToWorld(i,j,k))(0) +

                                         (v->getVolume(i,j+1,k) - smoke->getVolume(i,j,k))*
                                          solid->getVelocity(v->indexToWorld(i,j+1,k))(1) -
                                         (v->getVolume(i,j,k) - smoke->getVolume(i,j,k))*
                                          solid->getVelocity(v->indexToWorld(i,j,k))(1) +

                                         (w->getVolume(i,j,k+1) - smoke->getVolume(i,j,k))*
                                          solid->getVelocity(w->indexToWorld(i,j,k+1))(2) -
                                         (w->getVolume(i,j,k) - smoke->getVolume(i,j,k))*
                                          solid->getVelocity(w->indexToWorld(i,j,k))(2));
                }
    };

    // DEBUG!!!!!!!
    // Test if A is valid (the 'apply' method raises if a failure bool
    // which was set in CholeskyPrecondMatrix is false)
    // CholeskyPrecondMatrix precond(*A);
    // precond.apply(*r, *r);

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

                setVolume(i,j,k, 8.0);

                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix0+iy0*2*nwidth+iz0*4*nwidth*nheight]);
                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix0+iy1*2*nwidth+iz0*4*nwidth*nheight]);
                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix0+iy0*2*nwidth+iz1*4*nwidth*nheight]);
                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix0+iy1*2*nwidth+iz1*4*nwidth*nheight]);
                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix1+iy0*2*nwidth+iz0*4*nwidth*nheight]);
                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix1+iy1*2*nwidth+iz0*4*nwidth*nheight]);
                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix1+iy0*2*nwidth+iz1*4*nwidth*nheight]);
                setVolume(i,j,k, getVolume(i,j,k) - supervol[ix1+iy1*2*nwidth+iz1*4*nwidth*nheight]);

                setVolume(i,j,k, getVolume(i,j,k)/8.0);

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
    double* supervol = new double[ncells*8];
    double pointisinside;

    //Perform 2x2x2 supersampling over entire grid, so each ijk cell will be sampled 8 times
    for (int k=0; k<2*ndepth; k++)
        for (int j=0; j<2*nheight; j++)
            for (int i=0; i<2*nwidth; i++){
                pointisinside = 0;
                if (solid != nullptr)
                    if (solid->PointIsInside((0.25+i*0.5)*cellwidth, (0.25+j*0.5)*cellwidth, (0.25+k*0.5)*cellwidth))
                        pointisinside = 1;
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
                    u->setQuantity(i,j,k, u->getQuantity(i,j,k) - scale*(getPressure(i,j,k) - getPressure(i-1,j,k)));
                else
                    u->setQuantity(i,j,k, 0);
            }
    
    for (int k = 0; k < v->getzSamples(); k++)
        for (int j = 1; j < v->getySamples()-1; j++)
            for (int i = 0; i < v->getxSamples(); i++) {
                
                if (v->getVolume(i,j,k) > 0)
                    v->setQuantity(i,j,k, v->getQuantity(i,j,k) - scale*(getPressure(i,j,k) - getPressure(i,j-1,k)));
                else
                    v->setQuantity(i,j,k, 0);
            }
    
    for (int k = 1; k < w->getzSamples()-1; k++)
        for (int j = 0; j < w->getySamples(); j++)
            for (int i = 0; i < w->getxSamples(); i++) {
                
                if (w->getVolume(i,j,k) > 0)
                    w->setQuantity(i,j,k, w->getQuantity(i,j,k) - scale*(getPressure(i,j,k) - getPressure(i,j,k-1)));
                else
                    w->setQuantity(i,j,k, 0);
            }

    // temperature->ExtrapolateVelocity();
    // smoke->ExtrapolateVelocity();
    u->ExtrapolateVelocity();
    v->ExtrapolateVelocity();
    w->ExtrapolateVelocity();

    if (solid != nullptr)
        SetSolidBoundaries();
    SetWallBoundaries();
}

void FluidGrid::SetWallBoundaries(){
    //Boundary conditions for the walls of the box
    //Set the normal component to the wall to zero
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++){
            u->setQuantity(0,j,k, 0.0);
            u->setQuantity(nwidth,j,k, 0.0);
        }
    for (int k = 0; k < ndepth; k++)
        for (int i = 0; i < nwidth; i++){
            v->setQuantity(i,0,k, 0.0);
            v->setQuantity(i,nheight,k, 0.0);
        }
    for (int j = 0; j < nheight; j++)
        for (int i = 0; i < nwidth; i++){
            w->setQuantity(i,j,0, 0.0);
            w->setQuantity(i,j,ndepth, 0.0);
        }
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
openvdb::Vec3d FluidGrid::getVelocity(openvdb::Vec3d xyz){

    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];
    openvdb::Vec3d temp;

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

    openvdb::Vec3d vfluid;
    openvdb::Vec3d vsolid;
    openvdb::Vec3d normal;
    
    //For cells which have zero volume, set the velocity samples bordering to have zero relative velocity in the
    //direction of the normal to the solid boundary
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++)
            for (int i = 0; i < nwidth; i++) {
                if (smoke->getVolume(i,j,k) == 0)
                {
                    vsolid = solid->getVelocity(u->indexToWorld(i,j,k));
                    vfluid = getVelocity(u->indexToWorld(i,j,k));
                    normal = solid->getNormal(u->indexToWorld(i,j,k));
                    
                    u->setBuffer(i,j,k, (vfluid - (vfluid - vsolid).dot(normal)*normal)(0));
                    
                    vsolid = solid->getVelocity(u->indexToWorld(i+1,j,k));
                    vfluid = getVelocity(u->indexToWorld(i+1,j,k));
                    normal = solid->getNormal(u->indexToWorld(i+1,j,k));
                    
                    u->setBuffer(i+1,j,k, (vfluid - (vfluid - vsolid).dot(normal)*normal)(0));
                    
                    vsolid = solid->getVelocity(v->indexToWorld(i,j,k));
                    vfluid = getVelocity(v->indexToWorld(i,j,k));
                    normal = solid->getNormal(v->indexToWorld(i,j,k));
                    
                    v->setBuffer(i,j,k, (vfluid - (vfluid - vsolid).dot(normal)*normal)(1));
                    
                    vsolid = solid->getVelocity(v->indexToWorld(i,j+1,k));
                    vfluid = getVelocity(v->indexToWorld(i,j+1,k));
                    normal = solid->getNormal(v->indexToWorld(i,j+1,k));
                    
                    v->setBuffer(i,j+1,k, (vfluid - (vfluid - vsolid).dot(normal)*normal)(1));
                    
                    vsolid = solid->getVelocity(w->indexToWorld(i,j,k));
                    vfluid = getVelocity(w->indexToWorld(i,j,k));
                    normal = solid->getNormal(w->indexToWorld(i,j,k));
                    
                    w->setBuffer(i,j,k, (vfluid - (vfluid - vsolid).dot(normal)*normal)(2));
                    
                    vsolid = solid->getVelocity(w->indexToWorld(i,j,k+1));
                    vfluid = getVelocity(w->indexToWorld(i,j,k+1));
                    normal = solid->getNormal(w->indexToWorld(i,j,k+1));
                    
                    w->setBuffer(i,j,k+1, (vfluid - (vfluid - vsolid).dot(normal)*normal)(2));
                }
            }
    u->swap();
    v->swap();
    w->swap();
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
 * @param density_value Smoke density value for the emitter.
 * @param temperature_value Temperature value for the emitter.
 */
void FluidGrid::addSmoke(double x, double y, double z, double wh, double h, double l, double density_value, double temperature_value) {

    double x0 = int(x*nwidth);
    double x1 = int((x+wh)*nwidth);
    double y0 = int(y*nheight);
    double y1 = int((y+h)*nheight);
    double z0 = int(z*ndepth);
    double z1 = int((z+l)*ndepth);

    temperature->addEmitter(x0, y0, z0, x1, y1, z1, temperature_value);;
    smoke->addEmitter(x0, y0, z0, x1, y1, z1, density_value);

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
                    setQuantity(x,y,z, value);
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

    auto xform = getTransform();
    for (int k=0; k<zsamples; k++)
        for (int j=0; j<ysamples; j++)
            for (int i=0; i<xsamples; i++){

                // use index to world here for x,y,z
                openvdb::Vec3d xyz = xform->indexToWorld(openvdb::Coord(i,j,k));
                x = xyz[0];
                y = xyz[1];
                z = xyz[2];

                uvel = u->InterpolateLinear(x,y,z);
                vvel = v->InterpolateLinear(x,y,z);
                wvel = w->InterpolateLinear(x,y,z);

                //Runge-Kutta 2
                xmid = x - 0.5*timestep*uvel;
                ymid = y - 0.5*timestep*vvel;
                zmid = z - 0.5*timestep*wvel;

                uvel = u->InterpolateLinear(xmid,ymid,zmid);
                vvel = v->InterpolateLinear(xmid,ymid,zmid);
                wvel = w->InterpolateLinear(xmid,ymid,zmid);

                x -= timestep*uvel;
                y -= timestep*vvel;
                z -= timestep*wvel;

                if (gridsolid != nullptr){
                //Are we inside the solid?
                    if (gridsolid->PointIsInside(x, y, z))
                        gridsolid->ClosestSurface(x, y, z);
                }

                setBuffer(i,j,k, InterpolateCubic(x,y,z));
    }  
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
 * @brief Solves the pressure Poisson equation and updates the velocity field of the fluid grid.
 *
 * This function first solves the pressure Poisson equation Ap = r, where A is the pressure coefficient matrix, 
 * p is the pressure vector, and r is the divergence vector. The solver used is the Preconditioned Conjugate Gradient (PCG) 
 * method.
 * The solution of the equation provides the pressure field that ensures the fluid is incompressible.
 * After the pressure field is calculated, the function calls the updateVelocities() function to update the velocity fields 
 * (u, v, w) based on the newly calculated pressure field. This ensures the fluid motion follows the Navier-Stokes equations.
 */
void FluidGrid::Project(){

    //DEBUG
    auto pressure_before(*pressure);

    CholeskyPrecondMatrix precond(*A);

    auto state = openvdb::math::pcg::terminationDefaults<double>();
    state.iterations = MAXITER; // Maximum number of iterations
    state.relativeError = TOL; // Tolerance for relative error

    openvdb::util::NullInterrupter interrupter;
    std::cout << "Solving..." << std::endl;
    auto result = openvdb::math::pcg::solve(*A, *r, *pressure, precond, interrupter, state);

    if (!result.success)
        std::cout << "WARNING: Solver failed!" << std::endl;

    if (pressure->eq(pressure_before))
        std::cout << "No change to pressure" << std::endl;

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
    double x,y,z;

    auto xform = v->getTransform();
    for (int k = 0; k < v->getzSamples(); k++)
        for (int j = 0; j < v->getySamples(); j++)
            for (int i = 0; i < v->getxSamples(); i++){

                // use index to world here for x,y,z
                openvdb::Vec3d xyz = xform->indexToWorld(openvdb::Coord(i,j,k));
                x = xyz[0];
                y = xyz[1];
                z = xyz[2];
    
                v->setQuantity(i,j,k, v->getQuantity(i,j,k) + getDeltaT() * (alpha*smoke->InterpolateLinear(x,y,z) - beta*(temperature->InterpolateLinear(x,y,z)-AMBIENT_TEMPERATURE)) * GRAVITY);
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

    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();
    setDeltaT();
    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "setDeltaT time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    Advect();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "Advect time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    AddForces();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "addForces time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    addSmoke(0.25, 0.1, 0.45, 0.1, 0.05, 0.1, 0.5, 450);
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "addSmoke time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    CalculateVolumes();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "CalculateVolumes time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    BuildLinearSystem();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "BuildLinearSystem time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    Project();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "Project time (seconds) = " << elapsed.count() * 1e-9 << std::endl;

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
                // The following line will make the solid visible by filling solid-occupied cells
                densityAccessor.setValue(ijk, smoke->getQuantity(i,j,k) + 1.0 - smoke->getVolume(i,j,k));

                // uncommenting this line will make the solid invisible
                // densityAccessor.setValue(ijk, smoke->getQuantity(i,j,k));
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
FluidQuantity::FluidQuantity(FluidGrid* parent, double ofx, double ofy, double ofz, double value)
{
    offsetx = ofx;
    offsety = ofy;
    offsetz = ofz;

    parentgrid = parent;

    if (offsetx == 0){
        xsamples = parent->getWidth() + 1;
    }
    else{
        xsamples = parent->getWidth();
    }
    if (offsety == 0){
        ysamples = parent->getHeight() + 1;
    }
    else{
        ysamples = parent->getHeight();
    }
    if (offsetz == 0){
        zsamples = parent->getDepth() + 1;
    }
    else{
        zsamples = parent->getDepth();
    }
    numsamples = xsamples*ysamples*zsamples;

    // Create OpenVDB grids with given dimensions
    quantity = openvdb::FloatGrid::create(value);
    buffer = openvdb::FloatGrid::create(value);
    volume = openvdb::FloatGrid::create(0);

    q_access = new openvdb::FloatGrid::Accessor(quantity->getAccessor());
    b_access = new openvdb::FloatGrid::Accessor(buffer->getAccessor());
    v_access = new openvdb::FloatGrid::Accessor(volume->getAccessor());

    // Adjust grid transform to match desired grid spacing and offset
    openvdb::math::Transform::Ptr q_transform = openvdb::math::Transform::createLinearTransform(parent->getCellwidth());
    q_transform->preTranslate(openvdb::Vec3d(offsetx, offsety, offsetz));
    quantity->setTransform(q_transform);
    buffer->setTransform(q_transform);
    volume->setTransform(q_transform);
}

/* Destuctor for FluidQuantity */
FluidQuantity::~FluidQuantity() {
    // No explicit cleanup required for OpenVDB grids
    delete q_access;
    delete b_access;
    delete v_access;
}

openvdb::Vec3d FluidQuantity::indexToWorld(int i, int j, int k){
    auto xform = getTransform();
    openvdb::Vec3d xyz = xform->indexToWorld(openvdb::Coord(i, j, k));
    return xyz;
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
double FluidQuantity::InterpolateLinear(double x, double y, double z) const
{
    // Create a linear interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> interpolator(*quantity);

    auto xyz = openvdb::Vec3d(x, y, z);
    auto bbox = parentgrid->getBBox();

    // Clamp position to the bbox
    for (int i = 0; i < 3; i++) {
        xyz[i] = std::max(bbox.min()[i], std::min(xyz[i], bbox.max()[i]));
    }

    // Interpolate the value at the given position
    double result = interpolator.wsSample(xyz);

    return result;
}

double FluidQuantity::InterpolateCubic(double x, double y, double z) const
{
    // Create a cubic interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::QuadraticSampler> interpolator(*quantity);

    auto xyz = openvdb::Vec3d(x, y, z);
    auto bbox = parentgrid->getBBox();

    // Clamp position to the bbox
    for (int i = 0; i < 3; i++) {
        xyz[i] = std::max(bbox.min()[i], std::min(xyz[i], bbox.max()[i]));
    }

    // Interpolate the value at the given position
    double result = interpolator.wsSample(xyz);

    return result;
}

double FluidQuantity::sum() const {
    double total = 0;
    int count = 0;
    for (openvdb::FloatGrid::ValueOnCIter iter = quantity->cbeginValueOn(); iter; ++iter) {
        total += *iter;
        count += 1;
    }
    total += (numsamples - count)*quantity->background(); 
    return total;
}

double FluidQuantity::max() const {
    double maxValue = 0;
    for (auto iter = quantity->cbeginValueOn(); iter; ++iter) {
        if (*iter > maxValue) maxValue = *iter;
    }
    return fmax(maxValue, quantity->background());
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
                            setQuantity(i,j,k, sum /(double)count);
                            valid[i + j*xsamples + k*xsamples*ysamples] = true;
                        }
                    }
        }
    }
    delete[] valid;
    delete[] old_valid;
};

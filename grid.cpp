#include <iostream>
#include "grid.h"
#include "maths.h"
#include "interp.h"



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


void FluidGrid::setDeltaT(){

    double maxvel = std::max(w->max(),std::max(u->max(),v->max()));
    std::cout << "max velocity before project = " << maxvel << std::endl;

    if (maxvel==0)
        deltaT = framedeltaT;
    else
        deltaT = CFLnumber*cellwidth / maxvel;

    if (currtime + 1.01*deltaT >= nextframetime){
        deltaT = nextframetime - currtime;
        writetocache = true;
        nextframetime += framedeltaT;
    }

    std::cout << "deltaT=" <<  deltaT << "currtime=" <<  currtime << "writetocache=" << writetocache << "framenumber=" << framenumber << std::endl;
}


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


void FluidQuantity::Advect(){

    openvdb::Vec3d xyz, xyzmid, velocity;

    double timestep = parentgrid->getDeltaT();
    double cellwidth = parentgrid->getCellwidth();

    FluidQuantity* u = parentgrid->getU();
    FluidQuantity* v = parentgrid->getV();
    FluidQuantity* w = parentgrid->getW();

    auto xform = getTransform();
    for (int k=0; k<zsamples; k++)
        for (int j=0; j<ysamples; j++)
            for (int i=0; i<xsamples; i++){

                // use index to world here for x,y,z
                xyz = xform->indexToWorld(openvdb::Coord(i,j,k));
                velocity = openvdb::Vec3d(
                    u->InterpolateLinear(xyz),
                    v->InterpolateLinear(xyz),
                    w->InterpolateLinear(xyz)
                );

                //Runge-Kutta 2
                xyzmid = xyz - 0.5*timestep*velocity;

                velocity = openvdb::Vec3d(
                    u->InterpolateLinear(xyzmid),
                    v->InterpolateLinear(xyzmid),
                    w->InterpolateLinear(xyzmid)
                );

                xyz -= timestep*velocity;

                // double sample = InterpolateLinear(xyz);
                // double sample = InterpolateQuadratic(xyz);
                double sample = InterpolateCubic(xyz);

                setBuffer(i,j,k, sample);  
    }
}


void FluidGrid::AddForces(){
    float alpha = 1;
    float beta = 1/AMBIENT_TEMPERATURE;
    openvdb::Vec3d xyz;

    auto xform = v->getTransform();
    for (int k = 0; k < v->getzSamples(); k++)
        for (int j = 1; j < v->getySamples()-1; j++)
            for (int i = 0; i < v->getxSamples(); i++){

                // use index to world here for x,y,z
                xyz = xform->indexToWorld(openvdb::Coord(i,j,k));
                v->setQuantity(i,j,k, v->getQuantity(i,j,k) + getDeltaT() * (alpha*smoke->InterpolateLinear(xyz) - beta*(temperature->InterpolateLinear(xyz)-AMBIENT_TEMPERATURE)) * GRAVITY);
            }
}


void FluidGrid::addSmoke(double x, double y, double z, double wh, double h, double l, double density_value, double temperature_value) {

    double x0 = int((x-wh/2)*nwidth);
    double x1 = int((x+wh/2)*nwidth);
    double y0 = int(y*nheight);
    double y1 = int((y+h)*nheight);
    double z0 = int((z-l/2)*ndepth);
    double z1 = int((z+l/2)*ndepth);

    temperature->addEmitter(x0, y0, z0, x1, y1, z1, temperature_value);;
    smoke->addEmitter(x0, y0, z0, x1, y1, z1, density_value);
}


void FluidQuantity::addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, double value) {
    for (int z = std::max(z0, 0); z < std::min(z1, zsamples); z++)
        for (int y = std::max(y0, 0); y < std::min(y1, ysamples); y++)
            for (int x = std::max(x0, 0); x < std::min(x1, xsamples); x++)
                setQuantity(x,y,z, value);
}


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


void FluidGrid::CalculateVolumes(){
    double* supervol = new double[ncells*8];
    double pointisinside;

    //Perform 2x2x2 supersampling over entire grid, so each ijk cell will be sampled 8 times
    for (int k=0; k<2*ndepth; k++)
        for (int j=0; j<2*nheight; j++)
            for (int i=0; i<2*nwidth; i++){
                pointisinside = 0;
                // if (solid != nullptr)
                //     if (solid->PointIsInside((0.25+i*0.5)*cellwidth, (0.25+j*0.5)*cellwidth, (0.25+k*0.5)*cellwidth))
                //         pointisinside = 1;
                supervol[i+j*2*nwidth+k*4*nwidth*nheight] = pointisinside;
            }
    u->CalculateVolumes(supervol);
    v->CalculateVolumes(supervol);
    w->CalculateVolumes(supervol);
    smoke->CalculateVolumes(supervol);

    delete[] supervol;
}


void FluidGrid::BuildLinearSystem(){
    double Ascale =  deltaT/(density*cellwidth*cellwidth);
    double rscale = 1.0/cellwidth;

    clearSparseMatrix(A);
    // pressure->fill(0);
    r->fill(0);

    //Setup up the matrix of pressure coefficients, weighted by volumes
    int ix = 0;
    int slice_size  = nwidth * nheight;
    int total_cells = nwidth * nheight * ndepth;
    for(int k=0; k<ndepth; k++)
        for(int j=0; j<nheight; j++)
            for(int i=0; i<nwidth; i++){

                ix = i + (j * nwidth) + (k * slice_size);
                addToA(ix, 0, Ascale*u->getVolume(i,j,k));
                addToA(ix, 0, Ascale*v->getVolume(i,j,k));
                addToA(ix, 0, Ascale*w->getVolume(i,j,k));
                addToA(ix, 0, Ascale*u->getVolume(i+1,j,k));
                addToA(ix, 0, Ascale*v->getVolume(i,j+1,k));
                addToA(ix, 0, Ascale*w->getVolume(i,j,k+1));

                // Applies to all cells EXCEPT the very last cell in the grid
                if (ix < ncells - 1) {
                    addToA(ix, 1, -Ascale*u->getVolume(i+1,j,k));
                }
                // Applies to all cells EXCEPT the last column of the last depth slice
                if (ix < ncells - nwidth) {
                    addToA(ix, 2, -Ascale*v->getVolume(i,j+1,k));
                }
                // Applies to all cells EXCEPT the last depth slice
                if (ix < ncells - slice_size) {
                    addToA(ix, 3, -Ascale*w->getVolume(i,j,k+1));
                }

                (*r)[ix] = -rscale*(u->getVolume(i+1,j,k)*u->getQuantity(i+1,j,k) - u->getVolume(i,j,k)*u->getQuantity(i,j,k) +
                                     v->getVolume(i,j+1,k)*v->getQuantity(i,j+1,k) - v->getVolume(i,j,k)*v->getQuantity(i,j,k) +
                                     w->getVolume(i,j,k+1)*w->getQuantity(i,j,k+1) - w->getVolume(i,j,k)*w->getQuantity(i,j,k));

    };

    if (!isMatrixSymmetric(*A, 1e-6)){
        std::cout << "A is not symmetric, this is a problem" << std::endl;
    }
    if (!rowsSumToZero(*A, 1e-6)){
        std::cout << "rows not summing to zero, or else no entries in row" << std::endl;
    }
}


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


template <typename SparseMatrixType>
void clearSparseMatrix(SparseMatrixType* matrix) {
    // decltype ensures we use the exact integer type the matrix expects (SizeType/Index32)
    // without needing to know its namespace.
    auto numRows = matrix->numRows();
    
    for (decltype(numRows) i = 0; i < numRows; ++i) {
        matrix->getRowEditor(i).clear();
    }
}


void FluidGrid::Project(){

    //DEBUG
    auto pressure_before(*pressure);

    CholeskyPrecondMatrix precond(*A);
    // openvdb::math::pcg::JacobiPreconditioner<SymmBandMatrix> precond(*A);

    auto state = openvdb::math::pcg::terminationDefaults<double>();
    state.iterations = MAXITER; // Maximum number of iterations
    state.relativeError = TOL; // Tolerance for relative error

    std::cout << "max divergence before solve = " <<     maxDivergence() << std::endl;

    openvdb::util::NullInterrupter interrupter;
    std::cout << "Solving..." << std::endl;
    auto result = openvdb::math::pcg::solve(*A, *r, *pressure, precond, interrupter, state);

    // bool result = CGSolver(*A, *r, *pressure, TOL, MAXITER);

    if (!result.success)
        std::cout << "WARNING: Solver failed!" << std::endl;
    // if (!result)
    //     std::cout << "WARNING: Solver failed!" << std::endl;

    if (pressure->eq(pressure_before))
        std::cout << "No change to pressure" << std::endl;

    updateVelocities();
    std::cout << "max divergence after solve = " <<     maxDivergence() << std::endl;

}


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
                
                if (v->getVolume(i,j,k) > 0){
                    v->setQuantity(i,j,k, v->getQuantity(i,j,k) - scale*(getPressure(i,j,k) - getPressure(i,j-1,k)));
                }
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

    auto begin = std::chrono::high_resolution_clock::now();
    setDeltaT();
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "setDeltaT time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    Advect();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "Advect time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    addSmoke(0.5, 0.1, 0.5, 0.25, 0.05, 0.25, 0.5, 450);
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "addSmoke time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    AddForces();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "addForces time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
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


void FluidGrid::WriteToCache(){
    // Create density grid and get its accessor
    openvdb::FloatGrid::Ptr densityGrid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor densityAccessor = densityGrid->getAccessor();

    // // Create temperature grid and get its accessor
    // openvdb::FloatGrid::Ptr temperatureGrid = openvdb::FloatGrid::create();
    // openvdb::FloatGrid::Accessor temperatureAccessor = temperatureGrid->getAccessor();

    // Define a coordinate with large signed indices.
    openvdb::Coord ijk;

    int& i = ijk[0];
    int& j = ijk[1];
    int& k = ijk[2];

    for (k=0; k<ndepth; k++)
        for (j=0; j<nheight; j++)
            for (i=0; i<nwidth; i++){
                densityAccessor.setValue(ijk, smoke->getQuantity(i,j,k));
                // temperatureAccessor.setValue(ijk, temperature->getQuantity(i,j,k));
            }

    char outfile[128];
    sprintf(outfile, "%sFrame%05d.vdb", outputpath.c_str(), framenumber++);

    std::cout << outfile << std::endl;

    // Create a VDB file object.
    openvdb::io::File file(outfile);

    densityGrid->setName("density");
    // temperatureGrid->setName("temperature");

    // Add the grid pointers to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(densityGrid);
    // grids.push_back(temperatureGrid);

    // Write out the contents of the container.
    file.write(grids);
    file.close();
}


double FluidQuantity::InterpolateLinear(openvdb::Vec3d xyz) const
{
    // Create a linear interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> interpolator(*quantity);

    auto bbox = parentgrid->getBBox();

    // Clamp position to the bbox
    for (int i = 0; i < 3; i++) {
        xyz[i] = std::max(bbox.min()[i], std::min(xyz[i], bbox.max()[i]));
    }

    // Interpolate the value at the given position
    double result = interpolator.wsSample(xyz);

    return result;
}


double FluidQuantity::InterpolateQuadratic(openvdb::Vec3d xyz) const
{
    // Create a cubic interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::QuadraticSampler> sampler(*quantity);

    auto bbox = parentgrid->getBBox();

    // Clamp position to the bbox
    for (int i = 0; i < 3; i++) {
        xyz[i] = std::max(bbox.min()[i], std::min(xyz[i], bbox.max()[i]));
    }

    // Interpolate the value at the given position
    double result = sampler.wsSample(xyz);

    return result;
}


double FluidQuantity::InterpolateCubic(openvdb::Vec3d xyz) const
{
    // Create a cubic interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, fluidsim::tools::CubicSampler> sampler(*quantity);

    auto bbox = parentgrid->getBBox();

    // Clamp position to the bbox
    for (int i = 0; i < 3; i++) {
        xyz[i] = std::max(bbox.min()[i], std::min(xyz[i], bbox.max()[i]));
    }

    // Interpolate the value at the given position
    double result = sampler.wsSample(xyz);

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


double FluidQuantity::sumVolumes() const {
    double total = 0;
    int count = 0;
    for (openvdb::FloatGrid::ValueOnCIter iter = volume->cbeginValueOn(); iter; ++iter) {
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


bool isMatrixSymmetric(const SymmBandMatrix& mat, SymmBandMatrix::ValueType eps) {
    
    // Iterate over every row
    for (int i = 0; i < mat.numRows(); ++i) {
        auto row = mat.getConstRow(i);
        
        // Iterate ONLY over the non-zero elements in this row
        for (auto it = row.cbegin(); it; ++it) {
            int j = it.column();
            SymmBandMatrix::ValueType val_ij = *it;
            
            // Fetch the mirrored value (A_ji)
            SymmBandMatrix::ValueType val_ji = mat.getValue(j, i);
            
            // Floating point comparison
            if (std::abs(val_ij - val_ji) > eps) {
                return false; 
            }
        }
    }
    
    return true;
}


bool rowsSumToZero(const SymmBandMatrix& mat, SymmBandMatrix::ValueType eps) {

    // Iterate over every row
    for (int i = 0; i < mat.numRows(); ++i) {
        auto row = mat.getConstRow(i);
        

        int count = 0;
        SymmBandMatrix::ValueType sum = 0;
        // Iterate ONLY over the non-zero elements in this row
        for (auto it = row.cbegin(); it; ++it) {
            sum += *it;
            count++;
        }

        if (!count) return false;

        if (std::abs(sum) > eps){
            return false;
        }
    }

    return true;
}


openvdb::Vec3d FluidQuantity::indexToWorld(int i, int j, int k){
    auto xform = getTransform();
    openvdb::Vec3d xyz = xform->indexToWorld(openvdb::Coord(i, j, k));
    return xyz;
}


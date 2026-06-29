#include <iostream>
#include "grid.h"
#include "interp.h"
#include <openvdb/tools/VolumeAdvect.h>



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

    r = new SymmBandMatrix::VectorType(ncells);
    pressure = new SymmBandMatrix::VectorType(ncells);
    A = new SymmBandMatrix(ncells);

    velocity = new FluidQuantity<openvdb::Vec3dGrid>(this, true, openvdb::Vec3d(0,0,0));
    smoke = new FluidQuantity<openvdb::FloatGrid>(this, false, 0);
    temperature = new FluidQuantity<openvdb::FloatGrid>(this, false, AMBIENT_TEMPERATURE);
}


template<class T>
FluidQuantity<T>::FluidQuantity(FluidGrid* parent, bool staggered, typename T::ValueType value) :   parentgrid(parent),
                                                                                                    staggered(staggered),
                                                                                                    offsetx(0.5),
                                                                                                    offsety(0.5),
                                                                                                    offsetz(0.5)
{
    using ValueType = typename T::ValueType;

    quantity = T::create(value);

    if (staggered){
        xsamples = parent->getWidth() + 1;
        ysamples = parent->getHeight() + 1;
        zsamples = parent->getDepth() + 1;
        // I don't know what effect, if any, this has
        quantity->setGridClass(openvdb::GRID_STAGGERED);
    }
    else{
        xsamples = parent->getWidth();
        ysamples = parent->getHeight();
        zsamples = parent->getDepth();
    }

    // Adjust grid transform to match desired grid spacing and offset
    openvdb::math::Transform::Ptr q_transform = openvdb::math::Transform::createLinearTransform(parent->getCellwidth());
    q_transform->preTranslate(openvdb::Vec3d(offsetx, offsety, offsetz));
    quantity->setTransform(q_transform);
}


void FluidGrid::setDeltaT(){

    double maxvel = std::max({velocity->max()[0],velocity->max()[1],velocity->max()[2]});
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
    velocity->Advect();

    // we need to enforce solid boundary conditions after advection
    SetWallBoundaries();
}


template <class T>
void FluidQuantity<T>::Advect()
{
    auto velocityGrid = parentgrid->getVelocity()->quantity;
    openvdb::tools::VolumeAdvection advector(*velocityGrid);
    typename T::Ptr outGrid = advector.template advect<T, fluidsim::tools::CubicSampler>(*quantity, parentgrid->getDeltaT());
    std::swap(quantity, outGrid); 
}


void FluidGrid::AddForces(){
    // this function only affects the y-component of the velocity field - x and z components are unaffected.

    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> tsampler(*temperature->quantity);
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> dsampler(*smoke->quantity);
    auto v_access = velocity->quantity->getAccessor();
    double alpha = 1;
    double beta = 1/AMBIENT_TEMPERATURE;
    for (int k = 0; k < velocity->getzSamples(); k++)
        for (int j = 1; j < velocity->getySamples()-1; j++)
            for (int i = 0; i < velocity->getxSamples(); i++){
                auto ijk = openvdb::Coord(i,j,k);
                auto result = getDeltaT() * (alpha*dsampler.isSample(ijk) - beta*(tsampler.isSample(ijk)-AMBIENT_TEMPERATURE)) * GRAVITY;
                v_access.setValue(ijk, v_access.getValue(ijk) + openvdb::Vec3d(0,result,0));
            }
}


void FluidGrid::addSmoke(double x, double y, double z, double wh, double h, double l, double density_value, double temperature_value) {

    double x0 = round((x-wh/2)*nwidth);
    double x1 = round((x+wh/2)*nwidth);
    double y0 = round(y*nheight);
    double y1 = round((y+h)*nheight);
    double z0 = round((z-l/2)*ndepth);
    double z1 = round((z+l/2)*ndepth);

    temperature->addEmitter(x0, y0, z0, x1, y1, z1, temperature_value);;
    smoke->addEmitter(x0, y0, z0, x1, y1, z1, density_value);
}


template <class T>
void FluidQuantity<T>::addEmitter(int x0, int y0, int z0, int x1, int y1, int z1, typename T::ValueType value) {
    auto q_access = quantity->getAccessor();
    for (int z = std::max(z0, 0); z < std::min(z1, zsamples); z++)
        for (int y = std::max(y0, 0); y < std::min(y1, ysamples); y++)
            for (int x = std::max(x0, 0); x < std::min(x1, xsamples); x++)
                q_access.setValue(openvdb::Coord(x,y,z), value);
}


void FluidGrid::updateVelocities() {

    auto v_access = velocity->quantity->getAccessor();
    double scale = deltaT/(density*cellwidth);
    // p.44 Bridson
    for (int k = 0; k < ndepth; k++)
        for (int j = 0; j < nheight; j++)
            for (int i = 0; i < nwidth; i++) {
                auto p = scale*getPressure(i,j,k);
                v_access.setValue(openvdb::Coord(i,j,k), v_access.getValue(openvdb::Coord(i,j,k)) - p);
                v_access.setValue(openvdb::Coord(i+1,j,k), v_access.getValue(openvdb::Coord(i+1,j,k)) + openvdb::Vec3d(p,0,0));
                v_access.setValue(openvdb::Coord(i,j+1,k), v_access.getValue(openvdb::Coord(i,j+1,k)) + openvdb::Vec3d(0,p,0));
                v_access.setValue(openvdb::Coord(i,j,k+1), v_access.getValue(openvdb::Coord(i,j,k+1)) + openvdb::Vec3d(0,0,p));
            }
}


void FluidGrid::SetWallBoundaries(){
    using namespace openvdb;
    auto v_access = velocity->quantity->getAccessor();
    //Boundary conditions for the walls of the box
    //Set the normal component to the wall to zero
    for (int k = 0; k < ndepth+1; k++)
        for (int j = 0; j < nheight+1; j++){
            Vec3d result(v_access.getValue(Coord(0,j,k)));
            v_access.setValue(Coord(0,j,k), Vec3d(0,result[1],result[2]));
            result = v_access.getValue(Coord(nwidth,j,k));
            v_access.setValue(Coord(nwidth,j,k), Vec3d(0,result[1],result[2]));
        }
    for (int k = 0; k < ndepth+1; k++)
        for (int i = 0; i < nwidth+1; i++){
            Vec3d result(v_access.getValue(Coord(i,0,k)));
            v_access.setValue(Coord(i,0,k), Vec3d(result[0],0,result[2]));
            result = v_access.getValue(Coord(i,nheight,k));
            v_access.setValue(Coord(i,nheight,k), Vec3d(result[0],0,result[2]));
        }
    for (int j = 0; j < nheight+1; j++)
        for (int i = 0; i < nwidth+1; i++){
            Vec3d result(v_access.getValue(Coord(i,j,0)));
            v_access.setValue(Coord(i,j,0), Vec3d(result[0],result[1],0));
            result = v_access.getValue(Coord(i,j,ndepth));
            v_access.setValue(Coord(i,j,ndepth), Vec3d(result[0],result[1],0));
        }
}


void FluidGrid::BuildLinearSystem(){
    double Ascale =  deltaT/(density*cellwidth*cellwidth);
    double rscale = 1.0/cellwidth;

    auto v_access = velocity->quantity->getAccessor();
    clearSparseMatrix(A);
    r->fill(0);

    //Setup up the matrix of pressure coefficients, weighted by volumes
    int ix = 0;
    int slice_size  = nwidth * nheight;
    int total_cells = nwidth * nheight * ndepth;
    for(int k=0; k<ndepth; k++)
        for(int j=0; j<nheight; j++)
            for(int i=0; i<nwidth; i++){

                ix = i + (j * nwidth) + (k * slice_size);
                if (i < nwidth-1){
                    addToA(ix, 0, Ascale);
                    addToA(ix + 1, 0, Ascale);
                    setA(ix, 1, -Ascale);
                }
                if (j < nheight-1){
                    addToA(ix, 0, Ascale);
                    addToA(ix + nwidth, 0, Ascale);
                    setA(ix, 2, -Ascale);
                }
                if (k < ndepth-1){
                    addToA(ix, 0, Ascale);
                    addToA(ix + slice_size, 0, Ascale);
                    setA(ix, 3, -Ascale);
                }

                (*r)[ix] = -rscale*(v_access.getValue(openvdb::Coord(i+1,j,k))[0] - v_access.getValue(openvdb::Coord(i,j,k))[0] +
                                    v_access.getValue(openvdb::Coord(i,j+1,k))[1] - v_access.getValue(openvdb::Coord(i,j,k))[1] +
                                    v_access.getValue(openvdb::Coord(i,j,k+1))[2] - v_access.getValue(openvdb::Coord(i,j,k))[2]);
            }
}


void FluidGrid::Project(){

    //DEBUG
    auto pressure_before(*pressure);

    CholeskyPrecondMatrix precond(*A);

    auto state = openvdb::math::pcg::terminationDefaults<double>();
    state.iterations = MAXITER; // Maximum number of iterations
    state.relativeError = TOL; // Tolerance for relative error

    std::cout << "max divergence before solve = " <<     maxDivergence() << std::endl;

    openvdb::util::NullInterrupter interrupter;
    std::cout << "Solving..." << std::endl;
    auto result = openvdb::math::pcg::solve(*A, *r, *pressure, precond, interrupter, state);

    if (!result.success)
        std::cout << "WARNING: Solver failed!" << std::endl;

    if (pressure->eq(pressure_before))
        std::cout << "No change to pressure" << std::endl;

    std::cout << "max divergence after solve = " <<     maxDivergence() << std::endl;
}


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
    AddForces();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "addForces time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
    begin = std::chrono::high_resolution_clock::now();
    addSmoke(0.5, 0.1, 0.5, 0.25, 0.05, 0.25, 0.5, 450);
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "addSmoke time (seconds) = " << elapsed.count() * 1e-9 << std::endl;
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
    updateVelocities();
    SetWallBoundaries();

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

    // Create u-velocity grid and get its accessor
    openvdb::FloatGrid::Ptr uvelocityGrid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor uvelocityAccessor = uvelocityGrid->getAccessor();

    // Create v-velocity grid and get its accessor
    openvdb::FloatGrid::Ptr vvelocityGrid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor vvelocityAccessor = vvelocityGrid->getAccessor();

    // Create w-velocity grid and get its accessor
    openvdb::FloatGrid::Ptr wvelocityGrid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor wvelocityAccessor = wvelocityGrid->getAccessor();


    // Define a coordinate with large signed indices.
    openvdb::Coord ijk;

    int& i = ijk[0];
    int& j = ijk[1];
    int& k = ijk[2];

    auto s_access = smoke->quantity->getAccessor();
    for (k=0; k<ndepth; k++)
        for (j=0; j<nheight; j++)
            for (i=0; i<nwidth; i++){
                densityAccessor.setValue(ijk, s_access.getValue(ijk));
                // uvelocityAccessor.setValue(ijk, fabs(velocity->getQuantity(i,j,k)[0]));
                // vvelocityAccessor.setValue(ijk, fabs(velocity->getQuantity(i,j,k)[1]));
                // wvelocityAccessor.setValue(ijk, fabs(velocity->getQuantity(i,j,k)[2]));

            }

    char outfile[128];
    sprintf(outfile, "%sFrame%05d.vdb", outputpath.c_str(), framenumber++);

    std::cout << outfile << std::endl;

    // Create a VDB file object.
    openvdb::io::File file(outfile);

    densityGrid->setName("density");
    uvelocityGrid->setName("uvelocity");
    vvelocityGrid->setName("vvelocity");
    wvelocityGrid->setName("wvelocity");


    // Add the grid pointers to a container.
    openvdb::GridPtrVec grids;
    grids.push_back(densityGrid);
    // grids.push_back(uvelocityGrid);
    // grids.push_back(vvelocityGrid);
    // grids.push_back(wvelocityGrid);

    // Write out the contents of the container.
    file.write(grids);
    file.close();
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


template <class T>
typename T::ValueType FluidQuantity<T>::max() const {
    typename T::ValueType maxValue{0};
    for (auto iter = quantity->cbeginValueOn(); iter; ++iter) {
        if (*iter > maxValue) maxValue = *iter;
    }
    return std::max(maxValue, quantity->background());
}


// https://stackoverflow.com/questions/28354752/template-vs-template-without-brackets-whats-the-difference
// Instantiations
template class FluidQuantity<openvdb::FloatGrid>;
template class FluidQuantity<openvdb::Vec3dGrid>;

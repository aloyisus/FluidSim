#include "grid.h"
#include <gtest/gtest.h>
#include "maths.h"
#include <numeric>
#include "interp.h"


TEST(GridTests, QuantityTests) {
    std::string filepath("dummy");
    FluidGrid* solver = new FluidGrid(4, 4, 4, 0.2, 1, filepath);

    // define a voxel-centred quantity
    FluidQuantity<openvdb::FloatGrid>* u = new FluidQuantity<openvdb::FloatGrid>(solver, false, 7.0);
    
    EXPECT_EQ(solver->getCellwidth(), 0.25);

    u->setQuantity(0,0,0, 3.0);
    u->setQuantity(1,1,1, 9.0);
    u->setQuantity(3,3,3, 27.0);

    EXPECT_EQ(u->getQuantity(1,0,0), 7);
    EXPECT_EQ(u->getQuantity(0,0,0), 3);
    EXPECT_EQ(u->getQuantity(1,1,1), 9);
    EXPECT_EQ(u->getQuantity(3,3,3), 27);

    u->swap();

    EXPECT_EQ(u->getQuantity(1,0,0), 7);
    EXPECT_EQ(u->getQuantity(0,0,0), 7);
    EXPECT_EQ(u->getQuantity(1,1,1), 7);
    EXPECT_EQ(u->getQuantity(3,3,3), 7);

    u->setQuantity(0,0,0, 4.0);
    u->setQuantity(1,1,1, 10.0);
    u->setQuantity(3,3,3, 28.0);

    u->swap();

    EXPECT_EQ(u->getQuantity(1,0,0), 7);
    EXPECT_EQ(u->getQuantity(0,0,0), 3);
    EXPECT_EQ(u->getQuantity(1,1,1), 9);
    EXPECT_EQ(u->getQuantity(3,3,3), 27);  

    u->swap();

    EXPECT_EQ(u->getQuantity(1,0,0), 7);
    EXPECT_EQ(u->getQuantity(0,0,0), 4);
    EXPECT_EQ(u->getQuantity(1,1,1), 10);
    EXPECT_EQ(u->getQuantity(3,3,3), 28);      

    // define a voxel-centred quantity
    FluidQuantity<openvdb::FloatGrid>* v = new FluidQuantity<openvdb::FloatGrid>(solver, false, 0);
    for (int k=0; k<4; k++)
        for (int j=0; j<4; j++)
            for (int i=0; i<4; i++){
                float value = i + 4*j + 16*k;
                v->setQuantity(i,j,k, value);
            }
    EXPECT_EQ(v->getQuantity(0,0,0), 0);
    EXPECT_EQ(v->getQuantity(3,3,3), 63);

    // InterpolateLinear takes world-space position as input
    EXPECT_EQ(v->InterpolateLinear(openvdb::Vec3d(0.125,0.125,0.125)), 0);
    EXPECT_EQ(v->InterpolateLinear(openvdb::Vec3d(0.875,0.875,0.875)), 63);
    EXPECT_EQ(v->InterpolateLinear(openvdb::Vec3d(0.625,0.625,0.625)), 42);
    EXPECT_EQ(v->InterpolateLinear(openvdb::Vec3d(0.375,0.625,0.875)), 57);
    EXPECT_EQ(v->InterpolateLinear(openvdb::Vec3d(0.5,0.5,0.5)), 31.5);

    std::cout << v->InterpolateCubic(openvdb::Vec3d(0.125,0.125,0.125)) << std::endl;
    std::cout << v->InterpolateCubic(openvdb::Vec3d(0.875,0.875,0.875)) << std::endl;
    std::cout << v->InterpolateCubic(openvdb::Vec3d(0.625,0.625,0.625)) << std::endl;

    EXPECT_EQ(v->InterpolateCubic(openvdb::Vec3d(0.125,0.125,0.125)), 0);
    EXPECT_EQ(v->InterpolateCubic(openvdb::Vec3d(0.875,0.875,0.875)), 63);
    EXPECT_EQ(v->InterpolateCubic(openvdb::Vec3d(0.625,0.625,0.625)), 42);
    EXPECT_EQ(v->InterpolateCubic(openvdb::Vec3d(0.375,0.625,0.875)), 57);
    EXPECT_EQ(v->InterpolateCubic(openvdb::Vec3d(0.5,0.5,0.5)), 31.5);

}


TEST(OpenVDBTests, FloatGrid) {

    openvdb::FloatGrid::Ptr u = openvdb::FloatGrid::create(5);
    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(/* spacing */ 1.0);
    u->setTransform(transform);
    openvdb::FloatGrid::Accessor u_access = u->getAccessor();

    // Define the bounding box in voxel coordinates
    //bbox(openvdb::Coord(0, 0, 0), openvdb::Coord(10, 10, 10));
    openvdb::CoordBBox bbox = openvdb::CoordBBox(openvdb::Coord(0, 0, 0), openvdb::Coord(10, 10, 10));

    u_access.setValue(openvdb::Coord(0, 0, 0), 3);
    u_access.setValue(openvdb::Coord(3, 3, 3), 27);
    u_access.setValue(openvdb::Coord(5, 0, 2), 45);
    u_access.setValue(openvdb::Coord(9, 8, 7), 99);   

    float sum = 0;
    for (int x = bbox.min().x(); x < bbox.max().x(); ++x) {
        for (int y = bbox.min().y(); y < bbox.max().y(); ++y) {
            for (int z = bbox.min().z(); z < bbox.max().z(); ++z) {
                openvdb::Coord coord(x, y, z);
                float value = u_access.getValue(coord);
                sum += value;
            }
        }
    }

    // 996 * 5 (the background value) + 3 + 27 + 45 + 99 = 5154
    EXPECT_EQ(sum, 5154);

    //Matrix for pressure coefficients A
    auto A = SymmBandMatrix(64);
    A.setValue(0, 0, 1);
    A.setValue(0, 1, 5);
    A.setValue(1, 0, 13);
    A.setValue(1, 1, 2);
    A.setValue(2, 2, 3);
    A.setValue(63, 63, 4);

    EXPECT_EQ(A.getValue(0, 0), 1);
    EXPECT_EQ(A.getValue(0, 1), 5);
    EXPECT_EQ(A.getValue(1, 0), 13);
    EXPECT_EQ(A.getValue(1, 1), 2);
    EXPECT_EQ(A.getValue(2, 2), 3);
    EXPECT_EQ(A.getValue(63, 63), 4);

    A.scale(0);
    EXPECT_EQ(A.getValue(0, 0), 0);
    EXPECT_EQ(A.getValue(0, 1), 0);
    EXPECT_EQ(A.getValue(1, 1), 0);
    EXPECT_EQ(A.getValue(2, 2), 0);
    EXPECT_EQ(A.getValue(63, 63), 0);

    // define a voxel-centred quantity and populate the grid
    openvdb::FloatGrid::Ptr w = openvdb::FloatGrid::create(0);
    transform = openvdb::math::Transform::createLinearTransform(/* spacing */ 1.0);
    transform->preTranslate(openvdb::Vec3d(0.5, 0.5, 0.5));
    w->setTransform(transform);
    openvdb::FloatGrid::Accessor w_access = w->getAccessor();
    for (int k=0; k<4; k++)
        for (int j=0; j<4; j++)
            for (int i=0; i<4; i++){
                float value = i + 4*j + 16*k;
                w_access.setValue(openvdb::Coord(i, j, k), value);
            }
    EXPECT_EQ(w_access.getValue(openvdb::Coord(0, 0, 0)), 0);
    EXPECT_EQ(w_access.getValue(openvdb::Coord(3, 3, 3)), 63);

    // Create a trilinear interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> interpolator(*w);

    EXPECT_EQ(interpolator.wsSample(openvdb::Vec3d(0.5,0.5,0.5)), 0);
    EXPECT_EQ(interpolator.wsSample(openvdb::Vec3d(3.5,3.5,3.5)), 63);
    EXPECT_EQ(interpolator.wsSample(openvdb::Vec3d(2.5,2.5,2.5)), 42);
    EXPECT_EQ(interpolator.wsSample(openvdb::Vec3d(1.5,2.5,3.5)), 57);
    EXPECT_EQ(interpolator.wsSample(openvdb::Vec3d(2,2,2)), 31.5);

    // Create a cubic interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::QuadraticSampler> cubic_interpolator(*w);

    EXPECT_EQ(cubic_interpolator.wsSample(openvdb::Vec3d(0.5,0.5,0.5)), 0);
    EXPECT_EQ(cubic_interpolator.wsSample(openvdb::Vec3d(3.5,3.5,3.5)), 63);
    EXPECT_EQ(cubic_interpolator.wsSample(openvdb::Vec3d(2.5,2.5,2.5)), 42);
    EXPECT_EQ(cubic_interpolator.wsSample(openvdb::Vec3d(1.5,2.5,3.5)), 57);
    EXPECT_EQ(cubic_interpolator.wsSample(openvdb::Vec3d(2,2,2)), 31.5);


    // Try to understand transform
    auto q = openvdb::FloatGrid::create(0);
    auto q_access = openvdb::FloatGrid::Accessor(q->getAccessor());

    // Adjust grid transform to match desired grid spacing and offset
    openvdb::math::Transform::Ptr q_transform = openvdb::math::Transform::createLinearTransform(0.33);

    q_transform->preTranslate(openvdb::Vec3d(0.5, 0.5, 0.5));
    q->setTransform(q_transform);

    // in world-space each cell is 0.33 wide, and co-ordinates map to the centre of the cells
    auto ijk = openvdb::Coord(1, 1, 1);
    openvdb::Vec3d worldSpacePoint = q_transform->indexToWorld(ijk);

    // the co-ordinate 1,1,1 therefore maps to (1 + 0.5) * 0.33 = 0.495
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0.495,0.495,0.495));

    // can get cell size from the transform
    EXPECT_EQ(q_transform->voxelSize(), openvdb::Vec3d(0.33,0.33,0.33));

    ijk = openvdb::Coord(1, 1, 1);
    auto lmn = openvdb::Coord(1, 1, 1);
    EXPECT_EQ(ijk, lmn);

    auto xyz = openvdb::Vec3d(10,11,12);
    EXPECT_EQ(xyz, openvdb::Vec3d(10,11,12));   

    openvdb::Vec3d centrepos;
    centrepos(0) = 10;
    centrepos(1) = 11;
    centrepos(2) = 77.9;

    openvdb::Vec3d pos(5,6,2);

    EXPECT_EQ(centrepos - pos, openvdb::Vec3d(5, 5, 75.9));

}


TEST(TestConjGradient, testJacobi)
{
    using namespace openvdb;

    typedef math::pcg::SparseStencilMatrix<double, 7> MatrixType;

    const math::pcg::SizeType rows = 5;

    MatrixType A(rows);
    A.setValue(0, 0, 24.0);
    A.setValue(0, 2,  6.0);
    A.setValue(1, 1,  8.0);
    A.setValue(1, 2,  2.0);
    A.setValue(2, 0,  6.0);
    A.setValue(2, 1,  2.0);
    A.setValue(2, 2,  8.0);
    A.setValue(2, 3, -6.0);
    A.setValue(2, 4,  2.0);
    A.setValue(3, 2, -6.0);
    A.setValue(3, 3, 24.0);
    A.setValue(4, 2,  2.0);
    A.setValue(4, 4,  8.0);

    EXPECT_TRUE(A.isFinite());

    MatrixType::VectorType
        x(rows, 0.0),
        b(rows, 1.0),
        expected(rows);

    expected[0] = 0.0104167;
    expected[1] = 0.09375;
    expected[2] = 0.125;
    expected[3] = 0.0729167;
    expected[4] = 0.09375;

    math::pcg::JacobiPreconditioner<MatrixType> precond(A);

    // Solve A * x = b for x.
    math::pcg::State result = math::pcg::solve(
        A, b, x, precond, math::pcg::terminationDefaults<double>());

    EXPECT_TRUE(result.success);
    EXPECT_TRUE(result.iterations <= 20);
    EXPECT_TRUE(x.eq(expected, 1.0e-5));
}


TEST(TestConjGradient, testIncompleteCholesky)
{
    using namespace openvdb;

    typedef math::pcg::SparseStencilMatrix<double, 7> MatrixType;
    typedef math::pcg::IncompleteCholeskyPreconditioner<MatrixType> CholeskyPrecond;

    const math::pcg::SizeType rows = 5;

    MatrixType A(5);
    A.setValue(0, 0, 24.0);
    A.setValue(0, 2,  6.0);
    A.setValue(1, 1,  8.0);
    A.setValue(1, 2,  2.0);
    A.setValue(2, 0,  6.0);
    A.setValue(2, 1,  2.0);
    A.setValue(2, 2,  8.0);
    A.setValue(2, 3, -6.0);
    A.setValue(2, 4,  2.0);
    A.setValue(3, 2, -6.0);
    A.setValue(3, 3, 24.0);
    A.setValue(4, 2,  2.0);
    A.setValue(4, 4,  8.0);

    EXPECT_TRUE(A.isFinite());

    CholeskyPrecond precond(A);
    {
        const CholeskyPrecond::TriangularMatrix lower = precond.lowerMatrix();

        CholeskyPrecond::TriangularMatrix expected(5);
        expected.setValue(0, 0,  4.89898);
        expected.setValue(1, 1,  2.82843);
        expected.setValue(2, 0,  1.22474);
        expected.setValue(2, 1,  0.707107);
        expected.setValue(2, 2,  2.44949);
        expected.setValue(3, 2, -2.44949);
        expected.setValue(3, 3,  4.24264);
        expected.setValue(4, 2,  0.816497);
        expected.setValue(4, 4,  2.70801);

        EXPECT_TRUE(lower.eq(expected, 1.0e-5));
    }
    {
        const CholeskyPrecond::TriangularMatrix upper = precond.upperMatrix();

        CholeskyPrecond::TriangularMatrix expected(5);
        {
            expected.setValue(0, 0,  4.89898);
            expected.setValue(0, 2,  1.22474);
            expected.setValue(1, 1,  2.82843);
            expected.setValue(1, 2,  0.707107);
            expected.setValue(2, 2,  2.44949);
            expected.setValue(2, 3, -2.44949);
            expected.setValue(2, 4,  0.816497);
            expected.setValue(3, 3,  4.24264);
            expected.setValue(4, 4,  2.70801);
        }

        EXPECT_TRUE(upper.eq(expected, 1.0e-5));
    }

    MatrixType::VectorType
        x(rows, 0.0),
        b(rows, 1.0),
        expected(rows);

    expected[0] = 0.0104167;
    expected[1] = 0.09375;
    expected[2] = 0.125;
    expected[3] = 0.0729167;
    expected[4] = 0.09375;

    // Solve A * x = b for x.
    math::pcg::State result = math::pcg::solve(
        A, b, x, precond, math::pcg::terminationDefaults<double>());

    EXPECT_TRUE(result.success);
    EXPECT_TRUE(result.iterations <= 20);
    EXPECT_TRUE(x.eq(expected, 1.0e-5));
}


TEST(TestConjGradient, testIncompleteCholesky2)
{
    using namespace openvdb;

    std::string filepath("dummy");
    const math::pcg::SizeType rows = 5;
    FluidGrid solver = FluidGrid(rows, 1, 1, 0.2, 1, filepath);

    solver.r->fill(1);

    solver.A->setValue(0, 0, 24.0);
    solver.A->setValue(0, 2,  6.0);
    solver.A->setValue(1, 1,  8.0);
    solver.A->setValue(1, 2,  2.0);
    solver.A->setValue(2, 0,  6.0);
    solver.A->setValue(2, 1,  2.0);
    solver.A->setValue(2, 2,  8.0);
    solver.A->setValue(2, 3, -6.0);
    solver.A->setValue(2, 4,  2.0);
    solver.A->setValue(3, 2, -6.0);
    solver.A->setValue(3, 3, 24.0);
    solver.A->setValue(4, 2,  2.0);
    solver.A->setValue(4, 4,  8.0);

    EXPECT_TRUE(solver.A->isFinite());

    CholeskyPrecondMatrix precond(*(solver.A));
    {
        const CholeskyPrecondMatrix::TriangularMatrix lower = precond.lowerMatrix();

        CholeskyPrecondMatrix::TriangularMatrix expected(5);
        expected.setValue(0, 0,  4.89898);
        expected.setValue(1, 1,  2.82843);
        expected.setValue(2, 0,  1.22474);
        expected.setValue(2, 1,  0.707107);
        expected.setValue(2, 2,  2.44949);
        expected.setValue(3, 2, -2.44949);
        expected.setValue(3, 3,  4.24264);
        expected.setValue(4, 2,  0.816497);
        expected.setValue(4, 4,  2.70801);

        EXPECT_TRUE(lower.eq(expected, 1.0e-5));
    }
    {
        const CholeskyPrecondMatrix::TriangularMatrix upper = precond.upperMatrix();

        CholeskyPrecondMatrix::TriangularMatrix expected(5);
        {
            expected.setValue(0, 0,  4.89898);
            expected.setValue(0, 2,  1.22474);
            expected.setValue(1, 1,  2.82843);
            expected.setValue(1, 2,  0.707107);
            expected.setValue(2, 2,  2.44949);
            expected.setValue(2, 3, -2.44949);
            expected.setValue(2, 4,  0.816497);
            expected.setValue(3, 3,  4.24264);
            expected.setValue(4, 4,  2.70801);
        }

        EXPECT_TRUE(upper.eq(expected, 1.0e-5));
    }

    SymmBandMatrix::VectorType
        expected(rows);

    expected[0] = 0.0104167;
    expected[1] = 0.09375;
    expected[2] = 0.125;
    expected[3] = 0.0729167;
    expected[4] = 0.09375;

    // Solve A * x = b for x.
    math::pcg::State result = math::pcg::solve(*(solver.A), *(solver.r), *(solver.pressure), precond, math::pcg::terminationDefaults<double>());

    EXPECT_TRUE(result.success);
    EXPECT_TRUE(result.iterations <= 20);
    EXPECT_TRUE(solver.pressure->eq(expected, 1.0e-5));
}


// TEST(TransformTests, indexToWorld)
// {
//     using namespace openvdb;

//     std::string filepath("dummy");
//     FluidGrid solver = FluidGrid(500, 2000, 1000, 0.2, 1, filepath);

//     auto xformu = solver.velocity->getTransform();
//     auto xformsmoke = solver.smoke->getTransform();

//     openvdb::Vec3d xyzu = xformu->indexToWorld(openvdb::Coord(50,50,50));
//     EXPECT_TRUE(xyzu[0] == 0.100);
//     EXPECT_TRUE(xyzu[1] == 0.101);
//     EXPECT_TRUE(xyzu[2] == 0.101);



//     // smoke is a cell-centred quantity, so offsets are (0.5,0.5,0.5)
//     openvdb::Vec3d xyzsmoke = xformsmoke->indexToWorld(openvdb::Coord(50,50,50));
//     EXPECT_TRUE(xyzsmoke[0] == 0.101);
//     EXPECT_TRUE(xyzsmoke[1] == 0.101);
//     EXPECT_TRUE(xyzsmoke[2] == 0.101);

//     auto cellwidth = solver.getCellwidth();
//     EXPECT_TRUE(cellwidth == 0.002);
// }

TEST(SymmBandMatrix, clear_sparse_matrix)
{
    //Matrix for pressure coefficients A
    auto A = SymmBandMatrix(64);
    SymmBandMatrix *pA = &A;
    A.setValue(0, 0, 1);
    A.setValue(0, 1, 5);
    A.setValue(1, 1, 2);
    A.setValue(2, 2, 3);
    A.setValue(63, 63, 4);

    EXPECT_EQ(A.getValue(0, 0), 1);
    EXPECT_EQ(A.getValue(0, 1), 5);
    EXPECT_EQ(A.getValue(1, 1), 2);
    EXPECT_EQ(A.getValue(2, 2), 3);
    EXPECT_EQ(A.getValue(63, 63), 4);

    clearSparseMatrix(pA);
    EXPECT_EQ(A.getValue(0, 0), 0);
    EXPECT_EQ(A.getValue(0, 1), 0);
    EXPECT_EQ(A.getValue(1, 1), 0);
    EXPECT_EQ(A.getValue(2, 2), 0);
    EXPECT_EQ(A.getValue(63, 63), 0);
}

TEST(CGSolver, instantiate_solver){
    SymmBandMatrix A = SymmBandMatrix(64);
    SymmBandMatrix::VectorType c(64,3);
    SymmBandMatrix::VectorType p(64,7);
    double tol;
    int maxIter;

    std::vector<double> a{1,2,3,4,5,6,7,8};
    std::vector<double> b{2,3,4,5,6,7,8,9};    
    double prod = std::inner_product(a.begin(),a.end(),b.begin(),0);
    EXPECT_EQ(prod, 240);

    SymmBandMatrix::VectorType result = c + p;
    EXPECT_EQ(result[63], 10);

    result = 5*(c+p);
    EXPECT_EQ(result[0], 50);

    result[0] = 999;
    double maxi = max(result);
    EXPECT_EQ(maxi, 999);


    // CGSolver(A, b, p, tol, maxIter);
}

TEST(Interpolation, gridsolver){

    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);

    openvdb::FloatGrid::Accessor grid_accessor = grid->getAccessor();
    // first cell, one sample is 66 opposite diagonal opposite sample is 1200
    grid_accessor.setValue(openvdb::Coord(0, 0, 0), 66);
    grid_accessor.setValue(openvdb::Coord(1, 1, 1), 1200);
    EXPECT_EQ(grid_accessor.getValue(openvdb::Coord(0, 0, 0)), 66);
    EXPECT_EQ(grid_accessor.getValue(openvdb::Coord(1, 1, 1)), 1200);


    // https://www.openvdb.org/documentation/doxygen/transformsAndMaps.html
    // create a linear transformation which scales i,j,k by 0.5
    auto linearTransform = openvdb::math::Transform::createLinearTransform(0.5);

    // Compute the location in world space that is the image of (0,0,0).
    // The result will be (0, 0, 0). There is no offset.
    openvdb::Coord ijk(0,0,0);
    openvdb::Vec3d worldSpacePoint = linearTransform->indexToWorld(ijk);
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0,0,0));

    // Compute the location in world space that is the image of (1,1,1).
    // The result will be (0.5, 0.5, 0.5). There is no offset.
    ijk = openvdb::Coord(1,1,1);
    worldSpacePoint = linearTransform->indexToWorld(ijk);
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0.5,0.5,0.5));

    // Compute the location in world space that is the image of (1,2,3).
    // The result will be (0.5, 1, 1.5). There is no offset.
    ijk = openvdb::Coord(1,2,3);
    worldSpacePoint = linearTransform->indexToWorld(ijk);
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0.5,1,1.5));

    // Compute the location in index space that is the pre-image of (0.5, 1, 1.5).
    // The result will be (1.0, 2.0, 3.0).
    openvdb::Vec3d indexSpacePoint = linearTransform->worldToIndex(worldSpacePoint);
    EXPECT_EQ(indexSpacePoint, openvdb::Vec3d(1,2,3));





    // create a linear transformation which scales i,j,k by 0.5
    // Also offset the voxels so the centre of a voxel in world space is (0,0.25,0.25) using preTranslate
    auto linearTransform_pretranslate = openvdb::math::Transform::createLinearTransform(0.5);
    linearTransform_pretranslate->preTranslate(openvdb::Vec3d(0, 0.5, 0.5));

    // Compute the location in world space that is the image of (0,0,0).
    ijk = openvdb::Coord(0,0,0);
    worldSpacePoint = linearTransform_pretranslate->indexToWorld(ijk);
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0,0.25,0.25));

    // Compute the location in world space that is the image of (1,2,3).
    // The result will be (0.5, 1.25, 1.75).
    ijk = openvdb::Coord(1,2,3);
    worldSpacePoint = linearTransform_pretranslate->indexToWorld(ijk);
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0.5,1.25,1.75));






    // create a linear transformation which scales i,j,k by 0.5
    // Also offset the voxels BUT THIS TIME USE POSTTRANSLATE METHOD
    // so the centre of (0,0,0) voxel in world space is (0,0.5,0.5)
    // note that the voxel size is still 0.5 x 0.5 x 0.5
    auto linearTransform_posttranslate = openvdb::math::Transform::createLinearTransform(0.5);
    linearTransform_posttranslate->postTranslate(openvdb::Vec3d(0, 0.5, 0.5));

    // Compute the location in world space that is the image of (0,0,0).
    ijk = openvdb::Coord(0,0,0);
    worldSpacePoint = linearTransform_posttranslate->indexToWorld(ijk);
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0,0.5,0.5));

    // Compute the location in world space that is the image of (1,2,3).
    // The result will be (1, 1.5, 2).
    ijk = openvdb::Coord(1,2,3);
    worldSpacePoint = linearTransform_posttranslate->indexToWorld(ijk);
    EXPECT_EQ(worldSpacePoint, openvdb::Vec3d(0.5,1.5,2));








    
    // associate the transformation with the grid
    grid->setTransform(linearTransform_pretranslate);

    // // Create a trilinear sampler
    // openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> box_sampler(*grid);

    // EXPECT_EQ(box_sampler.wsSample(openvdb::Vec3d(0,0,0)), 0);
    // EXPECT_EQ(box_sampler.wsSample(openvdb::Vec3d(0.25,0.25,0.25)), 8.25);
    // EXPECT_EQ(box_sampler.wsSample(openvdb::Vec3d(0,0.25,0.25)), 16.5);
    // EXPECT_EQ(box_sampler.wsSample(openvdb::Vec3d(0.5,0.75,0.75)), 300);
    // EXPECT_EQ(box_sampler.wsSample(openvdb::Vec3d(0.5,0.75,0.5)), 0);

    // // std::cout << box_sampler.wsSample(openvdb::Vec3d(0.5,0.5,0.5)) << std::endl;
    // // std::cout << box_sampler.wsSample(openvdb::Vec3d(0,0,0)) << std::endl;
    // // std::cout << box_sampler.wsSample(openvdb::Vec3d(0.25,0.25,0.25)) << std::endl;
    // // std::cout << box_sampler.wsSample(openvdb::Vec3d(0,0.25,0.25)) << std::endl;
    // // std::cout << box_sampler.wsSample(openvdb::Vec3d(0.5,0.75,0.75)) << std::endl;
    // // std::cout << box_sampler.wsSample(openvdb::Vec3d(0.5,0.75,0.5)) << std::endl;

    // // Create a cubic interpolator
    // openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::QuadraticSampler> cubic_sampler(*grid);

    // EXPECT_EQ(cubic_sampler.wsSample(openvdb::Vec3d(0.5,0.5,0.5)), 0);
    // EXPECT_EQ(cubic_sampler.wsSample(openvdb::Vec3d(0.25,0.25,0.25)), 6.9609375);
    // EXPECT_EQ(cubic_sampler.wsSample(openvdb::Vec3d(0,0.25,0.25)), 9.28125);
    // EXPECT_EQ(cubic_sampler.wsSample(openvdb::Vec3d(0.5,0.75,0.75)), 168.75);
    // EXPECT_EQ(cubic_sampler.wsSample(openvdb::Vec3d(0.5,0.75,0.5)), 0);
}



TEST(Interpolation, cubic){
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);

    auto transform = openvdb::math::Transform::createLinearTransform(1.0);
    transform->preTranslate(openvdb::Vec3d(0.5, 0.5, 0.5));
    grid->setTransform(transform);

    openvdb::FloatGrid::Accessor grid_accessor = grid->getAccessor();
    // Create a cubic interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::QuadraticSampler> cubic_interpolator(*grid);

    // i = 0
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ 1   1   1   │
    //   1   │ 1   1   1   │
    //   0   │ 1   1   1   │
    //       └─────────────┘
    // i = 1 (empty)
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ ·   ·   ·   │
    //   1   │ ·   ·   ·   │
    //   0   │ ·   ·   ·   │
    //       └─────────────┘
    // i = 2 (empty)
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ ·   ·   ·   │
    //   1   │ ·   ·   ·   │
    //   0   │ ·   ·   ·   │
    //       └─────────────┘

    grid_accessor.setValue(openvdb::Coord(0, 2, 2), 1);
    grid_accessor.setValue(openvdb::Coord(0, 2, 1), 1);
    grid_accessor.setValue(openvdb::Coord(0, 2, 0), 1);
    grid_accessor.setValue(openvdb::Coord(0, 1, 2), 1);
    grid_accessor.setValue(openvdb::Coord(0, 1, 1), 1);
    grid_accessor.setValue(openvdb::Coord(0, 1, 0), 1);
    grid_accessor.setValue(openvdb::Coord(0, 0, 2), 1);
    grid_accessor.setValue(openvdb::Coord(0, 0, 1), 1);
    grid_accessor.setValue(openvdb::Coord(0, 0, 0), 1);
    grid_accessor.setValue(openvdb::Coord(1, 2, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 0), 0);
    grid_accessor.setValue(openvdb::Coord(2, 2, 2), 0);
    grid_accessor.setValue(openvdb::Coord(2, 2, 1), 0);
    grid_accessor.setValue(openvdb::Coord(2, 2, 0), 0);
    grid_accessor.setValue(openvdb::Coord(2, 1, 2), 0);
    grid_accessor.setValue(openvdb::Coord(2, 1, 1), 0);
    grid_accessor.setValue(openvdb::Coord(2, 1, 0), 0);
    grid_accessor.setValue(openvdb::Coord(2, 0, 2), 0);
    grid_accessor.setValue(openvdb::Coord(2, 0, 1), 0);
    grid_accessor.setValue(openvdb::Coord(2, 0, 0), 0);

    // Centre of the 3 x 3 block of voxels defined above is (1.5,1.5,1.5) in world space

    // Note that interpolated value can be negative, even tbough none of the grid values are
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.5, 1.6, 1.5)),  1,       1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.6, 1.6, 1.5)),  0.99000, 1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.7, 1.6, 1.5)),  0.95999, 1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.8, 1.6, 1.5)),  0.91000, 1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.9, 1.6, 1.5)),  0.83999, 1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1,   1.6, 1.5)),  0.75,    1e-5);    
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.1, 1.6, 1.5)),  0.63999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.2, 1.6, 1.5)),  0.50999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.3, 1.6, 1.5)),  0.36000, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.4, 1.6, 1.5)),  0.18999, 1e-5);       
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.5, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.6, 1.6, 1.5)), -0.04500, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.7, 1.6, 1.5)), -0.07999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.8, 1.6, 1.5)), -0.10499, 1e-5);   
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.9, 1.6, 1.5)), -0.11999, 1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2,   1.6, 1.5)), -0.125,   1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.1, 1.6, 1.5)), -0.11999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.2, 1.6, 1.5)), -0.10499, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.3, 1.6, 1.5)), -0.07999, 1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.4, 1.6, 1.5)), -0.04500, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.5, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.6, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.7, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.8, 1.6, 1.5)),  0,       1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.9, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(3,   1.6, 1.5)),  0,       1e-5);
}


TEST(Interpolation, cubic2){
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);

    auto transform = openvdb::math::Transform::createLinearTransform(1.0);
    transform->preTranslate(openvdb::Vec3d(0.5, 0.5, 0.5));
    grid->setTransform(transform);

    openvdb::FloatGrid::Accessor grid_accessor = grid->getAccessor();
    // Create a cubic interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::QuadraticSampler> cubic_interpolator(*grid);

    // i = 0
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ ·   ·   ·   │
    //   1   │ ·   ·   ·   │
    //   0   │ ·   ·   ·   │
    //       └─────────────┘
    // i = 1 (empty)
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ ·   ·   ·   │
    //   1   │ ·   ·   ·   │
    //   0   │ ·   ·   ·   │
    //       └─────────────┘
    // i = 2 (empty)
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ 1   1   1   │
    //   1   │ 1   1   1   │
    //   0   │ 1   1   1   │
    //       └─────────────┘

    grid_accessor.setValue(openvdb::Coord(0, 2, 2), 0);
    grid_accessor.setValue(openvdb::Coord(0, 2, 1), 0);
    grid_accessor.setValue(openvdb::Coord(0, 2, 0), 0);
    grid_accessor.setValue(openvdb::Coord(0, 1, 2), 0);
    grid_accessor.setValue(openvdb::Coord(0, 1, 1), 0);
    grid_accessor.setValue(openvdb::Coord(0, 1, 0), 0);
    grid_accessor.setValue(openvdb::Coord(0, 0, 2), 0);
    grid_accessor.setValue(openvdb::Coord(0, 0, 1), 0);
    grid_accessor.setValue(openvdb::Coord(0, 0, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 0), 0);
    grid_accessor.setValue(openvdb::Coord(2, 2, 2), 1);
    grid_accessor.setValue(openvdb::Coord(2, 2, 1), 1);
    grid_accessor.setValue(openvdb::Coord(2, 2, 0), 1);
    grid_accessor.setValue(openvdb::Coord(2, 1, 2), 1);
    grid_accessor.setValue(openvdb::Coord(2, 1, 1), 1);
    grid_accessor.setValue(openvdb::Coord(2, 1, 0), 1);
    grid_accessor.setValue(openvdb::Coord(2, 0, 2), 1);
    grid_accessor.setValue(openvdb::Coord(2, 0, 1), 1);
    grid_accessor.setValue(openvdb::Coord(2, 0, 0), 1);

    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.5, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.6, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.7, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.8, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.9, 1.6, 1.5)),  0,       1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1,   1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.1, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.2, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.3, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.4, 1.6, 1.5)),  0,       1e-5);      
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.5, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.6, 1.6, 1.5)),  0.05499, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.7, 1.6, 1.5)),  0.11999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.8, 1.6, 1.5)),  0.19499, 1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.9, 1.6, 1.5)),  0.28000, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2,   1.6, 1.5)),  0.375,   1e-5);   
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.1, 1.6, 1.5)),  0.47999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.2, 1.6, 1.5)),  0.59500, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.3, 1.6, 1.5)),  0.72000, 1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.4, 1.6, 1.5)),  0.85500, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.5, 1.6, 1.5)),  1,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.6, 1.6, 1.5)),  0.99000, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.7, 1.6, 1.5)),  0.95999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.8, 1.6, 1.5)),  0.91000, 1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.9, 1.6, 1.5)),  0.83999, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(3,   1.6, 1.5)),  0.75,    1e-5);
}


TEST(Interpolation, monotonic_cubic){
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0);

    auto transform = openvdb::math::Transform::createLinearTransform(1.0);
    transform->preTranslate(openvdb::Vec3d(0.5, 0.5, 0.5));
    grid->setTransform(transform);

    openvdb::FloatGrid::Accessor grid_accessor = grid->getAccessor();
    // Create a cubic interpolator
    openvdb::tools::GridSampler<openvdb::FloatGrid, fluidsim::tools::CubicSampler> cubic_interpolator(*grid);

    // i = 0
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ ·   ·   ·   │
    //   1   │ ·   ·   ·   │
    //   0   │ ·   ·   ·   │
    //       └─────────────┘
    // i = 1 (empty)
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ ·   ·   ·   │
    //   1   │ ·   ·   ·   │
    //   0   │ ·   ·   ·   │
    //       └─────────────┘
    // i = 2 (empty)
    // k →
    // j ↓     0   1   2
    //       ┌─────────────┐
    //   2   │ 1   1   1   │
    //   1   │ 1   1   1   │
    //   0   │ 1   1   1   │
    //       └─────────────┘

    grid_accessor.setValue(openvdb::Coord(0, 2, 2), 0);
    grid_accessor.setValue(openvdb::Coord(0, 2, 1), 0);
    grid_accessor.setValue(openvdb::Coord(0, 2, 0), 0);
    grid_accessor.setValue(openvdb::Coord(0, 1, 2), 0);
    grid_accessor.setValue(openvdb::Coord(0, 1, 1), 0);
    grid_accessor.setValue(openvdb::Coord(0, 1, 0), 0);
    grid_accessor.setValue(openvdb::Coord(0, 0, 2), 0);
    grid_accessor.setValue(openvdb::Coord(0, 0, 1), 0);
    grid_accessor.setValue(openvdb::Coord(0, 0, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 2, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 1, 0), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 2), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 1), 0);
    grid_accessor.setValue(openvdb::Coord(1, 0, 0), 0);
    grid_accessor.setValue(openvdb::Coord(2, 2, 2), 1);
    grid_accessor.setValue(openvdb::Coord(2, 2, 1), 1);
    grid_accessor.setValue(openvdb::Coord(2, 2, 0), 1);
    grid_accessor.setValue(openvdb::Coord(2, 1, 2), 1);
    grid_accessor.setValue(openvdb::Coord(2, 1, 1), 1);
    grid_accessor.setValue(openvdb::Coord(2, 1, 0), 1);
    grid_accessor.setValue(openvdb::Coord(2, 0, 2), 1);
    grid_accessor.setValue(openvdb::Coord(2, 0, 1), 1);
    grid_accessor.setValue(openvdb::Coord(2, 0, 0), 1);

    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.5, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.6, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.7, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.8, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(0.9, 1.6, 1.5)),  0,       1e-5); 
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1,   1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.1, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.2, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.3, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.4, 1.6, 1.5)),  0,       1e-5);      
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.5, 1.6, 1.5)),  0,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.6, 1.6, 1.5)),  0.02800, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.7, 1.6, 1.5)),  0.10400, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.8, 1.6, 1.5)),  0.21600, 1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(1.9, 1.6, 1.5)),  0.35199, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2,   1.6, 1.5)),  0.5,     1e-5);   
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.1, 1.6, 1.5)),  0.64800, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.2, 1.6, 1.5)),  0.78399, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.3, 1.6, 1.5)),  0.89600, 1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.4, 1.6, 1.5)),  0.97200, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.5, 1.6, 1.5)),  1,       1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.6, 1.6, 1.5)),  0.97200, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.7, 1.6, 1.5)),  0.89600, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.8, 1.6, 1.5)),  0.78399, 1e-5);  
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(2.9, 1.6, 1.5)),  0.64800, 1e-5);
    EXPECT_NEAR(cubic_interpolator.wsSample(openvdb::Vec3d(3,   1.6, 1.5)),  0.5,     1e-5);
}


TEST(Grid, velocity){

    // We store u,v and w components of velocity in the middle of the cell faces, like this (diagram 2d only)
    //              v
    //      o-------o-------o
    //      |       |       |
    // u -->|       |       |<-- u
    //      |       |       |
    //      o-------o-------o
    //              v


    using namespace openvdb;

    Vec3dGrid::Ptr velocity = Vec3dGrid::create(Vec3d(0,0,0));
    Vec3dGrid::Accessor v_access = velocity->getAccessor();

    v_access.setValue(Coord(0,0,0), Vec3d(1,2,3));
    v_access.setValue(Coord(0,1,0), Vec3d(1,2,3));
    v_access.setValue(Coord(0,0,1), Vec3d(1,2,3));
    v_access.setValue(Coord(0,1,1), Vec3d(1,2,3)); 
    v_access.setValue(Coord(1,0,0), Vec3d(3,4,5));
    v_access.setValue(Coord(1,1,0), Vec3d(3,4,5));
    v_access.setValue(Coord(1,0,1), Vec3d(3,4,5));
    v_access.setValue(Coord(1,1,1), Vec3d(3,4,5));        

    EXPECT_EQ(v_access.getValue(Coord(0,0,0)), Vec3d(1,2,3));
    EXPECT_EQ(v_access.getValue(Coord(1,0,0)), Vec3d(3,4,5));
    EXPECT_EQ(v_access.getValue(Coord(0,1,0)), Vec3d(1,2,3));
    EXPECT_EQ(v_access.getValue(Coord(0,0,1)), Vec3d(1,2,3));

    auto transform = openvdb::math::Transform::createLinearTransform(1.0);
    // pretend the velocity values are co-located in the centre of the voxels
    transform->preTranslate(Vec3d(0.5, 0.5, 0.5));
    velocity->setTransform(transform);
    // using the staggered sampler offsets each component by 0.5 making it a staggered grid
    tools::GridSampler<Vec3dGrid, tools::StaggeredBoxSampler> sps(*velocity);

    EXPECT_EQ(sps.wsSample(Vec3d(1, 1, 1)), Vec3d(3,3,4));

}

TEST(General, max){

    using namespace openvdb;

    // https://stackoverflow.com/questions/2837854/initializing-an-object-to-all-zeroes
    // https://en.cppreference.com/cpp/language/value_initialization
    Vec3d a{};
    EXPECT_EQ(a, Vec3d(0,0,0));

    EXPECT_EQ(std::max(Vec3d(1,2,3), Vec3d(3,2,1)), Vec3d(3,2,1));

    Vec3d b(1,2,3);
    std::cout << std::max({b[0],b[1],b[2]});

}


TEST(General, add_vector_scalar){
    openvdb::Vec3d a(1,1,1);

    // vector + scalar works when I would expect it to complain
    EXPECT_EQ(a + 3, openvdb::Vec3d(4,4,4));
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
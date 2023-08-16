#include "grid.h"
#include <gtest/gtest.h>


TEST(GridTests, QuantityTests) {
    std::string filepath("dummy");
    FluidGrid* solver = new FluidGrid(4, 4, 4, 0.2, 1, filepath);

    // define a voxel-centred quantity
    FluidQuantity* u = new FluidQuantity(solver, 0.5, 0.5, 0.5, 7.0);
    
    EXPECT_EQ(solver->getSize(), 64);
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

    u->CopyQuantitytoBuffer();
    u->swap();

    EXPECT_EQ(u->getQuantity(1,0,0), 7);
    EXPECT_EQ(u->getQuantity(0,0,0), 4);
    EXPECT_EQ(u->getQuantity(1,1,1), 10);
    EXPECT_EQ(u->getQuantity(3,3,3), 28);

    // 4 * 4 * 4 * 7 - 3 * 7 + 4 + 10 + 28 = 469
    EXPECT_EQ(u->sum(), 469);
    EXPECT_EQ(u->max(), 28);
    u->swap();
    EXPECT_EQ(u->sum(), 469);
    EXPECT_EQ(u->max(), 28);


    // define a voxel-centred quantity
    FluidQuantity* v = new FluidQuantity(solver, 0.5, 0.5, 0.5, 0);
    for (int k=0; k<4; k++)
        for (int j=0; j<4; j++)
            for (int i=0; i<4; i++){
                float value = i + 4*j + 16*k;
                v->setQuantity(i,j,k, value);
            }
    EXPECT_EQ(v->getQuantity(0,0,0), 0);
    EXPECT_EQ(v->getQuantity(3,3,3), 63);


    // InterpolateLinear takes world-space position as input
    EXPECT_EQ(v->InterpolateLinear(0.125,0.125,0.125), 0);
    EXPECT_EQ(v->InterpolateLinear(0.875,0.875,0.875), 63);
    EXPECT_EQ(v->InterpolateLinear(0.625,0.625,0.625), 42);
    EXPECT_EQ(v->InterpolateLinear(0.375,0.625,0.875), 57);
    EXPECT_EQ(v->InterpolateLinear(0.5,0.5,0.5), 31.5);

    EXPECT_EQ(v->InterpolateCubic(0.125,0.125,0.125), 0);
    EXPECT_EQ(v->InterpolateCubic(0.875,0.875,0.875), 63);
    EXPECT_EQ(v->InterpolateCubic(0.625,0.625,0.625), 42);
    EXPECT_EQ(v->InterpolateCubic(0.375,0.625,0.875), 57);
    EXPECT_EQ(v->InterpolateCubic(0.5,0.5,0.5), 31.5);

    FluidGrid* solver2 = new FluidGrid(80, 80, 100, 0.02, 0.1, filepath);
    solver2->setSolid(new Cuboid(solver, 0.5, 0.5, 0.5, 0.7, 0.1, 0.33, 0.0, 0.0));

    EXPECT_EQ(solver2->getBBox().min(), openvdb::Vec3d(0,0,0));
    EXPECT_EQ(solver2->getBBox().max(), openvdb::Vec3d(1,1,1.25));
    EXPECT_EQ(solver2->getSolid()->PointIsInside(0.5,0.5,0.5), true);
    EXPECT_EQ(solver2->getSolid()->PointIsInside(0.0,0.0,0.0), false);
    EXPECT_EQ(solver2->getSolid()->PointIsInside(0.8,0.52,0.65), true);
    
    EXPECT_EQ(solver2->getU()->indexToWorld(32,32,40), openvdb::Vec3d(0.5,0.5,0.6));



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
    A.setValue(1, 1, 2);
    A.setValue(2, 2, 3);
    A.setValue(63, 63, 4);

    EXPECT_EQ(A.getValue(0, 0), 1);
    EXPECT_EQ(A.getValue(0, 1), 5);
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

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
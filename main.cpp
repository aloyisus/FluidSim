//Copyright (c) 2014 John Kelly
//
//
//
//


#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "grid.h"


int main(int argc, char* argv[]) {

    // Declare the supported options.
    po::options_description desc("Allowed options");
    std::string filepath;
    int width;
    int height;
    int depth;
    float stoptime;
    double density;
    desc.add_options()
        ("help", "produce help message")
        ("output", po::value<std::string>(&filepath), "path to vdb output")
        ("width", po::value<int>(&width)->default_value(80), "Width of the volume. Default is 80.")
        ("height", po::value<int>(&height)->default_value(120), "Height of the volume. Default is 120.")
        ("depth", po::value<int>(&depth)->default_value(80), "Depth of the volume. Default is 80.")
        ("stoptime", po::value<float>(&stoptime)->default_value(2.5), "Time for the simulation to run for.")
        ("density", po::value<double>(&density)->default_value(0.1), "Density of the smoke.")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    if (vm.count("output")) {
        std::cout << "Output path was set to " << filepath << ".\n";
    } else {
        std::cout << "Output path was not set.\n";
        return 1;
    }

    double timestep = 0.02;

    FluidGrid* solver = new FluidGrid(width, height, depth, timestep, density, filepath);
    
    //solver->setSolid(new Cuboid(solver, w/2-10,h/2+10,d/2-10,64,8,32,0.0,0.25));

    solver->setSolid(new Sphere(solver, width/2-10,height/2,depth/2,16));
    
    while (solver->getCurrtime() < stoptime) {
        solver->Update();
    }
    
    return 0;
}
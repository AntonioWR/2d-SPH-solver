// main cpp file

// Writen by Yuan Dong Wang Ruan
// latest update: 24/03/2021

// include the neccesary header file
#include "SPH.h"
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <mpi.h>
using namespace std;

// An alias to reduce typing
namespace po = boost::program_options;

//contents in the bracket to allow processes to communicate with each other.
int main(int argc, char** argv)
{

    // Specify the options we want to make available to the user
    po::options_description opts("Options to solve the SPH problems: ");
    opts.add_options()

        // Options available
        ("ic-dam-break", "Use dam-break initial condition.")
        ("ic-block-drop", "Use block-drop initial condition.")
        ("ic-droplet", "Use droplet initial condition.")
        ("ic-one-particle", "Use one particle validation case initial condition.")
        ("ic-two-particles", "Use two particles validation case initial condition.")
        ("ic-four-particles","Use four particles validation case initial condition.")
        ("dt", po::value<double>()->default_value(0.0001),"Time-step to use.")
        ("T", po::value<double>()->default_value(2), "Total integration time.")
        ("h",po::value<double>()->default_value(0.01),"Radius of influence of each particle.")
        ("help", "Print help message.");

    // Tell Boost to parse the command-line arguments using the list of
    // possible options and generate a map (vm) containing the options and
    // values actually specified by the user.
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    // Check if the user gave the "--help" option and print the usage.
    if(vm.count("help")) {
        cout << "Select the desired initial condition." << endl
             << "Input dt, T, or h if required to modify, otherwise default value will be used." << endl;
        cout << opts << endl;
        return 0;
    }

    // Extract the selection of initial condition with boolean type.
    const bool dam_ic = vm.count("ic-dam-break");
    const bool block_ic = vm.count("ic-block-drop");
    const bool droplet_ic = vm.count("ic-droplet");
    const bool one_ic = vm.count("ic-one-particle");
    const bool two_ic = vm.count("ic-two-particles");
    const bool four_ic = vm.count("ic-four-particles");
    // Extract the values given to other parameters using the appropriate
    // data type.
    const double dt = vm["dt"].as<double>();
    const double T = vm["T"].as<double>();
    const double h = vm["h"].as<double>();

    // Make sure only one initial condition is input otherwise the program will be finished.
    if(dam_ic + block_ic + droplet_ic + one_ic + two_ic + four_ic > 1) {
        cout << "Too many input initial conditions!" << endl;
        return 1;
    } else if(dam_ic + block_ic + droplet_ic + one_ic + two_ic + four_ic < 1) {
        cout << "Please input an initial condition!" << endl;
        return 1;
    }

    // Find out which initial condition was chosen, the variable (dam_ic,block_ic....) was of boolean type
    // and are casted into int, then a loop was used to find out which one is true and is equal to 1.
    // default selection will be 0 to elliminate the compiler warning here, but will be reassigned later.
    int ic_selected=0; 
    int ic_list[6] = { dam_ic, block_ic, droplet_ic, one_ic, two_ic, four_ic };
    for(int i = 0; i < 6; i++) {
        if(ic_list[i] == 1)
            ic_selected = i;
    }

    // initialise MPI, generate local rank number and total number of processes,
    // check if the communicator is functioning.
    MPI_Init(&argc, &argv);
    int rank, no_p, retval_rank, retval_size;
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &no_p);
    if(retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
        cout << " Invalid communicator" << endl;
        return 1;
    }

    // Initialize a single SPH class object to solve the SPh problem with class member functions.
    SPH myparticles;
    
    // APIs:
    // call setPara function to pass in required parameters which is input by the user or default value.
    myparticles.setPara(h, dt, T);
    // call setIC function to pass in the initial condition selected by the user.
    myparticles.setIC(ic_selected);
    // call setMpi function to pass in the local rank number and number of processes to class SPH.
    myparticles.setMpi(rank,no_p);

    // solve the sph problem by the algorithm.
    myparticles.SPHsolver();

    // terminate the parallel execution environment
    MPI_Finalize();
    return 0;
}

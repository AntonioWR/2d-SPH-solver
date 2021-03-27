// SPH class header file

// Writen by Yuan Dong Wang Ruan
// latest update: 24/03/2021

#ifndef SPH_H
#define SPH_H

// define class SPH
class SPH
{
    // The below parameters will be private to class as we want to encapsualte the data
    // such that only member function of SPH can modyfy them.
private:
    // Declare the static parameters which are identical to every SPH object
    static double k;        // gas constant
    static double rho_zero; // resting density
    static double mu;       // viscosity
    static double g;        // gravitational acceleration
    static double e;        // coefficient of restitution

    // Declare the parameters which will be passed in via the setPara member function
    double h;  // radius of influence
    double dt; // timestep
    double T;  // total run time
    int N;     // number of particles
    
    // Declare the parameters which will be passed in via the setMpi member function
    int rank;//local rank number
    int no_p;//total number of process.

    double* aLocation; // declare a pointer which the coordinate of particles will be stored

public:
    // define constructor
    SPH() = default;              // default constructor
    SPH(const SPH& rhs) = delete; // copy constructor is deleted as there's only a single object
    SPH(SPH&& rhs) = delete;      // move constructor is deleted as there's only a single object

    // destructor is not used as momory allocated from the heap will be deallocated explicitly
    // at the end of the member functions
    ~SPH() {};

    // the member functions are public so that they can be called outside the class.
    // declare function "setPara" to pass in the h,dt and T parameters
    void setPara(const double& ph, const double& pdt, const double& pT);

    // declare function "setIC" to pass in the locations of particles to aLocation and number of particles to N
    void setIC(int ic_selected);
    
    // declare function "setMpi" to pass in the local rank number and total number of process
    void setMpi(const int& prank, const int& pno_p);

    // declare function to solve the SQH problem by algorithm
    void SPHsolver();
};



#endif
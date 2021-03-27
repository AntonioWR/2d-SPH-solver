// SPH class cpp file

// Writen by Yuan Dong Wang Ruan
// latest update: 24/03/2021

#include "SPH.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>

#define PI 3.14159265359
using namespace std;

// initialize static variables
double SPH::k = 2000.0;
double SPH::rho_zero = 1000.0;
double SPH::mu = 1.0;
double SPH::g = 9.81;
double SPH::e = 0.5;

// define member functions
// member function setPara to pass in required parameters which is input by the user or default value.
void SPH::setPara(const double& ph, const double& pdt, const double& pT)
{
    h = ph;
    dt = pdt;
    T = pT;
}

// member function setIC to pass in the initial condition selected by the user.
void SPH::setIC(int ic_selected)
{
    // using a switch statement to set up the initial position of particles and
    // return the number of particles N to the class.
    switch(ic_selected) {

        // case 0-the dam break ic
    case 0: {
        double leftmost = 0.0 + h;              // left boundary
        double bottommost = 0.0 + h;            // left boundary
        double xrange = (0.2 - h) - leftmost;   // horizontal length between left and right boundary
        double yrange = (0.2 - h) - bottommost; // vertical length between top and bottom boundary

        int no_of_column = (xrange / h + 1.1); // suppose to be plus 1 but plus 1.1 here to prevent truncation error.
        int no_of_row = (yrange / h + 1.1);    // suppose to be plus 1 but plus 1.1 here to prevent truncation error.
        N = no_of_column * no_of_row;          // total number of particles
        aLocation = new double[2 * N];         // array of location of particles is allocated from the heap

        // initialize the array from the left bottom to top right, from left to right then from bottom to top.
        // i % no_of_column == 0 indicates this particular particles is at the left boundary, if it's not,
        // the x position is equal to the previous one plus the separation between particles.
        for(int i = 0; i < N; ++i) {
            if(i % no_of_column == 0) {
                aLocation[i] = leftmost;
            } else
                aLocation[i] = aLocation[i - 1] + h;
        }

        // Similarly, the y position is added into the array after the x position, in a row major array, it can be seen
        // as e.g. x0  x1  x2  x3  x4  y0  y1  y2  y3  y4
        for(int i = 0; i < N; ++i) {
            // the bottom row's y position is set to zero.
            if(i == 0) {
                aLocation[N + i] = bottommost;
            } else {
                // if at left most boundary, y is increased by h.
                if(i % no_of_column == 0) {
                    aLocation[N + i] = aLocation[N + i - 1] + h;
                } else
                    // within the same row, the y position is the same.
                    aLocation[N + i] = aLocation[N + i - 1];
            }
        }

        // add some noice(maximum of 0.1h)to the position allocation.
        srand(time(0));
        for(int i = 0; i < 2 * N; ++i) {
            aLocation[i] += (double)rand() / RAND_MAX * 0.1 * h;
        }
        break;
    }

        // case 1-the block drop ic
    case 1: {
        double leftmost = 0.1 + h;              // left boundary
        double bottommost = 0.3 + h;            // left boundary
        double xrange = (0.3 - h) - leftmost;   // horizontal length between left and right boundary
        double yrange = (0.6 - h) - bottommost; // vertical length between top and bottom boundary
        int no_of_column = (xrange / h + 1.1);  // suppose to be plus 1 but plus 1.1 here to prevent truncation error
        int no_of_row = (yrange / h + 1.1);     // suppose to be plus 1 but plus 1.1 here to prevent truncation error
        N = no_of_column * no_of_row;           // total number of particles
        aLocation = new double[2 * N];          // array of location of particles is allocated from the heap

        // Similar to case 0, the x/y positions of particles are put into the array aLocation
        for(int i = 0; i < N; ++i) {
            if(i % (no_of_column) == 0) {
                aLocation[i] = leftmost;
            } else
                aLocation[i] = aLocation[i - 1] + h;
        }

        for(int i = 0; i < N; ++i) {
            if(i == 0) {
                aLocation[N + i] = bottommost;
            } else {
                if(i % no_of_column == 0) {
                    aLocation[N + i] = aLocation[N + i - 1] + h;
                } else
                    aLocation[N + i] = aLocation[N + i - 1];
            }
        }

        // add some noice(maximum of 0.1h)to the position allocation.
        srand(time(0));
        for(int i = 0; i < 2 * N; ++i) {
            aLocation[i] += (double)rand() / RAND_MAX * 0.1 * h;
        }

        break;
    }

    // case 2-the droplet ic
    case 2: {
        // using Vogel's model for sunflower seed arrangement
        N = 200; // number of particles is chosen to be this value such
        // that the inter-particle spacing is approximately equal to h.
        aLocation = new double[2 * N];
        double* r = new double[N];
        double* theta = new double[N];

        double gr = (1 + sqrt(5)) / 2; // golden ratio
        double c = 0.1 / sqrt(N);      // factor to make sure the radius is within 0.1.

        // theta and r follows Vogel's model.
        for(int k = 0; k < N; ++k) {
            theta[k] = 2 * PI * k / pow(gr, 2);
            r[k] = c * sqrt(k);
        }

        // transfer from polar cordinates to Cartesian cordinates
        double x, y;
        for(int i = 0; i < N; ++i) {
            x = r[i] * cos(theta[i]) + 0.5;
            y = r[i] * sin(theta[i]) + 0.7;
            aLocation[i] = x;
            aLocation[N + i] = y;
        }
        break;
    }

        // case 3-single particle ic
    case 3: {
        N = 1;
        aLocation = new double[N * 2];

        // explicitly assign x and y value of this particle.
        aLocation[0] = 0.5;
        aLocation[1] = 0.5;

        break;
    }

        // case 4-two particles ic
    case 4: {
        N = 2;
        aLocation = new double[N * 2];

        // explicitly assign x and y value of these particle.
        double aLocation_list[] = { 0.5, 0.5, 0.5, h };
        for(int i = 0; i < N * 2; ++i) {
            aLocation[i] = aLocation_list[i];
        }
        break;
    }

        // case 5-four particles ic
    case 5: {
        N = 4;
        aLocation = new double[N * 2];

        // explicitly assign x and y value of these particle.
        double aLocation_list[] = { 0.505, 0.515, 0.51, 0.5, 0.5, 0.5, 0.45, 0.45 };
        for(int i = 0; i < N * 2; ++i) {
            aLocation[i] = aLocation_list[i];
        }
        break;
    }
    }
}

// member function setMpi to pass in the local rank number and number of processes to class SPH.
void SPH::setMpi(const int& prank, const int& pno_p)
{
    rank = prank;
    no_p = pno_p;
}

// member function SPHsolver where the algorithm is implemented.
void SPH::SPHsolver()
{
    // allocation of task to each process
    int mpi_quotient = N / no_p;
    int mpi_remainder = N % no_p;
    // calculate the local share of number of particles, if task is not evenly divisible, some processes with take one
    // more task than the others.
    if(rank < mpi_remainder) {
        mpi_quotient += 1;
    }

    double m = 1.0;   // the m is set to be one inthe first iteration and will be updated later.
    int size = N * N; // size of array which holds relation between particle i and j.

    // dynamically allocate arrays from the heap, will be deallocated an the end of the SPHsolver.
    // declare some variables.
    // They are outside the time loop as the declaration/allocation should be done only once.
    double* phi_d = new double[size];

    double xdiff, ydiff, q_temp;

    double* grad_phi_px = new double[size];
    double* grad_phi_py = new double[size];
    double* lap_phi_v = new double[size];
    double* vx = new double[N];
    double* vy = new double[N];
    double* vx_diff = new double[size];
    double* vy_diff = new double[size];

    double* rho_i = new double[N];
    double* loc_rho_i = new double[N];
    double* res_rho_i = new double[N];

    double* p_i = new double[N];

    double mp2r;
    int loc_p;
    double* loc_F_p_ix = new double[N];
    double* loc_F_p_iy = new double[N];
    double* F_p_ix = new double[N];
    double* F_p_iy = new double[N];

    double lmmr;
    double* loc_F_v_ix = new double[N];
    double* loc_F_v_iy = new double[N];
    double* F_v_ix = new double[N];
    double* F_v_iy = new double[N];

    double* F_g_iy = new double[N];

    double* ax = new double[N];
    double* ay = new double[N];

    // initialize ofstreeam object which will be used to output energy and position.
    ofstream vPosout("position.txt", ios::out | ios::trunc);
    ofstream vEPout("energy.txt", ios::out | ios::trunc);

    // initialize the boudary condition at where the particles will be bounced back.
    double right_bc = 1 - h;
    double left_bc = 0 + h;
    double top_bc = 1 - h;
    double bottom_bc = 0 + h;

    // fill in the diagonal of matrix phi_d where q=0, this is done outside the time loop as these value
    // on the diagonal are constant. and does not depend on the timestep.
    for(int i = 0; i < N; ++i) {
        phi_d[i * N + i] = 4 / (PI * pow(h, 2));
    }

    // initilize the loop to solve the SPH problems with time increase by dt after every loop.
    double t = 0.0;
    while(t < T) {

        for(int i = 0; i < N - 1; ++i) {
            for(int j = i + 1; j < N; ++j) {
                // xdiff and ydiff are intermediate value which will be used to calculate other inter-partilcles
                // relation.
                // xdiff is the x-component of vector r_ij
                xdiff = aLocation[i] - aLocation[j];
                // ydiff is the y-component of vector r_i
                ydiff = aLocation[i + N] - aLocation[j + N];

                // check wheather the difference in x coordinate and y coordinates are within the radius of influence,
                // if exceed, there's no need to calculate the inter-particle displacement. can save some computational
                // time
                if((abs(xdiff) < h) && (abs(ydiff) < h)) {
                    // q_temp is the displacement between particle i and j normalised by h.
                    q_temp = (double)sqrt(pow(abs(xdiff), 2) + pow(abs(ydiff), 2)) / h;

                    // fill in the matrix used in the density,pressure force and viscous force
                    // the below matrixs are in the format of a upper triangular, storing the corresponding relation
                    // between particle i and j. There is no need to calculate again the relation from particle j to
                    // i as it is related with the one from i to j.
                    //              j =0      j=1       j=2
                    // i=0             x         x           x
                    // i=1             0         x           x
                    // i=2             0         0           x
                    if(q_temp < 1) {
                        phi_d[i * N + j] = 4 / (PI * pow(h, 2)) * pow((1 - pow(q_temp, 2)), 3);
                        grad_phi_px[i * N + j] = -30 / (PI * pow(h, 3)) * pow(1 - q_temp, 2) / q_temp * xdiff;
                        grad_phi_py[i * N + j] = -30 / (PI * pow(h, 3)) * pow(1 - q_temp, 2) / q_temp * ydiff;
                        lap_phi_v[i * N + j] = 40 / (PI * pow(h, 4)) * (1 - q_temp);
                        vx_diff[i * N + j] = vx[i] - vx[j];
                        vy_diff[i * N + j] = vy[i] - vy[j];
                    }
                }
            }
        }

        // complete the lower triangular of phi_d matrix. phi_d_ij=phi_d_ji
        for(int i = 0; i < N - 1; ++i) {
            for(int j = i + 1; j < N; ++j) {
                phi_d[j * N + i] = phi_d[i * N + j];
            }
        }

        // evaluate the rho_i by summing up each row of the phi_d matrix *m
        // sum_rho_i is used to evaluate the m in the first iteration.
        double sum_rho_i = 0;
        for(int i = 0; i < N; ++i) {
            // initialise loc_rho_i element i at the start of each iteration as we want to sum up the contribution from
            // each p.
            loc_rho_i[i] = 0;
            // for rank that gets one more task if particles not evenly divisible, rank < mpi_remainder.
            if(rank < mpi_remainder) {
                for(int p = 0; p < mpi_quotient; ++p) {
                    loc_rho_i[i] +=
                        phi_d[i * N + rank * mpi_quotient + p] * m; // the index of phi_d should be the global index.
                }
            } else {
                for(int p = 0; p < mpi_quotient; ++p) {
                    loc_rho_i[i] += phi_d[i * N + rank * mpi_quotient + p + mpi_remainder] * m;
                }
            }
        }
        // sum up the local sum and broadcast to every processes for further calculation.
        MPI_Allreduce(loc_rho_i, rho_i, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // compute resiprocal of rho_i to avoid float point devision to save time.
        for(int i = 0; i < N; ++i) {
            res_rho_i[i] = 1.0 / rho_i[i];
            sum_rho_i += rho_i[i];
        }

        // calculate pressue of particle i by ideal gas law.
        for(int i = 0; i < N; ++i) {
            p_i[i] = k * (rho_i[i] - rho_zero);
        }

        // complete the grad_phi_p matrix, grad_phi_p_ij=-grad_phi_p_ji as displacement vector is involved.
        for(int i = 0; i < N - 1; ++i) {
            for(int j = i + 1; j < N; ++j) {
                grad_phi_px[j * N + i] = -grad_phi_px[i * N + j];
                grad_phi_py[j * N + i] = -grad_phi_py[i * N + j];
            }
        }

        // calculate pressure force acting on particle i.
        for(int i = 0; i < N; ++i) {
            // initialise F_p_ix element i at the start of each iteration as we want to sum up the contribution from
            // each p.
            loc_F_p_ix[i] = 0;
            loc_F_p_iy[i] = 0;
            // similar idea to rho_i cases, for rank that gets one more task if particles not evenly divisible, rank <
            // mpi_remainder.
            if(rank < mpi_remainder) {
                for(int p = 0; p < mpi_quotient; ++p) {
                    loc_p = rank * mpi_quotient + p;
                    // mp2r-precompute some part of the expression to save time.
                    mp2r = 0.5 * m * (p_i[i] + p_i[loc_p]) * res_rho_i[loc_p];
                    loc_F_p_ix[i] -= grad_phi_px[i * N + loc_p] * mp2r;
                    loc_F_p_iy[i] -= grad_phi_py[i * N + loc_p] * mp2r;
                }
            }

            else {
                for(int p = 0; p < mpi_quotient; ++p) {
                    loc_p = rank * mpi_quotient + p + mpi_remainder;
                    // mp2r-precompute some part of the expression to save time.
                    mp2r = 0.5 * m * (p_i[i] + p_i[loc_p]) * res_rho_i[loc_p];
                    loc_F_p_ix[i] -= grad_phi_px[i * N + loc_p] * mp2r;
                    loc_F_p_iy[i] -= grad_phi_py[i * N + loc_p] * mp2r;
                }
            }
        }
        // sum up the local sum and broadcast for further calculation.
        MPI_Allreduce(loc_F_p_ix, F_p_ix, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(loc_F_p_iy, F_p_iy, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // complete the lap_phi_v,vx_diff ,vy_diff matrix
        // vx_diff and vy_diff is equal to -vx_diff and -vy_diff respectively as they are vectors.
        for(int i = 0; i < N - 1; ++i) {
            for(int j = i + 1; j < N; ++j) {
                lap_phi_v[j * N + i] = lap_phi_v[i * N + j];
                vx_diff[j * N + i] = -vx_diff[i * N + j];
                vy_diff[j * N + i] = -vy_diff[i * N + j];
            }
        }

        // calculate viscous force acting on particle i.
        for(int i = 0; i < N; ++i) {
            // initialise F_p_ix element i at the start of each iteration as we want to sum up the F_p_ix for each p.
            loc_F_v_ix[i] = 0;
            loc_F_v_iy[i] = 0;
            // for rank that gets one more task if particles not evenly divisible, rank < mpi_remainder.
            if(rank < mpi_remainder) {
                for(int p = 0; p < mpi_quotient; ++p) {
                    loc_p = rank * mpi_quotient + p;
                    // lmmr-precompute some part of the expression to save time.
                    lmmr = lap_phi_v[i * N + loc_p] * m * mu * res_rho_i[loc_p];
                    loc_F_v_ix[i] -= vx_diff[i * N + loc_p] * lmmr;
                    loc_F_v_iy[i] -= vy_diff[i * N + loc_p] * lmmr;
                }
            } else {
                for(int p = 0; p < mpi_quotient; ++p) {
                    loc_p = rank * mpi_quotient + p + mpi_remainder;
                    // lmmr-precompute some part of the expression to save time.
                    lmmr = lap_phi_v[i * N + loc_p] * m * mu * res_rho_i[loc_p];
                    loc_F_v_ix[i] -= vx_diff[i * N + loc_p] * lmmr;
                    loc_F_v_iy[i] -= vy_diff[i * N + loc_p] * lmmr;
                }
            }
        }
        // sum up the local sum and broadcast for further calculation.
        MPI_Allreduce(loc_F_v_ix, F_v_ix, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(loc_F_v_iy, F_v_iy, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // calculate the gravity force for particle i.
        for(int i = 0; i < N; ++i) {
            F_g_iy[i] = -rho_i[i] * g;
        }

        // summarising all of the three forces acting on the particles and calculate acceleration vector
        for(int i = 0; i < N; ++i) {
            ax[i] = (F_p_ix[i] + F_v_ix[i]) * res_rho_i[i];
            ay[i] = (F_p_iy[i] + F_v_iy[i] + F_g_iy[i]) * res_rho_i[i];
        }

        // apply a forward euler scheme to update the velocity at next time step, because velocity is calculated
        // at half-steps, use only half of the time to calculate the velocity at the first time step
        for(int i = 0; i < N; ++i) {
            if(t == 0) {
                vx[i] += ax[i] * dt * 0.5;
                vy[i] += ay[i] * dt * 0.5;
            } else {
                vx[i] += ax[i] * dt;
                vy[i] += ay[i] * dt;
            }
        }

        // calculate the position of particles at next time step
        for(int i = 0; i < N; ++i) {
            aLocation[i] += vx[i] * dt;
            aLocation[N + i] += vy[i] * dt;

            // Apply boundary conditions, if the center of particle reaches any boundary+-h depends on whcih boundary
            // it's touching. The position and velocity is updated to reverse the motion of particles..
            // for left and right boundaries.
            if(aLocation[i] > right_bc) {
                aLocation[i] = 1 - h;
                vx[i] = -e * vx[i];
            } else if(aLocation[i] < left_bc) {
                aLocation[i] = 0 + h;
                vx[i] = -e * vx[i];
            }

            // for upper and lower boundaries.
            if(aLocation[N + i] > top_bc) {
                aLocation[N + i] = 1 - h;
                vy[i] = -e * vy[i];
            } else if(aLocation[N + i] < bottom_bc) {
                aLocation[N + i] = 0 + h;
                vy[i] = -e * vy[i];
            }
        }

        // calculate kinetic energy and potential energy
        double E_k = 0;
        double E_p = 0;
        double T_E = 0;

        for(int i = 0; i < N; ++i) {
            E_k += 0.5 * m * (pow(vx[i], 2) + pow(vy[i], 2));
            E_p += m * g * aLocation[N + i];
        }
        T_E = E_k + E_p;

        // output energy.txt with necesary data.
        vEPout << t << " " << E_k << " " << E_p << " " << T_E << endl;

        // scale m so the density is equal to the reference density.
        if(t == 0) {
            m = N * rho_zero / (sum_rho_i);
        }

        // update time
        t = t + dt;
    }

    // output location at the final time.
    for(int i = 0; i < N; ++i) {
        vPosout << aLocation[i] << " " << aLocation[N + i] << endl;
    }

    // deallocate the memory for arrays used in the SPHsolver function.
    delete[] aLocation;
    delete[] phi_d;
    delete[] grad_phi_px;
    delete[] grad_phi_py;
    delete[] lap_phi_v;
    delete[] vx;
    delete[] vy;
    delete[] vx_diff;
    delete[] vy_diff;
    delete[] rho_i;
    delete[] loc_rho_i;
    delete[] res_rho_i;
    delete[] p_i;
    delete[] loc_F_p_ix;
    delete[] loc_F_p_iy;
    delete[] F_p_ix;
    delete[] F_p_iy;
    delete[] loc_F_v_ix;
    delete[] loc_F_v_iy;
    delete[] F_v_ix;
    delete[] F_v_iy;
    delete[] F_g_iy;
    delete[] ax;
    delete[] ay;

    vEPout.close();
    vPosout.close();
}
/*
   md3d_argon.cxx is a basic Molecular Dynamics
   simulation of 3d box with argon
   Copyrigth 2008-2017 by Edmanuel Torres A.
   eetorres@gmail.com
   https://github.com/eetorres/mdargon
*/

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

using namespace std;

typedef long double real;

const int N = 1000;     // number of particles
real r[N][3];           // positions
real v[N][3];           // velocities
real a[N][3];           // acelerations
real dt = 0.001;        // time step

real vMax = 0.8;        // maximun velocity compontet
real rho  = 0.95;       // particle density
real lattice_contant;   // lattice constant from the box side
real L;                 // box size from the density
real L_2;

ofstream argon_xyz;
ofstream argon_energy;
string base_name_xyz = "argon_";

void initialize(void);
void computeAccelerations(void);
real computePotential(real);
real computeForce(real);
void velocityVerlet(real);
real instantaneousTemperature(void);
real potentialEnergy(void);
real kineticEnergy(void);

inline real pbc(real);

// Initial position in the cubic box
void initialize() {
    // time step
    dt = 0.005;
    // number of particles in each coordinate direction
    int n = int(ceil(pow(real(N), (real)(1.0/3) )));
    // Compute the dimension of the box side
    L = pow((real)real(N)/rho,(real)1.0/3.0); 
    L_2 = L/2;
    // lattice parameter
    lattice_contant = L / n;
    // RAmedon seed
    srand((unsigned)time(0));
    // Particle counter
    int p = 0; 
    // initial coordinates
    for (int x = 0; x < n; x++){
        for (int y = 0; y < n; y++){
            for (int z = 0; z < n; z++) {
                if (p < N) {
                    // offset displacement from the perfect cubic lattice
                    real shift = 0.1* (2.0*rand()/real(RAND_MAX) - 1.0);
                    r[p][0] = (x + 0.5 + shift) * lattice_contant - L_2;
                    r[p][1] = (y + 0.5 + shift) * lattice_contant - L_2;
                    r[p][2] = (z + 0.5 + shift) * lattice_contant - L_2;
                }
                ++p;
            }
        }
    }
    // initial velocities
    for (int p = 0; p < N; p++){
        for (int i = 0; i < 3; i++){
            v[p][i] = vMax * (2.0*rand() / real(RAND_MAX) - 1.0);
        }
    }
    // initial accelerations
    computeAccelerations();
}

// Potential energy of a particle i due the particle j to a distance r
real computePotential(real r2){
    real ir2 = 1.0/r2;
    // Lennard-Jones (LJ) potential for argon
    return(4.0 * (pow(ir2, 6) - pow(ir2, 3)) );
}

// Force between particles separated by a distance r
real computeForce(real r2){
    real ir = sqrt(1.0/r2);
    // Gradien of  the LJ potential
    return(24 * (2 * pow(ir, 13) - pow(ir, 7)));
}

// Compute the acelerations
void computeAccelerations() {
    real rij[3], f;
    for (int i = 0; i < N; i++)          // set acelerations to zero
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;

    for (int i = 0; i < N-1; i++){           // loop over all the particles with index (i)
        for (int j = i+1; j < N; j++) {      // sub-loop over the particles with index (j)
            real rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k]; // relatice distance of i relative to j
                rij[k] = pbc(rij[k]);
                rSqd += rij[k] * rij[k];
            }
            f = computeForce(rSqd);
            for (int k = 0; k < 3; k++) {
                 a[i][k] += rij[k] * f;    // components of the forece on i due to j
                 a[j][k] -= rij[k] * f;    // Second Newton law (force on j)
            }
        }
    }
    //for (int i = 0; i < N-1; i++){        // loop over all the particles i
          //cout<<fabs(a[i][0])+fabs(a[i][1])+fabs(a[i][2])<<endl;
    //}
}

// Time integration of the equatrion of motion using the Verlet algorithm
void velocityVerlet(real dt) {
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt; // New positions
            r[i][k] = pbc(r[i][k]);
            v[i][k] += 0.5 * a[i][k] * dt;                     // First step of the velocities
        }
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            v[i][k] += 0.5 * a[i][k] * dt;                     // Second step of the velocities
}

// Compute the instant temperature
real instantTemperature() {
    real sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}

// Compute the potential energy
real potentialEnergy() {
    real rij[3];
    real u_p=0;
    for (int i = 0; i < N-1; i++)           // loop over all the particles i
        for (int j = i+1; j < N; j++) {     // sub-loop over all the particles j
            real rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k]; // distance of i relative to j
                rSqd += rij[k] * rij[k];
            }
            u_p += computePotential(rSqd);
        }
    return u_p;
}

// Compute kinetic energy
real kineticEnergy() {
    real sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum;
}

// Apply periodic boundary conditions (PBC) to the box
inline real pbc(real _x){
    real _x_pbc;
    // Check if the coordinate is within the box limits
    // otherwise apply PBC
    if(_x >= L_2){
      _x_pbc = (_x - L);
    }else if(_x < -L_2){
      _x_pbc = (_x + L);
    }else{ return _x;} // do not change if it is inside
    return _x_pbc;
}

// end


#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <eigen3/Eigen/Dense>

//systems
#include "systems/systemFactory.h"
#include "systems/system.h"

//utils
#include "utils/particle.h"


using namespace Eigen;
using namespace std;
using namespace ParticleUtils;
using namespace Systems;

//config variables
const int PRINT_PRECISION = 15; //number of significative digits
ofstream out("out_gen.txt");  //output file
double h =800; //timestep
const int TIME_OFFSET = 0; //start the integrator from a specific timestep
const int MAX_T = 1E6; //number of integration steps
const int DEBUG_TIMESTEP_MULT = 10000; // 0=DISABLE, print to screen the timestep
const bool EXIT_ON_ERROR = 0; //exit if the error is too high
const double ERROR_THRESHOLD = 0.1; //error threshold
const int PRINT_MULTIPLE = 8; //print out every n timesteps
const int PRINT_OFFSET = 1; //start printing from a specific step
const double ORBIT_NORMALIZE = 50;  //normalize the timesteps 

const int DIM = 4;

Matrix<double,DIM,1> q0(0.050000000000000, 0.00000000000000 ,0.000000000000000 ,0.000390000000000);

int main(int argc, char* argv[]){
  
    cout.setf(std::ios::scientific);
    cout.precision(PRINT_PRECISION);
    out.setf(std::ios::fixed);
    out.precision(PRINT_PRECISION);
    out.setf(ios::showpoint);

    System<DIM>* system = systemFactory<DIM>("aaaaaa");
    Particle<DIM> particle = Particle<DIM>(system);
    particle.initialize(INIT_HAMILTONIAN);

    cout << "time step: " << h << endl;
    cout << "Initialization: " << endl;
    cout << "q0:\t" << particle.q0.transpose() << endl;
    cout << "p0:\t" << particle.p0.transpose() << endl;
  
    // ******
    //MAIN LOOP
    // ******
    for (int t=TIME_OFFSET + 1;t<=MAX_T+TIME_OFFSET;t++){

        if ((DEBUG_TIMESTEP_MULT>0) && (t%DEBUG_TIMESTEP_MULT==0)) cout << "Timestep " << t << endl;

        particle.StepForward();

        if (t==1) {
          cout << "q1:\t" << particle.q1.transpose() << endl;
          cout << "p1:\t" << particle.p1.transpose() << endl;
        }

        //EXIT IF THE ERROR IS TOO HIGH
        if ((EXIT_ON_ERROR) && (abs(particle.Eerr0)>ERROR_THRESHOLD)){
          cout << "Timestep: " << t << "\tError too high! Exiting." << endl;
          break;
        }

        //PRINT TO FILE
        if (((t+PRINT_OFFSET)%PRINT_MULTIPLE)==0){
          out << (t-1) << " " << (t-1)/ORBIT_NORMALIZE << " " << particle.Eerr0
          << " " << particle.q0.transpose() << " " << particle.p0.transpose() 
          << " " << (particle.p0-particle.mom0).transpose();
          out << endl;
        }

    }
  
    return 0;
}
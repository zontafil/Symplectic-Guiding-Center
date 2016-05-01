#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "utils/particle.h"
#include "utils/config.h"

//config variables
const int PRINT_PRECISION = 15; //number of significative digits
ofstream out("out_gen.txt");  //output file
const int DEBUG_TIMESTEP_MULT = 10000; // 0=DISABLE, print to screen the timestep
const int PRINT_MULTIPLE = 8; //print out every n timesteps
const int PRINT_OFFSET = 1; //start printing from a specific step

using namespace std;
using namespace Particles;

void setPrintPrecision(int print_precision);

int main(int argc, char* argv[]){


    setPrintPrecision(PRINT_PRECISION);

    //configuration and create new particle
    Config::Config* config = new Config::Config();
    Particle<Config::DIM> particle = Particle<Config::DIM>(config);

    cout << "time step: " << config->h << endl;
    cout << "Initialization: " << endl;
    cout << "q0:\t" << particle.get_q0().transpose() << endl;
    cout << "p0:\t" << particle.get_p0().transpose() << endl;
  
    // ******
    //MAIN LOOP
    // ******
    for (int t=config->time_offset + 1;t<=config->max_t+config->time_offset;t++){

        if ((DEBUG_TIMESTEP_MULT>0) && (t%DEBUG_TIMESTEP_MULT==0)) cout << "Timestep " << t << endl;
        particle.StepForward();

        if (t==1) {
          cout << "q1:\t" << particle.get_q1().transpose() << endl;
          cout << "p1:\t" << particle.get_p1().transpose() << endl;
        }

        //EXIT IF THE ERROR IS TOO HIGH
        if ((config->exit_on_error) && (abs(particle.Eerr0)>config->error_threshold)){
          cout << "Timestep: " << t << "\tError too high! Exiting." << endl;
          break;
        }

        //PRINT TO FILE
        if (((t+PRINT_OFFSET)%PRINT_MULTIPLE)==0){
            out 
                << (t-1) << " " 
                << (t-1)/config->orbit_normalize << " " 
                << particle.Eerr0 << " " 
                << particle.get_q0().transpose() << " " 
                << particle.get_p0().transpose() << " " 
                << (particle.get_p0()-particle.mom0).transpose();
            out << endl;
        }

    }
  
    return 0;
}

void setPrintPrecision(int print_precision){
    cout.setf(std::ios::scientific);
    cout.precision(print_precision);
    out.setf(std::ios::fixed);
    out.precision(print_precision);
    out.setf(ios::showpoint);
}
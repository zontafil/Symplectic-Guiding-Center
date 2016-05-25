//integrator for numerical simulation of physical systems
// and in particular for the symplectic integration of the non canonical guiding center problem.

//see configInterface.h ,config.h and the readme for configuration

//the format of the output file is:
// timestep (1 column)
// orbit (1 column)
// z (DIM columns)
// conserved quantities (specific for the system. The first column should always be the energy error.)


#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <stdexcept>
#include <vector>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#define BOOST_LOG_DYN_LINK 1

#include "utils/particle.h"
#include "config.h"

using namespace std;
using namespace Particles;

void setPrintPrecision(int print_precision, ofstream& out);
void printToFile(int t,Config::Config* config, Particle<Config::DIM>* particle, ofstream& out);
void setDebugLevel(Config::Config* config);

int main(int argc, char* argv[]){

    //configuration 
    Config::Config* config = new Config::Config();

    //open output file
    ofstream out(config->outFile.c_str()); 
    if (!out.is_open()) throw invalid_argument("Can't open output file for write");

    //set debug level and print precision
    setDebugLevel(config);
    setPrintPrecision(config->print_precision,out);

    //create a particle
    Particle<Config::DIM>* particle = new Particle<Config::DIM>(config);

    cout << "time step: " << config->h << endl;
    cout << "Initialization: " << endl;
    cout << "z_init:\t" << particle->z0.transpose() << endl;
    cout << "z0:\t" << particle->z1.transpose() << endl;
    cout << "conserved qnts: ";

    for (unsigned int i=0; i<particle->conservedQuantities1->size(); i++) cout << particle->conservedQuantities1->at(i) << " ";
    cout << endl;        

    printToFile(0,config,particle, out);
  
    // ******
    //MAIN LOOP
    // ******
    for (int t=config->time_offset + 1;t<=config->max_t+config->time_offset;t++){

        //PRINT TO SCREEN EVERY N STEPS
        if ((config->print_timestep_mult>0) && (t%config->print_timestep_mult==0)) cout << "Timestep " << t << endl;

        particle->StepForward();

        //PRINT TO SCREEN FIRST STEP
        if (t==1) cout << "z1:\t" << particle->z1.transpose() << endl;

        //EXIT IF THE ERROR IS TOO HIGH (ASSUMING ENERGY IS THE FIRST CONSERVED QUANTITY)
        if (config->exit_on_error){
            vector<double>* conservedQuantitiesError = particle->conservedQuantities_err1;
            if ((conservedQuantitiesError->size()>0) && (abs(conservedQuantitiesError->at(0))>config->error_threshold)){
              cout << "Timestep: " << t << "\tError too high! Exiting." << endl;
              break;
            }
        } 

        //PRINT TO FILE
        if (((t+config->print_timestep_offset)%config->file_timestep_mult)==0) printToFile(t,config,particle, out);
    }
  
    return 0;
}

void setPrintPrecision(int print_precision, ofstream& out){
    cout.setf(std::ios::scientific);
    cout.precision(print_precision);
    out.setf(std::ios::fixed);
    out.precision(print_precision);
    out.setf(ios::showpoint);
}

void printToFile(int t,Config::Config* config, Particle<Config::DIM>* particle, ofstream& out){
    out 
        << (t) << " " 
        << (t)/config->orbit_normalize << " " 
        << particle->z1.transpose() << " ";

    vector<double>* conservedQuantitiesError = particle->conservedQuantities_err1;
    for (unsigned int i=0; i< conservedQuantitiesError->size(); i++){
        out << conservedQuantitiesError->at(i) << " ";
    }
    out << endl;

}

void setDebugLevel(Config::Config* config){
   boost::log::add_console_log(std::cout, boost::log::keywords::format = ">> %Message%");
   if (config->debug_level=="debug")  boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::debug);
   else boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::info);
}

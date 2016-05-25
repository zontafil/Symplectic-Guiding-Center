//particle class 

#ifndef PARTICLE_H
#define PARTICLE_H

#include <eigen3/Eigen/Dense>
#include "../integrators/integratorFactory.h"
#include "particleUtils.h"
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace Eigen;
using namespace Integrators;
using namespace Particles;
using namespace std;

namespace Particles{

	template <int DIM> class Particle
	{
		private:
			//initialize the particle, set initial conditions
			void initialize(Config::Config* config);

			//get conserved quantities for the current step (i.e. energy, momenta)
			void buildConservedQuantities();
				
		public:
			Particle(Config::Config* config);
			~Particle(){};

			Matrix<double,DIM,1> z0,z1;

			//conserved quantities at step n,n+1 and difference with step 0
			vector<double> *conservedQuantities_init, *conservedQuantities0, *conservedQuantities1, *conservedQuantities_err0, *conservedQuantities_err1;

			//advance the particle by 1 time step
			void StepForward();

			Integrator<DIM>* integrator;
			Integrator<DIM>* auxiliaryIntegrator; //integrator used for initialization (INIT_LAGRANGIAN --> compute z1 from z0)

			//position getters and setters
			PhaseSpacePoints<DIM> get_z();
			void set_z(PhaseSpacePoints<DIM> z);

			//time step
			double h;
	};

	template <int DIM> Particle<DIM>::Particle(Config::Config* config){
		//build the integrator object and initialize the particle
		integrator = integratorFactory<DIM>(config->integrator,config);

		h = config->h;

		initialize(config);
	}


	template <int DIM> PhaseSpacePoints<DIM> Particle<DIM>::get_z(){
		PhaseSpacePoints<DIM> z;
		z.z0 = z0;
		z.z1 = z1;
		return z;
	};
	template <int DIM> void Particle<DIM>::set_z(PhaseSpacePoints<DIM> z){
		z0 = z.z0;
		z1 = z.z1;
	}



	template <int DIM> void Particle<DIM>::StepForward(){
		//shift values
		z0 = z1;
		conservedQuantities0 = conservedQuantities1;
		conservedQuantities_err0 = conservedQuantities_err1;

		//compute next time step
		z1 = integrator->StepForward(z0, h);

		//compute conserved quantitiesEerr0
		buildConservedQuantities();
	}

	template <int DIM> void Particle<DIM>::initialize(Config::Config* config){

		set_z(config->z);

		//initialize the particle
		initializationType init = config->initialization_type;
		if ((init==INIT_HAMILTONIAN) || (init==INIT_MANUAL) || (init==INIT_MANUAL_MULTISTEP)){
			set_z(integrator->initialize(get_z(),init, h, config));
		}
		else if (init==INIT_LAGRANGIAN){

			//initialize the particle with the help of an auxiliary integrator,
			//in case of the initial conditions are not sufficient for the discrete problem.
			//see PDF, paragraph 6.3

			//TODO: move this code to integrator classes

			//compute q1 from q0 using an auxiliary integrator
			if (config->init_steps<=0) throw invalid_argument("init_steps must be > 0");
			auxiliaryIntegrator = explicitIntegratorFactory<DIM>(config->auxiliary_integrator,config);
			z1 = z0;
			for (int i=0;i<config->init_steps;i++){
			  z1 = auxiliaryIntegrator->StepForward(z1, h/double(config->init_steps));
			}

			//compute z1 from (q0,q1) using main integrator (i.e. legendre right of symplectic integrators)
			set_z(integrator->initialize(get_z(),INIT_MANUAL_MULTISTEP, h, config));
		}
		else throw invalid_argument("Invalid initializatin type");

		//compute initial conserved quantities
		conservedQuantities_init = integrator->system->getConservedQuantities(get_z());
		conservedQuantities_err1 = integrator->system->getConservedQuantities(get_z());
		buildConservedQuantities();
	}

	template <int DIM> void Particle<DIM>::buildConservedQuantities(){
		//get conserved quantities for timestep 1, save and compute errors
		conservedQuantities1 = integrator->system->getConservedQuantities(get_z());
		for (unsigned int i=0;i<conservedQuantities1->size(); i++){
			conservedQuantities_err1->at(i)= ( conservedQuantities1->at(i)- conservedQuantities_init->at(i)) / conservedQuantities_init->at(i);
		}
	}

}

#endif
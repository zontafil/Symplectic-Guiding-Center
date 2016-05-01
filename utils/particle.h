#ifndef PARTICLE_H
#define PARTICLE_H

#include <eigen3/Eigen/Dense>
#include "../integrators/integratorFactory.h"
#include "particleUtils.h"
#include <iostream>

using namespace Eigen;
using namespace Integrators;
using namespace Particles;
using namespace std;

namespace Particles{

	template <int DIM> class Particle
	{
		public:
			Particle(Config::Config* config);
			~Particle(){};

			Matrix<double,DIM,1> q0,q1,p0,p1;
			PositionMomentumPoint<DIM> z0,z1;
			PositionPoints<DIM> q;
			Matrix<double,DIM,1> mom_initial,mom0,mom1;

			double E_initial,E0,E1,Eerr0,Eerr1;

			void initialize(initializationType init_type);
			void StepForward();

			Integrator<DIM>* integrator;
	};


	template <int DIM> Particle<DIM>::Particle(Config::Config* config){
		integrator = integratorFactory<DIM>(config->integrator,config);
	}
	template <int DIM> void Particle<DIM>::StepForward(){
		//shift values
		p0 = p1;
		q0 = q1;
		z0 = z1;
		E0 = E1;
		Eerr0 = Eerr1;
		mom0 = mom1;

		z1 = integrator->StepForward(z0);
		q1 = z1.q;
		p1 = z1.p;
		q.q0 = q0;
		q.q1 = q1;

		E1 = integrator->system->Hamiltonian(z1);
		Eerr1 = ( E1 - E_initial ) / E_initial;
		mom1 = integrator->system->momentum(q);
	}

	template <int DIM> void Particle<DIM>::initialize(initializationType init_type){
		q0 = z0.q;
		p0 = z0.p;
		if (init_type==INIT_HAMILTONIAN){
			q.q0 = q0;
			p0 = integrator->system->momentum(q);
			z0.p = p0;
			z0.q = q0;
			E0 = integrator->system->Hamiltonian(z0);
			E_initial = E0;
			mom0 = integrator->system->momentum(q);

			q.q1 = q0;
			z1.p = p0;
			z1.q = q0;
			E1 = E0;
			mom1 = mom0;
		}
	}

}

#endif
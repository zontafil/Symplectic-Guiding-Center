#ifndef PARTICLE_H
#define PARTICLE_H

#include <eigen3/Eigen/Dense>
#include "../integrators/integratorFactory.h"
#include "particleUtils.h"
#include <iostream>
#include <stdexcept>

using namespace Eigen;
using namespace Integrators;
using namespace Particles;

namespace Particles{

	template <int DIM> class Particle
	{
		private:
			Matrix<double,DIM,1> q0,q1,p0,p1;
			void initialize(Config::Config* config);
				
		public:
			Particle(Config::Config* config);
			~Particle(){};

			Matrix<double,DIM,1> mom_initial,mom0,mom1;
			double E_initial,E0,E1,Eerr0,Eerr1;

			void StepForward();

			Integrator<DIM>* integrator;
			Integrator<DIM>* auxiliaryIntegrator; //integrator used for initialization

			Matrix<double,DIM,1> get_q0(){return q0;};
			Matrix<double,DIM,1> get_q1(){return q1;};
			Matrix<double,DIM,1> get_p0(){return p0;};
			Matrix<double,DIM,1> get_p1(){return p1;};
			PositionMomentumPoint<DIM> get_z0();
			PositionMomentumPoint<DIM> get_z1();
			PositionPoints<DIM> get_q();
			PositionMomentumTwoPoints<DIM> get_z();

			void set_z0(PositionMomentumPoint<DIM> z0);
			void set_z1(PositionMomentumPoint<DIM> z1);
			void set_z(PositionMomentumTwoPoints<DIM> z);

			double h;
	};

	template <int DIM> Particle<DIM>::Particle(Config::Config* config){
		integrator = integratorFactory<DIM>(config->integrator,config);

		h = config->h;

		initialize(config);
	}


	template <int DIM> PositionMomentumPoint<DIM> Particle<DIM>::get_z0(){
		PositionMomentumPoint<DIM> z0;
		z0.q = q0;
		z0.p = p0;
		return z0;
	};
	template <int DIM> PositionMomentumPoint<DIM> Particle<DIM>::get_z1(){
		PositionMomentumPoint<DIM> z1;
		z1.q = q1;
		z1.p = p1;
		return z1;
	};
	template <int DIM> PositionPoints<DIM> Particle<DIM>::get_q(){
		PositionPoints<DIM> q;
		q.q0 = q0;
		q.q1 = q1;
		return q;
	}
	template <int DIM> PositionMomentumTwoPoints<DIM> Particle<DIM>::get_z(){
		PositionMomentumTwoPoints<DIM> z;
		z.q0 = q0;
		z.q1 = q1;
		z.p0 = p0;
		z.p1 = p1;
		return z;
	}
	template <int DIM> void Particle<DIM>::set_z0(PositionMomentumPoint<DIM> z0){
		q0 = z0.q;
		p0 = z0.p;
	}
	template <int DIM> void Particle<DIM>::set_z(PositionMomentumTwoPoints<DIM> z){
		q0 = z.q0;
		p0 = z.p0;
		q1 = z.q1;
		p1 = z.p1;
	}
	template <int DIM> void Particle<DIM>::set_z1(PositionMomentumPoint<DIM> z1){
		q1 = z1.q;
		p1 = z1.p;
	}

	template <int DIM> void Particle<DIM>::StepForward(){
		//shift values
		p0 = p1;
		q0 = q1;
		E0 = E1;
		Eerr0 = Eerr1;
		mom0 = mom1;

		PositionMomentumPoint<DIM> z1 = integrator->StepForward(get_z0(), h);
		set_z1(z1);
		E1 = integrator->system->Hamiltonian(z1);
		Eerr1 = ( E1 - E_initial ) / E_initial;
		mom1 = integrator->system->momentum(get_q());
	}

	template <int DIM> void Particle<DIM>::initialize(Config::Config* config){

		set_z(config->z);

		initializationType init = config->initialization_type;
		if ((init==INIT_HAMILTONIAN) || (init==INIT_MANUAL_POSITION_MOMENTUM) || (init==INIT_MANUAL_MULTISTEP)){
			set_z(integrator->initialize(get_z(),init, h));
		}
		else if (init==INIT_LAGRANGIAN){

			//compute q1 from q0 using an auxiliary integrator
			if (config->init_steps<=0) throw invalid_argument("init_steps must be > 0");
			auxiliaryIntegrator = integratorFactory<DIM>(config->auxiliary_integrator,config);

			set_z1(get_z0());
			for (int i=0;i<config->init_steps;i++){
			  set_z1(auxiliaryIntegrator->StepForward(get_z1(),h/double(config->init_steps)));
			}

			//compute z1 from (q0,q1) using main integrator (i.e. legendre right of symplectic integrators)
			set_z(integrator->initialize(get_z(),INIT_MANUAL_MULTISTEP, h));
		}
		else throw invalid_argument("Invalid initializatin type");

		//compute initial energy, momentum
		E1 = integrator->system->Hamiltonian(get_z1());
		E_initial = E1;
		mom1 = integrator->system->momentum(get_q());
		
	}

}

#endif
#ifndef PARTICLE_H
#define PARTICLE_H

#include <eigen3/Eigen/Dense>
#include "../integrators/integratorFactory.h"
#include "particleUtils.h"
#include <iostream>

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

			Matrix<double,DIM,1> get_q0(){return q0;};
			Matrix<double,DIM,1> get_q1(){return q1;};
			Matrix<double,DIM,1> get_p0(){return p0;};
			Matrix<double,DIM,1> get_p1(){return p1;};
			PositionMomentumPoint<DIM> get_z0();
			PositionMomentumPoint<DIM> get_z1();
			PositionPoints<DIM> get_q();

			void set_z0(PositionMomentumPoint<DIM> z0);
			void set_z1(PositionMomentumPoint<DIM> z1);
	};

	template <int DIM> Particle<DIM>::Particle(Config::Config* config){
		integrator = integratorFactory<DIM>(config->integrator,config);

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
	template <int DIM> void Particle<DIM>::set_z0(PositionMomentumPoint<DIM> z0){
		q0 = z0.q;
		p0 = z0.p;
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

		PositionMomentumPoint<DIM> z1 = integrator->StepForward(get_z0());
		set_z1(z1);
		E1 = integrator->system->Hamiltonian(z1);
		Eerr1 = ( E1 - E_initial ) / E_initial;
		mom1 = integrator->system->momentum(get_q());
	}

	template <int DIM> void Particle<DIM>::initialize(Config::Config* config){
		q0 = config->z0.q;
		if (config->initialization_type==INIT_HAMILTONIAN){
			q1 = q0;
			p1 = integrator->system->momentum(get_q());
			E1 = integrator->system->Hamiltonian(get_z1());
			E_initial = E1;
			mom1 = integrator->system->momentum(get_q());
		}
	}

}

#endif
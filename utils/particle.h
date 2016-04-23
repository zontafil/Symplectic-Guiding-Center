#ifndef PARTICLE_H
#define PARTICLE_H

#include <eigen3/Eigen/Dense>
#include "../systems/system.h"
#include "particleUtils.h"

using namespace Eigen;
using namespace Systems;
using namespace ParticleUtils;

typedef enum {
	INIT_HAMILTONIAN,
	INIT_LAGRANGIAN,
	INIT_MANUAL_POSITION_MOMENTUM,
	INIT_MANUAL_MULTISTEP
} initializationType;


template <int DIM> class Particle
{
	public:
		Particle();
		Particle(System<DIM>* system);
		Particle(Matrix<double,DIM,1> q0);
		~Particle();

		Matrix<double,DIM,1> q0,q1,p0,p1;
		PositionMomentumPoint<DIM> z0,z1;
		Matrix<double,DIM,1> mom_initial,mom0,mom1;

		double E_initial,E0,E1,Eerr0,Eerr1;

		void initialize(initializationType init_type);
		void StepForward();
	private:
		System<DIM> *_system;
};


template <int DIM> Particle<DIM>::Particle(System<DIM>* system){
	this->_system = system;	
}
// template <int DIM> Particle<DIM>::Particle(Matrix<double,DIM,1> q0){
// 	this->q0 = q0;
// }
template <int DIM> Particle<DIM>::~Particle(){}
template <int DIM> void Particle<DIM>::StepForward(){
	//shift q,p values
	this->p0 = this->p1;
	this->q0 = this->q1;
	this->z0 = this->z1;

	this->z1 = this->_system->StepForward(this->z0);
	this->q1 = this->z1.q;
	this->p1 = this->z1.p;
}

template <int DIM> void Particle<DIM>::initialize(initializationType init_type){
}


#endif
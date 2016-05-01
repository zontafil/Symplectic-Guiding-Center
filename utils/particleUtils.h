#ifndef PARTICLEUTILS_H
#define PARTICLEUTILS_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace Particles{

	typedef enum {
		INIT_HAMILTONIAN,
		INIT_LAGRANGIAN,
		INIT_MANUAL_POSITION_MOMENTUM,
		INIT_MANUAL_MULTISTEP
	} initializationType;

	template <int DIM> struct PositionPoints{
		Matrix<double,DIM,1> q0;
		Matrix<double,DIM,1> q1;
	};
	template <int DIM> struct PositionMomentumPoint{
		Matrix<double,DIM,1> q;
		Matrix<double,DIM,1> p;
	};
	template <int DIM> struct MomentumPoints{
		Matrix<double,DIM,1> p0;
		Matrix<double,DIM,1> p1;
	};

};

#endif
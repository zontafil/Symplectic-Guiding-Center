#ifndef PARTICLEUTILS_H
#define PARTICLEUTILS_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace Particles{

	//initialization type.
	//see paragraph 6.3 of PDF.

	//MULTISTEP: compute initial conditions (i.e. q0,p0 from q0,q1)
	typedef enum {
		INIT_HAMILTONIAN,
		INIT_LAGRANGIAN,
		INIT_MANUAL,
		INIT_MANUAL_MULTISTEP
	} initializationType;

	template <int DIM> struct PhaseSpacePoints{
		Matrix<double,DIM,1> z0;
		Matrix<double,DIM,1> z1;
	};

	//variables for canonical hamiltonian systems
	template <int DIM> struct PositionPoints{
		Matrix<double,DIM/2,1> q0;
		Matrix<double,DIM/2,1> q1;
	};
	template <int DIM> struct PositionMomentumPoint{
		Matrix<double,DIM/2,1> q;
		Matrix<double,DIM/2,1> p;
	};
	template <int DIM> struct PositionMomentumTwoPoints{
		Matrix<double,DIM/2,1> q0;
		Matrix<double,DIM/2,1> p0;
		Matrix<double,DIM/2,1> q1;
		Matrix<double,DIM/2,1> p1;
	};
	template <int DIM> struct MomentumPoints{
		Matrix<double,DIM/2,1> p0;
		Matrix<double,DIM/2,1> p1;
	};

};

#endif
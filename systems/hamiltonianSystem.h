//Canonical Hamiltonian system.

// DIM must be even

//Must implement an Hamiltonian H(q,p)
//		and a discretization of the continuous momentum p(q,dq)
//			--> useful with degenerate systems, when p = p(q), to compute hamiltonian init
//			--> unnecessary when hamiltonian init is not used	

#ifndef HAMILTONIANSYSTEM_H
#define HAMILTONIANSYSTEM_H
#include "../utils/particleUtils.h"
#include "../config.h"
#include <vector>

using namespace Particles;
using namespace std;

namespace Systems{
	template <int DIM> class HamiltonianSystem : public System<DIM>{
		public:
			HamiltonianSystem(Config::Config* config): System<DIM>(config), conservedQuantities(1+DIM/2){};
			~HamiltonianSystem(){};

			//hamiltonian, lagrangian and discretization of momentum (p(q0,q1))
			virtual double Hamiltonian(PositionMomentumPoint<DIM> z) = 0;
			virtual Matrix<double,DIM/2,1> momentum(PositionPoints<DIM> q) = 0; //this is useful for degenerate systems, when p=p(q)
			virtual double Lagrangian(Matrix<double,DIM/2,1> q, Matrix<double,DIM/2,1> v) = 0;

			const int conservedQuantities;
			vector<double>* getConservedQuantities(PhaseSpacePoints<DIM> z);

	};

	template <int DIM> vector<double>* HamiltonianSystem<DIM>::getConservedQuantities(PhaseSpacePoints<DIM> z){

		vector<double> *ret = new vector<double>(conservedQuantities);

		PositionMomentumPoint<DIM> z1;
		z1.q = z.z1.head(DIM/2);
		z1.p = z.z1.tail(DIM/2);
		PositionPoints<DIM> q;
		q.q0 = z.z1.head(DIM/2);

		//first conserved quantity is the hamiltonian
		ret->at(0) = Hamiltonian(z1);

		Matrix<double,DIM/2,1> mom = momentum(q);
		
		for (int i=0;i<DIM/2;i++){
			//compute (continuous momentum) - p1
			ret->at(i+1) = mom(i) - z.z1(i+DIM/2);
		}	

		return ret;
	}

}

#endif

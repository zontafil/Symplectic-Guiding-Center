#ifndef SYSTEM_H
#define SYSTEM_H
#include "../utils/particleUtils.h"
#include "../utils/config.h"

using namespace Particles;

namespace Systems{
	template <int DIM> class System{
		public:
			// System();
			System(Config::Config* config){};
			~System(){};

			virtual double Hamiltonian(PositionMomentumPoint<DIM> z) = 0;
			virtual Matrix<double,DIM,1> momentum(PositionPoints<DIM> q) = 0;
	};

}

#endif

#ifndef SYSTEM_H
#define SYSTEM_H
#include "../utils/particleUtils.h"
#include "../config.h"
#include <vector>

using namespace Particles;
using namespace std;

namespace Systems{
	template <int DIM> class System{
		public:
			System(Config::Config* config): conservedQuantities(0){};
			~System(){};

			const int conservedQuantities;

			virtual Matrix<double,DIM,1> f_eq_motion(Matrix<double,DIM,1> z) = 0;
			virtual vector<double>* getConservedQuantities(PhaseSpacePoints<DIM> z) = 0;
	};

}

#endif

## Symplectic Guiding Center

Integrator written in c++ for numerical simulations of physical systems, focused on the variational symplectic integration of the guiding center.

Other interators/system are implemented/possible.

#### warning

This code is not optimized and very, very, very slow. Don't try or think to use it in a production environment.

### Documentation

For a theoretical background and an explanation of the different integrators, please read the thesis (doc/thesis.pdf) or the comments inside each class.

### Prerequisites

* Boost and Boost.log
* gsl
* eigen3

### Compile and run

Make sure you have a valid configuration (config.h) file (see the following section)

      make
      ./bin/symplectic

if you don't have make, just use:

      g++ -DBOOST_LOG_DYN_LINK -Wall -o bin/symplectic symplectic.cpp -lgsl -lgslcblas -lboost_log -lboost_log_setup -lpthread -lboost_thread

### Configuration

See utils/configInterface.h for the structure of a valid configuration.

A valid config class must be put in ./config.h. To start, just copy an example configuration in examples/ to ./config.h.

The main configurable entities are Integrator (config variables "integrator", "auxiliary_integrator", "first_guess_integrator"), System (config variable "system") and EMField (config variable "emField").

#### Valid integrators:

* RK4 (explicit)
* VariationalMidpoint (implicit)
* SymplecticExplicit1 (explicit, system must be "GuidingCenter", DIM must be 8)
* SymplecticExplicit2 (explicit, system must be "GuidingCenter", DIM must be 8)
* SymplecticExplicit3 (explicit, system must be "GuidingCenter", DIM must be 8)
* SymplecticExplicit4 (explicit, system must be "GuidingCenter", DIM must be 8)
* SymplecticImplicit1 (implicit, system must be "GuidingCenter", DIM must be 8)
* symplecticSemiexplicitQin (implicit, system must be "GuidingCenter", DIM must be 8)
* SymplecticImplicit3D (implicit, system must be "GuidingCenter", DIM must be 6)

#### Valid systems

* GuidingCenter

#### Valid EM fields

* Tokamak
* ForceFree
* TwoDimField
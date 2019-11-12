#ifndef SPLINEFIELD_H
#define SPLINEFIELD_H

#include "guidingfield.h"
#include <stdexcept>
#include "./eqdskReader/eqdskReader.h"
#include "./eqdskReader/eqdskPsiInterp.h"
#include "ascot5-spline/common_spline.h"

namespace EMFields{

	template <int DIM> class SplineField: public GuidingFieldConfiguration<DIM>
	{
		public:
			SplineField(Config::Config* config);
			~SplineField(){};

			//compute the field from q
			GuidingField compute(Matrix<double,DIM/2,1> q);
	};

	template <int DIM> SplineField<DIM>::SplineField(Config::Config* config): GuidingFieldConfiguration<DIM>(config) {
	};

    template <int DIM> GuidingField SplineField<DIM>::compute(Matrix<double,DIM/2,1> q){
        eqdsk eqdsk_neq = readEqdskFile("./data/neqdsk_66832");

        interp2D_data* interp2Dc = eqdskPsiInterp(eqdsk_neq);
        interp1D_data* interp1Dc = eqdskFpolInterp(eqdsk_neq);

        // TODO
        // BdB_rz Bdata = evalBrz(2, 2, interp2Dc, interp1Dc, eqdsk_neq);

        GuidingField ret;
        return ret;
    }




}

#endif
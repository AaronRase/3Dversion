#ifndef TEMPLAT_LATTICE_IO_HDF5_XIFUNCTIONS_H
#define TEMPLAT_LATTICE_IO_HDF5_XIFUNCTIONS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2020


#include "TempLat/util/prettytostring.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/parameters/parameterparser.h"

#include "TempLat/lattice/IO/HDF5/helpers/hdf5file.h"
/////////////////////////////////////////////////////////////////////////////////
#include "TempLat/util/getcpptypename.h"

namespace TempLat {


    /** \brief A class which implements saving in pure HDF5.
     *
     *
     * Unit test: make test-filesaverhdf5
     **/

    

    class XiFunctions {
    public:
        /* Put public methods here. These should change very little over time. */

        template<typename R>
        double data(R r, std::vector<ptrdiff_t> coords){ //used to store an entity directly to a dataset, using it's own name.
            ConfirmSpace::apply(r,r.getToolBox()->mLayouts.getConfigSpaceLayout(), SpaceStateInterface::SpaceType::Configuration);
            GhostsHunter::apply(r);
            return Data(r,r.getToolBox()->mNDimensions-1,coords);
        }


        //To save our fields, we use the fact that the last dimension is not parallelised.
        //We iterate over the first N-1 dimensions, and for each of these we save the whole
        //last dimension to file.
        template<typename R>
        double Data(R r, int dim,  std::vector<ptrdiff_t> crds)
        {
		
			auto toolBox = r.getToolBox();
		/*	std::vector<ptrdiff_t> coords;
			//for(size_t i = crds.size()-1; i > 0; i--) coords.emplace_back( crds[i]);
			for(size_t i = 0; i <crds.size()-1 ; i++) coords.emplace_back( crds[i]);
			double endpoint = crds[crds.size()-1];
			
			coords.emplace_back(0);
			ptrdiff_t offset = r.getJumps().getTotalOffsetFromSpatialCoordinates(coords);

			sayMPI << " coords: " << coords;
			sayMPI << " offset: " << offset;*/

			ptrdiff_t offset = r.getJumps().getTotalOffsetFromSpatialCoordinates(crds);
			
			return r.get(offset);
			           
        }

    //std::vector<std::vector<double>> alldata;
    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */

    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };



} /* TempLat */

#ifdef TEMPLATTEST
#include "TempLat/lattice/IO/HDF5/xifunctions_test.h"
#endif



#endif

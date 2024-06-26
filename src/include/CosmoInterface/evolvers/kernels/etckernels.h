#ifndef COSMOINTERFACE_EVOLVERS_KERNELS_ETCKERNELS_H
#define COSMOINTERFACE_EVOLVERS_KERNELS_ETCKERNELS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Adrien Florio, Nicolás Layza,  Year: 2022

#include "TempLat/util/tdd/tdd.h"
#include "CosmoInterface/definitions/potential.h"
#include "CosmoInterface/definitions/PITensor.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"

namespace TempLat {


    /** \brief A class which stores the kernel for the GWs fields.
     *
     * 
     * Unit test: make test-gwskernels
     **/


    class ETCKernels {
    public:
        /* Put public methods here. These should change very little over time. */
        ETCKernels() = delete;

        template <class Model, int N>
        static auto get(Model& model, Tag<N> n){
        	return (PITensor::totalTensor(model,n));
        }
		
		
    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */



    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };



} /* FCN */

#ifdef TEMPLATTEST
#include "CosmoInterface/evolvers/kernels/etckernels_test.h"
#endif


#endif

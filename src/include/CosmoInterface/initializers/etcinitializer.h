#ifndef COSMOINTERFACE_INITIALIZERS_ETCINITIALIZER_H
#define COSMOINTERFACE_INITIALIZERS_ETCINITIALIZER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Adrien Florio, Nicolás Layza,  Year: 2022

#include "TempLat/util/rangeiteration/for_in_range.h"
#include "CosmoInterface/initializers/fluctuationsgenerator.h"

namespace TempLat {


    /** \brief A class which is used to initialize the gravitational waves fields.
     *
     *
     * Unit test: make test-scalarinitializer
     **/

    class ETCInitializer {
    public:
        /* Put public methods here. These should change very little over time. */
        ETCInitializer() = delete;

        template<class Model>
        static void initializeETC(Model& model)
        {

            ForLoop(i, 0, Model::NGWs - 1, (*model.ETC)(i)= 0.;);

        }


    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */

    };

    struct ETCInitializerTester{
      #ifdef TEMPLATTEST
              static inline void Test(TDDAssertion& tdd);
      #endif
    };


} /* TempLat */

#ifdef TEMPLATTEST
#include "CosmoInterface/initializers/etcinitializer_test.h"
#endif


#endif

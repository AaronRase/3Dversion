#ifndef TEMPLAT_LATTICE_MEASUREMENTS_RADIALPROJECTORETC_H
#define TEMPLAT_LATTICE_MEASUREMENTS_RADIALPROJECTORETC_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Adrien Florio, Nicolás Layza,  Year: 2022

#include <algorithm>

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionresult.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglequantity.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialbincomputer.h"
#include "TempLat/lattice/algebra/operators/squareroot.h"


#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getvalue.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/algebra/spacestateinterface.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"

#include "CosmoInterface/definitions/ETCProjector.h"

// This class is based on the TempLat/lattice/algebra/measuringtools/radialprojector.h class.
// It is a reimplementation of the radial projection which allows to compute the projected
// GW spectrum point by point, circumventing the exprission templates for efficiency.


namespace TempLat {


    template <typename T = double>
    class RadialProjectorETC {
    public:

        typedef typename GetGetReturnType<T>::type vType;
        typedef typename GetFloatType<vType>::type sType;

        static constexpr bool isComplexValued = GetGetReturnType<T>::isComplex;
        typedef typename RadialProjectionResult<sType>::floatType floatType;
        typedef RadialProjectionResult<sType> resultType;

        template<class Model>
        RadialProjectorETC(Model& model, bool pUseCentralBinValues, int& PRJType, int& PSType, int& PSVersion) :
        mToolBox(model.getToolBox()),
        mUseBinCentralValues(pUseCentralBinValues),
        mPRJType(PRJType),
        mPSType(PSType),
        mPSVersion(PSVersion),
        etcProjector(mPRJType)
        {
        }

        template<class Model>
        RadialProjectionResult<sType> measure(Model& model, ptrdiff_t nLinearBins = -1, sType customRange = -1, bool excludeOrigin = true) {
            if ( nLinearBins < 0 ) {
                ptrdiff_t nGrid = getLayout().getGlobalSizes()[0];
                ptrdiff_t nDim =  getLayout().getGlobalSizes().size();
                nLinearBins = std::pow(nGrid, std::max((ptrdiff_t)1, nDim - 1));
            }


            RadialProjectionResult<sType> baseWorkSpace(nLinearBins, mUseBinCentralValues,  true);

            sType minValue = excludeOrigin ? 1.0 : 0.0;

            RadialProjectionResult<sType> myResult = computeETCSpectrum(model, makeBinComputer(nLinearBins, minValue, customRange), baseWorkSpace, excludeOrigin);

            myResult.finalize(mToolBox->mGroup.getBaseComm());

            return myResult;
        }

    private:
        std::shared_ptr<MemoryToolBox> mToolBox;
        bool mUseBinCentralValues;
        int mPRJType;
        int mPSType;
        int mPSVersion;
        ETCProjector etcProjector;


        const auto& getLayout()
        {
          return  mToolBox->mLayouts.getFourierSpaceLayout();
        }


        template <class Model, typename BINCOMPUTETYPE>
        auto computeETCSpectrum(Model& model, BINCOMPUTETYPE binComputer, RadialProjectionResult<sType> baseWorkSpace, bool excludeOrigin ) {

            auto it = mToolBox->itP();
            //Trying to fix multiplicity for ETC
            auto layout = getLayout();
            HermitianRedundancy quality;
        
           for(it.begin();it.end();++it)
            {
            	 //bool amIRoot = amIRoot(model.getToolBox());
            	
                if((not excludeOrigin) or (not it.isAtOrigin())) {
                quality=layout.getHermitianPartners()->qualify(it.getVec());
                 
                if ( quality != HermitianRedundancy::negativePartner ) {
                    
                    sType r = rFromCoords(it.getVec());

                    ptrdiff_t bin = binComputer(r);
                    
                    /* don't over-weight the real-valued entries: only one float value, only half the weight. */
                    floatType weight = quality == HermitianRedundancy::realValued ? 0.5 : 1;
                    
                    if (mPSVersion != 3) baseWorkSpace.add(bin, etcProjector.getProjectionatpoint(model, it), r, weight);
                    else {
                        WaveNumber ntilde(model.getOneField().getToolBox());
                        if (mPSType == 2) baseWorkSpace.add(bin, etcProjector.getProjectionatpoint(model, it) * GetValue::get(pow<3>(ntilde.norm()), it()), r, weight);
                        else baseWorkSpace.add(bin, etcProjector.getProjectionatpoint(model, it) * GetValue::get(ntilde.norm(), it()), r, weight);
                    	}
				    }
                }
            }
            binComputer.setCentralBinBounds(baseWorkSpace.getCentralBinBounds());

            return baseWorkSpace;
        }

        /** \brief Creates the lambda that maps the IterationCoordinates to a bin. */
        inline auto makeBinComputer(ptrdiff_t nLinearBins, sType minValue, sType customRange = -1) {
            auto layout = getLayout();

            auto rMax = customRange < 0 ? layout.getMaxRadius() : customRange;

            return RadialBinComputer(minValue, rMax, nLinearBins);
        }

        floatType rFromCoords(const std::vector<ptrdiff_t> & pCoord) {
            using namespace std;
            floatType r2 = 0;
            for ( auto&& it : pCoord ){
              r2 += it * it;
            //  say << it;
            }
            return sqrt(r2);
        }
    };

    template <class Model, typename T =double>
    RadialProjectorETC<T> projectRadiallyETC(Model& model, bool useBinCentralValues, int PRJType, int PSType, int PSVersion ) {
        return RadialProjectorETC<T>(model, useBinCentralValues, PRJType, PSType, PSVersion);
    }


    struct RadialProjectorETCTester {
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };
}

#ifdef TEMPLATTEST
#include "CosmoInterface/measurements/radialprojectorETC_test.h"
#endif


#endif

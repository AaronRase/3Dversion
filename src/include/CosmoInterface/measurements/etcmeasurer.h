#ifndef COSMOINTERFACE_MEASUREMENTS_ETCMEASURER_H
#define COSMOINTERFACE_MEASUREMENTS_ETCMEASURER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Adrien Florio, Nicolás Layza,  Year: 2022

#include "CosmoInterface/measurements/meansmeasurer.h"
#include "CosmoInterface/measurements/measurementsIO/spectrumsaver.h"
#include "CosmoInterface/measurements/measurementsIO/spectrumGWsaver.h"
#include "CosmoInterface/measurements/measurementsIO/filesmanager.h"
#include "CosmoInterface/measurements/measurementsIO/measurementssaver.h"
#include "TempLat/util/templatvector.h"

#include "CosmoInterface/definitions/ETCProjector.h"
// #include "CosmoInterface/definitions/checkTT.h"

#include "CosmoInterface/measurements/etcpowerspectrum.h"





namespace TempLat {


    /** \brief A class which contains standard measurements for the ETC.
     *
     *
     * Unit test: make test-etcsmeasurer
     **/

    template <typename T>
    class ETCMeasurer {
    public:
        /* Put public methods here. These should change very little over time. */
        template <typename Model>
        ETCMeasurer(Model& model, FilesManager& filesManager, const RunParameters<T>& par, bool append):
        PSType(par.powerSpectrumType)
 		{

            bool amIRoot = model.getToolBox()->amIRoot();
            // We create a file containing the spectra
        	spectraOut.emplace_back(
                			SpectrumGWSaver<T>(filesManager, "etc", amIRoot, append, par, !model.ETC)
                                    ); // File for spectra
            // standardOut.emplace_back( MeasurementsSaver<T>(filesManager, "TT_Test", amIRoot, append, getTTHeaders(model)) ); // File for TT_check

        }

        // template <typename Model, typename Check>
        // void measureStandard(Model& model, T t, Check& TestTransTrace) {
        //     standardOut(0).addAverage(t);
        //     standardOut(0).addAverage(TestTransTrace.checkTrans(model, 1_c));
        //     standardOut(0).addAverage(TestTransTrace.checkTrans(model, 2_c));
        //     standardOut(0).addAverage(TestTransTrace.checkTrans(model, 3_c));
        //     standardOut(0).addAverage(TestTransTrace.checkTrace(model));
        //     standardOut(0).save();
        // }

        template <typename Model>
        void measureSpectra(Model& model, T t, ETCPowerSpectrumMeasurer& ETCPSMeasurer) {
                   if (model.ETC != nullptr) spectraOut(0).save(t, model, ETCPSMeasurer.etcpowerSpectrum(model));


        }

     	// Returns string with the header of the strin parameters file
        // template <typename Model>
        // std::vector<std::string> getTTHeaders(Model& model) const
        // {
        //    std::vector<std::string> ret;
        //    ret.emplace_back("t");
        //    ret.emplace_back("lambda_1");
        //    ret.emplace_back("lambda_2");
        //    ret.emplace_back("lambda_3");
        //    ret.emplace_back("delta");

        //    return ret;
        // }



    private:

        /* Put all member variables and private methods here. These may change arbitrarily. */

        TempLatVector<SpectrumGWSaver<T>> spectraOut;
        TempLatVector<MeasurementsSaver<T>> standardOut;

        const int PSType;
    };


    class ETCMeasurerTester{
    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };


} /* FCN */

#ifdef TEMPLATTEST
#include "CosmoInterface/measurements/etcmeasurer_test.h"
#endif


#endif

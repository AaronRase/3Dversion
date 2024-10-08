#ifndef COSMOINTERFACE_MEASUREMENTS_SPECTRUMETCSAVER_H
#define COSMOINTERFACE_MEASUREMENTS_SPECTRUMETCSAVER_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Adrien Florio, Nicolás Layza,  Year: 2022

#include "CosmoInterface/measurements/measurementsIO/std/spectrumETCsaverstd.h"


namespace TempLat {

    /** \brief A class which saves spectra to files.
     *
     * Unit test: make test-spectrumsaver
     **/

    template<typename T>
    class SpectrumETCSaver {
    public:
        /* Put public methods here. These should change very little over time. */

        SpectrumETCSaver(FilesManager& fm, std::string fn, bool amIRoot, bool appendMode, const RunParameters<T>& rPar, bool dontCreate = false) :
        useHDF5(fm.getUseHDF5Spectra())
        {
            if(not dontCreate) {
                saverStd = std::make_shared<SpectrumETCSaverStd<T>>(fm, fn, amIRoot, appendMode, rPar);
            }
        }

        SpectrumETCSaver(FilesManager& fm, const Field<T>& fld, bool amIRoot, bool appendMode, const RunParameters<T>& rPar,  bool dontCreate = false) :
        useHDF5(fm.getUseHDF5Spectra())
        {
            saverStd = std::make_shared<SpectrumETCSaverStd<T>>(fm, fld, amIRoot, appendMode, rPar);
        }


        template<template<typename> class... Spectra, class Model>
        void save(T& t, Model& model, Spectra<T>... spectra){
            saverStd->save(std::vector<std::shared_ptr<RadialProjectionResult<T>>> {std::make_shared<RadialProjectionResult<T>>(spectra)...}, t, model);
        }


    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */

        bool useHDF5;
        std::shared_ptr<SpectrumETCSaverStd<T>> saverStd;
#ifdef HDF5
        std::shared_ptr<SpectrumSaverHDF5<T>> saverHDF5;
#endif

    };

    struct SpectrumETCSaverTester{
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };


} /* TempLat */

#ifdef TEMPLATTEST
#include "CosmoInterface/measurements/measurementsIO/spectrumETCsaver_test.h"
#endif


#endif

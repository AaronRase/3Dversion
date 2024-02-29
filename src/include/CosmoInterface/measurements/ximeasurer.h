#ifndef COSMOINTERFACE_MEASUREMENTS_XIMEASURER_H
#define COSMOINTERFACE_MEASUREMENTS_XIMEASURER_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "CosmoInterface/runparameters.h"
#include "CosmoInterface/measurements/meansmeasurer.h"
#include "CosmoInterface/measurements/measurementsIO/filesmanager.h"
#include "TempLat/util/templatvector.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/algebra/helpers/getngrid.h"
#include "CosmoInterface/definitions/energies.h"
#include "CosmoInterface/definitions/hubbleconstraint.h"
#include "CosmoInterface/definitions/gaugederivatives.h"
#include "CosmoInterface/definitions/fieldfunctionals.h"
//////////////////////////////////////////////////////////////////////////////////
#include "TempLat/util/getcpptypename.h"
#include "TempLat/lattice/algebra/helpers/getvalue.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "CosmoInterface/definitions/potential.h"
#include "TempLat/util/tdd/tdd.h"
#include "CosmoInterface/measurements/xifunctions.h"
#include <cstdlib> 
#include <time.h> 


namespace TempLat {

    /** \brief A class which contains measurements of energies and scale factor.
     *
     **/

    template <typename T>
    class XiMeasurer {
    public:
        /* Put public methods here. These should change very little over time. */
        template <typename Model>
        XiMeasurer(Model& model, FilesManager& filesManager, const RunParameters<T>& par, bool append) :
               amIRoot( model.getToolBox()->amIRoot()),
               N(par.N),
                xi(filesManager, "xi", amIRoot, append, getXiHeaders(model)) // Output file for volume-average energies.
                {
        }

        template <class Model>
        void measure(Model& model, T t, std::string ximethod)
        {
            xi.addAverage(t);  // add to file
            
            
            if(model.getToolBox()->mNDimensions == 3){
            if(ximethod == "1"){
                double Empi = rhoDW(model);
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Empi);
                xi.addAverage(Etotal);
            };
            if(ximethod == "2"){
                double Areampi = area(model);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Areampi);
                xi.addAverage(Areatotal);
            };
            if(ximethod == "3"){
                double linempi = line(model);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(linempi))/model.getToolBox()->getNProcesses();
                xi.addAverage(linetotal);
            }
            if(ximethod == "12"){
                std::vector<double> VVmpi = func2(model); 
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[0]);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[1]);
                xi.addAverage(Etotal);
                xi.addAverage(Areatotal);
            }
            if(ximethod == "13"){
                double Empi = rhoDW(model);
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Empi);
                double linempi = line(model);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(linempi))/model.getToolBox()->getNProcesses();
                xi.addAverage(Etotal);
                xi.addAverage(linetotal);
            }
            if(ximethod == "23"){
                double Areampi = area(model);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Areampi);
                double linempi = line(model);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(linempi))/model.getToolBox()->getNProcesses();
                xi.addAverage(Areatotal);
                xi.addAverage(linetotal);
            }
            if(ximethod == "123"){
                std::vector<double> VVmpi = all(model); 
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[0]);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[1]);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[2]))/model.getToolBox()->getNProcesses();
                xi.addAverage(Etotal);
                xi.addAverage(Areatotal);
                xi.addAverage(linetotal);
            }
            };

            if(model.getToolBox()->mNDimensions == 2){
            if(ximethod == "1"){
                double Empi = rhoDW2D(model);
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Empi);
                xi.addAverage(Etotal);
            };
            if(ximethod == "2"){
                double Areampi = area2D(model);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Areampi);
                xi.addAverage(Areatotal);
            };
            if(ximethod == "3"){
                double linempi = line2D(model);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(linempi))/model.getToolBox()->getNProcesses();
                xi.addAverage(linetotal);
            }
            if(ximethod == "12"){
                std::vector<double> VVmpi = func22D(model); 
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[0]);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[1]);
                xi.addAverage(Etotal);
                xi.addAverage(Areatotal);
            }
            if(ximethod == "13"){
                double Empi = rhoDW2D(model);
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Empi);
                double linempi = line2D(model);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(linempi))/model.getToolBox()->getNProcesses();
                xi.addAverage(Etotal);
                xi.addAverage(linetotal);
            }
            if(ximethod == "23"){
                double Areampi = area2D(model);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(Areampi);
                double linempi = line2D(model);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(linempi))/model.getToolBox()->getNProcesses();
                xi.addAverage(Areatotal);
                xi.addAverage(linetotal);
            }
            if(ximethod == "123"){
                std::vector<double> VVmpi = all2D(model); 
                double Etotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[0]);
                double Areatotal = model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[1]);
                double linetotal = (model.getToolBox()->mGroup.getBaseComm().computeAllSum(VVmpi[2]))/model.getToolBox()->getNProcesses();
                xi.addAverage(Etotal);
                xi.addAverage(Areatotal);
                xi.addAverage(linetotal);
            }
            };

            xi.save();
        };

    template<class Model>
    double rhoDW(Model& model)
    {
            auto toolBox = model.getToolBox();
            auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
            auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.

            bool test;
            auto phi = model.fldS(0_c);
            double width = 0.761594;
            double totE = 0;
           
            auto grad = Energies::gradientS(model,FieldFunctionals::grad2S(model,0_c));
            auto kin = Energies::kineticS(model,FieldFunctionals::pi2S(model,0_c));
            auto pot = Potential::potential(model);
            auto tot = grad + kin + pot;

            for(int nx = starts[0]; nx < starts[0]+sizes[0]; nx++){
                for(int ny = starts[1]; ny < starts[1]+sizes[1]; ny++){
                    for(int nz = starts[2]; nz < starts[2]+sizes[2]; nz++){
                        if(abs(phi(test,{nx,ny,nz})) < width){
                            totE += xifunction.data(tot, {nx,ny,nz});
                        };
                    };
                };
            };

            double rho = 1/(pow(N,3)) * totE;
            return rho;
    };


    template<class Model>
    double area(Model& model)
    {
        auto toolBox = model.getToolBox();
        auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
        auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.
        
        bool test;
        auto phi = model.fldS(0_c);
        double Area = 0;

        auto DX = GaugeDerivatives::forwardGradient(model, 0_c, 1_c);
        auto DY = GaugeDerivatives::forwardGradient(model, 0_c, 2_c);
        auto DZ = GaugeDerivatives::forwardGradient(model, 0_c, 3_c);
        auto grad2 = FieldFunctionals::grad2S(model, 0_c);


         for(int nx = starts[0]; nx < starts[0]+sizes[0]; nx++){
                for(int ny = starts[1]; ny < starts[1]+sizes[1]; ny++){
                    for(int nz = starts[2]; nz < starts[2]+sizes[2]; nz++){
                        
                        if((xifunction.data(DX,{nx,ny,nz})/phi(test,{nx,ny,nz}))*model.dx + 1 < 0){

                            Area += pow(model.dx, 2) * sqrt(xifunction.data(grad2,{nx,ny,nz}))/(abs(xifunction.data(DX,{nx,ny,nz})) + abs(xifunction.data(DY,{nx,ny,nz})) + abs(xifunction.data(DZ,{nx,ny,nz})));

                        };
                        if((xifunction.data(DY,{nx,ny,nz})/phi(test,{nx,ny,nz}))*model.dx + 1 < 0){
                            
                            Area += pow(model.dx, 2) * sqrt(xifunction.data(grad2,{nx,ny,nz}))/(abs(xifunction.data(DX,{nx,ny,nz})) + abs(xifunction.data(DY,{nx,ny,nz})) + abs(xifunction.data(DZ,{nx,ny,nz})));
                            
                        };
                        if((xifunction.data(DZ,{nx,ny,nz})/phi(test,{nx,ny,nz}))*model.dx + 1 < 0){
                            
                            Area += pow(model.dx, 2) * sqrt(xifunction.data(grad2,{nx,ny,nz}))/(abs(xifunction.data(DX,{nx,ny,nz})) + abs(xifunction.data(DY,{nx,ny,nz})) + abs(xifunction.data(DZ,{nx,ny,nz})));
                
                        };
                    };
                };
            };
        
        
        return Area;
    };

   /* template<class Model>
    double line(Model& model)
    {
        bool test;
        auto phi = model.fldS(0_c);
        double amountDW = 0;
        srand(time(0));
       for(int ii=0; ii<500; ii++){
        int xyzplane = rand() % 3;
        int g1 = (rand() % N);
        int g2 = (rand() % N);

        if(xyzplane == 0){
            for(int coord = 0; coord < N; coord++){
                if(((coord < (N-1)) && (phi(test, {coord, g1, g2})/phi(test, {coord+1, g1, g2}) < 0)) || ((coord == (N-1)) && (phi(test, {coord, g1, g2})/phi(test, {0, g1, g2}) < 0))){
                    amountDW++;
                };
            };
        } else if (xyzplane == 1){
            for(int coord = 0; coord < N; coord++){
                if(((coord < (N-1)) && (phi(test, {g1, coord, g2})/phi(test, {g1, coord+1, g2}) < 0)) || ((coord == (N-1)) && (phi(test, {g1, coord, g2})/phi(test, {g1, 0, g2}) < 0))){
                    amountDW++;
                };
            };
        } else{
            for(int coord = 0; coord < N; coord++){
                if(((coord < (N-1)) && (phi(test, {g1, g2, coord})/phi(test, {g1, g2, coord+1}) < 0)) || ((coord == (N-1)) && (phi(test, {g1, g2, coord})/phi(test, {g1, g2, 0}) < 0))){
                    amountDW++;
                };
            };
        }
       };
        

        return amountDW/500;
    };*/

    template<class Model>
    double line(Model& model)
    {
        auto toolBox = model.getToolBox();
        auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
        auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.

        bool test;
        auto phi = model.fldS(0_c);
        double amountDW = 0;
        srand(time(0));
        
       for(int ii=0; ii<500; ii++){
            int g1 = starts[0] + (rand() % sizes[0]);
            int g2 = starts[1] + (rand() % sizes[1]);

            for(int coord = 0; coord < N; coord++){
                if(((coord < (N-1)) && (phi(test, {g1, g2, coord})/phi(test, {g1, g2, coord+1}) < 0)) || ((coord == (N-1)) && (phi(test, {g1, g2, coord})/phi(test, {g1, g2, 0}) < 0))){
                    amountDW++;
                };
            };
        }
        

        return amountDW/500;
    };

    template<class Model>
    std::vector<double> func2(Model& model)
    {
            auto toolBox = model.getToolBox();
            auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
            auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.
            
            bool test;
            auto phi = model.fldS(0_c);
            double width = 0.761594;
            double totE = 0;
            double Area = 0;

            auto DX = GaugeDerivatives::forwardGradient(model, 0_c, 1_c);
            auto DY = GaugeDerivatives::forwardGradient(model, 0_c, 2_c);
            auto DZ = GaugeDerivatives::forwardGradient(model, 0_c, 3_c);
            auto grad2 = FieldFunctionals::grad2S(model, 0_c);
           
            auto grad = Energies::gradientS(model,FieldFunctionals::grad2S(model,0_c));
            auto kin = Energies::kineticS(model,FieldFunctionals::pi2S(model,0_c));
            auto pot = Potential::potential(model);
            auto tot = grad + kin + pot;

            for(int nx = starts[0]; nx < starts[0]+sizes[0]; nx++){
                for(int ny = starts[1]; ny < starts[1]+sizes[1]; ny++){
                    for(int nz = starts[2]; nz < starts[2]+sizes[2]; nz++){
                        if(abs(phi(test,{nx,ny,nz})) < width){
                            totE += xifunction.data(tot, {nx,ny,nz});
                        };
                        if((xifunction.data(DX,{nx,ny,nz})/phi(test,{nx,ny,nz}))*model.dx + 1 < 0){
                            
                            Area += pow(model.dx, 2) * sqrt(xifunction.data(grad2,{nx,ny,nz}))/(abs(xifunction.data(DX,{nx,ny,nz})) + abs(xifunction.data(DY,{nx,ny,nz})) + abs(xifunction.data(DZ,{nx,ny,nz})));

                        };
                        if((xifunction.data(DY,{nx,ny,nz})/phi(test,{nx,ny,nz}))*model.dx + 1 < 0){
                            
                            Area += pow(model.dx, 2) * sqrt(xifunction.data(grad2,{nx,ny,nz}))/(abs(xifunction.data(DX,{nx,ny,nz})) + abs(xifunction.data(DY,{nx,ny,nz})) + abs(xifunction.data(DZ,{nx,ny,nz})));;
                            
                        };
                        if((xifunction.data(DZ,{nx,ny,nz})/phi(test,{nx,ny,nz}))*model.dx + 1 < 0){
                            
                            Area += pow(model.dx, 2) * sqrt(xifunction.data(grad2,{nx,ny,nz}))/(abs(xifunction.data(DX,{nx,ny,nz})) + abs(xifunction.data(DY,{nx,ny,nz})) + abs(xifunction.data(DZ,{nx,ny,nz})));;
                
                        };
                    };
                };
            };

            double rho = 1/(pow(N,3)) * totE;
            return {rho, Area};
    };

   

    template<class Model>
    std::vector<double> all(Model& model)
    {
        std::vector<double> VECTOR;
        VECTOR = func2(model);
        VECTOR.emplace_back(line(model));
        return VECTOR;
    }


    template<class Model>
    double rhoDW2D(Model& model)
    {
            
            auto toolBox = model.getToolBox();
            auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
            auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.
            
            bool test;
            auto phi = model.fldS(0_c);
            
            double width = 0.761594;
            double totE = 0;
           
            auto grad = Energies::gradientS(model,FieldFunctionals::grad2S(model,0_c));
            auto kin = Energies::kineticS(model,FieldFunctionals::pi2S(model,0_c));
            auto pot = Potential::potential(model);
            auto tot = grad + kin + pot;

            
            for(int nx = starts[0]; nx < starts[0]+sizes[0]; nx++){
                for(int ny = starts[1]; ny < starts[1]+sizes[1]; ny++){
                    double val = abs(phi(test,{nx,ny}));
                         if( val < width){
                          totE += xifunction.data(tot, {nx,ny});
                    };
                };
            };
            

            double rho = 1/(pow(N,2)) * totE;
            return rho;
    };


    template<class Model>
    double area2D(Model& model)
    {
        auto toolBox = model.getToolBox();
        auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
        auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.
        
        bool test;
        auto phi = model.fldS(0_c);
        double Area = 0;

        auto DX = GaugeDerivatives::forwardGradient(model, 0_c, 1_c);
        auto DY = GaugeDerivatives::forwardGradient(model, 0_c, 2_c);
        auto grad2 = FieldFunctionals::grad2S(model, 0_c);
        

         for(int nx = starts[0]; nx < starts[0]+sizes[0]; nx++){
                for(int ny = starts[0]; ny < starts[0]+sizes[0]; ny++){
                        
                        if((xifunction.data(DX,{nx,ny})/phi(test,{nx,ny}))*model.dx + 1 < 0){

                            Area += model.dx * sqrt(xifunction.data(grad2,{nx,ny}))/(abs(xifunction.data(DX,{nx,ny})) + abs(xifunction.data(DY,{nx,ny})) );

                        };
                        if((xifunction.data(DY,{nx,ny})/phi(test,{nx,ny}))*model.dx + 1 < 0){
                        
                            Area += model.dx * sqrt(xifunction.data(grad2,{nx,ny}))/(abs(xifunction.data(DX,{nx,ny})) + abs(xifunction.data(DY,{nx,ny})));
                            
                        };
                    
                };
            };
        
        
        return Area;
    };

    /* template<class Model>
    double line2D(Model& model)
    {
        
        bool test;
        auto phi = model.fldS(0_c);
        double amountDW = 0;
        srand(time(0));
       for(int ii=0; ii<500; ii++){
        int xyzplane = rand() % 2;
        int g1 = (rand() % N);
        
        if(xyzplane == 0){
            for(int coord = 0; coord < N; coord++){
                if(((coord < (N-1)) && (phi(test, {coord, g1})/phi(test, {coord+1, g1}) < 0)) || ((coord == (N-1)) && (phi(test, {coord, g1})/phi(test, {0, g1}) < 0))){
                    amountDW++;
                };
            };
        } else {
            for(int coord = 0; coord < N; coord++){
                if(((coord < (N-1)) && (phi(test, {g1, coord})/phi(test, {g1, coord+1}) < 0)) || ((coord == (N-1)) && (phi(test, {g1, coord})/phi(test, {g1, 0}) < 0))){
                    amountDW++;
                };
            };
        } 
       };
        

        return amountDW/500;
    };*/

    template<class Model>
    double line2D(Model& model)
    {
        auto toolBox = model.getToolBox();
        auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
        auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.

        bool test;
        auto phi = model.fldS(0_c);
        double amountDW = 0;
        srand(time(0));
       for(int ii=0; ii<50; ii++){
           int g1 = starts[0] + (rand() % sizes[0]);
       
           for(int coord = 0; coord < N; coord++){
                if(((coord < (N-1)) && (phi(test, {g1, coord})/phi(test, {g1, coord+1}) < 0)) || ((coord == (N-1)) && (phi(test, {g1, coord})/phi(test, {g1, 0}) < 0))){
                    amountDW++;
                };
            } 
       };
        

        return amountDW/50;
    };


    template<class Model>
    std::vector<double> func22D(Model& model)
    {
            auto toolBox = model.getToolBox();
            auto starts =  toolBox->mLayouts.getConfigSpaceStarts(); //Local mpi offset.
            auto sizes =  toolBox->mLayouts.getConfigSpaceSizes(); //Local mpi sizes.
            
            bool test;
            auto phi = model.fldS(0_c);
            double width = 0.761594;
            double totE = 0;
            double Area = 0;

            auto DX = GaugeDerivatives::forwardGradient(model, 0_c, 1_c);
            auto DY = GaugeDerivatives::forwardGradient(model, 0_c, 2_c);
            auto grad2 = FieldFunctionals::grad2S(model, 0_c);
           
            auto grad = Energies::gradientS(model,FieldFunctionals::grad2S(model,0_c));
            auto kin = Energies::kineticS(model,FieldFunctionals::pi2S(model,0_c));
            auto pot = Potential::potential(model);
            auto tot = grad + kin + pot;

            for(int nx = starts[0]; nx < starts[0]+sizes[0]; nx++){
                for(int ny = starts[1]; ny < starts[1]+sizes[1]; ny++){
                        if(abs(phi(test,{nx,ny})) < width){
                            totE += xifunction.data(tot, {nx,ny});
                        };
                        if((xifunction.data(DX,{nx,ny})/phi(test,{nx,ny}))*model.dx + 1 < 0){
                            
                            Area += model.dx * sqrt(xifunction.data(grad2,{nx,ny}))/(abs(xifunction.data(DX,{nx,ny})) + abs(xifunction.data(DY,{nx,ny})) );

                        };
                        if((xifunction.data(DY,{nx,ny})/phi(test,{nx,ny}))*model.dx + 1 < 0){
                            
                            Area += model.dx * sqrt(xifunction.data(grad2,{nx,ny}))/(abs(xifunction.data(DX,{nx,ny})) + abs(xifunction.data(DY,{nx,ny})) );
                            
                        };
                };
            };

            double rho = 1/(pow(N,2)) * totE;
            return {rho, Area};
    };

   

    template<class Model>
    std::vector<double> all2D(Model& model)
    {
        std::vector<double> VECTOR;
        VECTOR = func22D(model);
        VECTOR.emplace_back(line2D(model));

        return VECTOR;
    }
        
       

    private:
    
     	// Returns string with the header of the energies file.
        template <typename Model>
        std::vector<std::string> getXiHeaders(Model& model) const 
        {
            std::vector<std::string> ret;
            ret.emplace_back("t");

            ret.emplace_back("xi");

            return ret;
        }



        /* Put all member variables and private methods here. These may change arbitrarily. */
        const bool amIRoot;
        const int N;
        MeasurementsSaver<T> xi;
        XiFunctions xifunction;



    };

    class XiMeasurerTester{
    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };



} /* TempLat */

#ifdef TEMPLATTEST
#include "CosmoInterface/measurements/energiesmeasurer_test.h"
#endif


#endif

/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include <fstream>
/* This file is generated from the cellml files provided at github.com/chaste/cellml */
#include "beeler_reuter_model_1977Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  void TestTusscherSimulation()
  {
#ifdef CHASTE_CVODE
      
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    boost::shared_ptr<AbstractCvodeCell> p_model(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	
    const double period = 1000;
    const double start_time = p_regular_stim->GetStartTime();
    const double duration   = p_regular_stim->GetDuration();
      
    p_regular_stim->SetPeriod(period);
    p_model->SetTolerances(1e-12,1e-12);
    p_model->SetMaxSteps(1e5);
    //unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
      
    double sampling_timestep = 1;
    unsigned int paces  = 10000;
    OdeSolution current_solution;
    std::ofstream apd_file;
    std::string username = std::string(getenv("USER"));
    
    //apd_file.open("/tmp/"+username+"/apd90plot.ssv");
    /*Set the output to be as precise as possible */
    //apd_file.precision(17);
    /*Print variable names on the first line*/
    /*	for(unsigned int i = 0; i < state_variable_names.size(); i++){
	variables_file << state_variable_names[i] << " ";
	}
	variables_file << "\n";
	for(int i=0; i < steps; i++){
	if(i>0){
	for(unsigned int j=0; j < state_variables->size(); j++){
	variables_file << (*state_variables)[j] << " ";
	}
	variables_file << "\n";
	}
    */
    for(unsigned int i = 0; i < paces; i++){
      /*Set the initial values to be the terminal values of the last solution*/
      current_solution = p_model->Compute(0, start_time, sampling_timestep);
      current_solution = p_model->Compute(start_time, start_time + duration, sampling_timestep);
      current_solution = p_model->Compute(start_time+duration, period, sampling_timestep);
    }
    /*Calculate and output apd to file*/
    /*std::vector<double> voltages = current_solution.GetVariableAtIndex(voltage_index);
    CellProperties cell_props = CellProperties(voltages, current_solution.rGetTimes());
    double apd = cell_props.GetLastActionPotentialDuration(90);
    apd_file << apd << " ";
    apd_file << "\n";
    apd_file.close();
	
    current_solution.WriteToFile("TestCvodeCells","final_trace","ms");
    */	
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

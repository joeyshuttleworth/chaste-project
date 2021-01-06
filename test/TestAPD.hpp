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
#include "SimulationTools.hpp"
#include <boost/filesystem.hpp>
#include <fstream>

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  void TestTusscherSimulation()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;    
    boost::shared_ptr<AbstractCvodeCell> p_model(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
    const double period = 500;
   
    const double duration   = p_regular_stim->GetDuration();
    const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
    std::ofstream output_file;
    std::cout << "Testing model: " + model_name + "\n";

    p_regular_stim->SetPeriod(period);
    p_model->SetTolerances(1e-12, 1e-12);
    p_model->SetMaxSteps(1e5);
    p_model->SetMaxTimestep(1000);
    p_regular_stim->SetStartTime(0);
    
    unsigned int paces  = 1000;
    OdeSolution current_solution;
    std::vector<std::vector<double>> state_variables;
    std::string username = std::string(getenv("USER"));
      
    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;

    boost::filesystem::create_directory("/tmp/"+username);
    boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
    output_file.open("/tmp/"+username+"/"+model_name+"/apdplot.dat");

    TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat"), 0);

    std::vector<double> initial_conditions = MakeStdVec(p_model->rGetStateVariables());
    
    /*Set the output to be as precise as possible */
    output_file.precision(18);
    
    TS_ASSERT_EQUALS(output_file.is_open(), true);
    
    const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();
   
    unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
    std::vector<double> times;
    std::vector<double> sampling_timesteps = {1, 0.1, 0.01, 0.001};
    
    for(unsigned int i = 0; i < sampling_timesteps.size(); i++){
      p_model->SetStateVariables(initial_conditions);
      for(unsigned int j = 0; j < paces; j++){
	current_solution = p_model->Compute(0, duration, sampling_timesteps[i]);
	state_variables = current_solution.rGetSolutions();
	times = current_solution.rGetTimes();
	current_solution = p_model->Compute(duration, period, sampling_timesteps[i]);
	std::vector<std::vector<double>> current_state_variables = current_solution.rGetSolutions();
	std::vector<double> current_times = current_solution.rGetTimes();
	state_variables.insert(state_variables.end(), current_state_variables.begin(), current_state_variables.end());
	times.insert(times.end(), current_times.begin(), current_times.end());
	const std::vector<double> voltages = GetNthVariable(state_variables, voltage_index);
	CellProperties cell_props = CellProperties(voltages, times); 
	double current_apd90 = cell_props.GetLastActionPotentialDuration(90);
	output_file << current_apd90 << " ";
      }
      output_file << "\n";
    }
    output_file.close();
  }
#else
  std::cout << "Cvode is not enabled.\n";
#endif
};
  

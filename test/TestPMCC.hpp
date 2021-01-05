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

#include <boost/filesystem.hpp>
#include <fstream>
#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include "SimulationTools.hpp"
#include <limits>
/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
private:
  const unsigned int buffer_size = 100;
 
  std::vector<std::vector<double>> LogDifferences(std::vector<double> &values){
    std::vector<double> x_vals, y_vals;
    x_vals.reserve(buffer_size);
    y_vals.reserve(buffer_size);
    for(unsigned int i=0; i<values.size()-1; i++){
      double tmp  = abs(values[i] - values[i+1]);
      if(tmp != 0){
	tmp = log(tmp);
	y_vals.push_back(tmp);
	x_vals.push_back(i);
      }
    }
    return std::vector<std::vector<double>>({x_vals, y_vals});
  }
  
public:
  void TestTusscherSimulation()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
    
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    
    std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);
    
    for(unsigned int i = 0; i < 4; i++){
      double period = 1000;
      if(i<4)
	period = 500;

      boost::shared_ptr<AbstractCvodeCell> p_model = models[i];
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
      
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);

      const double duration   = p_regular_stim->GetDuration();

      p_regular_stim->SetPeriod(period);
      p_regular_stim->SetStartTime(0);
      p_model->SetTolerances(1e-8, 1e-8);
      p_model->SetMaxSteps(1e5);
      
      const unsigned int paces = 5000;
      OdeSolution current_solution;
      std::ofstream output_file, mrms_file;
      
      const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();

      if(period == 500){
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
      }
      else{
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth2Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
      }

      std::cout << "Testing " << model_name << " with period " << period << "\n";
      p_model->SetMaxTimestep(1000);
      if(period==500){
	output_file.open("/tmp/"+username+"/"+model_name+"/pmcc1Hz2Hz.dat");
	mrms_file.open("/tmp/"+username+"/"+model_name+"/mrms1Hz2Hz.dat");
      }
      else{
	output_file.open("/tmp/"+username+"/"+model_name+"/pmcc2Hz1Hz.dat");
	mrms_file.open("/tmp/"+username+"/"+model_name+"/mrms2Hz1Hz.dat");
      }
      TS_ASSERT_EQUALS(output_file.is_open(), true);
      
      /*Set the output to be as precise as possible */
      output_file.precision(18);

      std::vector<std::string> state_names = p_model->GetSystemInformation()->rGetStateVariableNames();

      std::cout << "Number of state variable names: " << state_names.size() << "\n";
      
      boost::circular_buffer<double> mrms_values(buffer_size);
      boost::circular_buffer<std::vector<double>> state_values(buffer_size);            
    
      /*Run the simulation*/

      std::vector<double> current_state_variables = p_model->GetStdVecStateVariables(), previous_state_variables;
      std::cout << "Number of state variables : " << current_state_variables.size() << "\n";
      for(unsigned int i = 0; i < paces; i++){
	previous_state_variables = current_state_variables;
	p_model->SolveAndUpdateState(0, duration);
	p_model->SolveAndUpdateState(duration, period);
	current_state_variables = p_model->GetStdVecStateVariables();

	state_values.push_back(current_state_variables);
	double mrms_val = mrms(current_state_variables, previous_state_variables);
	/*If values is full, this will remove the oldest value so that values will only contain up to buffer_size points*/
	mrms_values.push_back(mrms_val);
	mrms_file << mrms_val  << "\n";

	if(i>buffer_size && i % 10 == 0){
	  output_file << CalculatePMCC(mrms_values) << " ";
	  // std::vector<std::vector<double>> vals;
	  // for(unsigned int j = 0; j < current_state_variables.size(); j++){
	  //   std::vector<double> state_trace = cGetNthVariable(state_values, j);
	  //   vals = LogDifferences(state_trace);
	  //   output_file << CalculatePMCC(vals) << " ";
	  //
	  output_file << "\n";
	}
	  
	  
	//	  std::vector<double> regression_values = FitExponential(vals[0], vals[1]);
	  //std::cout  << "Regression values: "<< regression_values[0] << ", " << regression_values[1] << "\n";

      }
      output_file.close();
    }
    
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

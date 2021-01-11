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
#include "decker_2009Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
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

    for(unsigned int i = 0; i < models.size(); i++){
      boost::shared_ptr<AbstractCvodeCell> p_model = models[i];
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
      
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      std::string username = std::string(getenv("USER"));
    
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
    
      const double duration   = p_regular_stim->GetDuration();
      double period = 1000;
      if(i < 4){
	period = 500;
      }
    
      p_regular_stim->SetPeriod(period);
      p_regular_stim->SetStartTime(0);
      p_model->SetTolerances(1e-8, 1e-8);
      p_model->SetMaxSteps(1e6);
    
      const unsigned int paces = 10000;
      std::ofstream output_file;
      if(period==500)
	output_file.open("/tmp/"+username+"/"+model_name+"/StateVariablesPlot1Hz2Hz.dat");
      else
	output_file.open("/tmp/"+username+"/"+model_name+"/StateVariablesPlot2Hz1Hz.dat");
      std::vector<double> current_states;
    
      if(period == 500){
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
      }
      else{
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth2Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
      }
    
      std::vector<std::string> names = p_model->GetSystemInformation()->rGetStateVariableNames();
    
      for(auto i = names.begin(); i != names.end(); i++){
	output_file << *i << " ";
      }
      output_file << "\n";
      

      current_states = p_model->GetStdVecStateVariables();

      
      for(auto i = current_states.begin(); i != current_states.end(); i++){
	output_file << *i << " ";
      }

      output_file << "\n";

      for(unsigned int i = 0; i < paces; i++){
	p_model->SolveAndUpdateState(0, duration);
	p_model->SolveAndUpdateState(duration, period);
	current_states = p_model->GetStdVecStateVariables();
	for(auto i = current_states.begin(); i != current_states.end(); i++){
	  output_file << *i << " ";
	}
	output_file << "\n";
      }
      output_file.close();
    }
#else
      std::cout << "Cvode is not enabled.\n";
#endif
    } 
  };
  

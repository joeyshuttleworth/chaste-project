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
#include "Simulation.hpp"
#include <boost/filesystem.hpp>
#include <fstream>

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"

/*Output the total number of paces to reach limiting behaviour by SmartSimulation over all models for different choices of buffer_size and extrapolation_constant*/

class TestBenchmark : public CxxTest::TestSuite
{
private:
  const double threshold = 1.8e-07;
  const unsigned int paces = 5000;
  std::ofstream output_file;
  std::string username;

  const std::vector<unsigned int> buffer_sizes = {50}; //{25, 50, 100, 150, 200, 300 ,400};
  const std::vector<double>       extrapolation_constants = {0.9};
  const std::vector<double> apds = {231.87874670168091, 186.378828929472235, 228.391097821924944, 186.568452124915893, 268.928719750840457, 212.495013520340706, 268.49004986811957, 211.93350538338558};
public:
  unsigned int RunModel(boost::shared_ptr<AbstractCvodeCell> p_model, double period, unsigned int buffer_size, double extrapolation_constant, unsigned int index){
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      std::cout << "Testing " << model_name << " with period " << period << "\n";
      const double duration = p_model->UseCellMLDefaultStimulus()->GetDuration();      
      std::string input_path;

      if(period == 500)
	input_path = "/home/"+username+"/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat";
      else if(period == 1000)
	input_path = "/home/"+username+"/code/chaste-project-data/"+model_name+"/GroundTruth2Hz/final_state_variables.dat";
      
      unsigned int j;
      /*Run the simulations*/
      SmartSimulation simulation(p_model, period, input_path);
      simulation.Initialise(buffer_size, extrapolation_constant);
      for(j = 0; j < paces; j++){
	if(simulation.RunPace())
	  std::cout << "Model " << model_name << " period " << period << " Extrapolation method finished after " << j << " paces \n";
	
	if(simulation.is_finished()){
	  break;
	}	
      }

      /*Check that the methods have converged to the same place*/
      std::vector<double> vec1 = simulation.GetStateVariables();
      output_file << model_name << " " << period << " ";
      WriteStatesToFile(vec1, output_file);
      output_file << CalculateAPD(p_model, period, duration, 90) << "\n";
      double apd = CalculateAPD(p_model, period, duration, 90);
      double apd_error = CalculateAPD(p_model, period, duration, 90) - apds[index];
      std::cout << "apd error " << apds[index] << " " << apd_error << " " <<apd <<  "\n";
      TS_ASSERT(abs(CalculateAPD(p_model, period, duration, 90) - apds[index]) < 0.1);
      return j;
  }
  
  void TestMain(){
#ifdef CHASTE_CVODE
    username = std::string(getenv("USER"));
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
    
    boost::filesystem::create_directory("/tmp/"+username);
    output_file.open("/tmp/"+username+"/BenchmarkStates.dat");
    output_file.precision(18);
    std::ofstream f_results;
    f_results.open("/tmp/joey/BenchmarkResults.dat");

    for(unsigned int i = 0; i < extrapolation_constants.size(); i++){
      f_results << extrapolation_constants[i] << "\t";
    }

    f_results << "\n";
    
    for(unsigned int j = 0; j < buffer_sizes.size(); j++){
      f_results << buffer_sizes[j] << " ";
      for(unsigned int k = 0; k < extrapolation_constants.size(); k++){
	unsigned int benchmark = 0;
	for(unsigned int i = 0; i < 8; i++){
	  double period = 1000;
	  if(i<4)
	    period = 500;
	  benchmark += RunModel(models[i], period, buffer_sizes[j], extrapolation_constants[k], i);
	}
	std::cout << "Score is: " << benchmark << "\n";
	f_results << benchmark << "\t";
	TS_ASSERT(output_file.is_open());
      }
        f_results <<"\n";
    }
    output_file << "\n";
    output_file.close();
    f_results.close();
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

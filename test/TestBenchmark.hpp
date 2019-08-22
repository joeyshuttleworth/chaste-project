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
  const double thresholds[4] = {0.01, 0.01, 0.01, 0.01}; 
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
    
    for(unsigned int i = 0; i < models.size(); i++){
      double period = 1000;
      if(i<4)
	period = 500;

      boost::shared_ptr<AbstractCvodeCell> p_model = models[i];
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
      
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/GroundTruth2Hz");
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/GroundTruth1Hz");
      const double duration   = p_regular_stim->GetDuration();

      p_regular_stim->SetPeriod(period);
      p_regular_stim->SetStartTime(0);
      p_model->SetTolerances(1e-12, 1e-12);
      p_model->SetMaxSteps(1e5);
      
      double sampling_timestep = 0.01;
      const unsigned int paces  = 10000;
      OdeSolution current_solution;
      std::ofstream output_file;
      const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();
      std::vector<bool> threshold_met = {false, false, false, false};
      
      std::cout << "Testing " << model_name << " with period " << period << "\n";
      p_model->SetMaxTimestep(1000);
      if(period==500){
	output_file.open("/tmp/"+username+"/"+model_name+"/1Hz2HzBenchmark.dat");
      }
      else{
	output_file.open("/tmp/"+username+"/"+model_name+"/2Hz1HzBenchmark.dat");
      }
      TS_ASSERT_EQUALS(output_file.is_open(), true);

      /*Set the output to be as precise as possible */
    output_file.precision(18);

    p_model->SetMinimalReset(true);
    /*Run the simulation*/
    std::vector<std::vector<double>> current_states, previous_states;
    for(unsigned int i = 0; i < paces - 1; i++){
      OdeSolution solution = p_model->Compute(0, duration, sampling_timestep);
      if(i>0)
	previous_states = current_states;
      current_states = solution.rGetSolutions();
      solution = p_model->Compute(duration, period, sampling_timestep);
      current_states.insert(current_states.end(), solution.rGetSolutions().begin(), solution.rGetSolutions().end());
      if(i>0){
	if(!threshold_met[0] && TwoNorm(current_states.back(), previous_states.back()) < thresholds[0]){
	  output_file << model_name << "with period " << period << "took " << i << " for the two norm to less than " << thresholds[0]; 
	}
	if(!threshold_met[1] && mrms(current_states.back(), previous_states.back()) < thresholds[1]){
	  output_file << model_name << "with period " << period << "took " << i << " for the mrms to less than " << thresholds[1]; 
	}
	if(!threshold_met[2] && TwoNormTrace(current_states, previous_states) < thresholds[2]){
	  output_file << model_name << "with period " << period << "took " << i << " for the pace two norm to less than " << thresholds[2]; 
	}
	if(!threshold_met[3] && mrmsTrace(current_states, previous_states) < thresholds[3]){
	  output_file << model_name << "with period " << period << "took " << i << " for the pace mrms to be less than " << thresholds[3]; 
	}
      }
    }
    
    std::vector<std::vector<double>> final_trace;
    std::vector<double> times;

    /*Store the solution to the final pace and create std::vectors of the voltages and the times*/
    current_solution = p_model->Compute(0, duration, sampling_timestep);
    times = current_solution.rGetTimes();
    final_trace = current_solution.rGetSolutions();
    current_solution = p_model->Compute(duration, period, sampling_timestep);
    times.insert(times.end(), ++current_solution.rGetTimes().begin(), current_solution.rGetTimes().end());
    final_trace.insert(final_trace.end(), ++current_solution.rGetSolutions().begin(), current_solution.rGetSolutions().end());
	
    }
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

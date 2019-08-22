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
#include "luo_rudy_1994Cvode.hpp"
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

    //models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    //models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    //models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));

    std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);

    const double sampling_timestep = 0.01;
    const unsigned int paces  = 10000;
    OdeSolution current_solution;
    std::ofstream errors_file;
    
    for(unsigned int i=0; i < 2; i++){
      boost::shared_ptr<AbstractCvodeCell> p_model = models[0];
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
      const double duration   = p_regular_stim->GetDuration();
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      std::cout << "Testing model: " + model_name + "\n";

      double period = 1000;
      if(i < 1)
	period = 500;
      
      p_regular_stim->SetPeriod(period);
      p_model->SetTolerances(1e-12, 1e-12);
      p_model->SetMaxSteps(1e5);
      p_model->SetMaxTimestep(1000);
      p_regular_stim->SetStartTime(0);
    

      std::vector<std::vector<double>> state_variables;

      std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
    
      const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();

      std::string errors_file_path;
	
      if(period == 500){
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
	errors_file_path = "/tmp/"+username+"/"+model_name+"/1Hz2Hzerrors.dat";
      }
      else{
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth2Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
	errors_file_path = "/tmp/"+username+"/"+model_name+"/2Hz1Hzerrors.dat";
      }
      errors_file.open(errors_file_path);
      TS_ASSERT_EQUALS(errors_file.is_open(), true);

      errors_file.precision(18);
      
      errors_file << "APD 2-Norm MRMS 2-Norm-Trace MRMS-Trace ";

      std::vector<std::string> names = p_model->GetSystemInformation()->rGetStateVariableNames();

      for(unsigned int i = 0; i < names.size(); i++){
	errors_file << names[i] << " ";
      }
      errors_file << "\n";
      
      unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
      std::vector<double> times;
      std::ifstream apd_file;
      std::vector<std::vector<double>> current_state_variables, previous_state_variables;
      for(unsigned int i = 0; i < paces; i++){
	previous_state_variables = current_state_variables;
	current_solution = p_model->Compute(0, duration, sampling_timestep);
	times = current_solution.rGetTimes();
	current_state_variables = current_solution.rGetSolutions();
	current_solution = p_model->Compute(duration, period, sampling_timestep);
	state_variables = current_solution.rGetSolutions();
	current_state_variables.insert(current_state_variables.end(), ++state_variables.begin(), state_variables.end()); 
	times.insert(times.end(), ++current_solution.rGetTimes().begin(), current_solution.rGetTimes().end());
	if(i>0){
	  std::vector<double> voltages = GetNthVariable(current_state_variables, voltage_index);
	  CellProperties cell_props = CellProperties(voltages, times);
	  double apd = cell_props.GetLastActionPotentialDuration(90);
	  errors_file << apd << " ";
	  errors_file << TwoNorm(current_state_variables.back(), previous_state_variables.back()) << " ";
	  errors_file << mrms(current_state_variables.back(),  previous_state_variables.back()) << " ";
	  errors_file << TwoNormTrace(current_state_variables, previous_state_variables) << " ";
	  errors_file << mrmsTrace(current_state_variables, previous_state_variables) << " ";
	  /* for(unsigned int k = 0; k < current_state_variables.back().size(); k++){
	    errors_file << current_state_variables.back()[k] << " ";
	    }*/
	}
	errors_file << "\n";
      }
      errors_file.close();
    }
#else
      std::cout << "Cvode is not enabled.\n";
#endif
    }
};

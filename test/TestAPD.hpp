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
    boost::shared_ptr<AbstractCvodeCell> p_model(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
    const double period = 1000;
   
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
      
    const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();
   
    unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
    std::vector<double> times;
    std::vector<double> sampling_timesteps = {1, 0.1, 0.01, 0.001};
    
    for(unsigned int i = 0; i < sampling_timesteps.size(); i++){
      for(unsigned int j = 0; j < paces; j++){
	current_solution = p_model->Compute(0, duration);
	state_variables = current_solution.rGetSolutions();
	times = current_solution.rGetTimes();
	current_solution = p_model->Compute(duration, period);
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
  

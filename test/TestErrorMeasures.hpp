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
    
    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;

    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));

    for(unsigned int i=0; i < models.size(); i++){
      boost::shared_ptr<AbstractCvodeCell> p_model = models[i];
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
      const double period = 500;
      const double duration   = p_regular_stim->GetDuration();
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      std::cout << "Testing model: " + model_name + "\n";
      p_regular_stim->SetPeriod(period);
      p_model->SetTolerances(1e-12, 1e-12);
      p_model->SetMaxSteps(1e5);
      p_model->SetMaxTimestep(1000);
      p_regular_stim->SetStartTime(0);
    
      double sampling_timestep = 0.01;
      unsigned int paces  = 100;
      OdeSolution current_solution;
      std::ofstream errors_file;
      std::vector<std::vector<double>> state_variables;
      std::vector<std::vector<double> >  *current_state_variables = NULL, *previous_state_variables= NULL; 
      std::string username = std::string(getenv("USER"));

      std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
    
      const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();
      TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat"), 0);
      std::string errors_file_path = "/tmp/"+username+"/"+model_name+"/errors.dat";
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);

      errors_file.open(errors_file_path);
      errors_file << "2-Norm MRMS 2-Norm-Trace MRMS-Trace\n ";
    
      double current_apd90=0, previous_apd90=0;
      unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
      std::vector<double> times;
      std::ifstream apd_file;
      double true_apd90;

      apd_file.open("/home/joey/code/chaste-project-data/" + model_name + "/GroundTruth2Hz/apd90.dat");
      TS_ASSERT_EQUALS(apd_file.is_open(), true);
      apd_file >> true_apd90;

      for(unsigned int i = 0; i < paces; i++){
	if(previous_state_variables)
	  delete(previous_state_variables);
	previous_state_variables = current_state_variables;
	current_solution = p_model->Compute(0, duration, 1000);
	if(i==0)
	  times = current_solution.rGetTimes();
	current_state_variables = new std::vector<std::vector<double> >(current_solution.rGetSolutions());
	current_solution = p_model->Compute(duration + sampling_timestep, period, sampling_timestep);
	state_variables = current_solution.rGetSolutions();
	current_state_variables->insert(current_state_variables->end(), state_variables.begin(), state_variables.end()); 
	if(i>0){
	  errors_file << TwoNorm(current_state_variables->back(), previous_state_variables->back()) << " ";
	  errors_file << mrms(current_state_variables->back(),  previous_state_variables->back()) << " ";
	  errors_file << TwoNormTrace(*current_state_variables, *previous_state_variables) << " ";
	  errors_file << mrmsTrace(*current_state_variables, *previous_state_variables) << " ";
	  previous_apd90 = current_apd90;
	}
	else
	  times.insert(times.end(), current_solution.rGetTimes().begin(), current_solution.rGetTimes().end());
	const std::vector<double> voltages = GetNthVariable(current_state_variables, voltage_index);
	CellProperties cell_props = CellProperties(voltages, times); 
	current_apd90 = cell_props.GetLastActionPotentialDuration(90);
	if(i>0){
	  errors_file << abs(current_apd90 - previous_apd90) << " " << current_apd90 - true_apd90 << "\n";
	}
      }
      errors_file.close();
    }
#else
      std::cout << "Cvode is not enabled.\n";
#endif
    }
};

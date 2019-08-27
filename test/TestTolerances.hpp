#include <cxxtest/TestSuite.h>
#include <sstream>
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

#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  std::string username;
  void TestTolerances()
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
    username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);
    std::vector<double> tolerances = {1e-4, 1e-6, 1e-8, 1e-10};

    for(unsigned int i = 0; i < tolerances.size(); i++){
      unsigned int j = 1;
	if(j<4)
	  PrintErrors(models[j], tolerances[i], 1000);
	else
	  PrintErrors(models[j], tolerances[i], 500);   
    }
  
  void PrintErrors(boost::shared_ptr<AbstractCvodeCell> p_model, double tolerance, const double period){

    boost::shared_ptr<RegularStimulus> p_regular_stim  = p_model->UseCellMLDefaultStimulus();
    p_regular_stim->SetPeriod(period);
    p_model->SetTolerances(tolerance, tolerance);
    p_model->SetMaxSteps(1e5);
    p_model->SetMaxTimestep(1000);
    p_regular_stim->SetStartTime(0);
     
    const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();
    const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
    const unsigned int paces = 10000;
    std::stringstream file_path;
	
      if(period==500){
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
	file_path << "/tmp/"+username+"/"+model_name+"/1Hz2Hz-st-" << std::scientific << tolerance << ".dat";
      }
      else{
	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth2Hz/final_state_variables.dat"), 0);
	boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);      
        file_path << "/tmp/"+username+"/"+model_name+"/2Hz1Hz-st-" << std::scientific << tolerance << ".dat";
      }
      std::ofstream errors_file;
      errors_file.open(file_path.str());
      TS_ASSERT_EQUALS(errors_file.is_open(), true);

      errors_file.precision(18);
      
      errors_file << "APD90 2-Norm MRMS 2-Norm-Trace MRMS-Trace";
      errors_file << "\n";
      
      //      unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
      std::vector<double> times;
      std::ifstream apd_file;
      std::vector<std::vector<double>> current_state_variables, previous_state_variables;
      double duration = p_regular_stim->GetDuration();
      for(unsigned int i = 0; i < paces; i++){
	previous_state_variables = current_state_variables;
	if(i>0 && i % 10 == 1){
	  std::vector<double> previous_state_variables = p_model->GetStdVecStateVariables();
	  p_model->SolveAndUpdateState(0, duration);
	  p_model->SolveAndUpdateState(duration, period);
	  std::vector<double> current_state_variables = p_model->GetStdVecStateVariables();
	  errors_file << CalculateAPD(p_model, period, duration, 90.0) << " ";     
	  errors_file << TwoNorm(current_state_variables, previous_state_variables) << " ";
	  errors_file << mrms(current_state_variables,  previous_state_variables) << " ";
	  errors_file << "\n";
	}

	else{
	  p_model->SolveAndUpdateState(0, duration);
	  p_model->SolveAndUpdateState(duration, period);
	}
      }
      errors_file.close();
  }
#else
      std::cout << "Cvode is not enabled.\n";
#endif
};


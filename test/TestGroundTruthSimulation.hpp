#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
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

    std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);
    
    for(unsigned int i = 0; i < 8; i++){
      double period = 1000;
      if(i<4)
	period = 500;

      boost::shared_ptr<AbstractCvodeCell> p_model = models[i%4];
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();

      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/GroundTruth2Hz");
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/GroundTruth1Hz");
      const double start_time = p_regular_stim->GetStartTime();
      const double duration   = p_regular_stim->GetDuration();
      
      p_regular_stim->SetPeriod(period);
      p_regular_stim->SetStartTime(0);
      p_model->SetTolerances(1e-12, 1e-12);
      p_model->SetMaxSteps(1e5);
      
      double sampling_timestep = 0.01;
      unsigned int paces  = 10000;
      OdeSolution current_solution;
      std::ofstream apd_file;
      std::ofstream variables_file;
      const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();
      
      p_model->SetMaxTimestep(1000);
      if(i<4){
	apd_file.open("/tmp/"+username+"/"+model_name+"/GroundTruth2Hz/apd90.dat");
	variables_file.open("/tmp/"+username+"/"+model_name+"/GroundTruth2Hz/final_state_variables.dat");
      }
      else{
	apd_file.open("/tmp/"+username+"/"+model_name+"/GroundTruth1Hz/apd90.dat");
	variables_file.open("/tmp/"+username+"/"+model_name+"/GroundTruth1Hz/final_state_variables.dat");
      }
    /*Set the output to be as precise as possible */
    apd_file.precision(18);
    variables_file.precision(18);
    /*Run the simulation*/
    for(unsigned int i = 0; i < paces - 1; i++){
      /*Set the initial values to be the terminal values of the last solution*/
      p_model->SolveAndUpdateState(0, start_time);
      p_model->SolveAndUpdateState(start_time, start_time + duration);
      p_model->SolveAndUpdateState(start_time + duration, period);
    }
    
    std::vector<std::vector<double>> final_trace;
    std::vector<double> voltages;
    std::vector<double> times;

    /*Store the solution to the final pace and create std::vectors of the voltages and the times*/
    current_solution = p_model->Compute(0, duration, sampling_timestep);
    times = current_solution.rGetTimes();
    final_trace = current_solution.rGetSolutions();
    current_solution = p_model->Compute(duration, period, sampling_timestep);
    times.insert(times.end(), ++current_solution.rGetTimes().begin(), current_solution.rGetTimes().end());
    final_trace.insert(final_trace.end(), ++current_solution.rGetSolutions().begin(), current_solution.rGetSolutions().end());
	
    unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
		
    voltages.reserve(final_trace.size());
    for(unsigned int i = 0; i < final_trace.size(); i++){
      voltages.push_back(final_trace[i][voltage_index]);
    }

    /*Output terminal state variables to file*/
    for(unsigned int i = 0; i < final_trace[0].size(); i++){
      variables_file << final_trace.back()[i] << " ";
    }
    variables_file << "\n";
	
    /*Calculate and output apd to file*/
    CellProperties cell_props = CellProperties(voltages, times); 
    double apd = cell_props.GetLastActionPotentialDuration(90);
    apd_file << apd << " ";
    variables_file.close();
    apd_file << "\n";
    apd_file.close();	
    }
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

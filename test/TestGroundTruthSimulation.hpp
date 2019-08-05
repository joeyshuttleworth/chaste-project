#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include <fstream>

/* This file is generated from the cellml files provided at github.com/chaste/cellml */
#include "beeler_reuter_model_1977Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
    void TestTusscherSimulation()
    {
#ifdef CHASTE_CVODE
      boost::shared_ptr<RegularStimulus> p_stimulus;
      boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
      boost::shared_ptr<AbstractCvodeCell> p_model(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus));
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	
	const double period = 1000;
	const double start_time = p_regular_stim->GetStartTime();
	const double duration   = p_regular_stim->GetDuration();
	
	p_regular_stim->SetPeriod(period);
     	p_model->SetTolerances(1e-12, 1e-12);
	p_model->SetMaxSteps(1e5);
	
        double sampling_timestep = 0.001;
	unsigned int paces  = 10000;
	OdeSolution current_solution;
	std::ofstream apd_file;
	std::ofstream variables_file;
	std::string username = std::string(getenv("USER"));
	const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();

	p_model->SetMaxTimestep(1000);
	
	apd_file.open("/tmp/"+username+"/apd90.dat");
	variables_file.open("/tmp/"+username+"/final_state_variables.dat");
	/*Set the output to be as precise as possible */
        apd_file.precision(17);
	variables_file.precision(17);
	/*Run the simulation*/
	for(unsigned int i = 0; i < paces - 1; i++){
	  /*Set the initial values to be the terminal values of the last solution*/
	  p_model->SolveAndUpdateState(0, start_time);
	  p_model->SolveAndUpdateState(start_time, start_time + duration);
	  p_model->SolveAndUpdateState(start_time+duration, period);
	}
	std::vector<std::vector<double>> final_trace;
	std::vector<double> voltages;
	std::vector<double> times;

	/*Store the solution to the final pace and create std::vectors of the voltages and the times*/
	current_solution = p_model->Compute(0, start_time, sampling_timestep);
	final_trace.insert(final_trace.end(), current_solution.rGetSolutions().begin(), current_solution.rGetSolutions().end());
	times = current_solution.rGetTimes();
	current_solution = p_model->Compute(start_time, start_time + duration, sampling_timestep);
	times.insert(times.end(), current_solution.rGetTimes().begin(), current_solution.rGetTimes().end());
	final_trace.insert(final_trace.end(), current_solution.rGetSolutions().begin(), current_solution.rGetSolutions().end());	
	current_solution = p_model->Compute(start_time+duration, period, sampling_timestep);
	times.insert(times.end(), current_solution.rGetTimes().begin(), current_solution.rGetTimes().end());
	final_trace.insert(final_trace.end(), current_solution.rGetSolutions().begin(), current_solution.rGetSolutions().end());
	
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
	
  	current_solution.WriteToFile("TestCvodeCells","final_trace","ms");
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

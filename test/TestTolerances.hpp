#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "TenTusscher2004EpiCvode.hpp"
#include "FakePetscSetup.hpp"
#include <fstream>
class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
    void TestTusscherSimulation()
    {
#ifdef CHASTE_CVODE
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellTenTusscher2004EpiFromCellMLCvode(p_solver, p_stimulus));
	boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	
	const double period = 1000;
	const std::vector<double> tolerances = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9};
        p_regular_stim->SetPeriod(period);

	double max_timestep = p_regular_stim->GetDuration()/2;

        p_model->SetMaxTimestep(max_timestep);
	//unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane__V");
	
        double sampling_timestep = max_timestep;
	int steps = 1000;
	OdeSolution current_solution;
	std::ofstream tolerances_file;
	std::string username = std::string(getenv("USER"));
	const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames(); 
	
	tolerances_file.open("/tmp/"+username+"/tolerances_and_tolerances.ssv");
	/*Set cout to be as precise as possible */
	tolerances_file.precision(17);
	/*Print variable names on the first line*/
	//	tolerances_file << "APD90 ";
	for(unsigned int i = 0; i < state_variable_names.size(); i++){
	  tolerances_file << state_variable_names[i] << " ";
	}
	tolerances_file << "\n";

	std::vector<double> state_variables;
	std::vector<double> state_variables2;
 	for(unsigned int k = 0; k < tolerances.size(); k++){
	  std::cout << "Testing with tolerance " << tolerances[k] << "\n";
	  p_model->SetTolerances(tolerances[k], tolerances[k]); 
	  for(int i=0; i < steps; i++){
	    double start_time = i*period;
	    double end_time   = start_time + period;
		
	    /*Set the initial values to be the termial values of the last solution*/
	    if(current_solution){
	      state_variables = current_solution->rGetSolutions()[current_solution->GetNumberOfTimeSteps()-1];
	      state_variables2 = p_model->rGetVariables[current_solution->GetNumberOfTimeSteps()-1];
	      for(unsigned int i = 0; i < state_variabvles_names.size(); i++){
		std::cout << state_variables[i] << " vs " << state_variables2[i] << "\n"; 
	      }
	      p_model->SetStateVariables(state_variables);
	    }	      
	    current_solution = p_model->Compute(start_time, end_time, sampling_timestep);
	  }
	  //std::vector<double> voltages = current_solution->GetVariableAtIndex(voltage_index);
	  //CellProperties cell_props(voltages, current_solution->rGetTimes());
	  //double apd = cell_props.GetLastActionPotentialDuration(90);
	  tolerances_file << tolerances[k] << " ";
	  for(unsigned int j=0; j < state_variables.size(); j++){
	    tolerances_file << state_variables[j] << " ";
	    //tolerances_file << apd << " " << state_variables[j] << " ";
	  }
	  tolerances_file << "\n";
	}
	  tolerances_file.close();	
#else
	std::cout << "Cvode is not enabled.\n";
#endif
    }
};

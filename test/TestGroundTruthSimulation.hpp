#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "TenTusscher2006EpiCvode.hpp"
#include "FakePetscSetup.hpp"
#include <fstream>
class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
    void TestShannonSimulation()
    {
#ifdef CHASTE_CVODE
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellTenTusscher2006EpiFromCellMLCvode(p_solver, p_stimulus));
	boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	
	const double period = 1000;
        p_regular_stim->SetPeriod(period);
     	p_model->SetTolerances(1e-7,1e-7);

	double max_timestep = p_regular_stim->GetDuration()/2;

        p_model->SetMaxTimestep(max_timestep);
	unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane__V");
	
        double sampling_timestep = max_timestep;
	int steps = 10000;
	OdeSolution *current_solution = NULL;
	std::ofstream apd_file;
	std::ofstream variables_file;
	std::string username = std::string(getenv("USER"));
	const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames(); 
		
	apd_file.open("/tmp/"+username+"/apd90plot.ssv");
	variables_file.open("/tmp/"+username+"/state_variables.ssv");
	/*Set cout to be as precise as possible */
        apd_file.precision(17);
	variables_file.precision(17);
	/*Print variable names on the first line*/
	for(unsigned int i = 0; i < state_variable_names.size(); i++){
	  variables_file << state_variable_names[i] << " ";
	}
	variables_file << "\n";
	for(int i=0; i < steps; i++){
	  double start_time = i*period;
	  double end_time   = start_time + period;
	  double apd;
	  std::vector<double> state_variables;
		
	  /*Set the initial values to be the termial values of the last solution*/
	  if(current_solution){
	    state_variables = current_solution->rGetSolutions()[current_solution->GetNumberOfTimeSteps()-1];
	    std::vector<double> voltages = current_solution->GetVariableAtIndex(voltage_index);
	    CellProperties cell_props(voltages, current_solution->rGetTimes());
	    
	    apd = cell_props.GetLastActionPotentialDuration(90);
	    p_model->SetStateVariables(state_variables);
	  }
	  current_solution = new OdeSolution;
	  *current_solution = p_model->Compute(start_time, end_time, sampling_timestep);
	  apd_file << apd << " ";
	  for(unsigned int j=0; j < state_variables.size(); j++){
	    variables_file << state_variables[j] << " ";
	  }
	  variables_file << "\n";
	}
	variables_file.close();
	apd_file << "\n";
	apd_file.close();
	//        solution.WriteToFile("TestCvodeCells","TenTusscher2006EpiGroundTruth", "ms");
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

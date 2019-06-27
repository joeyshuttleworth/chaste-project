#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "TenTusscher2006EpiCvode.hpp"
#include "FakePetscSetup.hpp"

double ValueAtTime(std::vector<double> times,
		   std::vector<std::vector<double>> values,
		   double time,
		   int    variable_index){
  unsigned int i;
  for(i = 0; i < values.size();i++){
    if(times[i] > time)
      break;
  }
  if(i == 0){
    return values[0][variable_index];
  }
  else{
    //interpolate between i - 1 and i + 1
    return values[i-1][variable_index] + (times[i] - time) * (values[i+1][variable_index] - values[i][variable_index]) / (times[i] - times[i-1]);
  }
}

class TestUndsteadyNorms : public CxxTest::TestSuite
{
public:
    void TestSimulation()
    {
#ifdef CHASTE_CVODE
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellTenTusscher2006EpiFromCellMLCvode(p_solver, p_stimulus));

        boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	const double period = 500.0;
        p_regular_stim->SetPeriod(period);
     	p_model->SetTolerances(1e-9,1e-9);

	double max_timestep = p_regular_stim->GetDuration()/2;

        p_model->SetMaxTimestep(max_timestep);

        double steps = 1000;
	
	p_model->ForceUseOfNumericalJacobian();

	std::vector<double> change_between_paces;
	const std::vector<std::vector<double>> values = solution.rGetSolutions();
        TS_ASSERT_EQUALS(values.size()==0, false);
	unsigned int number_of_variables = values[0].size();
	std::vector<double> changes;
	OdeSolution *current_solution, *last_solution;
	std::ofstream changes_file;
	std::string username = std::string(getenv("USER"));
	
	changes_file.open("/tmp/"+username+"/changes.ssv");
	
	for(int i=0; i < steps; i++){
	  double start_time = i*period;
	  double end_time   = start_time + period;

	  std::vector<std::vector<double>> *current_state_variables;
	  std::vector<std::vector<double>> *last_state_variables;
	  
	  if(last_solution)
	    delete(last_solution);
	  last_solution = current_solution;
	  /*Set the initial values to be the termial values of the last solution*/
	  if(current_solution){
	    last_state_variables = current_state_variables;
	    p_model->SetStateVariables(current_state_variables);
	    current_solution  = new OdeSolution;
	    *current_solution = p_model->Compute(start_time, end_time, max_timestep); 
	    *current_state_variables = current_solution->rGetSolutions();
	  }
	  if(i > 1){
	    /*Calculate the "change" between the state variables over the last two beats using the Euclidean norm*/
	    double sum = 0;
	    std::vector<double> last_times = last_solution->rGetTimes();
	    std::vector<double> current_times = current_solution->rGetTimes();
	    for(unsigned int j = 0; j < number_of_variables; j++){
	      for(unsigned int k = 0; k < current_state_variables->size(); k++){
		sum = sum + pow(current_state_variables[k][j] - ValueAtTime(times, last_state_variables, current_times[k], j), 2);
	      }
	    }
	  }
	}
#else
        std::cout << "Cvode is not enabled.\n";
#endif
	}
};

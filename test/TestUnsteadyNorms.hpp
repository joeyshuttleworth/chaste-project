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

/*Uses linear interpolation to find the value at a certain time */
double ValueAtTime(std::vector<double> times,
		   std::vector<std::vector<double>> values,
		   double time,
		   int    variable_index){
  unsigned int i;
  for(i = 0; i < values.size() - 1; i++){
    if(times[i] > time)
      break;
  }
  if(i == 0){
    return values[0][variable_index];
  }
  else{
    /*Interpolate between i - 1 and i */
    return values[i-1][variable_index] + (time - times[i]) * (values[i][variable_index] - values[i-1][variable_index])/ (times[i] - times[i-1]);
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
	const unsigned int number_of_variables = p_model->GetNumberOfStateVariables();
	  
        p_regular_stim->SetPeriod(period);
     	p_model->SetTolerances(1e-9,1e-9);

	double max_timestep = p_regular_stim->GetDuration()/2;

        p_model->SetMaxTimestep(max_timestep);

        double steps = 100;
	
	OdeSolution *current_solution = NULL, *last_solution = NULL;
	std::ofstream changes_file;
	std::string username = std::string(getenv("USER"));
	
	changes_file.open("/tmp/"+username+"/changes.ssv");
	
	for(int i=0; i < steps; i++){
	  double start_time = i*period;
	  double end_time   = start_time + period;
	  
	  if(last_solution)
	    delete(last_solution);
	  if(current_solution)
	    last_solution = current_solution;

	  current_solution  = new OdeSolution;
	  *current_solution = p_model->Compute(start_time, end_time, max_timestep); 

	  /*Set the initial values to be the terminal values of the last solution*/
	  p_model->SetStateVariables(current_solution->rGetSolutions().back());

	  if(i>1){
	    /*Calculate the change between the state variables over the last two beats using the Euclidean norm*/
	    double sum_2 = 0;
	    double sum_mrms = 0;
	    std::vector<double> last_times = last_solution->rGetTimes();
	    std::vector<double> current_times = current_solution->rGetTimes();
	    for(unsigned int j = 0; j < number_of_variables; j++){
	      for(unsigned int k = 0; k < current_solution->GetNumberOfTimeSteps(); k++){
		sum_2 = sum_2 + pow(current_solution->rGetSolutions()[k][j] - ValueAtTime(last_times, last_solution->rGetSolutions(), current_times[k], j), 2);
	      }
	    }
	    changes_file << sqrt(sum_2) << " " << sqrt(sum_mrms / number_of_variables) << "\n";
	  }
	}
#else
        std::cout << "Cvode is not enabled.\n";
#endif
	}
};

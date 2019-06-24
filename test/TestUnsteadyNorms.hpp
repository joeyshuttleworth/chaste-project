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
	
        p_regular_stim->SetPeriod(500.0);
     	p_model->SetTolerances(1e-9,1e-9);

	double max_timestep = 1;

        p_model->SetMaxTimestep(max_timestep);

        double sampling_timestep = max_timestep;
        double start_time = 0.0;
        double end_time = 1000000.0;
	p_model->ForceUseOfNumericalJacobian();
        OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

        solution.WriteToFile("TestCvodeCells", "TenTusscher2004Epi2Hz", "ms");

	int number_of_timesteps = solution.GetNumberOfTimeSteps();
	int pace = 1, sum = 0;
	std::vector<double> change_between_paces;
	const std::vector<double> times = solution.rGetTimes();
	const std::vector<std::vector<double>> values = solution.rGetSolutions();
        TS_ASSERT_EQUALS(values.size()==0, false);
	int number_of_variables = values[0].size();
	std::vector<double> changes;
	
	for(int i = 0; i <= number_of_timesteps; i++){
	  if(times[i] > pace * 500 || i == number_of_timesteps){
	    /*We are now in the next pace. If this is not the first pace, 
	      calculate the "change" between this pace and the previous one.  */
	    if(pace > 1){
	      changes.push_back(sqrt(sum));
	      std::cout << sqrt(sum) << " ";
	    }
	    /*Update variables for the next pace*/
	    sum = 0;
	    pace++;
	  }
	  if(pace > 1){
	    for(int j = 0; j < number_of_variables; j++)
	      sum = sum + pow(values[i][j] - ValueAtTime(times, values, times[i] - 500, j), 2);
	  }
	}
	std::cout <<"\n";
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

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
        p_regular_stim->SetPeriod(period);
     	p_model->SetTolerances(1e-8,1e-8);

	double max_timestep = p_regular_stim->GetDuration()/2;

        p_model->SetMaxTimestep(max_timestep);
	
        double sampling_timestep = max_timestep;
	int steps = 100;
	OdeSolution current_solution;
	
	for(int i=0; i < steps; i++){
	  double start_time = p_regular_stim->GetStartTime();
	  if(start_time>0){
	    p_model->Compute(0, start_time, sampling_timestep);
	  }
	  time_1 = i*period + start_time;
	  time_2 = i*period + start_time + duration;
	  p_model->Compute(start_time, start_time + duration, sampling_timestep);
	  time_1 = i*period + start_time + duration;
	  time_2 = (i + 1) * period;
	  current_solution = p_model->Compute(start_time + duration, period, sampling_timestep);
	  std::cout << current_solution->rGetSolution()[0][0] << "\n";
	}
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

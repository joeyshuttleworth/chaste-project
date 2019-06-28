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
    void TestTusscherSimulation()
    {
#ifdef CHASTE_CVODE
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellTenTusscher2006EpiFromCellMLCvode(p_solver, p_stimulus));
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
	  double time_1 = i*period;
	  double time_2   = start_time + period;
	  double start_time = p_regular_stim->GetStartTime();
	  if(start_time>0){
	    time_1 = i*period;
	    time_2 = i*period + start_time;
	    p_model->SetStateVariables(p_model->Compute(time_1, time_2, samplig_timestep)->rGetSolutions().back());
	  }
	  time_1 = i*period + start_time;
	  time_2 = i*period + start_time + duration;
	  p_model->SetStateVaraibles(p_model->Compute(time_1, time_2, sampling_timestep)->rGetSolutions().back());
	  time_1 = i*period + start_time + duration;
	  time_2 = (i + 1) * period;
	  p_model->SetStateVaraibles(p_model->Compute(time_1, time_2, sampling_timestep)->rGetSolutions().back());
	}
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include <fstream>

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
        p_regular_stim->SetPeriod(period);
     	p_model->SetTolerances(1e-12,1e-12);

	double max_timestep = p_regular_stim->GetDuration();

        p_model->SetMaxTimestep(max_timestep);
	p_model->SetMaxSteps(1e5);
	
	int steps = 10000;
	OdeSolution current_solution;
	
	for(int i=0; i < steps; i++){
	  p_model->SolveAndUpdateState(0, period);
	}
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

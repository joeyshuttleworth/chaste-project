#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "TenTusscher2006EpiCvode.hpp"
#include "FakePetscSetup.hpp"

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
	
        p_regular_stim->SetPeriod(1000.0);
     	p_model->SetTolerances(1e-9,1e-9);

	double max_timestep = 1;

	std::cout << "Stimulus duration is: " << p_regular_stim->GetDuration() << "\n";
        p_model->SetMaxTimestep(max_timestep);

        double sampling_timestep = max_timestep;
        double start_time = 0.0;
        double end_time = 10000000.0;
	p_model->ForceUseOfNumericalJacobian();
        OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

        solution.WriteToFile("TestCvodeCells","TenTusscher2006EpiGroundTruth", "ms");

#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

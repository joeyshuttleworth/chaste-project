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
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
    void TestTusscherSimulation()
    {
#ifdef CHASTE_CVODE
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus));
	boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	std::ofstream output_file;

	output_file.open("/tmp/joey/OnePart.dat");
	const double period = 1000;
        p_regular_stim->SetPeriod(period);
	p_regular_stim->SetStartTime(0);
	p_model->SetTolerances(1e-5,1e-5);
	
	double max_timestep = p_regular_stim->GetDuration();
	double sampling_timestep = 0.5;
	
        p_model->SetMaxTimestep(max_timestep);
	p_model->SetMaxSteps(1e5);
	
	int paces = 1000;
	OdeSolution solution;
   
	for(int i=0; i < paces; i++){
	  solution = p_model->Compute(0, period, sampling_timestep);
	  std::vector<double> states = solution.rGetSolutions().back();
	  for(unsigned int j = 0; j < states.size(); j++){
	    output_file << states[j] << " ";
	  }
	  output_file << "\n";
	}

	solution.WriteToFile("ohara_rudy", "OnePartTrace", "ms");
	output_file.close();
       
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include <boost/filesystem.hpp>
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
	
	int paces = 1000;
	OdeSolution current_solution;
	const double period = 1000;
	double sampling_timestep = 0.5;
	std::ofstream output_file;
	std::ofstream trace_file;

	output_file.open("/tmp/joey/TwoPart.dat");
	trace_file.open("/tmp/joey/TwoPartTrace.dat");
	
	p_regular_stim->SetPeriod(period);
	p_regular_stim->SetStartTime(0);
     	p_model->SetTolerances(1e-5,1e-5);
	p_model->SetMaxSteps(1e5);
	
	p_regular_stim->SetStartTime(0);
	double duration = p_regular_stim->GetDuration();
	
	for(int i=0; i < paces - 1; i++){
	 p_model->Compute(0, duration, sampling_timestep);
	 OdeSolution solution = p_model->Compute(duration, period, sampling_timestep);	  
	 std::vector<double> states = solution.rGetSolutions().back();
	 for(unsigned int j = 0; j < states.size(); j++){
	   output_file << states[j] << " ";
	 }
	 output_file << "\n";
	}
	
	OdeSolution solution = p_model->Compute(0, duration, sampling_timestep);
	std::vector<std::vector<double>> states = solution.rGetSolutions();
	std::vector<double> times = solution.rGetTimes();
	  
	for(unsigned int i = 0; i < states.size(); i++){
	  trace_file << times[i] << " ";
	  for(unsigned int j = 0; j < states.size(); j++){
	    trace_file << states[i][j] << " ";
	  }
	  trace_file << "\n";
	}
	solution = p_model->Compute(duration, period, sampling_timestep);
	states = solution.rGetSolutions();
	times = solution.rGetTimes();
	for(unsigned int i = 1; i < states.size(); i++){
	  trace_file << times[i] << " ";
	  for(unsigned int j = 0; j < states[0].size(); j++){
	    trace_file << states[i][j] << " ";
	  }
	  trace_file << "\n";
	}
	
	trace_file.close();
	//solution.WriteToFile("1PartVs2Part", "2Part.dat", "ms");
#else 
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

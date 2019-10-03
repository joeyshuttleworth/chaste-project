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

#include "decker_2009Cvode.hpp"
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
      
      const int paces = 10000;
      const double period = 1000;
      
      std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
      
      models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
      models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
      models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
      models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));

      for(unsigned int i = 0; i < 4; i++){
	boost::shared_ptr<AbstractCvodeCell> p_model = models[i];
	boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	
	p_regular_stim->SetPeriod(period);
	p_regular_stim->SetStartTime(0);
     	p_model->SetTolerances(1e-12,1e-12);
	p_model->SetMaxSteps(1e5);
	
	p_regular_stim->SetStartTime(0);
	//	double duration = p_regular_stim->GetDuration();
	
	for(int i=0; i < paces; i++){
	  p_model->Solve(0, period, 1000);
	}
      }
#else 
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

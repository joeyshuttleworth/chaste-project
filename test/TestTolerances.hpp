#include <cxxtest/TestSuite.h>
#include <sstream>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include "SimulationTools.hpp"
#include <boost/filesystem.hpp>
#include <fstream>

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"
#include "iyer_2004.hpp"
#include "ToRORd_dynCl_endo.hpp"

// Models modified to have analytic voltage:
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  std::string username;
  void TestTolerances(){
    const std::vector<double> tolerances={1e-8, 1e-6, 1e-04, 1e-10, 1e-12};
    const std::vector<double> periods={1000};
    const std::vector<double> IKrBlocks={0};
    const std::string filename_suffix = "test_tolerances";


    for(auto tolerance : tolerances){

      boost::shared_ptr<RegularStimulus> p_stimulus;
      boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
      for(auto IKrBlock : IKrBlocks){
        for(auto period : periods){
            std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
            models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
            models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
            models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
            models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
            models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus)));
            models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus)));
            models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellhund_rudy_2004FromCellMLCvode(p_solver, p_stimulus)));
            for(auto model : models){
              compare_error_measures(model, period, IKrBlock, tolerance, filename_suffix);
          }
        }
      }
    }
  }
};

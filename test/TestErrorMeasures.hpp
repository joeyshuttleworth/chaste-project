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


class TestErrorMeasures : public CxxTest::TestSuite
{
public:
  std::string username;
  const int default_max_paces = 10000;
  void TestErrors(){
    int paces = get_paces();
    paces = paces==INT_UNSET?default_max_paces:paces;
    const std::vector<double> periods= get_periods(); {500, 750, 1000, 1250};
    const std::vector<double> IKrBlocks= get_IKr_blocks();

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models = get_models();

    const std::string filename_suffix = "error_measures";
    const double tolerance = 1e-10;

    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    for(auto model : models){
      for(auto period : periods){
        for(auto IKrBlock : IKrBlocks){
          compare_error_measures(paces, model, period, IKrBlock, tolerance, filename_suffix);
        }
      }
    }
  }
};

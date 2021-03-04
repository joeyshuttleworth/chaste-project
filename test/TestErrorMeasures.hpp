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
  void TestErrors(){
    const std::vector<double> periods={500, 750, 1000, 1250};
    const std::vector<double> IKrBlocks={0, 0.25, 0.5};
    const std::string filename_suffix = "error_measures";
    const double tolerance = 1e-10;

    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellhund_rudy_2004FromCellMLCvode(p_solver, p_stimulus)));

    for(auto model : models){
      for(auto period : periods){
        for(auto IKrBlock : IKrBlocks){
          compare_error_measures(model, period, IKrBlock, tolerance, filename_suffix);
        }
      }
    }
  }
};

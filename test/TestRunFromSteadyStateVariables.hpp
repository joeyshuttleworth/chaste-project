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

// Analytic models
#include "decker_2009_analytic_voltageCvode.hpp"
#include "hund_rudy_2004_analytic_voltageCvode.hpp"
#include "iyer_2004_analytic_voltageCvode.hpp"
#include "ohara_rudy_2011_epi_analytic_voltageCvode.hpp"
#include "ohara_rudy_cipa_2017_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2006_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2004_epi_analytic_voltageCvode.hpp"
#include "ToRORd_dyn_chloride_epi_analytic_voltageCvode.hpp"



class TestRunFromSteadyState : public CxxTest::TestSuite
{
public:
  std::string username;
  const int default_max_paces = 10000;
  void TestRun(){
    int paces = get_max_paces();

    double IKrBlock = 0.5;
    std::string input_filename ="";
    const std::string option = "--input_filename";

    if(CommandLineArguments::Instance()->OptionExists(option)){
      input_filename = CommandLineArguments::Instance()->GetStringCorrespondingToOption(option);
    }

    if(input_filename==""){
      EXCEPTION("No input filename given");
    }


    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    auto model = boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_epi_analytic_voltageFromCellMLCvode(p_solver, p_stimulus));
    Simulation sim(model, 10000, input_filename, 1e-8, 1e-8);
    sim.SetIKrBlock(IKrBlock);
    // sim.RunPaces(1000);

    std::cout << "apd is " << sim.GetApd(90) << "\n";

  }
};

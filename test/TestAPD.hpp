#include <cxxtest/TestSuite.h>
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
#include <iomanip>
#include "Simulation.hpp"

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

/*  Test to see how APD90s vary when the model has reached a steady state. Two
    different tolerances are used (1e-12, 1e-12) and (1e-8, 1e-8). The former is
    what the steady state was calculated with and the latter is more likely
    to be used in practice.

    The results should show that these APD90s are similar apart from maybe one
    or two paces immediately after the tolerances have been modified.

    This test outputs these sets of APD90s as two separate files.
 */

class TestAPD : public CxxTest::TestSuite
{
  const unsigned int paces=10;
public:
  void TestRunAPDSimulation()
  {
#ifdef CHASTE_CVODE

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    boost::shared_ptr<RegularStimulus> p_stimulus;
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellhund_rudy_2004FromCellMLCvode(p_solver, p_stimulus)));

    std::vector<double> periods = {1000, 500, 750, 1250};
    std::vector<double> IKrBlocks = {0, 0.25, 0.5};

    for(auto model : models){
      for(auto period : periods){
        for(auto IKrBlock : IKrBlocks)
          RunModel(model, period, IKrBlock);
      }
    }
  }

#else
  std::cout << "Cvode is not enabled.\n";
#endif
  void RunModel(boost::shared_ptr<AbstractCvodeCell> model, const double period, const double IKrBlock){
    std::stringstream dirname;
    const std::string CHASTE_TEST_OUTPUT = std::string(getenv("CHASTE_TEST_OUTPUT"));
    const std::string model_name = model->GetSystemInformation()->GetSystemName();
    std::cout << "Testing " << model_name << " with period " << int(period) <<" and IKrBlock "<< IKrBlock << "\n";

    const std::string GKrParameterName = "membrane_rapid_delayed_rectifier_potassium_current_conductance";
    const double default_GKr = model->GetParameter(GKrParameterName);
    model->SetParameter(GKrParameterName, default_GKr*(1-IKrBlock));

    dirname << model_name+"_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";
    boost::filesystem::create_directories(dirname.str());


    Simulation simulation(model, period, CHASTE_TEST_OUTPUT + "/" + dirname.str() + "final_states.dat");
    simulation.SetTerminateOnConvergence(false);

    std::ofstream output_file(CHASTE_TEST_OUTPUT + "/" + dirname.str() + "apds_using_groundtruth_1e-8.dat");
    std::ofstream output_file2(CHASTE_TEST_OUTPUT + "/" + dirname.str() + "apds_using_groundtruth_1e-12.dat");

    output_file << std::setprecision(20);
    output_file2 << std::setprecision(20);
    std::cout   << std::setprecision(20);

    // Calculate Ground Truth APD
    simulation.SetTolerances(1e-12, 1e-12);

    std::vector<double> last_state_vars, state_vars = model->GetStdVecStateVariables();
    for(unsigned int i=0; i<paces; i++){
      last_state_vars = state_vars;
      const double reference_apd = simulation.GetApd(90, false);
      simulation.RunPace();
      state_vars = simulation.GetModel()->GetStdVecStateVariables();
      std::cout << "mrms = " << mrms(state_vars, last_state_vars) << "\n";
      output_file2 << reference_apd << "\n";
    }
    simulation.SetTolerances(1e-8, 1e-8);

    for(unsigned int i = 0; i < paces; i++){
      simulation.RunPace();
      const double apd = simulation.GetApd(90, false);
      output_file << apd << "\n";
    }

    output_file.close();
    output_file2.close();
    model->SetParameter(GKrParameterName, default_GKr);
  }
};

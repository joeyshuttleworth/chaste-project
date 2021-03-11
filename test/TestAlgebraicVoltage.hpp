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
#include <algorithm>

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
#include "ToRORd_dynCl_epi_analytic_voltageCvode.hpp"
#include "iyer_2004_analytic_voltageCvode.hpp"


/* Run the models under different scenarios and output:
   - All variables over the final pace
   - The APD90
   - Terminal state variables
   to separate files.
 */

class TestGroundTruthSimulation : public CxxTest::TestSuite
{

  const std::string GKrParameterName = "membrane_rapid_delayed_rectifier_potassium_current_conductance_scaling_factor";
public:
  const int paces = 10000;
  void TestRunSimulation()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    const std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directories("/home/" + username + "/testoutput/");

    std::vector<double> periods = {1000, 500};
    std::vector<double> IKrBlocks = {0.25};

    std::vector<boost::shared_ptr<AbstractCvodeCell>> algebraic_models;
    std::vector<boost::shared_ptr<AbstractCvodeCell>> original_models;

    algebraic_models.push_back(boost::make_shared<CellToRORd_dynCl_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Celliyer_2004_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));


    original_models.push_back(boost::make_shared<CellToRORd_dynCl_epiFromCellMLCvode>(p_solver, p_stimulus));
    original_models.push_back(boost::make_shared<Celliyer_2004FromCellMLCvode>(p_solver, p_stimulus));

    for(unsigned i : indicies(algebraic_models)){
      const N_Vector initial_states = model->GetStateVariables();
      for(double period : periods){
        for(double IKrBlock : IKrBlocks){
          Simulation sim(analytic_models[i], period);
          Simulation sim_o(original_models[i], period);
          sim.SetIKrBlock(IKrBlock);
          sim.RunPaces(100);
          sim_0.RunPaces(100);

          std::vector<double> original_vars = original_models[i]->GetStdVecStateVariables();
          original_vars.pop_front();

          const double error = mrms(original_vars, algebraic_models[i]->GetStdVecStateVariables());
          original_models[i]->SetStateVariables(initial_states);
          algebraic_models[i]->SetStateVariables(initial_states);
        }
      }
    }

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

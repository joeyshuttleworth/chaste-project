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
#include "iyer_2004Cvode.hpp"
#include "ToRORd_dynCl_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"

// Models modified to have analytic voltage:
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"
#include "ToRORd_dyn_chloride_epi_analytic_voltageCvode.hpp"
#include "iyer_2004_analytic_voltageCvode.hpp"
#include "hund_rudy_2004_analytic_voltageCvode.hpp"

class TestAlgebraicVoltage : public CxxTest::TestSuite
{

  const std::string GKrParameterName = "membrane_rapid_delayed_rectifier_potassium_current_conductance_scaling_factor";
public:
  const int paces = 100;
  void TestRunSimulation()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    std::vector<double> periods = {1000};
    std::vector<double> IKrBlocks = {0, 0.25};

    std::vector<boost::shared_ptr<AbstractCvodeCell>> algebraic_models;
    std::vector<boost::shared_ptr<AbstractCvodeCell>> original_models;

    algebraic_models.push_back(boost::make_shared<CellToRORd_dyn_chloride_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Celliyer_2004_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Cellhund_rudy_2004_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));


    original_models.push_back(boost::make_shared<CellToRORd_dynCl_epiFromCellMLCvode>(p_solver, p_stimulus));
    original_models.push_back(boost::make_shared<Celliyer_2004FromCellMLCvode>(p_solver, p_stimulus));
    original_models.push_back(boost::make_shared<Cellhund_rudy_2004FromCellMLCvode>(p_solver, p_stimulus));

    for(unsigned int i = 0; i < algebraic_models.size(); i++){
      std::vector<double> original_initial_states = original_models[i]->GetStdVecStateVariables();
      std::vector<double> algebraic_initial_states = algebraic_models[i]->GetStdVecStateVariables();

      std::cout << "sizes " << original_initial_states.size() << " " << algebraic_initial_states.size() << "\n";
      for(double period : periods){
        for(double IKrBlock : IKrBlocks){
          const std::string model_name = original_models[i]->GetSystemInformation()->GetSystemName();
          std::cout << "Testing " << model_name << " with period " << int(period) <<" and IKrBlock "<< IKrBlock << "\n";

          std::vector<double> original_vars = original_initial_states;
          original_vars.erase(original_vars.begin());
          algebraic_models[i]->SetStateVariables(algebraic_initial_states);
          original_models[i]->SetStateVariables(original_initial_states);

          const double error1 = mrms(original_vars, algebraic_models[i]->GetStdVecStateVariables());
          std::cout << "initial error is " << error1 << "\n";


          Simulation sim(algebraic_models[i], period);
          Simulation sim_o(original_models[i], period);

          sim.SetTerminateOnConvergence(false);
          sim_o.SetTerminateOnConvergence(false);
          // sim.SetIKrBlock(IKrBlock);
          // sim_o.SetIKrBlock(IKrBlock);

          std::stringstream filename_base;
          filename_base << model_name+"_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block";

          /*  Write initial paces to file */
          sim.WritePaceToFile("TestAlgebraicVoltage", filename_base.str() + "_analytic_initial");
          sim_o.WritePaceToFile("TestAlgebraicVoltage", filename_base.str() + "_initial");

          sim.RunPaces(paces);
          sim_o.RunPaces(paces);

          /*  Write end paces to file */
          sim.WritePaceToFile("TestAlgebraicVoltage", filename_base.str() + "_analytic_final");
          sim_o.WritePaceToFile("TestAlgebraicVoltage", filename_base.str() + "_final");

          original_vars = original_models[i]->GetStdVecStateVariables();
          original_vars.erase(original_vars.begin());

          const double error2 = mrms(original_vars, algebraic_models[i]->GetStdVecStateVariables());

          std::cout << "error is " << error2 << "\n";
          assert(error1<1);

          /* Print APD90s */
          std::cout << "APD90s are " << sim.GetApd(90,false) << " " << sim_o.GetApd(90, false) << "\n";

          sim_o.SetStateVariables(original_initial_states);
          sim.SetStateVariables(algebraic_initial_states);
        }
      }
    }
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

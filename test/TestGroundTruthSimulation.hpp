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
#include "CommandLineArguments.hpp"

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
#include "ToRORd_dynCl_endo_analytic_voltageCvode.hpp"


/* Run the models under different scenarios and output:
   - All variables over the final pace
   - The APD90
   - Terminal state variables
   to separate files.
 */

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  const int paces = 10000;
  void TestRunSimulation()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    const std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directories("/home/" + username + "/testoutput/");

    std::vector<double> periods = {1000, 500, 750, 1250};
    std::vector<double> IKrBlocks = {0, 0.25, 0.5};

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;

    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellhund_rudy_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::make_shared<ToRORd_dynCl_endo_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));

    /* If "--models" option is provided, use only those models that are specified */
    if(CommandLineArguments::Instance()->OptionExists("-models")){
      std::vector<std::string> model_names;
      model_names = CommandLineArguments::Instance()->GetStringsCorrespondingToOption("-models");
        TS_ASSERT(model_names.size()>0);
    }
    std::vector<boost::shared_ptr<AbstractCvodeCell>> new_models;
      for(auto name : model_names){
        /* Find all models with the name that has been provided */
        std::list<boost::shared_ptr<AbstractCvodeCell>>> found_models = std::find_if(models.begin(), models.end(), [&](boost::shared_ptr<AbstractCvodeCell>m)->{return m->GetSystemInformation()->GetSystemName()==name;});
        TS_ASSERT(found_models.size()==1);
        new_models.push_back(found_models.front());
      }
      /* new_models contains all of the models which have been specified - use this instead of models */
      models = new_models;
  }

    for(auto model : models){
      const N_Vector initial_states = model->GetStateVariables();
      for(double period : periods){
        for(double IKrBlock : IKrBlocks){
          ComputeGroundTruth(model, period, IKrBlock);
          model->SetStateVariables(initial_states);
        }
      }
    }

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
#ifdef CHASTE_CVODE
  void ComputeGroundTruth(boost::shared_ptr<AbstractCvodeCell> model, double period, double IKrBlock){
    const std::string username = std::string(getenv("USER"));
    const std::string CHASTE_TEST_OUTPUT = std::string(getenv("CHASTE_TEST_OUTPUT"));
    const std::string model_name = model->GetSystemInformation()->GetSystemName();
    std::cout << "Testing " << model_name << " with period " << int(period) <<" and IKrBlock "<< IKrBlock << "\n";
    std::stringstream dirname;
    dirname << model_name+"_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";
    std::cout << dirname.str() << std::endl;
    boost::filesystem::path dir(dirname.str());
    dir = boost::filesystem::path(CHASTE_TEST_OUTPUT) / dir;

    // Initialise simulation with fine tolerances
    Simulation simulation(model, period, "", 1e-12, 1e-12);

    // Turn off convergence criteria
    simulation.SetTerminateOnConvergence(false);

    // Set Gkr if it exists
    std::vector<std::string> param_names = model->GetSystemInformation()->rGetParameterNames();
    double default_GKr = DOUBLE_UNSET;
    bool set_GKr = false;

    const std::string GKrParameterName = "membrane_rapid_delayed_rectifier_potassium_current_conductance";
    default_GKr = model->GetParameter(GKrParameterName);
    model->SetParameter(GKrParameterName, default_GKr*(1-IKrBlock));

    try{
      // Run the simulation for a large number of paces
      simulation.RunPaces(paces);
    }
    catch(const Exception &ex){
      std::cout << "caught an exception after " << simulation.GetPaces() << " paces\n";
      throw(ex);
    }

    // Output the final pace
    const std::string pace_filename = "final_pace";

    simulation.WritePaceToFile(dirname.str(), pace_filename);

    // Output the APD90 of the final pace
    const std::string apd_filename = "final_apd90.dat";
    std::ofstream apd_file(dirname.str() + apd_filename);
    apd_file << simulation.GetApd(90) << "\n";
    apd_file.close();

    // Print final mrms
    std::cout << "final mrms is " << simulation.GetMrms(false) << "\n";

    simulation.WriteStatesToFile(dir, "final_states.dat");
    model->SetParameter(GKrParameterName, default_GKr);
  }
#endif
};

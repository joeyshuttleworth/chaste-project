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

<<<<<<< HEAD

=======
#include "Simulation.hpp"

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"

>>>>>>> e40a0a2dda788944b67c844ec4db5252a0b39e60
/* Run the models under different scenarios and output:
   - All variables over the final pace
   - The APD90
   - Terminal state variables
   to separate files.
 */

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  const int default_max_paces = 25000;
  void TestRunSimulation()
  {
#ifdef CHASTE_CVODE
    int paces = get_max_paces();
    // Use default_max_paces if get_max_paces returns DOUBLE_UNSET
    paces = paces==DOUBLE_UNSET?default_max_paces:paces;

    const std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directories("/home/" + username + "/testoutput/");

    auto models = get_models();

    // List the names of the models we're testing

    std::cout << "Testing models:\n";
    for(auto model : models){
      std::cout << model->GetSystemInformation()->GetSystemName() << "\n";
    }

    auto IKrBlocks = get_IKr_blocks();
    auto periods = get_periods();

    for(auto model : models){
      const N_Vector initial_states = model->GetStateVariables();
      for(double period : periods){
        for(double IKrBlock : IKrBlocks){
          ComputeGroundTruth(paces, model, period, IKrBlock);
          model->SetStateVariables(initial_states);
        }
      }
    }

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
#ifdef CHASTE_CVODE
  void ComputeGroundTruth(int paces, boost::shared_ptr<AbstractCvodeCell> model, double period, double IKrBlock){
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
    simulation.SetIKrBlock(IKrBlock);

    // Turn off convergence criteria
    simulation.SetTerminateOnConvergence(false);

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
  }
#endif
};

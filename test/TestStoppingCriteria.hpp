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

class TestStoppingCriteria : public CxxTest::TestSuite
{
public:
  const int default_max_paces = 200000;
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
          Simulation sim(model, period, "", 1e-8, 1e-8);
          const std::string model_name = model->GetSystemInformation()->GetSystemName();
          sim.SetIKrBlock(IKrBlock);
          sim.RunPaces(paces);

          const int paces_taken = sim.GetPaces();

          if(sim.HasTerminated())
            std::cout << model_name << " with period " << period << " and IKrBlock " << IKrBlock << "finished after " << paces_taken << " paces\n";

          TS_ASSERT(sim.HasTerminated());
          model->SetStateVariables(initial_states);
        }
      }
    }

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
#ifdef CHASTE_CVODE
#endif
};

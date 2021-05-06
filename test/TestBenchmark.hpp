#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include "Simulation.hpp"
#include <boost/filesystem.hpp>
#include <fstream>
#include "SmartSimulation.hpp"
#include "SimulationTools.hpp"

class TestBenchmark : public CxxTest::TestSuite
{
private:
  std::ofstream output_file;
  int baseline_score = 0;

  std::vector<int> buffer_sizes = {250, 25, 50, 100, 500, 750, 1000, 2000};
  std::vector<double> extrapolation_constants ={1, 0.1, 0.5, 0.75, 0.9, 1.05, 1.1, 1.25};
  const unsigned int max_paces = 100000;
  // Get max jumps parameter
  int max_jumps;
public:

  void TestBenchmarkRun(){
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    std::vector<double> periods = get_periods();
    std::vector<double> IKrBlocks = get_IKr_blocks();

    auto models = get_models("algebraic");

    max_jumps = get_max_jumps();
    max_jumps = max_jumps==INT_UNSET?2000:max_jumps;

    const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));

    std::cout << "using extrapolation_constants: ";

    for(auto e_c : extrapolation_constants)
      std::cout << e_c << " ";

    std::cout << "and buffer sizes: ";
    for(auto bs : buffer_sizes)
      std::cout << bs << " ";

    std::cout << "and suffix " << suffix << "\n";

    boost::filesystem::path directory = (boost::filesystem::path(test_dir) / (boost::filesystem::path("TestBenchmark"+suffix)));
    // Remake the empty directoy
    std::cout << "outputting to directory: " << directory.string() << "\n";
    boost::filesystem::create_directories(directory.string());

    for(auto model : models){
      // Construct and destruct a Simulation to prevent any discrepencies
      // Without this there seems to be some small differences in the output for
      // small numbers of max_paces. Probably due to some internal Cvode settings
      {
        Simulation sim(model);
      }
      const std::string model_name = model->GetSystemInformation()->GetSystemName();
      const std::string filepath = (directory / boost::filesystem::path(model_name + "_results.dat")).string();
      // std::cout << "filepath: " << filepath << "\n";

      std::ofstream output_file(filepath);

        output_file << "what_modified model_name buffer_size extrapolation_constant ic_period ic_block period IKrBlock score jumps_used APD90 last_mrms reference_mrms reference_trace_mrms reference_2_norm reference_trace_2_norm ";

        for(auto var_name : model->GetSystemInformation()->rGetStateVariableNames())
          output_file << var_name << " ";

        output_file << "\n";

        for(auto buffer_size : buffer_sizes){
          for(auto extrapolation_constant : extrapolation_constants){
            OutputScores(model, periods, IKrBlocks, extrapolation_constant, buffer_size, output_file);
          }
        }
        // Next, get the score with no extrapolation
        OutputScores(model, periods, IKrBlocks, 0, 100, output_file);
      }
  }

  void OutputScores(boost::shared_ptr<AbstractCvodeCell> model, std::vector<double> periods, std::vector<double> IKrBlocks, double extrapolation_constant, unsigned int buffer_size, std::ofstream& output_file){
    /*  Output a scores for each setting. Test 3 times, first varying both IKrBlock and pacing then once with each individually.  */
    output_file.precision(20);
    const std::string model_name = model->GetSystemInformation()->GetSystemName();
    for(auto period : periods){
      for(auto IKrBlock : IKrBlocks){
        OutputScore("period_IKrBlock", model, period, IKrBlock, extrapolation_constant, buffer_size, output_file);
        OutputScore("period", model, period, IKrBlock, extrapolation_constant, buffer_size, output_file);
        OutputScore("IKrBlock", model, period, IKrBlock, extrapolation_constant, buffer_size, output_file);
      }
    }
  }
};

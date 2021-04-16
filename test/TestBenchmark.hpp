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

class TestBenchmark : public CxxTest::TestSuite
{
private:
  std::ofstream output_file;
  std::string username;
  int baseline_score = 0;

  const std::vector<unsigned int> buffer_sizes = {100, 250, 500, 750, 1000, 1250, 1500, 2000};//{25, 50, 100, 150, 200, 300 ,400};
  const std::vector<double> extrapolation_constants ={0.25, 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5};
  const unsigned int max_paces = 100000;
public:

  void TestBenchmarkRun(){
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    std::vector<double> periods = get_periods();
    std::vector<double> IKrBlocks = get_IKr_blocks();

    auto models = get_analytic_models();

    const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));


    for(auto model : models){
        const std::string model_name = model->GetSystemInformation()->GetSystemName();
        std::string filepath = (boost::filesystem::path(test_dir) / boost::filesystem::path("TestBenchmark/" + model_name + "_results.dat")).string();
        boost::filesystem::create_directories(filepath);
        std::ofstream output_file(filepath);

        output_file << "what_modified model_name buffer_size extrapolation_constant period IKrBlock score APD90";

        // First get the score with no extrapolation
        OutputScore(model, periods, IKrBlocks, 0, 100, output_file);

        for(auto buffer_size : buffer_sizes){
          for(auto extrapolation_constant : extrapolation_constants){
            OutputScore(model, periods, IKrBlocks, extrapolation_constant, buffer_size, output_file);
          }
        }
      }
  }

  void OutputScore(boost::shared_ptr<AbstractCvodeCell> model, std::vector<double> periods, std::vector<double> IKrBlocks, double extrapolation_constant, unsigned int buffer_size, std::ofstream& output_file){
    /*  Output a scores for each setting. Test 3 times, first varying both IKrBlock and pacing then once with each individually.  */
    const std::string model_name = model->GetSystemInformation()->GetSystemName();
    for(auto period : periods){
      for(auto IKrBlock : IKrBlocks){
        // Change both IKrBlock and pacing frequency
        double ic_period = period==1000?500:1000;
        double ic_block  = IKrBlock==0?0.5:0;

        // Print the performance of this model with these settings
        output_file << "period_IKrBlock " << model_name << " " << buffer_size << " " << extrapolation_constant << " " << period << " " << IKrBlock << " ";
        output_file << RunModel(model, ic_period, ic_block, period, IKrBlock, buffer_size, extrapolation_constant) << " ";
        // Output APD90
        output_file << Simulation(model).GetApd(90) << "\n";

        // Change both only pacing frequency
        ic_period = period==1000?500:1000;
        ic_block  = IKrBlock;

        // Print the performance of this model with these settings
        output_file << "period " << model_name << " " << buffer_size << " " << extrapolation_constant << " " << period << " " << IKrBlock << " ";
        output_file << RunModel(model, ic_period, ic_block, period, IKrBlock, buffer_size, extrapolation_constant) << " ";
        // Output APD90
        output_file << Simulation(model).GetApd(90) << "\n";

        // Change only IKrBlock
        ic_period = period;
        ic_block  = IKrBlock==0?0.5:0;

        // Print the performance of this model with these settings
        output_file << "IKrBlock " << model_name << " " << buffer_size << " " << extrapolation_constant << " " << period << " " << IKrBlock << " ";
        output_file << RunModel(model, ic_period, ic_block, period, IKrBlock, buffer_size, extrapolation_constant) << " ";
        // Output APD90
        output_file << Simulation(model, period, "" ).GetApd(90) << "\n";

      }
    }
  }


  unsigned int RunModel(boost::shared_ptr<AbstractCvodeCell> model, double ic_period, double ic_IKrBlock, double period, double IKrBlock, unsigned int buffer_size, double extrapolation_constant){
    const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));

    const std::string model_name = model->GetSystemInformation()->GetSystemName();

    std::cout << "For model " << model_name << "\n";

    std::stringstream dirname;
    dirname << "/" << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";

    std::stringstream input_dirname_ss;
    input_dirname_ss << model_name+"_" << std::to_string(int(ic_period)) << "ms_" << int(100*ic_IKrBlock)<<"_percent_block/";
    const std::string input_dirname = (test_dir / boost::filesystem::path(input_dirname_ss.str())).string();

    const std::string input_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("final_states.dat")).string();

    SmartSimulation smart_simulation(model, period, input_path, 1e-8, 1e-8);
    smart_simulation.SetIKrBlock(IKrBlock);
    smart_simulation.RunPaces(max_paces);

    unsigned int paces = smart_simulation.GetPaces();
    std::cout << "took " << paces << " paces\n";
    TS_ASSERT(paces+2<max_paces);

    return paces;
  }
};

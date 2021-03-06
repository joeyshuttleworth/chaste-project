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
  const double threshold = 1.8e-07;
  std::ofstream output_file;
  std::string username;
  int baseline_score = 0;

  const std::vector<unsigned int> buffer_sizes = {100};//{25, 50, 100, 150, 200, 300 ,400};
  const std::vector<double> extrapolation_constants ={1}; //{0.5, 0.75, 0.9, 1, 1.1};
  const unsigned int max_paces = 2000;
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

        output_file << "model_name buffer_size extrapolation_constant period IKrBlock score";

        for(auto buffer_size : buffer_sizes){
          for(auto extrapolation_constant : extrapolation_constants){
            for(auto period : periods){
              for(auto IKrBlock : IKrBlocks){
                output_file << model_name << " " << buffer_size << " " << extrapolation_constant << " " << period << " " << IKrBlock << " ";
                output_file << RunModel(model, period, IKrBlock, buffer_size, extrapolation_constant) << "\n";
              }
            }
          }
        }
      }
  }

  unsigned int RunModel(boost::shared_ptr<AbstractCvodeCell> model, double period, double IKrBlock, unsigned int buffer_size, double extrapolation_constant){
    const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));
    const double default_GKr = model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance");
    model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", (1 - IKrBlock)*default_GKr);

    const std::string model_name = model->GetSystemInformation()->GetSystemName();

    std::cout << "For model " << model_name << "\n";
    model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", (1-IKrBlock)*default_GKr);

    const double starting_period = period==1000?500:1000;
    const double starting_block  = period==0?0.5:0;

    std::stringstream dirname;
    dirname << "/" << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";

    std::stringstream input_dirname_ss;
    input_dirname_ss << model_name+"_" << std::to_string(int(starting_period)) << "ms_" << int(100*starting_block)<<"_percent_block/";
    const std::string input_dirname = (test_dir / boost::filesystem::path(input_dirname_ss.str())).string();

    const std::string input_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("final_states.dat")).string();

    SmartSimulation smart_simulation(model, period, input_path, 1e-8, 1e-8);
    smart_simulation.SetThreshold(threshold);
    smart_simulation.RunPaces(max_paces);

    unsigned int paces = smart_simulation.GetPaces();
    std::cout << "took " << paces << " paces\n";
    // TS_ASSERT(paces+2<max_paces);

    model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr);
    return paces;
  }
};

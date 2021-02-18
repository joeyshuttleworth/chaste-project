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

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"

/*Output the total number of paces to reach limiting behaviour by SmartSimulation over all models for different choices of buffer_size and extrapolation_constant*/

class TestBenchmark : public CxxTest::TestSuite
{
private:
  const double threshold = 1.8e-07;
  const unsigned int paces = 5000;
  std::ofstream output_file;
  std::string username;

  const std::vector<unsigned int> buffer_sizes = {25, 50, 100, 150, 200, 300 ,400};
  const std::vector<double> extrapolation_constants = {0.5, 0.75, 0.9, 1, 1.1};
  const int max_paces = 10000;
public:

  void TestBenchmarkRun(){
    const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));
    std::ofstream output_file((boost::filesystem::path(test_dir) / boost::filesystem::path("TestBenchark/results.dat")).string());
    output_file << "brute force method takes " << RunModels(0,0) << "\n";

    for(auto buffer_size : buffer_sizes){
      for(auto extrapolation_constant : extrapolation_constants){
        const unsigned int score = RunModels(buffer_size, extrapolation_constant);
        output_file << buffer_size << ", \t" << extrapolation_constant << " \t" << score << "\n";
      }
    }
    output_file.close();
  }


  unsigned int RunModels(unsigned int buffer_size, double extrapolation_constant){
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    std::vector<double> periods = {500, 750, 1000, 1250};
    std::vector<double> IKrBlocks = {0, 0.25, 0.5, 0.75};
    unsigned int score=0;
    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus)));

    for(auto model : models){
      for(auto period : periods){
        for(auto IKrBlock : IKrBlocks){
          score += RunModel(model, period, IKrBlock, buffer_size, extrapolation_constant);
        }
      }
    }
    return score;
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
    std::cout << "Testing model: " << model_name << " with period " << period << "ms and IKrBlock " << IKrBlock << "\n";

    const std::string input_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("final_states.dat")).string();

    SmartSimulation smart_simulation(model, period, input_path, 1e-8, 1e-8);

    smart_simulation.RunPaces(max_paces);

    unsigned int paces = smart_simulation.GetPaces();
    assert(paces<max_paces-10);

    model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr);
    return paces;
  }
};

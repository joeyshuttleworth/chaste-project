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
  int baseline_score = 0;

  const std::vector<unsigned int> buffer_sizes = {250, 50, 100, 500, 750, 1000, 2000};
  const std::vector<double> extrapolation_constants ={1, 0.1, 0.5, 0.75, 0.9, 1.1, 1.25};
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

    boost::filesystem::path directory = (boost::filesystem::path(test_dir) / boost::filesystem::path("TestBenchmark/"));
    // Remake the empty directoy
    boost::filesystem::create_directories(directory);

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

        output_file << "what_modified model_name buffer_size extrapolation_constant ic_period ic_block period IKrBlock score jumps_used APD90 last_mrms reference_mrms reference_trace_mrms reference_2_norm reference_trace_2_norm\n";

        for(auto buffer_size : buffer_sizes){
          for(auto extrapolation_constant : extrapolation_constants){
            OutputScores(model, periods, IKrBlocks, extrapolation_constant, buffer_size, output_file);
          }
        }
        // Next, get the score with no extrapolation
        OutputScores(model, periods, IKrBlocks, 0, 100, output_file);
      }
  }

  void OutputScore(std::string to_change, boost::shared_ptr<AbstractCvodeCell> model, double period, double IKrBlock, double extrapolation_constant, unsigned int buffer_size, std::ofstream& output_file){

    const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));
    const std::string model_name = model->GetSystemInformation()->GetSystemName();
    std::stringstream reference_path;
    reference_path << test_dir.string() << "/" << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/final_states.dat";

    std::vector<double> reference_states = LoadStatesFromFile(reference_path.str());

    double ic_period;
    double ic_block;
    if(to_change=="period_IKrBlock"){
      // Change both IKrBlock and pacing frequency
      ic_period = period==1000?500:1000;
      ic_block  = IKrBlock==0?0.5:0;
    }
    else if(to_change=="period"){
      ic_period = period==1000?500:1000;
      ic_block  = IKrBlock;
    }
    else if(to_change=="IKrBlock"){
      ic_period = period;
      ic_block  = IKrBlock==0?0.5:0;
    }
    else{
      EXCEPTION("string argument doesn't match any option");
    }

    // Print the performance of this model with these settings
    output_file << "period_IKrBlock " << model_name << " " << buffer_size << " " << extrapolation_constant << " " << ic_period << " " << ic_block << " " << period << " " << IKrBlock << " ";

    try{
      auto sim = RunModel(model, ic_period, ic_block, period, IKrBlock, buffer_size, extrapolation_constant);
      output_file << sim->GetPaces() << " " << sim->GetNumberOfJumps() << " ";
      // Output APD90
      output_file << Simulation(model).GetApd(90) << " " << sim->GetMRMS() << " ";

      //Compute error measures with the reference solution
      auto states = sim->GetStateVariables();
      const double two_norm = TwoNorm(reference_states, states);
      const double reference_mrms = mrms(reference_states, states);

      auto pace = sim->GetPace().rGetSolutions();

      // Compute reference pace
      Simulation reference_sim(model, period, reference_path.str());
      auto reference_pace = reference_sim.GetPace().rGetSolutions();

      double reference_trace_two_norm = TwoNormTrace(reference_pace, pace);
      double reference_trace_mrms = mrmsTrace(reference_pace, pace);

      output_file << reference_mrms << " " << reference_trace_mrms << " " << two_norm << " " << reference_trace_two_norm << "\n";


    }
    catch(Exception& e){
      std::cout << "Something went wrong when running the model! " << e.GetMessage() << "\n";;
      output_file << "-1 -1 -1 -1 -1 -1 -1 -1\n";
      // Output APD90
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

  std::shared_ptr<SmartSimulation> RunModel(boost::shared_ptr<AbstractCvodeCell> model, double ic_period, double ic_IKrBlock, double period, double IKrBlock, unsigned int buffer_size, double extrapolation_constant){

    const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));

    const std::string model_name = model->GetSystemInformation()->GetSystemName();

    std::cout << "For model " << model_name << " with  n = " << buffer_size << " and e_c = " << extrapolation_constant << "\n";

    std::stringstream dirname;
    dirname << "/" << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";

    std::stringstream input_dirname_ss;
    input_dirname_ss << model_name+"_" << std::to_string(int(ic_period)) << "ms_" << int(100*ic_IKrBlock)<<"_percent_block/";
    const std::string input_dirname = (test_dir / boost::filesystem::path(input_dirname_ss.str())).string();

    const std::string input_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("final_states.dat")).string();

    std::shared_ptr<SmartSimulation> smart_simulation = std::make_shared<SmartSimulation>(model, period, input_path, 1e-8, 1e-8, buffer_size, extrapolation_constant);
    smart_simulation->SetIKrBlock(IKrBlock);
    smart_simulation->SetMaxJumps(max_jumps);

    int paces_to_run = get_max_paces();
    paces_to_run = paces_to_run==INT_UNSET?max_paces:paces_to_run;

    smart_simulation->RunPaces(paces_to_run);

    unsigned int paces = smart_simulation->GetPaces();

    std::cout << "took " << paces << " paces\n";
    TS_ASSERT(paces+2<max_paces);

    return smart_simulation;
  }
};

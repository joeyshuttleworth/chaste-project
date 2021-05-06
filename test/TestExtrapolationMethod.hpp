#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include "Simulation.hpp"
#include "SmartSimulation.hpp"
#include "SimulationTools.hpp"
#include "CellProperties.hpp"

#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>

class TestExtrapolationMethod : public CxxTest::TestSuite
{
private:
  const unsigned int default_paces = 200000;

  bool CompareMethods(double e_c, int bs, int paces, boost::shared_ptr<AbstractCvodeCell> brute_force_model, boost::shared_ptr<AbstractCvodeCell> smart_model, double ic_period, double ic_block, double period, double IKrBlock){
    std::string username = std::string(getenv("USER"));
    std::string CHASTE_TEST_DIR = std::string(getenv("CHASTE_TEST_OUTPUT"));
    const std::string model_name = brute_force_model->GetSystemInformation()->GetSystemName();
    const std::string dirname  = "/home/"+username+"/testoutput/"+model_name+"/TestExtrapolationMethod";
    boost::filesystem::create_directories(dirname);

    std::cout << "Testing " << model_name << " " << ic_period << " " << ic_block << " " << period << " " << IKrBlock << "\n";

    std::stringstream input_dirname;
    input_dirname << CHASTE_TEST_DIR << "/" << model_name << "_" << std::to_string(int(ic_period)) << "ms_" << int(100*ic_block)<<"_percent_block/";

    const std::string input_file = input_dirname.str() + "final_states.dat";

    // Uses a method to extrapolate to the steady state
    SmartSimulation smart_simulation(smart_model, period, input_file, 1e-08, 1e-08, bs, e_c, "/home/chaste/testoutput/" + model_name + "/TestExtrapolationMethod");
    // Runs the model without using this method
    Simulation simulation(brute_force_model, period, input_file, 1e-08, 1e-08);

    simulation.SetIKrBlock(IKrBlock);;
    smart_simulation.SetIKrBlock(IKrBlock);

    // Setup directories for output
    std::cout << "-------------------------------\n\n\nTesting " << model_name  << "\n";

    //Open files to output states to at the start/end of each pace
    std::ofstream smart_output_file, brute_output_file;
    smart_output_file.open("/home/"+username+"/testoutput/"+model_name+"/TestExtrapolationMethod/smart.dat");
    brute_output_file.open("/home/"+username+"/testoutput/"+model_name+"/TestExtrapolationMethod/bruteforce.dat");
    smart_output_file << std::setprecision(20);
    brute_output_file << std::setprecision(20);

    std::vector<std::string> state_names = smart_model->rGetStateVariableNames();

    // Set up header line
    smart_output_file << "pace mrms ";
    brute_output_file << "pace mrms ";
    for(unsigned int i = 0; i < state_names.size(); i++){
      smart_output_file << state_names[i] << " ";
      brute_output_file << state_names[i] << " ";
    }
    smart_output_file << "\n";
    brute_output_file << "\n";

    // Calculate difference in APD90s using initial conditions
    {
      const double smart_apd = smart_simulation.GetApd(90);
      const double apd_difference = smart_apd - simulation.GetApd(90);
      std::cout << "Difference in APD90s using initial conditions " << apd_difference << "\n";
    }

    // Run the simulations until they finish
    bool brute_finished = false;
    bool smart_finished = false;
    for(int j = 0; j < paces; j++){
      if(!smart_finished){
        if(smart_simulation.RunPace()){
          std::cout << "Model " << model_name << " period " << period << " extrapolation method finished after " << j << " paces \n";
          smart_finished = true;
        }
        std::vector<double> state_vars = smart_model->GetStdVecStateVariables();
        smart_output_file << j << " ";
        smart_output_file << smart_simulation.GetMRMS() << " ";
        for(unsigned int i = 0; i < state_vars.size(); i++){
          smart_output_file << state_vars[i] << " ";
        }
        smart_output_file << "\n";
      }

      if(!brute_finished){
        if(simulation.RunPace()){
          std::cout << "Model " << model_name << " period " << period << " brute force method finished after " << j << " paces \n";
          brute_finished = true;
        }
        std::vector<double> state_vars = brute_force_model->GetStdVecStateVariables();
        brute_output_file << j << " ";
        brute_output_file << simulation.GetMRMS() << " ";
        //Don't print membrane_voltage (usually the first state variable)
        for(unsigned int i = 0; i < state_vars.size(); i++){
          brute_output_file << state_vars[i] << " ";
        }
        brute_output_file << "\n";
      }
      if(smart_finished && brute_finished)
        break;
    }

    std::vector<double> brute_states = simulation.GetStateVariables();
    std::vector<double> smart_states = smart_simulation.GetStateVariables();

    if(brute_states.size() == smart_states.size()+1){
      brute_states.erase(brute_states.begin());
    }

    /*Check that the methods have converged to the same place*/
    // First calculate and output the mrms error between the solutions
    double mrms_difference = mrms(brute_states, smart_states);
    std::cout << "MRMS between solutions is " << mrms_difference << "\n";

    // Calculate difference in APD90s
    const double smart_apd = smart_simulation.GetApd(90);
    const double apd_difference = smart_apd - simulation.GetApd(90);
    std::cout << "Difference in APD90s " << apd_difference << "\n";

    // Compare smart apd with reference version
    std::stringstream apd_file_ss;
    const std::string CHASTE_TEST_OUTPUT = getenv("CHASTE_TEST_OUTPUT");
    boost::filesystem::path apd_filepath(CHASTE_TEST_OUTPUT);
    apd_file_ss << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block";
    apd_filepath = apd_filepath / apd_file_ss.str() / "final_apd90.dat";
    try{
      std::ifstream apd_file(apd_filepath.string());
      std::string line;
      std::getline(apd_file, line);
      const double reference_apd = std::stod(line);
      std::cout << "Difference between apd and reference apd " << smart_apd - reference_apd << "\n";
    }
    catch(const std::exception& e){
      std::cout << "Couldn't find reference voltage at " << apd_filepath.string() << " - ignoring\n" << e.what() << "\n---------------------------\n";
    }


    std::string smart_filename = "smart_final_pace.dat";
    std::string brute_filename = "brute_final_pace.dat";
    std::cout << "Outputting final paces as " << smart_filename << " and " << brute_filename << "\n";
    const std::string final_pace_dir = model_name + "/TestExtrapolationMethod/" + std::to_string(int(period));
    simulation.WritePaceToFile(final_pace_dir + "smart_pace", brute_filename);
    smart_simulation.WritePaceToFile(final_pace_dir + "brute_pace" , smart_filename);

    std::ofstream voltage_trace_file(dirname + "/voltage_traces.dat");
    voltage_trace_file << "smart_pace brute_pace\n";
    std::vector<double> smart_voltages = smart_simulation.GetVoltageTrace();
    std::vector<double> brute_voltages  =  simulation.GetVoltageTrace();

    for(unsigned int i = 0; i < smart_voltages.size(); i++){
      voltage_trace_file << smart_voltages[i] << " " << brute_voltages[i] << "\n";
    }

    TS_ASSERT_LESS_THAN(mrms_difference, 1e-3);
    TS_ASSERT(smart_finished && brute_finished);

    return brute_finished&&smart_finished;
  }


public:
  void TestExtrapolationRun()
  {
#ifdef CHASTE_CVODE
    int max_paces = get_max_paces();
    auto extrapolation_constants = get_extrapolation_constants();
    auto buffer_sizes = get_buffer_sizes();
    auto periods = get_periods();
    auto blocks = get_IKr_blocks();
    auto models = get_analytic_models();
    auto brute_models = get_analytic_models();

    std::cout << "Running each scenario for " << max_paces << " paces.\n";
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    /* Compare each model with a version modified for the extrapolation method */

    // TS_ASSERT(original_models.size()==algebraic_models.size());

    const double ic_period = 1000;
    const double ic_block  = 0;
    int i=0;
    for(auto model : models){
      auto brute_model = brute_models[i];
      for(auto e_c : extrapolation_constants){
        for(auto bs : buffer_sizes){
          for(auto period : periods){
            for(auto block : blocks){
              // Just change IKrBlock
              CompareMethods(e_c, bs, max_paces, brute_model, model, ic_period, ic_block,  period,  block);
            }
          }
        }
      }
    }
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

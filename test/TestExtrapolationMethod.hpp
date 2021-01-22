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
#include "CellProperties.hpp"

#include <boost/filesystem.hpp>
#include <fstream>

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"


class TestExtrapolationMethod : public CxxTest::TestSuite
{
private:
  const unsigned int buffer_size = 100;
  const double extrapolation_coefficient = 1;
  const unsigned int paces = 2000;

  bool CompareMethodsPeriod(boost::shared_ptr<AbstractCvodeCell> brute_force_model, boost::shared_ptr<AbstractCvodeCell> smart_model){
    std::string username = std::string(getenv("USER"));
    const std::string model_name = brute_force_model->GetSystemInformation()->GetSystemName();
    const std::string dirname  = "/home/"+username+"/testoutput/"+model_name+"/TestExtrapolationMethod";
    boost::filesystem::create_directories(dirname);
    const int period = 750;
    // Uses a method to extrapolate to the steady state
    SmartSimulation smart_simulation(smart_model, period, "", 1e-8, 1e-8, 200, 1, "/home/chaste/testoutput/" + model_name + "/TestExtrapolationMethod");
    // Runs the model without using this method
    Simulation      simulation(brute_force_model, period, "", 1e-8, 1e-8);

    smart_simulation.SetThreshold(0);
    simulation.SetThreshold(0);

    // Setup directories for output
    std::cout << "Testing " << model_name  << "\n";

    //Open files to output states to at the start/end of each pace
    std::ofstream smart_output_file, brute_output_file;
    smart_output_file.open("/home/"+username+"/testoutput/"+model_name+"/TestExtrapolationMethod/smart.dat");
    brute_output_file.open("/home/"+username+"/testoutput/"+model_name+"/TestExtrapolationMethod/bruteforce.dat");

    std::vector<std::string> state_names = smart_model->rGetStateVariableNames();

    // Set up header line
    smart_output_file << "pace ";
    brute_output_file << "pace ";
    for(unsigned int i = 0; i < state_names.size(); i++){
      smart_output_file << state_names[i] << " ";
      brute_output_file << state_names[i] << " ";
    }
    smart_output_file << "\n";
    brute_output_file << "\n";

    // Run the simulations until they finish
    bool brute_finished = false;
    bool smart_finished = false;
    for(unsigned int j = 0; j < 5000; j++){
      if(!smart_finished){
        if(smart_simulation.RunPace()){
          std::cout << "Model " << model_name << " period " << period << " extrapolation method finished after " << j << " paces \n";
          smart_finished = true;
        }
        std::vector<double> state_vars = smart_model->GetStdVecStateVariables();
        smart_output_file << j << " ";
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
        //Don't print membrane_voltage (usually the first state variable)
        for(unsigned int i = 1; i < state_vars.size(); i++){
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


    for(unsigned int i = 0; i < smart_states.size() && i < brute_states.size(); i++){
      std::cout << brute_states[i] << " " << smart_states[i] << "\n";
    }

    /*Check that the methods have converged to the same place*/
    // First calculate and output the mrms error between the solutions
    double mrms_difference = mrms(brute_states, smart_states);
    std::cout << "MRMS between solutions is " << mrms_difference << "\n";

    // Calculate difference in APD90s
    // Currently broken
    const double apd_difference = smart_simulation.GetApd(90) - simulation.GetApd(90);
    std::cout << "Difference in APD90s " << apd_difference << "\n";

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

  void CompareMethodsKrBlock(boost::shared_ptr<AbstractCvodeCell> brute_model, boost::shared_ptr<AbstractCvodeCell> smart_model){
    const std::string model_name = brute_model->GetSystemInformation()->GetSystemName();

    std::cout << "Testing " << model_name  << " with 50% block of IKr\n";

    //set Gkr parameter
    double default_GKr = brute_model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance");
    brute_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr*0.5);
    smart_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr*0.5);

    Simulation simulation(brute_model, 1000);
    SmartSimulation smart_simulation(smart_model, 1000);

    simulation.RunPaces(paces);
    smart_simulation.RunPaces(paces);

    std::vector<double> brute_states = simulation.GetStateVariables();
    std::vector<double> smart_states = smart_simulation.GetStateVariables();

    // drop voltage (assuming it's the first variable)
    if(brute_states.size() == smart_states.size()+1){
      assert(brute_model->GetSystemInformation()->rGetStateVariableNames()[0] == "membrane_voltage");
      brute_states.erase(brute_states.begin());
    }

    double mrms_difference = mrms(brute_states, smart_states);

    const double smart_apd = smart_simulation.GetApd(90);
    smart_simulation.SetStateVariables(brute_states);
    const double brute_apd = smart_simulation.GetApd(90);
    smart_simulation.SetStateVariables(smart_states);
    std::cout << "Difference in APD90s " << smart_apd - brute_apd << "\n";

    std::cout << "MRMS between solutions is " << mrms_difference << "\n";

    const std::string output_dir = model_name+"/TestExtrapolationMethod/IKrBlock";
    simulation.WritePaceToFile(output_dir+"/brute", "brute_pace");
    smart_simulation.WritePaceToFile(output_dir+"/smart", "smart_pace");

    TS_ASSERT_LESS_THAN(mrms_difference, 1e-3);
    TS_ASSERT(smart_simulation.IsFinished() && simulation.IsFinished());

    // reset GKr parameter
    brute_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr);
    smart_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr);
  }

public:
  void TestTusscherSimulation()
  {
#ifdef CHASTE_CVODE

    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    /* Compare each model with a version modified for the extrapolation method */

    // ohara_rudy_cipa_v1_2017 model
    boost::shared_ptr<AbstractCvodeCell> p_model1(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<AbstractCvodeCell> p_model2(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus));
    CompareMethodsPeriod(p_model1, p_model2);
    CompareMethodsKrBlock(p_model1, p_model2);

    // ten_tusscher_model_2006_epiFromCellMLCvode
    boost::shared_ptr<AbstractCvodeCell> p_model3(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<AbstractCvodeCell> p_model4(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus));
    CompareMethodsPeriod(p_model3, p_model4);
    CompareMethodsKrBlock(p_model3, p_model4);

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

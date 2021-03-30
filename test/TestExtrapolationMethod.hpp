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

// Non algebraic models
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_epiCvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"
#include "iyer_2004Cvode.hpp"
#include "ToRORd_dynCl_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"


// Analytic models
#include "decker_2009_analytic_voltageCvode.hpp"
#include "hund_rudy_2004_analytic_voltageCvode.hpp"
#include "iyer_2004_analytic_voltageCvode.hpp"
#include "ohara_rudy_2011_epi_analytic_voltageCvode.hpp"
#include "ohara_rudy_cipa_2017_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2006_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2004_epi_analytic_voltageCvode.hpp"
#include "ToRORd_dyn_chloride_epi_analytic_voltageCvode.hpp"

class TestExtrapolationMethod : public CxxTest::TestSuite
{
private:
  const unsigned int buffer_size = 100;
  const double extrapolation_coefficient = 1;
  const unsigned int paces = 2000;

  bool CompareMethodsPeriod(boost::shared_ptr<AbstractCvodeCell> brute_force_model, boost::shared_ptr<AbstractCvodeCell> smart_model){
    const double IKrBlock = 0;
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
    const double smart_apd = smart_simulation.GetApd(90);
    const double apd_difference = smart_apd - simulation.GetApd(90);
    std::cout << "Difference in APD90s " << apd_difference << "\n";

    // Compare smart apd with reference version
    std::stringstream apd_file_ss;
    const std::string CHASTE_TEST_OUTPUT = getenv("CHASTE_TEST_OUTPUT");
    apd_file_ss << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/apds_using_groundtruth.dat";
    std::ifstream apd_file(apd_file_ss.str());
    std::string line;
    std::getline(apd_file, line);
    const double reference_apd = std::stod(line);
    std::cout << "Difference between apd and reference apd " << smart_apd - reference_apd << "\n";

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

  void CompareMethodsIKrBlock(boost::shared_ptr<AbstractCvodeCell> brute_model, boost::shared_ptr<AbstractCvodeCell> smart_model){
    const std::string model_name = brute_model->GetSystemInformation()->GetSystemName();
    std::cout << "Testing " << model_name  << " with 50% block of IKr\n";

    Simulation simulation(brute_model, 1000);
    SmartSimulation smart_simulation(smart_model, 1000);

    simulation.SetIKrBlock(0.5);
    smart_simulation.SetIKrBlock(0.5);

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
  }

public:
  void TestTusscherSimulation()
  {
#ifdef CHASTE_CVODE

    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    /* Compare each model with a version modified for the extrapolation method */

    std::vector<boost::shared_ptr<AbstractCvodeCell>> original_models;
    original_models.push_back(boost::make_shared<CellToRORd_dynCl_epiFromCellMLCvode>(p_solver, p_stimulus));
    original_models.push_back(boost::make_shared<Celliyer_2004FromCellMLCvode>(p_solver, p_stimulus));
    original_models.push_back(boost::make_shared<Cellhund_rudy_2004FromCellMLCvode>(p_solver, p_stimulus));
    original_models.push_back(boost::make_shared<Celldecker_2009FromCellMLCvode>(p_solver, p_stimulus));
    original_models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    original_models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus)));
    original_models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_epiFromCellMLCvode(p_solver, p_stimulus)));
    original_models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus)));

    std::vector<boost::shared_ptr<AbstractCvodeCell>> algebraic_models;
    algebraic_models.push_back(boost::make_shared<CellToRORd_dyn_chloride_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Celliyer_2004_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Cellhund_rudy_2004_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Celldecker_2009_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Cellten_tusscher_2004_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Cellten_tusscher_2006_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Cellohara_rudy_2011_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    algebraic_models.push_back(boost::make_shared<Cellohara_rudy_cipa_2017_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));

    for(unsigned int i = 0; i < original_models.size(); i++){
      CompareMethodsPeriod(original_models[i], algebraic_models[i]);
      // CompareMethodsIKrBlock(original_models[i], algebraic_models[i]);
    }

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

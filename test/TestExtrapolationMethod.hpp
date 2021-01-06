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

#include <boost/filesystem.hpp>
#include <fstream>

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"


class TestGroundTruthSimulation : public CxxTest::TestSuite
{
private:
  const unsigned int buffer_size = 100;
  const double extrapolation_coefficient = 1;

  bool compareMethods(boost::shared_ptr<AbstractCvodeCell> brute_force_model, boost::shared_ptr<AbstractCvodeCell> smart_model){
    SmartSimulation smart_simulation(smart_model, 500);
    Simulation      simulation(brute_force_model, 500);

    std::string username = std::string(getenv("USER"));

    const std::string model_name = brute_force_model->GetSystemInformation()->GetSystemName();
    const std::string model_name2 = smart_model->GetSystemInformation()->GetSystemName();

    std::cout << "Testing " << model_name  << "\n";
    boost::filesystem::create_directory("/tmp/"+username);
    boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
    boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/TestExtrapolation");
    TS_ASSERT_EQUALS(model_name , model_name2);

    std::ofstream smart_output_file, brute_output_file;
    smart_output_file.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/smart.dat");
    brute_output_file.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/bruteforce.dat");

    std::vector<std::string> state_names = brute_force_model->rGetStateVariableNames();

    smart_output_file << "pace ";
    brute_output_file << "pace ";

    for(unsigned int i = 0; i < state_names.size(); i++){
      smart_output_file << state_names[i] << " ";
      brute_output_file << state_names[i] << " ";
    }
    smart_output_file << "\n";
    brute_output_file << "\n";

    /*Run the simulations*/
    bool brute_finished = false;
    bool smart_finished = false;
    const unsigned int paces  = 2000;
    for(unsigned int j = 0; j < paces; j++){
      if(!smart_finished){
        if(smart_simulation.RunPace()){
          std::cout << "Model " << model_name << " period " << 500 << " extrapolation method finished after " << j << " paces \n";
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
          std::cout << "Model " << model_name << " period " << 500 << " brute force method finished after " << j << " paces \n";
          brute_finished = true;
        }
        std::vector<double> state_vars = brute_force_model->GetStdVecStateVariables();
        brute_output_file << j << " ";
        for(unsigned int i = 0; i < state_vars.size(); i++){
          brute_output_file << state_vars[i] << " ";
        }
        brute_output_file << "\n";
      }
      if(smart_finished && brute_finished)
        break;
    }

    // /* Run the numerical voltage version */
    // numerical_voltage_simulation = Simulation(numerical_comparison_model, 500);
    // numerical_voltage_simulation.RunPaces(10000);

    std::vector<double> brute_states = brute_force_model->GetStdVecStateVariables();
    std::vector<double> smart_states = smart_model->GetStdVecStateVariables();

    for(unsigned int i = 0; i < brute_states.size(); i++){
      std::cout << brute_states[i] << "\t" << smart_states[i] << "\n";
    }

    /*Check that the methods have converged to the same place*/
    double mrms_difference = mrms(simulation.GetStateVariables(), smart_simulation.GetStateVariables());

    std::cout << "MRMS between solutions is " << mrms_difference << "\n";

    TS_ASSERT_LESS_THAN(mrms_difference, 1e-3);
    TS_ASSERT(smart_finished && brute_finished);

    return brute_finished&&smart_finished;
  }

public:
  void TestTusscherSimulation()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    boost::shared_ptr<AbstractCvodeCell> p_model1(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<AbstractCvodeCell> p_model2(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<AbstractCvodeCell> p_model3(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<AbstractCvodeCell> p_model4(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus));

    // compareMethods(p_model1, p_model2);
    compareMethods(p_model3, p_model4);

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

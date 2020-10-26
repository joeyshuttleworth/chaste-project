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

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
private:
  const unsigned int buffer_size = 100;
  const double e_c = 1;
public:
  void TestTusscherSimulation()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    boost::shared_ptr<Cellten_tusscher_model_2006_epiFromCellMLCvode> p_model3(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<Cellten_tusscher_model_2006_epiFromCellMLCvode> p_model4(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<Cellohara_rudy_cipa_v1_2017FromCellMLCvode> p_model1(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus));
    boost::shared_ptr<Cellohara_rudy_cipa_v1_2017FromCellMLCvode> p_model2(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus));

    std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);

    // int benchmark = 0;

    const std::string model_name = p_model1->GetSystemInformation()->GetSystemName();
    std::cout << "Testing " << model_name  << "\n";
    boost::filesystem::create_directory("/tmp/"+username);
    boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
    boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/TestExtrapolation");

    // const double duration = 2;

    const unsigned int paces  = 1500;
    SmartSimulation smart_simulation(p_model1, 500);
    Simulation      simulation(p_model2, 500);
    smart_simulation.Initialise(buffer_size, e_c);

    std::ofstream smart_output_file;
    std::ofstream brute_output_file;

    smart_output_file.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/smart.dat");
    brute_output_file.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/bruteforce.dat");

    TS_ASSERT_EQUALS(smart_output_file.is_open(), true);
    TS_ASSERT_EQUALS(brute_output_file.is_open(), true);

    /*Set the output to be very precise */
    smart_output_file.precision(18);
    brute_output_file.precision(18);

    smart_output_file << "membrane_voltage" << "\n";
    brute_output_file << "membrane_voltage" << "\n";

    // p_model->SetIntegrationConstant(initial_v);
    // p_model2->SetIntegrationConstant(initial_v);

    bool finished1 = false;
    bool finished2 = false;

    p_model1->SetIntegrationConstant(-85.0133852794686);
    p_model2->SetIntegrationConstant(-85.0133852794686);


    /*Run the simulations*/
    for(unsigned int j = 0; j < paces; j++){
      std::cout << "pace " << j << "\n";
      if(!finished1){
        if(smart_simulation.RunPace()){
          std::cout << "Model " << model_name << " period " << 500 << " extrapolation method finished after " << j << " paces \n";
          finished1 = true;
        }
        // else
        //   smart_output_file << p_model->CalculateAnalyticVoltage() << "\n";
      }

      if(!finished2){
        if(simulation.RunPace()){
          std::cout << "Model " << model_name << " period " << 500 << " brute force method finished after " << j << " paces \n";
          finished2 = true;
        }
        // else
        //   brute_output_file << p_model2->CalculateAnalyticVoltage() << "\n";
      }
      if(finished1 && finished2)
        break;
    }

    std::vector<double> states = p_model1->GetStdVecStateVariables();
    std::vector<std::string> state_names = p_model1->rGetStateVariableNames();

    smart_output_file.close();
    brute_output_file.close();

    /*Check that the methods have converged to the same place*/

    std::cout << "MRMS between solutions is " << mrms(simulation.GetStateVariables(), smart_simulation.GetStateVariables()) << "\n";
    // std::cout << "Pace MRMS between solutions is " << CalculatePaceMrms(p_model1, simulation.GetStateVariables(), smart_simulation.GetStateVariables(), 500, p_model->UseCellMLDefaultStimulus()) << "\n";

    // std::cout << "Difference in APDs " << CalculateAPD(p_model, 500, 2, 90) - CalculateAPD(p_model2, 500, 2, 90) << "\n";

#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

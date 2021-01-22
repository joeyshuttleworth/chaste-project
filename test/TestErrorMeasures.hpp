#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include "SimulationTools.hpp"
#include <boost/filesystem.hpp>
#include <fstream>
#include "Simulation.hpp"
#include "OutputFileHandler.hpp"

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "luo_rudy_1994Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"

// The modified "analytic_voltage" models
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"

class TestErrorMeasures : public CxxTest::TestSuite
{
public:
  void TestAllErrorMeasures()
  {
#ifdef CHASTE_CVODE
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;

    const std::vector<std::string> models_with_redudant_voltage = { "ohara_rudy_2011_endo",
                                                                    "decker_2009",
                                                                    "ten_tusscher_model_2004",
                                                                    "shannon_wang_puglisi_weber_epi"
    };

    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017_analyticFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epi_analyticFromCellMLCvode(p_solver, p_stimulus)));

    std::string username = std::string(getenv("USER"));

    const unsigned int paces  = 100;

    for(unsigned int i=0; i < models.size(); i++){
      std::vector<double> periods = {500, 1000};
      for(auto period : periods){
        const std::string model_name = models[i]->GetSystemInformation()->GetSystemName();
        boost::filesystem::create_directory("/home/" + username + "/testoutput/");
        // Use the 1000ms period version unless we're at 1000ms
        const double starting_period = period==1000?500:1000;
        const std::string output_dirname = "/home/" + username + "/testoutput/TestErrorMeasure_" + model_name + "_" + std::to_string(int(starting_period)) + "ms/";
        const std::string input_path = "/home/" + username + "/testoutput/TestErrorMeasure_" + model_name + "_" + std::to_string(int(period)) + "ms/final_pace.dat";
        boost::filesystem::create_directory(output_dirname);
        auto p_model = models[i];
        std::cout << "Testing model: " + model_name + "\n";

        Simulation simulation(models[i], period, input_path, 1e-12, 1e-12);

        std::vector<std::vector<double>> state_variables;

        const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();

        std::string errors_file_path = output_dirname + "error_measures.dat";

        std::ofstream errors_file(errors_file_path);
        TS_ASSERT_EQUALS(errors_file.is_open(), true);

        errors_file.precision(18);

        errors_file << "APD 2-Norm MRMS Trace-2-Norm Trace-MRMS \n";

        std::vector<std::string> names = p_model->GetSystemInformation()->rGetStateVariableNames();

        std::vector<double> times;
        std::vector<std::vector<double>> current_pace, previous_pace;
        for(unsigned int j = 0; j < paces; j++){
          OdeSolution current_solution = simulation.GetPace(true);
          std::vector<std::vector<double>> current_pace = current_solution.rGetSolutions();
          previous_pace = current_pace;
          times = current_solution.rGetTimes();

          // Can't do error measures if this is the first pace!
          if(j==0)
            continue;
          const std::vector<double> current_states = current_pace.back();
          const std::vector<double> previous_states = previous_pace.back();

          if(j%1000==0){
            std::cout << "pace " << j << std::endl;
          }
          if(i % 10==0){
            unsigned int starting_index = 0;
            if(models_with_redundant_voltage.find(model_name)!=models_with_redudant_voltage.end()){
              start_index = 1;
            }
            errors_file << simulation.GetApd(90, false, starting_index) << " ";
            errors_file << TwoNorm(current_states, previous_states, starting_index) << " ";
            errors_file << mrms(current_states,  previous_states, starting_index) << " ";
            errors_file << TwoNormTrace(current_pace, previous_pace, starting_index) << " ";
            errors_file << mrmsTrace(current_pace, previous_pace, starting_index) << " ";
          }
          //Print state variables
          for(unsigned int k = 0; k < current_states.size(); k++){
            errors_file << current_states[k] << " ";
          }
          errors_file << "\n";
        }
        errors_file.close();
      }
    }
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

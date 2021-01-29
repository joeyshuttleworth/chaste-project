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

    boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;

    const std::vector<std::string> models_with_redundant_voltage = { "ohara_rudy_2011_endo",
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

    const unsigned int paces  = 5000;

    std::vector<double> periods = {1000, 500};
    std::vector<double> IKrBlocks = {0, 0.5};

    for(auto model : models){
      for(auto period : periods){
        for(double IKrBlock : IKrBlocks){
          const std::string model_name = model->GetSystemInformation()->GetSystemName();
          const double starting_period = period==1000?500:1000;
          const double starting_block  = period==0?0.5:0;

          std::stringstream dirname;
          dirname << "/" << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";

          std::stringstream input_dirname_ss;
          input_dirname_ss << model_name+"_" << std::to_string(int(starting_period)) << "ms_" << int(100*starting_block)<<"_percent_block/";
          const std::string input_dirname = (test_dir / boost::filesystem::path(input_dirname_ss.str())).string();
          std::cout << "Testing model: " << model_name << " with period " << period << "ms and IKrBlock " << IKrBlock << "\n";

          const std::string input_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("final_states.dat")).string();

          Simulation simulation(model, period, input_path, 1e-12, 1e-12);

          std::vector<std::vector<double>> state_variables;

          const std::vector<std::string> state_variable_names = model->rGetStateVariableNames();

          const std::string errors_file_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("error_measures.dat")).string();

          std::cout << "outputting to " << errors_file_path << "\n";

          std::ofstream errors_file(errors_file_path);
          TS_ASSERT_EQUALS(errors_file.is_open(), true);

          errors_file.precision(18);

          errors_file << "APD 2-Norm MRMS Trace-2-Norm Trace-MRMS \n";

          std::vector<std::string> names = model->GetSystemInformation()->rGetStateVariableNames();

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
            if(j % 10==0){
              unsigned int starting_index = 0;
              if(find(models_with_redundant_voltage.begin(), models_with_redundant_voltage.end(), model_name)!=models_with_redundant_voltage.end()){
                starting_index = 1;
              }
              errors_file << simulation.GetApd(90, false) << " ";
              errors_file << TwoNorm(current_states, previous_states, starting_index) << " ";
              errors_file << mrms(current_states,  previous_states, starting_index) << " ";
              errors_file << TwoNormTrace(current_pace, previous_pace, starting_index) << " ";
              errors_file << mrmsTrace(current_pace, previous_pace, starting_index) << " ";
            //Print state variables
            for(unsigned int k = 0; k < current_states.size(); k++){
              errors_file << current_states[k] << " ";
            }
            errors_file << "\n";
            }
          }
          errors_file.close();
        }
      }
    }
#else
      std::cout << "Cvode is not enabled.\n";
#endif
    }
  };

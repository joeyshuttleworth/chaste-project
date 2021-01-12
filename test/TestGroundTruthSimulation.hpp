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

/* These header files are generated from the cellml files provided at github.com/chaste/cellml */

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017_analyticCvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_epi_analyticCvode.hpp"

/* Run the models under different scenarios and output:
   - All variables over the final pace
   - The APD90
   - Terminal state variables
   to separate files.
 */

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  void TestRunSimulation()
  {
#ifdef CHASTE_CVODE
    const int paces = 100;
    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;

    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus)));

    std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);

    for(unsigned int i = 0; i < models.size(); i++){
      const auto p_model = models[i];
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      std::vector<double> periods = {500, 1000};
      std::vector<double> initial_states = models[i]->GetStdVecStateVariables();
      for(unsigned int j =0; j < periods.size(); j++){
        std::cout << "Testing " << model_name << " with period " << int(periods[j]) << "\n";
        const std::string dirname = "/home/" + username + "/testoutput/TestGroundTruth_" + model_name  + std::to_string(int(periods[j])) + "ms/";
        boost::filesystem::create_directory(dirname);

        // Initialise simulation with fine tolerances
        Simulation simulation(models[i], periods[j], "", 1e-12, 1e-12);
        simulation.SetStateVariables(initial_states);

        // Turn off convergence criteria
        simulation.SetThreshold(0);

        // Run the simulation for a large number of paces
        simulation.RunPaces(paces);

        // Output the final pace
        const std::string pace_filename = "final_pace.dat";

        simulation.WritePaceToFile(dirname, pace_filename);

        // Output the APD90 of the final pace
        const std::string apd_filename = "final_apd90.dat";
        std::ofstream apd_file(dirname + apd_filename);
        apd_file << simulation.GetApd(90) << "\n";
        apd_file.close();
      }
   }
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

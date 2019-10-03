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

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
      
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Celldecker_2009FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));

    
    std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);

    std::ofstream output_file;
    int benchmark = 0;
    
    for(unsigned int i = 1; i < 8; i++){
      double period = 1000;
      if(i<4)
	period = 500;
      boost::shared_ptr<AbstractCvodeCell> p_model = models[i];
      
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      std::cout << "Testing " << model_name << " with period " << period << "\n";
      boost::filesystem::create_directory("/tmp/"+username);
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/TestExtrapolation");

      const double duration = 2;
      
      std::string input_path;
      if(period==500){
	//      	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat"), 0);
	input_path = "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth1Hz/final_state_variables.dat";
      }
      
      else{
	//	TS_ASSERT_EQUALS(LoadStatesFromFile(p_model, "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth2Hz/final_state_variables.dat"), 0);
	input_path = "/home/joey/code/chaste-project-data/"+model_name+"/GroundTruth2Hz/final_state_variables.dat";
      }
      const unsigned int paces  = 10000;
      
      SmartSimulation smart_simulation(models[i], period, input_path); 
      smart_simulation.Initialise(buffer_size, e_c);
      Simulation simulation(models[i+8], period, input_path);

      if(period==500)
	output_file.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/extrapolations1Hz2Hz.dat");
      else
	output_file.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/extrapolations2Hz1Hz.dat");

      TS_ASSERT_EQUALS(output_file.is_open(), true);

      /*Set the output to be very precise */
      output_file.precision(18);
      
      /*Run the simulations*/
      for(unsigned int j = 0; j < paces; j++){
	// if(simulation.RunPace())
	//   std::cout << "Model " << model_name << " period " << period << " Brute force method finished after " << j << " paces \n";
	if(smart_simulation.RunPace())
	  std::cout << "Model " << model_name << " period " << period << " Extrapolation method finished after " << j << " paces \n";

	// output_file << simulation.GetMrms() << " " << smart_simulation.GetMrms() << " ";

	// for(unsigned int k = 0; k < smart_simulation.GetStateVariables().size(); k++){
	//   output_file << smart_simulation.GetStateVariables()[k] << " ";
	// }

	
	// for(unsigned int k = 0; k < smart_simulation.GetStateVariables().size(); k++){
	//   output_file << simulation.GetStateVariables()[k] << " ";
	// }
	// output_file << "\n";
	if(simulation.is_finished()){
	  benchmark += j;
	  break;
	}
	
      }
      
      /*Check that the methods have converged to the same place*/
      std::ofstream tmp1, tmp2;
      tmp1.precision(18);
      tmp2.precision(18);

      std::cout << "Benchmark: " << benchmark << "\n";
      
      tmp1 << "\n";
      tmp1.close();
      if(period==500){
	tmp1.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/1Hz2Hzfinal_apd90.dat");
      }
      else{
	tmp1.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/2Hz1Hzfinal_apd90.dat");
      }
      tmp1 << CalculateAPD(p_model, period, duration, 90) << " " << CalculateAPD(models[i+8], period, duration, 90) << "\n";
      tmp1.close();
      if(period == 500){
	tmp1.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/1Hz2HzExtrapolationTrace.dat");
	tmp2.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/1Hz2HzBruteForceTrace.dat");
      }
      else{
	tmp1.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/2Hz1HzExtrapolationTrace.dat");
	tmp2.open("/tmp/"+username+"/"+model_name+"/TestExtrapolation/2Hz1HzBruteForceTrace.dat");
      }
      std::vector<double> vec1 = smart_simulation.GetStateVariables();
      std::vector<double> vec2 = simulation.GetStateVariables();
      std::vector<std::vector<double>> pace1 = GetPace(vec1, models[i], period, duration);
      std::vector<std::vector<double>> pace2 = GetPace(vec1, models[i], period, duration);
      for(unsigned int j = 0; j < pace1.size(); j++){
	WriteStatesToFile(pace1[j], tmp1);
	WriteStatesToFile(pace2[j], tmp2);
      }
      tmp1.close();
      tmp2.close();
	  
      output_file.close();
    }
#else
      std::cout << "Cvode is not enabled.\n";
#endif
    }
  };

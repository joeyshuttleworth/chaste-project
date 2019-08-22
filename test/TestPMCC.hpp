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
#include <boost/circular_buffer.hpp>

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

  double CalculatePMCC(std::vector<double> values){
    const unsigned int N = values.size();
    double sum_xy = 0, sum_x, sum_y = 0, sum_x2, sum_y2=0;

    sum_x  = N*N/2;
    sum_x2 = N*(N+1)*(2*N+1)/6;
    
    for(unsigned int i = 0; i < N; i++){
      sum_y  += values[i];
      sum_y2 += values[i]*values[i];
      sum_xy += values[i]*i;
    }

    double pmcc = (N*sum_xy - sum_x*sum_y)/pow((N*(sum_x2 - (sum_x*sum_x)))*(N*sum_y2 - sum_y*sum_y), 0.5);
    return pmcc;
  }
  
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
    std::string username = std::string(getenv("USER"));
    boost::filesystem::create_directory("/tmp/"+username);
    
    for(unsigned int i = 0; i < 2; i++){
      double period = 1000;
      if(i<1)
	period = 500;

      boost::shared_ptr<AbstractCvodeCell> p_model = models[i];
      boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
      
      const std::string model_name = p_model->GetSystemInformation()->GetSystemName();
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name);
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/GroundTruth2Hz");
      boost::filesystem::create_directory("/tmp/"+username+"/"+model_name+"/GroundTruth1Hz");
      const double duration   = p_regular_stim->GetDuration();

      p_regular_stim->SetPeriod(period);
      p_regular_stim->SetStartTime(0);
      p_model->SetTolerances(1e-12, 1e-12);
      p_model->SetMaxSteps(1e5);
      
      const unsigned int paces  = 10000;
      OdeSolution current_solution;
      std::ofstream output_file;
      const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames();

      std::cout << "Testing " << model_name << " with period " << period << "\n";
      p_model->SetMaxTimestep(1000);
      if(period==500){
	output_file.open("/tmp/"+username+"/"+model_name+"/1Hz2Hzpmcc.dat");
      }
      else{
	output_file.open("/tmp/"+username+"/"+model_name+"/2Hz1Hzpmcc.dat");
      }
      TS_ASSERT_EQUALS(output_file.is_open(), true);
      
      /*Set the output to be as precise as possible */
      output_file.precision(18);

      boost::circular_buffer<std::vector<double>> values(buffer_size);            
    
      /*Run the simulation*/

      unsigned int number_of_state_variables  = p_model->GetNumberOfStateVariables();
      
      for(unsigned int i = 0; i < paces - 1; i++){
	p_model->SolveAndUpdateState(0, duration);
	p_model->SolveAndUpdateState(duration, period);
	std::vector<double> state_variables = p_model->GetStdVecStateVariables();

	/*If values is full, this will remove the oldest value so that values will only contain up to 50 points*/  
	values.push_back(state_variables);

	if(i>0 && i%50 == 0){
	  for(unsigned int i = 0; i < number_of_state_variables; i++){
	    output_file << CalculatePMCC(values[i]) << " ";
	  }
	  output_file << "\n";
	}
      }
      output_file.close();
    }
#else
    std::cout << "Cvode is not enabled.\n";
#endif
  }
};

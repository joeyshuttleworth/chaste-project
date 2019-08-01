#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include <fstream>

#include "beeler_reuter_model_1977Cvode.hpp"


double twoNorm(std::vector<double> *A, std::vector<double> *B){
  double norm = 0;
  for(unsigned int i=0; i < A->size(); i++){
    double a = (*A)[i];
    double b = (*B)[i];
    norm += pow(a - b, 2);
  }
  return sqrt(norm);
}

double mrms(std::vector<double> *A, std::vector<double> *B){
  double norm = 0;
  
  for(unsigned int i=0; i < A->size(); i++){
    double a = (*A)[i];
    double b = (*B)[i];
    norm += pow((a - b)/(1 + abs(a)), 2);   
  }
  return sqrt(norm/A->size());
}

double twoNormTrace(OdeSolution *A, OdeSolution *B){
  double norm = 0;
  for(unsigned int i = 0; i < A->GetNumberOfTimeSteps(); i++){
    for(unsigned int j = 0; j < B->rGetSolutions()[0].size(); j++){
      norm = pow(A->rGetSolutions()[i][j] - B->rGetSolutions()[i][j],2);
    }
  }
  return sqrt(norm);  
}

double mrmsTrace(OdeSolution *A, OdeSolution *B){
  double norm = 0;
  for(unsigned int i = 0; i < A->GetNumberOfTimeSteps(); i++){
    for(unsigned int j = 0; j < A->rGetSolutions()[0].size(); j++){
      double a = A->rGetSolutions()[i][j];
      double b = B->rGetSolutions()[i][j];
      norm += pow((a - b)/(1+abs(a)), 2);
    }
  }
  return sqrt(norm/(A->GetNumberOfTimeSteps() * A->rGetSolutions().size()));  
}

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
    void TestTusscherSimulation()
    {
#ifdef CHASTE_CVODE
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus));
	boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	std::vector<std::vector<double>> solutions;
	std::vector<OdeSolution> finalTraces;
	const double period = 1000;
	const std::vector<double> tolerances = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12};
        p_regular_stim->SetPeriod(period);

	double max_timestep = p_regular_stim->GetDuration();

        p_model->SetMaxTimestep(max_timestep);
	p_model->SetMaxSteps(1e5);
	unsigned int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
	
        double sampling_timestep = max_timestep;
	int steps = 10000;
	OdeSolution current_solution;
	std::ofstream tolerances_file;
	std::ofstream norms_file;
	std::string username = std::string(getenv("USER"));
	const std::vector<std::string> state_variable_names = p_model->rGetStateVariableNames(); 
	norms_file.open("/tmp/"+username+"/norms.dat");
	tolerances_file.open("/tmp/"+username+"/tolerances.dat");
	/*Set cout to be as precise as possible */
	tolerances_file.precision(17);
	/*Print variable names on the first line*/
	tolerances_file << "APD90 ";
	for(unsigned int i = 0; i < state_variable_names.size(); i++){
	  tolerances_file << state_variable_names[i] << " ";
	}
	tolerances_file << "\n";

	std::vector<double> state_variables;
 	for(unsigned int k = 0; k < tolerances.size(); k++){
	  std::cout << "Testing with tolerance " << tolerances[k] << "\n";
	  p_model->SetTolerances(tolerances[k], tolerances[k]); 
	  for(int i=0; i < steps; i++){
	    current_solution = p_model->Compute(0, period, sampling_timestep);
	    state_variables = current_solution.rGetSolutions()[current_solution.GetNumberOfTimeSteps()-1];
	    p_model->SetStateVariables(state_variables);
	  }
	  
	  current_solution.WriteToFile("Trace" + std::to_string(k), "1e-" + std::to_string(k+3) + "-tolerance", "ms");
	  finalTraces.push_back(current_solution);
	  solutions.push_back(state_variables);
	  std::vector<double> voltages = current_solution.GetVariableAtIndex(voltage_index);
	  CellProperties cell_props(voltages, current_solution.rGetTimes());
	  double apd = cell_props.GetLastActionPotentialDuration(90);
	  tolerances_file << tolerances[k] <<apd  << " ";
	  for(unsigned int j=0; j < state_variables.size(); j++){
	    tolerances_file << state_variables[j] << " ";
	  }
	  tolerances_file << "\n";
	}
	std::cout << "Tolerance \t 2Norm  \t MRMS \t 2Norm over final trace\n";
	for(unsigned int i=0; i < solutions.size()-1; i++){
	  norms_file << tolerances[i] << "\t" << twoNorm(&solutions[i], &(solutions[solutions.size()-1])) << "\t" << mrms(&solutions.back(), &solutions[i]) << "\t" << twoNormTrace(&finalTraces[i], &finalTraces.back()) << "\t" << mrmsTrace(&finalTraces[i], &finalTraces.back()) << "\n";
	}
	tolerances_file.close();
	norms_file.close();
#else
	std::cout << "Cvode is not enabled.\n";
#endif
    }
};

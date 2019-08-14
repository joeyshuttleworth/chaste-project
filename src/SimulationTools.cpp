#include "SimulationTools.hpp"

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"

void RunSimulation(boost::shared_ptr<AbstractCvodeCell> p_model, unsigned int paces, unsigned int period, double tolerances){
  boost::shared_ptr<RegularStimulus> p_stimulus;
  boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
  
  p_stimulus = p_model->UseCellMLDefaultStimulus();
  p_stimulus->SetStartTime(0);
  for(unsigned int i = 0; i < paces - 1; i++){
    p_model->Solve(0, period, 1000);
  }
  OdeSolution solution = p_model->Compute(0, period, 0.1);
  solution.WriteToFile(p_model->GetSystemInformation()->GetSystemName(), "final_trace", "ms");
  return;
}

int LoadStatesFromFile(boost::shared_ptr<AbstractCvodeCell> p_model, std::string file_path){
  std::ifstream file_in;
  file_in.open(file_path);
  if(!file_in.is_open()){
    std::cout << "Couldn't open file! " + file_path + " \n";
    return -1;
  }
  std::string line;
  std::vector<std::string> state_variables_str;
  std::vector<double> state_variables;
  std::getline(file_in, line);
  boost::split(state_variables_str, line, boost::is_any_of(" "));
  for(unsigned int i = 0; i < p_model->GetNumberOfStateVariables(); i++){
    state_variables.push_back(std::stod(state_variables_str[i]));
  }
  p_model->SetStateVariables(state_variables);
  return 0;
}

std::vector<double> GetNthVariable(std::vector<std::vector<double>> states, unsigned int index){
  std::vector<double> vec;
  vec.reserve(states.size());
  for(unsigned int i = 0; i < states.size(); i++){
    vec.push_back(states[i][index]);
  }
  return vec;
}

double TwoNorm(std::vector<double> A, std::vector<double> B){
  double norm = 0;
  for(unsigned int i=0; i < A.size(); i++){
    double a = A[i];
    double b = B[i];
    norm += pow(a - b, 2);
  }
  return sqrt(norm);
}

double mrms(std::vector<double> A, std::vector<double> B){
  double norm = 0;
  
  for(unsigned int i=0; i < A.size(); i++){
    double a = A[i];
    double b = B[i];
    norm += pow((a - b)/(1 + abs(a)), 2);   
  }
  return sqrt(norm/A.size());
}

double TwoNormTrace(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){
  double norm = 0;
  for(unsigned int i = 0; i < A.size(); i++){
    for(unsigned int j = 0; j < A[0].size(); j++){
      norm += pow(A[i][j] - B[i][j], 2);
    }
  }
  return sqrt(norm);  
}

double mrmsTrace(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){
  double norm = 0;
  for(unsigned int i = 0; i < A.size(); i++){
    for(unsigned int j = 0; j < A[0].size(); j++){
      double a = A[i][j];
      double b = B[i][j];
      norm += pow((a - b)/(1+abs(a)), 2);
    }
  }
  return sqrt(norm/(A.size() * A[0].size()));  
}


#include "SimulationTools.hpp"

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
  file_in.close();
  p_model->SetStateVariables(state_variables);
  return 0;
}

std::vector<double> GetNthVariable(std::vector<std::vector<double>> states, unsigned int index){
  std::vector<double> vec;
  vec.reserve(states.size());
  for(auto i = states.begin(); i!=states.end(); i++){
    vec.push_back((*i)[index]);
  }
  return vec;
}

std::vector<double> cGetNthVariable(boost::circular_buffer<std::vector<double>> states, unsigned int index){
  std::vector<double> vec;
  vec.reserve(states.size());
  for(auto i = states.begin(); i != states.end(); i++){
    vec.push_back((*i)[index]);
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
  assert(A.size()==B.size());
  for(unsigned int i=0; i < A.size(); i++){
    double a = A[i];
    double b = B[i];
    norm += pow((a - b)/(1 + abs(a)), 2);
  }
  const double return_val = sqrt(norm/A.size());
  return return_val;
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

double CalculateAPD(boost::shared_ptr<AbstractCvodeCell> p_model, double period, double duration, double percentage){
  double apd;

  double sampling_timestep = 0.1;
  const std::vector<double> initial_conditions = p_model->GetStdVecStateVariables();
  const double rel_tol = p_model->GetRelativeTolerance();
  const double abs_tol = p_model->GetAbsoluteTolerance();

  p_model->SetMaxSteps(1e5);
  p_model->SetTolerances(1e-12, 1e-12);

  OdeSolution solution = p_model->Compute(0, duration, sampling_timestep);
  std::vector<std::vector<double>> state_variables = solution.rGetSolutions();
  std::vector<double> times = solution.rGetTimes();

  solution = p_model->Compute(duration, period, sampling_timestep);

  state_variables.insert(state_variables.end(), ++solution.rGetSolutions().begin(), solution.rGetSolutions().end());
  times.insert(times.end(), ++solution.rGetTimes().begin(), solution.rGetTimes().end());
  int voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");

  const std::vector<double> voltages = GetNthVariable(state_variables, voltage_index);
  CellProperties cell_props = CellProperties(voltages, times);
  
  apd = cell_props.GetLastActionPotentialDuration(percentage);

  p_model->SetTolerances(rel_tol, abs_tol);
  p_model->SetStateVariables(initial_conditions);
  return apd;
}

std::vector<std::vector<double>> GetPace(std::vector<double> initial_conditions, boost::shared_ptr<AbstractCvodeCell> p_model, double period, double duration){ 
  double sampling_timestep = 0.1;
  const std::vector<double> original_states = p_model->GetStdVecStateVariables();
  const double rel_tol = p_model->GetRelativeTolerance();
  const double abs_tol = p_model->GetAbsoluteTolerance();  
  
  p_model->SetMaxSteps(1e5);
  p_model->SetTolerances(1e-12, 1e-12);
  p_model->SetStateVariables(initial_conditions);

  OdeSolution solution = p_model->Compute(0, duration, sampling_timestep);
  std::vector<std::vector<double>> state_variables = solution.rGetSolutions();
  solution = p_model->Compute(duration, period, sampling_timestep);
  
  state_variables.insert(state_variables.end(), ++solution.rGetSolutions().begin(), solution.rGetSolutions().end());
  

  p_model->SetTolerances(rel_tol, abs_tol);
  p_model->SetStateVariables(original_states);
  return state_variables;
}


double CalculatePace2Norm(boost::shared_ptr<AbstractCvodeCell> p_model, std::vector<double> first_states, std::vector<double> second_states, double period, double duration){
  std::vector<std::vector<double>> A = GetPace(first_states, p_model, period, duration);
  std::vector<std::vector<double>> B = GetPace(second_states, p_model, period, duration);
  return TwoNormTrace(A, B);
}


double CalculatePaceMrms(boost::shared_ptr<AbstractCvodeCell> p_model, std::vector<double> first_states, std::vector<double> second_states, double period, double duration){
  std::vector<std::vector<double>> A = GetPace(first_states, p_model, period, duration);
    std::vector<std::vector<double>> B = GetPace(second_states, p_model, period, duration);
    return mrmsTrace(A, B);
}



double CalculatePMCC(std::vector<double> x, std::vector<double> y){
  const unsigned int N = x.size();
  if(x.size() <= 2){
    return -NAN;
  }
  double sum_x = 0, sum_x2 = 0, sum_y = 0, sum_y2 = 0, sum_xy = 0;

  for(unsigned int i = 0; i < N; i++){
    sum_x  += x[i];
    sum_x2 += x[i]*x[i];
    sum_y  += y[i];
    sum_y2 += y[i]*y[i];
    sum_xy += x[i]*y[i];
  }

  double pmcc = (N*sum_xy - sum_x*sum_y)/sqrt((N*sum_x2 - sum_x*sum_x)*(N*sum_y2 - sum_y*sum_y));
  assert(abs(pmcc <= 1.001) || pmcc == NAN);

  return pmcc;
}

void WriteStatesToFile(std::vector<double> states, std::ofstream &f_out){
  for(auto i = states.begin(); i!=states.end(); ++i){
    f_out << *i << " ";
  }
  f_out << "\n";
  return;
}

#include "SimulationTools.hpp"
#include <boost/filesystem.hpp>
#include "Simulation.hpp"

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
    EXCEPTION("Couldn't open file " + file_path);
    return -1;
  }
  std::string line;
  std::vector<std::string> state_variables_str;
  std::vector<double> state_variables;

  // Get second line
  std::getline(file_in, line);
  std::getline(file_in, line);

  boost::split(state_variables_str, line, boost::is_any_of(" "), boost::token_compress_on);
  for(auto state_var : state_variables_str){
    if(state_var!="")
      state_variables.push_back(std::stod(state_var));
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

double TwoNorm(std::vector<double> A, std::vector<double> B, unsigned int starting_index){
  double norm = 0;
  for(unsigned int i=starting_index; i < A.size(); i++){
    double a = A[i];
    double b = B[i];
    norm += pow(a - b, 2);
  }
  return sqrt(norm);
}

double mrms(std::vector<double> A, std::vector<double> B, unsigned int starting_index){
  double norm = 0;
  assert(A.size()==B.size());
  for(unsigned int i=starting_index; i < A.size(); i++){
    double a = A[i];
    double b = B[i];
    norm += pow((a - b)/(1 + abs(a)), 2);
  }
  const double return_val = sqrt(norm/A.size());
  return return_val;
}

double TwoNormTrace(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B, unsigned int starting_index){
  double norm = 0;
  for(unsigned int i = 0; i < A.size(); i++){
    for(unsigned int j = starting_index; j < A[0].size(); j++){
      norm += pow(A[i][j] - B[i][j], 2);
    }
  }
  return sqrt(norm);
}

double mrmsTrace(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B, unsigned int starting_index){
  double norm = 0;
  for(unsigned int i = 0; i < A.size(); i++){
    for(unsigned int j = starting_index; j < A[0].size(); j++){
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

  p_model->SetMaxSteps(1e5);

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

  p_model->SetStateVariables(initial_conditions);
  return apd;
}

std::vector<std::vector<double>> GetPace(std::vector<double> initial_conditions, boost::shared_ptr<AbstractCvodeCell> p_model, double period, double duration){
  double sampling_timestep = 0.1;
  const std::vector<double> original_states = p_model->GetStdVecStateVariables();

  p_model->SetMaxSteps(1e5);
  p_model->SetStateVariables(initial_conditions);

  OdeSolution solution = p_model->Compute(0, duration, sampling_timestep);
  std::vector<std::vector<double>> state_variables = solution.rGetSolutions();
  solution = p_model->Compute(duration, period, sampling_timestep);

  state_variables.insert(state_variables.end(), ++solution.rGetSolutions().begin(), solution.rGetSolutions().end());

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

void compare_error_measures(boost::shared_ptr<AbstractCvodeCell> model, double period, double IKrBlock, double tolerance, std::string filename_suffix){
  const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));
  const unsigned int paces  = 2500;

  const double default_GKr = model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance");
  model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr*(1- IKrBlock));

  const std::string model_name = model->GetSystemInformation()->GetSystemName();
  const unsigned int starting_index = 0;

  std::cout << "For model " << model_name << " using starting_index " << starting_index << "\n";
  model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", (1-IKrBlock)*default_GKr);

  const double starting_period = period==1000?500:1000;
  const double starting_block  = period==0?0.5:0;

  std::stringstream dirname;
  dirname << "/" << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";

  std::stringstream input_dirname_ss;
  input_dirname_ss << model_name+"_" << std::to_string(int(starting_period)) << "ms_" << int(100*starting_block)<<"_percent_block/";
  const std::string input_dirname = (test_dir / boost::filesystem::path(input_dirname_ss.str())).string();
  std::cout << "Testing model: " << model_name << " with period " << period << "ms and IKrBlock " << IKrBlock << "\n";

  const std::string input_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("final_states.dat")).string();

  Simulation simulation(model, period, input_path, tolerance, tolerance);
  model->SetTolerances(tolerance, tolerance);
  simulation.SetTerminateOnConvergence(false);

  std::vector<std::vector<double>> state_variables;

  const std::vector<std::string> state_variable_names = model->rGetStateVariableNames();

  std::stringstream error_file_name;
  error_file_name << filename_suffix << "_" << tolerance << ".dat";
  const std::string errors_file_path = (test_dir / boost::filesystem::path(dirname.str()) / boost::filesystem::path(error_file_name.str())).string();

  std::cout << "outputting to " << errors_file_path << "\n";

  std::ofstream errors_file(errors_file_path);
  if(!errors_file.is_open()){
    EXCEPTION("Failed to open file " + errors_file_path);
  }

  errors_file.precision(18);

  errors_file << "APD 2-Norm MRMS Trace-2-Norm Trace-MRMS ";

  std::vector<std::string> names = model->GetSystemInformation()->rGetStateVariableNames();
  for(std::string name : names){
    errors_file << name << " ";
  }
  errors_file << "\n";

  std::vector<double> times;
  for(unsigned int j = 0; j < paces; j++){
    // std::cout << "pace = " << j << "\n";
    OdeSolution current_solution = simulation.GetPace(1, false);
    const std::vector<std::vector<double>> previous_pace = current_solution.rGetSolutions();
    simulation.RunPace();
    current_solution = simulation.GetPace(1, false);
    const std::vector<std::vector<double>> current_pace = current_solution.rGetSolutions();
    times = current_solution.rGetTimes();

    const std::vector<double> current_states = current_pace.back();
    const std::vector<double> previous_states = previous_pace.back();

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

    simulation.RunPaces(9);
    j+=9;
  }
  errors_file.close();

  model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", default_GKr);
}

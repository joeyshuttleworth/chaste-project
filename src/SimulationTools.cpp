#include "SimulationTools.hpp"
#include <boost/filesystem.hpp>
#include "Simulation.hpp"

#include "CommandLineArguments.hpp"

// These header files are generated from the cellml files provided at github.com/chaste/cellml
#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_epiCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"
#include "iyer_2004Cvode.hpp"
#include "ToRORd_dynCl_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"

// Analytic models
#include "decker_2009_analytic_voltageCvode.hpp"
#include "hund_rudy_2004_analytic_voltageCvode.hpp"
#include "iyer_2004_analytic_voltageCvode.hpp"
#include "ohara_rudy_2011_epi_analytic_voltageCvode.hpp"
#include "ohara_rudy_cipa_2017_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2006_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2004_epi_analytic_voltageCvode.hpp"
#include "ToRORd_dyn_chloride_epi_analytic_voltageCvode.hpp"

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

std::vector<double> LoadStatesFromFile(std::string file_path){
  std::ifstream file_in;
  file_in.open(file_path);
  if(!file_in.is_open()){
    std::cout << "Couldn't open file! " + file_path + " \n";
    EXCEPTION("Couldn't open file " + file_path);
    return {};
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
  return state_variables;
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
  return sqrt(norm/(A.size()-starting_index));
}

double mrms(std::vector<double> A, std::vector<double> B, unsigned int starting_index){
  double norm = 0;
  assert(A.size()==B.size());
  for(unsigned int i=starting_index; i < A.size(); i++){
    double a = A[i];
    double b = B[i];
    norm += pow((a - b)/(1 + abs(a)), 2);
  }
  const double return_val = sqrt(norm/(A.size()-starting_index));
  return return_val;
}

double TwoNormTrace(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B, unsigned int starting_index){
  double norm = 0;
  for(unsigned int i = 0; i < A.size(); i++){
    for(unsigned int j = starting_index; j < A[0].size(); j++){
      norm += pow(A[i][j] - B[i][j], 2);
    }
  }
  return sqrt(norm/(A.size() * (A[0].size()-starting_index)));
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
  return sqrt(norm/(A.size() * (A[0].size()-starting_index)));
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

void compare_error_measures(int paces, boost::shared_ptr<AbstractCvodeCell> model, double period, double IKrBlock, double tolerance, std::string filename_prefix){
  /* Run the models and output each of the error measures after every paces.
     Also output the APD90 resulting fromcomputing a pace using the current
     state variables as the initial conditions in the model and integrating
     forward one pace with fine tolerances. This is very slow and involves a lot
     of work.
   */

  if(paces == INT_UNSET)
    EXCEPTION("Invalid number of paces");

  const boost::filesystem::path test_dir(getenv("CHASTE_TEST_OUTPUT"));

  const std::string model_name = model->GetSystemInformation()->GetSystemName();
  const unsigned int starting_index = 0;

  std::cout << "For model " << model_name << " using starting_index " << starting_index << "\n";

  const double starting_period = period==1000?500:1000;
  const double starting_block  = IKrBlock==0?0.5:0;

  std::stringstream dirname;
  dirname << "/" << model_name << "_" << std::to_string(int(period)) << "ms_" << int(100*IKrBlock)<<"_percent_block/";

  std::stringstream input_dirname_ss;
  input_dirname_ss << model_name+"_" << std::to_string(int(starting_period)) << "ms_" << int(100*starting_block)<<"_percent_block/";
  const std::string input_dirname = (test_dir / boost::filesystem::path(input_dirname_ss.str())).string();
  std::cout << "Testing model: " << model_name << " with period " << period << "ms and IKrBlock " << IKrBlock << " and tolerances " << tolerance << "\n";

  const std::string input_path = (test_dir / boost::filesystem::path(input_dirname_ss.str()) / boost::filesystem::path("final_states.dat")).string();

  Simulation simulation(model, period, input_path, tolerance, tolerance);
  simulation.SetIKrBlock(IKrBlock);

  simulation.SetTerminateOnConvergence(false);

  std::vector<std::vector<double>> state_variables;

  const std::vector<std::string> state_variable_names = model->rGetStateVariableNames();

  std::stringstream output_file_name;
  output_file_name << filename_prefix << "_" << tolerance << ".dat";
  const std::string output_file_path = (test_dir / boost::filesystem::path(dirname.str()) / boost::filesystem::path(output_file_name.str())).string();

  std::cout << "outputting to " << output_file_path << "\n";

  std::ofstream output_file(output_file_path);
  if(!output_file.is_open()){
    EXCEPTION("Failed to open file " + output_file_path);
  }

  output_file.precision(18);

  output_file << "APD 2-Norm MRMS Trace-2-Norm Trace-MRMS ";

  std::vector<std::string> names = model->GetSystemInformation()->rGetStateVariableNames();
  for(std::string name : names){
    output_file << name << " ";
  }
  output_file << "\n";

  output_file << simulation.GetApd(90) << " 0 0 0 0 ";

  auto current_states = simulation.GetStateVariables();

  for(auto state : current_states)
    output_file << state << " ";

  output_file << "\n";

  std::vector<double> times;
  for(int j = 0; j < paces; j++){
    OdeSolution current_solution = simulation.GetPace(1, false);
    const std::vector<std::vector<double>> previous_pace = current_solution.rGetSolutions();
    bool failed = false;
    try{
      simulation.RunPace();
    }
    catch(const Exception &exc){
      failed = true;
    }

    current_solution = simulation.GetPace(1, false);
    const std::vector<std::vector<double>> current_pace = current_solution.rGetSolutions();
    times = current_solution.rGetTimes();

    const std::vector<double> current_states = current_pace.front();
    const std::vector<double> previous_states = previous_pace.front();

    // Expensive so only compute every 50 paces
    if(j%50==0)
      output_file << simulation.GetApd(90, false) << " ";
    else
      output_file << NAN << " ";

    output_file << TwoNorm(current_states, previous_states, starting_index) << " ";
    output_file << mrms(current_states,  previous_states, starting_index) << " ";
    //Computationally expensive so only plot every 50
    if(j%50==0){
      output_file << TwoNormTrace(current_pace, previous_pace, starting_index) << " ";
      output_file << mrmsTrace(current_pace, previous_pace, starting_index) << " ";
    }
    else{
      output_file << NAN << " ";
      output_file << NAN << " ";
    }

    //Print state variables
    for(unsigned int k = 0; k < current_states.size(); k++){
      output_file << current_states[k] << " ";
    }
    output_file << "\n";
    if(failed){
      std::cout << "Terminated early after " << j  << " paces";
      break;
    }
  }
  output_file.close();
}

std::vector<boost::shared_ptr<AbstractCvodeCell>> get_models(const std::string& type){
  /* Get models using specified command line argument. If no argument is given use all models */
  std::vector<boost::shared_ptr<AbstractCvodeCell>> models;
  boost::shared_ptr<RegularStimulus> p_stimulus;
  boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

  if(type=="all"){
  // Original models with no algebraic voltage version
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellbeeler_reuter_model_1977FromCellMLCvode(p_solver, p_stimulus)));
  }

  if(type=="original"||type=="all"){
    // Original models which have a corresponding algebraic voltage version
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2004_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_2011_epiFromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::shared_ptr<AbstractCvodeCell>(new Cellohara_rudy_cipa_v1_2017FromCellMLCvode(p_solver, p_stimulus)));
    models.push_back(boost::make_shared<CellToRORd_dynCl_epiFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Celliyer_2004FromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Cellhund_rudy_2004FromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Celldecker_2009FromCellMLCvode>(p_solver, p_stimulus));
  }
  if(type=="algebraic" || type=="all"){
    // Algebraic voltage models
    models.push_back(boost::make_shared<Cellten_tusscher_2004_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Cellten_tusscher_2006_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Cellohara_rudy_2011_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Cellohara_rudy_cipa_2017_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<CellToRORd_dyn_chloride_epi_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Celliyer_2004_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Cellhund_rudy_2004_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
    models.push_back(boost::make_shared<Celldecker_2009_analytic_voltageFromCellMLCvode>(p_solver, p_stimulus));
 }

  /* If "--models" option is provided, use only those models that are specified */
  if(CommandLineArguments::Instance()->OptionExists("--models")){
    std::vector<std::string> model_names;
    model_names = CommandLineArguments::Instance()->GetStringsCorrespondingToOption("--models");
    assert(model_names.size()>0);
    std::vector<boost::shared_ptr<AbstractCvodeCell>> new_models;
    for(auto name : model_names){
      /* Find all models with the name that has been provided */
      auto found_model = std::find_if(models.begin(), models.end(), [&](boost::shared_ptr<AbstractCvodeCell>m)->bool {return m->GetSystemInformation()->GetSystemName()==name;});
      if(found_model!=models.end())
        new_models.push_back(*found_model);
    }
    /* new_models contains all of the models which have been specified - use this instead of models */
    models = new_models;
  }

  return models;
}

std::vector<double> get_IKr_blocks(){
  const std::string option = "--IKrBlocks";
  /* Get IKrBlocks using specified command line argument. If no argument is given default to defaults */
  std::vector<double> IKrBlocks = {0, 0.25, 0.5};

  /* If "--IKrBlocks" option is provided, use only those IKrBlocks that are specified */
  if(CommandLineArguments::Instance()->OptionExists(option)){
    IKrBlocks = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption(option);
  }

  return IKrBlocks;
}

std::vector<double> get_periods(){
  const std::string option = "--periods";
  /* Get periods using specified command line argument. If no argument is given, default to defaults */
  std::vector<double> periods={1000, 500, 750, 1250};

    if(CommandLineArguments::Instance()->OptionExists(option)){
    periods = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption(option);
  }

  return periods;
}

double get_max_paces(){
  double paces = INT_UNSET;
  const std::string option = "--paces";

  /* Get maximum number of paces to run using specified command line argument. If no argument is given, return DOUBLE_UNSET */

  if(CommandLineArguments::Instance()->OptionExists(option)){
    paces = CommandLineArguments::Instance()->GetIntCorrespondingToOption(option);
  }
  return paces;
}

std::vector<boost::shared_ptr<AbstractCvodeCell>> get_analytic_models(){
  return(get_models("algebraic"));
}

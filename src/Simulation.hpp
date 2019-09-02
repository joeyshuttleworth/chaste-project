#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <sstream>
#include <iostream>

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp" 
#include "SimulationTools.hpp"

class Simulation
{
private:
protected:
  bool finished;
  boost::shared_ptr<AbstractCvodeCell> p_model;
  unsigned int number_of_state_variables;
  std::vector<double> state_variables;
  std::ofstream output_file;
  double period = 1000;
  double TolAbs = 1e-8, TolRel = 1e-8;
  double sampling_timestep = 1;
  double current_mrms = NAN;
  boost::shared_ptr<RegularStimulus> p_stimulus;
  const double threshold = 1.5e-7;
public:
  Simulation(){
    return;
  }
  Simulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-8, double _tol_rel=1e-8) : p_model(_p_model), period(_period), TolAbs(_tol_abs), TolRel(_tol_rel){ 
    finished = false;
    p_stimulus = p_model->UseCellMLDefaultStimulus();
    p_stimulus->SetStartTime(0);
    p_stimulus->SetPeriod(period);
    p_model->SetMaxSteps(1e5);
    p_model->SetMaxTimestep(1000);
    p_model->SetTolerances(TolAbs, TolRel);
    number_of_state_variables = p_model->GetSystemInformation()->rGetStateVariableNames().size();
    state_variables = p_model->GetStdVecStateVariables();
    if(input_path.length()>=1){
      LoadStatesFromFile(p_model, input_path);
    }  
  }

  bool RunPaces(unsigned int paces){
    if(finished)
      return false;
    for(unsigned int i = 0; i < paces; i++){
      /*Solve in two parts*/
      std::vector<double> tmp_state_variables = p_model->GetStdVecStateVariables();
      p_model->SolveAndUpdateState(0, p_stimulus->GetDuration());
      p_model->SolveAndUpdateState(p_stimulus->GetDuration(), period);
      std::vector<double> new_state_variables = p_model->GetStdVecStateVariables();
      current_mrms = mrms(tmp_state_variables, new_state_variables);
      p_model->SetStateVariables(new_state_variables);
      if(current_mrms < threshold){
	finished = true;
	return true;
      }
    }
    return false;
  }

  void RunPace(){
    RunPaces(1);
    return;
  }

  void WriteToFile(std::string){
    return;
  }

  double GetMrms(){
     if(finished)
       return NAN;
     else
      return current_mrms;
  }
  
  bool is_finished(){
    return finished;
  }
  
  std::vector<double> GetStateVariables(){
    return p_model->GetStdVecStateVariables();
  }
};

class SmartSimulation : public Simulation{

private:
  const unsigned int buffer_size = 200;
  boost::circular_buffer<std::vector<double>>  states_buffer;
  boost::circular_buffer<double> mrms_buffer;
  unsigned int jumps = 0;
  const unsigned int max_jumps = 1;

  double ExtrapolateState(unsigned int state_index){
    
    /* Calculate the log absolute differences of the state and store these in y_vals. Store the corresponding x values in x_vals*/
    std::vector<double> y_vals, x_vals;
    std::vector<double> state = cGetNthVariable(states_buffer, state_index);
     // if(CalculatePMCC(state)<-0.9)  //return the unchanged value if there is no negative correlation (PMCC > -0.9 or PMCC = NAN) 
     //   return state.back();
    for(unsigned int i = 0; i < buffer_size; i++){
      double tmp = abs(state[i] - state[i+1]);
      if(tmp != 0){
	y_vals.push_back(log(tmp));
	x_vals.push_back(i);
      }
    }

    /*Compute the required sums*/

    double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_xy = 0;
    const unsigned int N = x_vals.size();

    if(N<=2){
      std::cout << "N = "<< N << "!\n";
      return state.back();
    }
    
    for(unsigned int i = 0; i < N; i++){
      sum_x += x_vals[i];
      sum_y += y_vals[i];
      sum_x2+= x_vals[i]*x_vals[i];
      sum_xy+= x_vals[i]*y_vals[i];
    }

    double tau  = -(N*sum_xy - sum_x*sum_y) / (N*sum_x2 - sum_x*sum_x);

    double intercept = (sum_y*sum_x2 - sum_x*sum_xy) / (N*sum_x2 - sum_x*sum_x);

    double alpha = exp(intercept-buffer_size*tau+tau);
    
    /*Is V(t) increasing or decreasing?*/    
    
    if(state.back() - state.front() < 0)
      alpha = - alpha;
    
    std::cout << "(Gradient, Intercept) = (" << -tau << "," << intercept << ")\n";

    double change_in_variable = alpha / (exp(tau)-1);

    double new_value = state.back() + change_in_variable;   //Being careful 
    std::cout << "Change in " << p_model->GetSystemInformation()->rGetStateVariableNames()[state_index] << " is: " << change_in_variable << "\n" << "New value is " << new_value << "\n";
    if(std::isfinite(new_value)) 
      return new_value;
    else
      return state.back();
  }
  
  void ExtrapolateStates(){
    if(jumps>=max_jumps)
      return;
    if(!states_buffer.full())
      return;
    double mrms_pmcc = CalculatePMCC(mrms_buffer);
    std::vector<double> new_state_variables;
    new_state_variables.reserve(number_of_state_variables);
    if(mrms_pmcc < -0.95){
      for(unsigned int i = 0; i < number_of_state_variables; i++){
	new_state_variables.push_back(ExtrapolateState(i));
      }
    std::cout << "Jumped to new variables\n";

    std::ofstream f_out;
    if(period == 500)
      f_out.open("/tmp/joey/"+p_model->GetSystemInformation()->GetSystemName() + "/TestBenchmark/1Hz2HzJump.dat");
    else
      f_out.open("/tmp/joey/"+p_model->GetSystemInformation()->GetSystemName() + "/TestBenchmark/2Hz1HzJump.dat");
    WriteStatesToFile(new_state_variables, f_out);
    f_out.close();
    
    jumps++;
    
    /*Write new variables to file&*/
    std::ofstream output;
    output.precision(18);
    output.open("/tmp/jump.dat");
    for(unsigned int i = 0; i < state_variables.size(); i++){
      output << state_variables[i] << " ";
    }
    output << "\n";
    state_variables = new_state_variables;
    p_model->SetStateVariables(new_state_variables);
    for(unsigned int i = 0; i < state_variables.size(); i++){
      output << state_variables[i] << " ";
    }
    output.close();
    }
    return;
  }
  
public:

  using Simulation::Simulation;
  
  bool RunPaces(unsigned int paces){
    if(is_finished())
      return false;
    for(unsigned int i = 0; i < paces; i++){
      /*Solve in two parts*/
      p_model->SolveAndUpdateState(0, p_stimulus->GetDuration());
      p_model->SolveAndUpdateState(p_stimulus->GetDuration(), period);
      std::vector<double> new_state_variables = p_model->GetStdVecStateVariables();
      states_buffer.push_back(new_state_variables);
      current_mrms = mrms(new_state_variables, state_variables);
      mrms_buffer.push_back(current_mrms);
      state_variables = new_state_variables;
      if(current_mrms < threshold){
  	finished = true;
  	return true;
      }
      ExtrapolateStates();
    }
    return false;
  }
  
  void Initialise(){
    states_buffer.set_capacity(buffer_size);
    mrms_buffer.set_capacity(buffer_size);
  }

};


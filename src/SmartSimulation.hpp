#ifndef SMARTSIMULATION_HPP
#define SMARTSIMULATION_HPP
#include "Simulation.hpp"


class SmartSimulation : public Simulation{

private:
  unsigned int buffer_size = 200;
  double  extrapolation_coefficient;
  boost::circular_buffer<std::vector<double>>  states_buffer;
  boost::circular_buffer<double> mrms_buffer;
  unsigned int jumps = 0;
  unsigned int max_jumps = 100;
  std::vector<double> safe_state_variables;
  unsigned int pace = 0;
  std::ofstream errors;
  bool ExtrapolateState(unsigned int state_index){
    /* Calculate the log absolute differences of the state and store these in y_vals. Store the corresponding x values in x_vals*/
    std::vector<double> state = cGetNthVariable(states_buffer, state_index);

    for(unsigned int i = 0; i < buffer_size - 1; i++){
      double tmp = abs(state[i] - state[i+1]);
      if(tmp != 0){
        y_vals.push_back(log(tmp));
        x_vals.push_back(i);
      }
    }
    const double pmcc = CalculatePMCC(x_vals, y_vals);

    /*Compute the required sums*/

    double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_xy = 0;
    const unsigned int N = x_vals.size();

    if(N<=2){
      return false;
    }

    for(unsigned int i = 0; i < N; i++){
      sum_x += x_vals[i];
      sum_y += y_vals[i];
      sum_x2+= x_vals[i]*x_vals[i];
      sum_xy+= x_vals[i]*y_vals[i];
    }

    const double beta  = (N*sum_xy - sum_x*sum_y) / (N*sum_x2 - sum_x*sum_x);

    const double alpha = (sum_y*sum_x2 - sum_x*sum_xy) / (N*sum_x2 - sum_x*sum_x);


    output_file << state_index << " " << beta << " " << alpha << "\n";

    if(pmcc>-0.75){  //return the unchanged value if there is no negative correlation (PMCC > -0.9 or PMCC = NAN) 
      //      std::cout <<  p_model->GetSystemInformation()->rGetStateVariableNames()[state_index]<< ": PMCC was " << pmcc << " Ignoring. \n";
      return false;
    }

    if(exp(alpha) < 100*TolRel*state.back()){
      //The difference will be about as small as solver tolerances so there is no point going any further
      //std::cout << p_model->GetSystemInformation()->rGetStateVariableNames()[state_index] << ": Alpha too small!\n";
      return false;
    }

    if(beta > 0){
      //The difference is increasing
      //std::cout << p_model->GetSystemInformation()->rGetStateVariableNames()[state_index] << ": Beta is positive\n";
      return false;
    }


    const double tau = -1/beta;
    double change_in_variable =  abs(extrapolation_coefficient * exp(alpha - buffer_size/tau + 1/tau) / (exp(1/tau) - 1));

    /*Is V(t) increasing or decreasing?*/

    if(state.back() - state.front() < 0)
      change_in_variable = - change_in_variable;

    // std::cout << "(Gradient, Intercept) = (" << beta << "," << alpha << ")\n";

    double new_value = state.back() + change_in_variable;

    // std::cout << "Change in " << p_model->GetSystemInformation()->rGetStateVariableNames()[state_index] << " is: " << change_in_variable << "\n" << "New value is " << new_value << "\n";
    //std::cout << "Old value: " << state.back() << "\n";
    if(std::isfinite(new_value)){
      state_variables[state_index] = new_value;
      return true;
    }
    else{
      return false;
    }
  }

  bool ExtrapolateStates(){
    if(jumps>=max_jumps)
      return false;
    if(!mrms_buffer.full())
      return false;
    double mrms_pmcc = CalculatePMCC(mrms_buffer);
    std::vector<double> new_state_variables;
    bool extrapolated = false;
    std::ofstream f_out;
    std::string model_name = p_model->GetSystemInformation()->GetSystemName();
    if(mrms_pmcc < -0.975){
      safe_state_variables = state_variables;
      if(period == 500)
        f_out.open("/tmp/joey/"+ model_name + "/1Hz2HzJump.dat");
      else
        f_out.open("/tmp/joey/"+p_model->GetSystemInformation()->GetSystemName() + "/2Hz1HzJump.dat");
      std::cout << "Extrapolating - start of buffer is " << pace - buffer_size + 1<< "\n";
      pace++;
      output_file.open("/tmp/joey/"+p_model->GetSystemInformation()->GetSystemName()+"/"+ std::to_string(int(period)) + "JumpParameters.dat");
      output_file << pace << " " << buffer_size << " " << extrapolation_coefficient << "\n";

      WriteStatesToFile(state_variables, f_out);

      std::ofstream f_buffer;
      f_buffer.open("/tmp/joey/"+model_name+"/Buffer.dat");

      for(unsigned int i = 0; i < number_of_state_variables; i++){
        if(ExtrapolateState(i))
          extrapolated = true;
        WriteStatesToFile(cGetNthVariable(states_buffer, i), f_buffer);
      }

      output_file.close();

      WriteStatesToFile(state_variables, f_out);

      f_buffer.close();

      f_out.close();

      /*Write new variables to file*/
      std::ofstream output;
      for(unsigned int i = 0; i < state_variables.size(); i++){
        output << state_variables[i] << " ";
      }
      for(unsigned int i = 0; i < state_variables.size(); i++){
        output << state_variables[i] << " ";
      }
      output.close();
      if(extrapolated){
        mrms_buffer.clear();
        states_buffer.clear();

        //	std::cout << "Jumped to new variables\n";
        jumps++;
      }
      p_model->SetStateVariables(state_variables);
      return extrapolated;
    }
    else
      return false;
  }

public:

  using Simulation::Simulation;

  bool RunPace(){
    bool extrapolated = false;
    if(is_finished())
      return false;
    extrapolated = ExtrapolateStates();
    if(!extrapolated){
      /*Solve in two parts*/
      try{
        p_model->SolveAndUpdateState(0, p_stimulus->GetDuration());
        p_model->SolveAndUpdateState(p_stimulus->GetDuration(), period);
        pace++;
      }
      catch(Exception &e){
        std::cout << "RunPaces failed - returning to old state_variables\n";
        state_variables = safe_state_variables;
        p_model->SetStateVariables(state_variables);
        std::ofstream errors;
        errors.open("/tmp/joey/ExtrapolationErrors.dat", std::fstream::app);
        errors << p_model->GetSystemInformation()->GetSystemName() << " " << period << " " << buffer_size << " " << extrapolation_coefficient << "\n \n \n";
        errors << e.GetMessage();
        errors << "\n\n\n\n";
        errors.close();
        mrms_buffer.clear();
        states_buffer.clear();
        max_jumps = 0;
        return false;
      }
      std::vector<double> new_state_variables = GetStateVariables();
      states_buffer.push_back(new_state_variables);
      current_mrms = mrms(new_state_variables, state_variables);
      mrms_buffer.push_back(current_mrms);
      state_variables = new_state_variables;
      if(current_mrms < threshold){
        finished = true;
        return true;
      }
    }
    else{
      states_buffer.clear();
      mrms_buffer.clear();
      current_mrms = 0;
    }
    return false;
  }

  void Initialise(unsigned int _buffer_size, double _extrapolation_constant){
    buffer_size = _buffer_size;
    states_buffer.set_capacity(buffer_size);
    mrms_buffer.set_capacity(buffer_size);
    states_buffer.push_back(GetStateVariables());
    extrapolation_coefficient = _extrapolation_constant;
    safe_state_variables = state_variables;
  }
};

#endif

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include "SimulationTools.hpp"

class Simulation
{
private:
protected:
  bool mFinished;
  boost::shared_ptr<AbstractCvodeCell> mpModel;
  unsigned int mNumberOfStateVariables;
  std::vector<double> mStateVariables;
  std::ofstream mOutputFile;
  double mPeriod = 1000;
  double mTolAbs;
  double mTolRel;
  double mSamplingTimestep = 1;
  double mCurrentMrms = NAN;
  const double mThreshold = 1.8e-07;
  boost::shared_ptr<RegularStimulus> mpStimulus;
public:
  Simulation(){
    return;
  }
  Simulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-7, double _tol_rel=1e-7) : mpModel(_p_model), mPeriod(_period), mTolAbs(_tol_abs), mTolRel(_tol_rel){
    mFinished = false;
    mpStimulus = mpModel->UseCellMLDefaultStimulus();
    mpStimulus->SetStartTime(0);
    //There's no need to be near the second stimulus because Solve is called for
    //each pace
    mpStimulus->SetPeriod(2*mPeriod);
    mpModel->SetMaxSteps(1e5);
    mpModel->SetMaxTimestep(1000);
    mpModel->SetTolerances(mTolAbs, mTolRel);
    mpModel->SetMinimalReset(false); //Not sure if this is needed
    mNumberOfStateVariables = mpModel->GetSystemInformation()->rGetStateVariableNames().size();
    mStateVariables = mpModel->GetStdVecStateVariables();
    if(input_path.length()>=1){
      LoadStatesFromFile(mpModel, input_path);
    }
  }

  bool RunPace(){
    if(mFinished)
      return false;
    /*Solve in two parts*/
    std::vector<double> tmp_state_variables = mpModel->GetStdVecStateVariables();
    mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
    mpStimulus->SetPeriod(mPeriod*2);
    mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
    mpStimulus->SetPeriod(mPeriod);
    std::vector<double> new_state_variables = mpModel->GetStdVecStateVariables();
    mCurrentMrms = mrms(tmp_state_variables, new_state_variables);
    mpModel->SetStateVariables(new_state_variables);
    // mpModel->SetVoltage(mpModel->CalculateAnalyticVoltage());
    if(mCurrentMrms < mThreshold){
      mFinished = true;
      return true;
    }
    return false;
  }

  /**Output a pace to file*/
  void WriteToFile(std::string filename){
    OdeSolution solution = mpModel->Compute(0, mpStimulus->GetDuration(), 1);
    solution.WriteToFile(filename, filename, "ms");
    solution = mpModel->Compute(mpStimulus->GetDuration(), mPeriod, 1);
    solution.WriteToFile(filename, filename, "ms");
    return;
  }

  double GetMrms(){
    if(mFinished)
      return NAN;
    else
      return mCurrentMrms;
  }
  bool is_mFinished(){
    return mFinished;
  }
  std::vector<double> GetStateVariables(){
    return mpModel->GetStdVecStateVariables();
  }
};

class SmartSimulation : public Simulation{

private:
  unsigned int buffer_size = 200;
  double  extrapolation_coefficient = 1;
  boost::circular_buffer<std::vector<double>>  states_buffer;
  boost::circular_buffer<double> mrms_buffer;
  unsigned int jumps = 0;
  unsigned int max_jumps = 100;
  std::vector<double> safe_state_variables;
  unsigned int pace = 0;
  std::ofstream errors;

  bool ExtrapolateState(unsigned int state_index){
    /* Calculate the log absolute differences of the state and store these in y_vals. Store the corresponding x values in x_vals*/
    std::vector<double> y_vals;
    std::vector<double> x_vals;

    y_vals.reserve(buffer_size);
    x_vals.reserve(buffer_size);

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


    mOutputFile << state_index << " " << beta << " " << alpha << "\n";

    if(pmcc>-0.9){  //Don't change the value if there is no negative correlation (PMCC > -0.9 or PMCC = NAN)
      //      std::cout <<  p_model->GetSystemInformation()->rGetStateVariableNames()[state_index]<< ": PMCC was " << pmcc << " Ignoring. \n";
      return false;
    }

    if(exp(alpha) < 100*mTolRel*state.back()){
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
      mStateVariables[state_index] = new_value;
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
    std::string model_name = mpModel->GetSystemInformation()->GetSystemName();
    if(mrms_pmcc < -0.99){
      safe_state_variables = mStateVariables;
      if(mPeriod == 500)
        f_out.open("/tmp/joey/"+ model_name + "/1Hz2HzJump.dat");
      else
        f_out.open("/tmp/joey/"+mpModel->GetSystemInformation()->GetSystemName() + "/2Hz1HzJump.dat");
      std::cout << "start of buffer " << pace - buffer_size + 1<< "\n";
      pace++;
      mOutputFile.open("/tmp/joey/"+mpModel->GetSystemInformation()->GetSystemName()+"/"+ std::to_string(int(mPeriod)) + "JumpParameters.dat");
      mOutputFile << pace << " " << buffer_size << " " << extrapolation_coefficient << "\n";

      WriteStatesToFile(mStateVariables, f_out);

      std::ofstream f_buffer;
      f_buffer.open("/tmp/joey/"+model_name+"/Buffer.dat");

      for(unsigned int i = 0; i < mNumberOfStateVariables; i++){
        if(ExtrapolateState(i))
          extrapolated = true;
        WriteStatesToFile(cGetNthVariable(states_buffer, i), f_buffer);
      }

      mOutputFile.close();

      WriteStatesToFile(mStateVariables, f_out);

      f_buffer.close();

      f_out.close();

      /*Write new variables to file*/
      std::ofstream output;
      for(unsigned int i = 0; i < mStateVariables.size(); i++){
        output << mStateVariables[i] << " ";
      }
      for(unsigned int i = 0; i < mStateVariables.size(); i++){
        output << mStateVariables[i] << " ";
      }
      output.close();
      if(extrapolated){
        mrms_buffer.clear();
        states_buffer.clear();

        std::cout << "Jumped to new variables\n";
        jumps++;
      }
      mpModel->SetStateVariables(mStateVariables);
      return extrapolated;
    }
    else
      return false;
  }

public:
  using Simulation::Simulation;

  bool RunPace(){
    bool extrapolated = false;
    extrapolated = ExtrapolateStates();
    if(!extrapolated){
      /*Solve in two parts*/
      try{
        mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
        mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
        pace++;
      }
      catch(Exception &e){
        if(jumps==0){
          std::cout << "Failed to run model\n";
          return true;
        }

        std::cout << "RunPaces failed - returning to old mStateVariables\n";
        mStateVariables = safe_state_variables;
        mpModel->SetStateVariables(mStateVariables);
        std::ofstream errors;
        errors.open("/tmp/joey/ExtrapolationErrors.dat", std::fstream::app);
        errors << mpModel->GetSystemInformation()->GetSystemName() << " " << mPeriod << " " << buffer_size << " " << extrapolation_coefficient << "\n \n \n";
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
      mCurrentMrms = mrms(new_state_variables, mStateVariables);
      mrms_buffer.push_back(mCurrentMrms);
      mStateVariables = new_state_variables;
      if(mCurrentMrms < mThreshold){
        mFinished = true;
        return true;
      }
    }
    else{
      states_buffer.clear();
      mrms_buffer.clear();
      mCurrentMrms = 0;
    }
    return false;
  }

  void Initialise(unsigned int _buffer_size, double _extrapolation_constant){
    buffer_size = _buffer_size;
    states_buffer.set_capacity(buffer_size);
    mrms_buffer.set_capacity(buffer_size);
    states_buffer.push_back(GetStateVariables());
    extrapolation_coefficient = _extrapolation_constant;
    safe_state_variables = mStateVariables;
  }
};

#endif

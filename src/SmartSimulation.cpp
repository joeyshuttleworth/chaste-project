#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "SmartSimulation.hpp"

bool SmartSimulation::ExtrapolateState(unsigned int state_index, bool& stop_extrapolation){
  /* Calculate the log absolute differences of the state and store these in y_vals. Store the corresponding x values in x_vals*/
  std::vector<double> y_vals;
  std::vector<double> x_vals;
  y_vals.reserve(mBufferSize);
  x_vals.reserve(mBufferSize);

  std::vector<double> state = cGetNthVariable(mStatesBuffer, state_index);
  for(unsigned int i = 0; i < mBufferSize - 1; i++){
    double tmp = abs(state[i] - state[i+1]);
    if(tmp != 0){
      y_vals.push_back(log(tmp));
      x_vals.push_back(i);
    }
  }
  // const double pmcc = CalculatePMCC(x_vals, y_vals);

  // if(pmcc>-0.8 && std::isfinite(pmcc)){  //keep the unchanged value if there is no negative correlation (PMCC > -0.9 or PMCC = NAN)
  //   std::cout <<  mpModel->GetSystemInformation()->rGetStateVariableNames()[state_index]<< ": PMCC was " << pmcc << " Ignoring. \n";
  //   return false;
  // }


  /*Compute the required sums*/

  double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_xy = 0, sum_y2 = 0;
  const unsigned int N = x_vals.size();

  if(N<=2){
    return false;
  }

  for(unsigned int i = 0; i < N; i++){
    sum_x += x_vals[i];
    sum_y += y_vals[i];
    sum_x2+= x_vals[i]*x_vals[i];
    sum_xy+= x_vals[i]*y_vals[i];
    sum_y2+= y_vals[i]*y_vals[i];
  }

  const double r2 = pow(sum_xy/N - sum_x*sum_y/(N*N), 2)/((sum_x2/N - sum_x*sum_x/(N*N))*(sum_y2/N - sum_y*sum_y/(N*N)));

  if(r2<0.5 && std::isfinite(r2)){
    // std::cout <<  mpModel->GetSystemInformation()->rGetStateVariableNames()[state_index]<< ": r^2 was " << r2 << " Ignoring. \n";
    return false;
  }

  const double beta  = (N*sum_xy - sum_x*sum_y) / (N*sum_x2 - sum_x*sum_x);

  const double alpha = (sum_y*sum_x2 - sum_x*sum_xy) / (N*sum_x2 - sum_x*sum_x);

  if(beta > 0){
    //The difference is increasing or
    // std::cout << mpModel->GetSystemInformation()->rGetStateVariableNames()[state_index] << ": Beta is positive\n";
    return false;
  }


  const double tau = -1/beta;

  double V_diff = exp(alpha - mBufferSize/tau + 1/tau) / (exp(1/tau) - 1);

  /*Is V(t) increasing or decreasing?*/

  if(state.back() - state.front() < 0)
    V_diff = -V_diff;

  const double change_in_variable =  mExtrapolationConstant * V_diff;

  const double new_value = state.back() + change_in_variable;

  // Check that actual residue at the start of the buffer isn't too big
  const double V_total_difference = exp(alpha)/(1 - exp(-1/tau));
  const double predicted_v0 = state.back() + V_diff*(1 - exp(mBufferSize/tau));
  const double check_val = std::abs((state.front() - predicted_v0)/ V_total_difference);

  if(check_val > 0.5){
    // std::cout << "Not extrapolating - first residual too high \t" << check_val << "\n";
    return false;
  }

  // Check timescale isn't too big
  if(tau > mBufferSize * 50){
    // std::cout << "timescale too long, ignoring: tau = \t" << tau << "\n";
    return false;
  }
  // std::cout << "Change in " << p_model->GetSystemInformation()->rGetStateVariableNames()[state_index] << " is: " << change_in_variable << "\n" << "New value is " << new_value << "\n";

  //std::cout << "Old value: " << state.back() << "\n"
if(std::isfinite(new_value)){
    mStateVariables[state_index] = new_value;
    const std::string var_name = mpModel->GetSystemInformation()->rGetStateVariableNames()[state_index];
    mOutputFile << state_index << " " << var_name << " " << beta << " " << alpha <<  " " << state.back()+change_in_variable << "\n";
    return true;
  }
  else{
    return false;
  }
}


bool SmartSimulation::RunPace(){
  bool extrapolated = false;

  if(mHasTerminated)
    return true;

  extrapolated = ExtrapolateStates();
  if(!extrapolated){
    try{
      mHasTerminated = Simulation::RunPace();
    }
    catch(Exception& e){
      std::cout << "Failed to integrate pace. Returning to last known safe state.\n";
      if(mSafeStateVariables.size()==0){
        EXCEPTION("Solver crashed and there's no safe state to return to.");
      }
      else{
        mStateVariables = mSafeStateVariables;
        mpModel->SetStateVariables(mStateVariables);
        mSafeStateVariables = {};
        ClearBuffers();
      }
    }
    /* If the extrapolation method has been and a good number of paces have
       passed. We should check the pace-to-pace MRMS error is lower than it was
       before the jump. If not, we may have jumped further away from the
       solution and should reset to the state before the extrapolation. This
       seems to be quite rare so print a warning and stop any further extrapolations */
    if(mLastExtrapolationPace + mBufferSize == mPace  && mPreviousMinimalMRMSs.size()>0){
      if(mPreviousMinimalMRMSs.back() > mMRMSBeforeExtrapolation){
          std::cout << "Warning: the extrapolation appears to have gone wrong. The pace-to-pace MRMS error is greater than it was before the extrapolation. Will not perform further extrapolations\n";
          mStateVariables = mSafeStateVariables;
          mpModel->SetStateVariables(mStateVariables);
          mSafeStateVariables = {};

          //  Ensures that no further extrapolation will take place
          mMaxJumps=0;
        }
    }
  }
  mStatesBuffer.push_back(mStateVariables);
  return mHasTerminated;
}

bool SmartSimulation::ExtrapolateStates(){
    if(mJumps>=mMaxJumps)
      return false;
    if(!mStatesBuffer.full())
      return false;
    // double mrms_pmcc = CalculatePMCC(mMRMSBuffer);
    bool extrapolated = false;
    std::string model_name = mpModel->GetSystemInformation()->GetSystemName();
    const std::string dir_name = mOutputDir;
    boost::filesystem::create_directory(dir_name);
    if(true){// if(mrms_pmcc < -0.90){
      mSafeStateVariables = mStateVariables;
      std::cout << "Extrapolating - start of buffer is " << mPace - mBufferSize + 1<< "\n";

      mOutputFile.open(dir_name + "/" + std::to_string(int(mPeriod)) + "JumpParameters" + std::to_string(mJumps) + ".dat");
      mOutputFile << std::setprecision(20);
      mOutputFile << mPace << " " << mBufferSize << " " << mExtrapolationConstant << "\n";

      bool stop_extrapolation = false;
      for(unsigned int i = 0; i < mStateVariables.size(); i++){
          if(ExtrapolateState(i, stop_extrapolation)){
          extrapolated = true;
        }
        if(stop_extrapolation)
          break;
      }

      if(stop_extrapolation){
        mpModel->SetStateVariables(mSafeStateVariables);
        extrapolated=false;
      }

      mOutputFile.close();

      if(extrapolated){
        mJumps++;
        mpModel->SetStateVariables(mStateVariables);
        std::cout << "Extrapolated \n";

        // Debugging
        // std::cout << "new state variables are:\n";
        // for(auto variable : mStateVariables){
        //   std::cout << variable << " ";
        // }
        // std::cout<<"\n";
        ClearBuffers();
        mPreviousMinimalMRMSsPMCC = DOUBLE_UNSET;
        mLastExtrapolationPace=mPace;
        mMRMSBeforeExtrapolation=mCurrentMRMS;
      }
      mStatesBuffer.clear();
     return extrapolated;
    }
    else
      return false;
  }

void SmartSimulation::ClearBuffers(){
  Simulation::ClearBuffers();
  mStatesBuffer.clear();
  return;
}

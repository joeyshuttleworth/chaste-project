/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include "SmartSimulation.hpp"

bool SmartSimulation::ExtrapolateState(unsigned int state_index){
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

  if(pmcc>-0.75){  //return the unchanged value if there is no negative correlation (PMCC > -0.9 or PMCC = NAN)
    // std::cout <<  mpModel->GetSystemInformation()->rGetStateVariableNames()[state_index]<< ": PMCC was " << pmcc << " Ignoring. \n";
    return false;
  }

  if(beta > 0){
    //The difference is increasing
    // std::cout << mpModel->GetSystemInformation()->rGetStateVariableNames()[state_index] << ": Beta is positive\n";
    return false;
  }


  const double tau = -1/beta;
  double change_in_variable =  abs(mExtrapolationConstant * exp(alpha - mBufferSize/tau + 1/tau) / (exp(1/tau) - 1));

  /*Is V(t) increasing or decreasing?*/

  if(state.back() - state.front() < 0)
    change_in_variable = - change_in_variable;

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


bool SmartSimulation::RunPace(){
  bool extrapolated = false;

  extrapolated = ExtrapolateStates();
  if(!extrapolated){
    /*Solve in two parts*/
    try{
      mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
      mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
    }
    catch(Exception &e){
      std::cout << "RunPace failed - returning to old mStateVariables\n";
      mStateVariables = mSafeStateVariables;
      mpModel->SetStateVariables(mStateVariables);
      std::ofstream errors;
      mMrmsBuffer.clear();
      mStatesBuffer.clear();
      // The solver has been crashed so don't do any more extrapolations.
      mMaxJumps=0;
      if(mMaxJumps==0){
        throw std::exception();
      }
      return false;
    }
    std::vector<double> new_state_variables = GetStateVariables();
    mStatesBuffer.push_back(new_state_variables);
    mCurrentMrms = mrms(new_state_variables, mStateVariables);
    mMrmsBuffer.push_back(mCurrentMrms);
    mStateVariables = new_state_variables;
    if(mCurrentMrms < mThreshold){
      mFinished = true;
      return true;
    }
  }
  else{
    mCurrentMrms = 0;
  }
  return false;
}

bool SmartSimulation::ExtrapolateStates(){
    if(mJumps>=mMaxJumps)
      return false;
    if(!mMrmsBuffer.full())
      return false;
    double mrms_pmcc = CalculatePMCC(mMrmsBuffer);
    bool extrapolated = false;
    std::string model_name = mpModel->GetSystemInformation()->GetSystemName();
    const std::string dir_name = mOutputDir + model_name;
    boost::filesystem::create_directory(dir_name);
    if(mrms_pmcc < -0.985){
      mSafeStateVariables = mStateVariables;
      std::cout << "Extrapolating - start of buffer is " << mPaces - mBufferSize + 1<< "\n";

      boost::filesystem::create_directory(dir_name + "/LastExtrapolationLog");
      mOutputFile.open(dir_name + "/" + std::to_string(int(mPeriod)) + "JumpParameters.dat");
      mOutputFile << mPaces << " " << mBufferSize << " " << mExtrapolationConstant << "\n";

      for(unsigned int i = 0; i < mStateVariables.size(); i++){
        if(ExtrapolateState(i)){
          extrapolated = true;
        }
      }

      mOutputFile.close();

      if(extrapolated){
        mJumps++;
        mpModel->SetStateVariables(mStateVariables);
      }
      mMrmsBuffer.clear();
      mStatesBuffer.clear();
      return extrapolated;
    }
    else
      return false;
  }


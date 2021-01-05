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

  /* Run paces until max_paces is exceeded or the model reaches a steady state */
  bool RunPaces(int max_paces){
    for(int i = 0; i <= max_paces; i++){
      if(RunPace())
        return true;
    }
    return false;
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

#endif

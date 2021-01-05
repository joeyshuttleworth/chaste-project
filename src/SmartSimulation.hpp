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

#ifndef SMART_SIMULATION_HPP
#define SMART_SIMULATION_HPP

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include "Simulation.hpp"


class SmartSimulation : public Simulation{
public:
  SmartSimulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-7, double _tol_rel=1e-7, int _buffer_size = 200, double _extrapolation_constant = 1){
    mBufferSize = _buffer_size;
    mExtrapolationConstant = _extrapolation_constant;
    mpModel = _p_model;

    mFinished = false;
    mpStimulus = mpModel->UseCellMLDefaultStimulus();
    mpStimulus->SetStartTime(0);
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

    Simulation(_p_model, _period, input_path, _tol_abs, _tol_rel);
  }


  bool RunPace();
private:
  unsigned int mBufferSize;
  double  mExtrapolationConstant;
  boost::circular_buffer<std::vector<double>>  mStatesBuffer;
  boost::circular_buffer<double> mMrmsBuffer;
  unsigned int mJumps = 0;
  unsigned int mMaxJumps = 100;
  std::vector<double> mSafeStateVariables;
  unsigned int pace = 0;
  std::ofstream errors;

  bool ExtrapolateState(unsigned int state_index);
  bool ExtrapolateStates();
};

#endif

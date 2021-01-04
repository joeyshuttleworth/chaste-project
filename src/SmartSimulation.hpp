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

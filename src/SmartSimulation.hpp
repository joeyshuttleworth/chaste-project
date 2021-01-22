#ifndef SMART_SIMULATION_HPP
#define SMART_SIMULATION_HPP

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include "Simulation.hpp"


class SmartSimulation : public Simulation{
public:
  SmartSimulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-8, double _tol_rel=1e-8, int _buffer_size = 200, double _extrapolation_constant = 1, std::string output_dir = ""){
    mBufferSize = _buffer_size;
    mExtrapolationConstant = _extrapolation_constant;
    mpModel = _p_model;

    mFinished = false;
    mpStimulus = mpModel->UseCellMLDefaultStimulus();
    mpStimulus->SetStartTime(0);
    mpStimulus->SetPeriod(2*mPeriod);
    mPeriod = _period;
    mpModel->SetMaxSteps(1e5);
    mpModel->SetMaxTimestep(1000);
    mpModel->SetTolerances(mTolAbs, mTolRel);
    mpModel->SetMinimalReset(false); //Not sure if this is needed
    mNumberOfStateVariables = mpModel->GetSystemInformation()->rGetStateVariableNames().size();
    mStateVariables = mpModel->GetStdVecStateVariables();
    mSafeStateVariables=mStateVariables;
    if(input_path.length()>=1){
      LoadStatesFromFile(mpModel, input_path);
    }

    mMrmsBuffer.set_capacity(mBufferSize);
    mStatesBuffer.set_capacity(mBufferSize);

    if(output_dir=="")
      mOutputDir = std::string(getenv("CHASTE_TEST_OUTPUT"));
    else
      mOutputDir = output_dir;

    boost::filesystem::create_directories(mOutputDir);

    Simulation(_p_model, _period, input_path, _tol_abs, _tol_rel);
  }

  bool RunPace();

  // Getters and Setters

  void SetExtrapolationConstant(double e_c){mExtrapolationConstant=e_c;}
  void SetBufferSize(unsigned int buffer_size){
    mBufferSize = buffer_size;
    mMrmsBuffer.set_capacity(mBufferSize);
    mStatesBuffer.set_capacity(mBufferSize);
  }

  void SetMaxJumps(unsigned int max_jumps){mMaxJumps = max_jumps;}
  void SetOutputDir(std::string dir){
    mOutputDir = dir;
    boost::filesystem::create_directories(dir);
  }
private:
  std::string mOutputDir;
  unsigned int mBufferSize;
  double  mExtrapolationConstant;
  boost::circular_buffer<std::vector<double>>  mStatesBuffer;
  boost::circular_buffer<double> mMrmsBuffer;
  unsigned int mJumps = 0;
  unsigned int mMaxJumps = 1;
  std::vector<double> mSafeStateVariables;
  unsigned int pace = 0;
  std::ofstream errors;

  bool ExtrapolateState(unsigned int state_index);
  bool ExtrapolateStates();
};

#endif

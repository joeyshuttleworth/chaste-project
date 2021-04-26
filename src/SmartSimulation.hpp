#ifndef SMART_SIMULATION_HPP
#define SMART_SIMULATION_HPP

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include "Simulation.hpp"


class SmartSimulation : public Simulation{
public:
  SmartSimulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-8, double _tol_rel=1e-8, int _buffer_size = 200, double _extrapolation_constant = 1, std::string output_dir = "") : Simulation(_p_model, _period, input_path, _tol_abs, _tol_rel), mBufferSize(_buffer_size), mExtrapolationConstant(_extrapolation_constant){

    mNumberOfStateVariables = mpModel->GetSystemInformation()->rGetStateVariableNames().size();
    mSafeStateVariables=mStateVariables;
    mStatesBuffer.set_capacity(mBufferSize);

    if(output_dir=="")
      mOutputDir = std::string(getenv("CHASTE_TEST_OUTPUT"));
    else
      mOutputDir = output_dir;

    boost::filesystem::create_directories(mOutputDir);
  }

  bool RunPace();

  // Getters and Setters

  int GetNumberOfJumps(){
    return mJumps;
  }

  void SetExtrapolationConstant(double e_c){mExtrapolationConstant=e_c;}
  void SetBufferSize(unsigned int buffer_size){
    mBufferSize = buffer_size;
    mMRMSBuffer.set_capacity(mBufferSize);
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
  boost::circular_buffer<std::vector<double>>  mStatesBuffer;
  double  mExtrapolationConstant;
  unsigned int mJumps = 0;
  unsigned int mMaxJumps = 1000;
  std::vector<double> mSafeStateVariables;
  std::ofstream errors;

  unsigned int mLastExtrapolationPace = 0;
  double mMRMSBeforeExtrapolation = DOUBLE_UNSET;

  bool ExtrapolateState(unsigned int state_index, bool& stop_extrapolation);
  bool ExtrapolateStates();

  void ClearBuffers();
};

#endif

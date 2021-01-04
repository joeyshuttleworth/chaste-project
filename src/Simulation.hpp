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

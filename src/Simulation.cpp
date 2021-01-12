#include "Simulation.hpp"

Simulation::Simulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path, double _tol_abs, double _tol_rel) : mpModel(_p_model), mPeriod(_period), mTolAbs(_tol_abs), mTolRel(_tol_rel){
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


bool Simulation::RunPaces(int max_paces){
  for(int i = 0; i <= max_paces; i++){
    if(RunPace())
      return true;
  }
  return false;
}

bool Simulation::RunPace(){
  if(mFinished)
    return false;
  /*Solve in two parts*/
  std::vector<double> tmp_state_variables = mpModel->GetStdVecStateVariables();
  // mpModel->SetForceReset(true);
  mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
  mpStimulus->SetPeriod(mPeriod*2);
  mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
  mpStimulus->SetPeriod(mPeriod);
  std::vector<double> new_state_variables = mpModel->GetStdVecStateVariables();
  mCurrentMrms = mrms(tmp_state_variables, new_state_variables);
  mpModel->SetStateVariables(new_state_variables);
  if(mCurrentMrms < mThreshold){
    mFinished = true;
    return true;
  }
  return false;
}

OdeSolution Simulation::WritePaceToFile(std::string dirname, std::string filename, double sampling_timestep){
  OdeSolution solution = mpModel->Compute(0, mPeriod, 1);
  solution.WriteToFile(dirname, filename, "ms");
  return solution;
}

double Simulation::GetMrms(){
  if(mFinished)
    return NAN;
  else
    return mCurrentMrms;
}

bool Simulation::IsFinished(){
  return mFinished;
}

OdeSolution Simulation::GetPace(bool update_vars){
  const OdeSolution solution = mpModel->Compute(0, mPeriod, mSamplingTimestep);
  if(update_vars)
    mStateVariables=mpModel->GetStdVecStateVariables();
  else
    mpModel->SetStateVariables(mStateVariables);
  return solution;
}

std::vector<double> Simulation::GetStateVariables(){
  return mpModel->GetStdVecStateVariables();
}

void Simulation::SetStateVariables(std::vector<double> states){
  mpModel->SetStateVariables(states);
}

void Simulation::SetThreshold(double threshold){
  /* The threshold must be non-negative */
  if(threshold<0){
    throw std::exception();
  }
  mThreshold = threshold;
}

double Simulation::GetApd(double percentage, bool update_vars){
  mpStimulus->SetPeriod(mPeriod);
  OdeSolution solution = mpModel->Compute(mPeriod-10, 2*mPeriod, 1);
  solution.CalculateDerivedQuantitiesAndParameters(mpModel.get());
  CellProperties cell_props(solution.GetAnyVariable("membrane_voltage"), solution.rGetTimes());
  if(update_vars)
    mStateVariables = mpModel->GetStdVecStateVariables();
  else
    mpModel->SetStateVariables(mStateVariables);
  return cell_props.GetLastActionPotentialDuration(90);
}

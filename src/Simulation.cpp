#include "Simulation.hpp"
#include <iomanip>
#include <algorithm>
#include "SimulationTools.hpp"

Simulation::Simulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path, double _tol_abs, double _tol_rel) : mpModel(_p_model), mPeriod(_period), mTolAbs(_tol_abs), mTolRel(_tol_rel){
  mHasTerminated = false;
  mpStimulus = mpModel->UseCellMLDefaultStimulus();
  mpStimulus->SetStartTime(0);
  //There's no need to be near the second stimulus because Solve is called for
  //each pace
  mpStimulus->SetPeriod(2*mPeriod);
  mpModel->SetMaxSteps(1e5);
  mpModel->SetMaxTimestep(1000);
  mpModel->SetTolerances(mTolAbs, mTolRel);
  mNumberOfStateVariables = mpModel->GetSystemInformation()->rGetStateVariableNames().size();

  if(input_path.length()>=1){
    mStateVariables=LoadStatesFromFile(input_path);
    mpModel->SetStateVariables(mStateVariables);
  }
  else
    mStateVariables = mpModel->GetStdVecStateVariables();

  SetTolerances(_tol_abs, _tol_rel);

  mPreviousMinimalMRMSs.set_capacity(mPreviousMinimalMRMSsWindowSize);
  mPreviousMaximalMRMSs.set_capacity(mPreviousMinimalMRMSsWindowSize);
  mMRMSBuffer.set_capacity(mMRMSBufferSize);

  // Initialise GKr parameter
  SetIKrBlock(0);
}

Simulation::~Simulation(){
  mpModel->ResetToInitialConditions();
  mpStimulus->SetPeriod(mPeriod);
  mpModel->SetForceReset(true);

  // Return IKr block to its original value
  if(mDefaultGKr!=DOUBLE_UNSET)
    SetIKrBlock(0);
}

bool Simulation::RunPaces(int max_paces){
  RunPace();
  mpModel->SetForceReset(false);
  for(int i = 1; i < max_paces; i++){
    mHasTerminated=RunPace();
    if(mHasTerminated)
      return true;
  }
  return false;
}

bool Simulation::RunPace(){
  mPace++;
  if(mHasTerminated)
    return mHasTerminated;

  /*Solve in two parts*/
  std::vector<double> tmp_state_variables = mpModel->GetStdVecStateVariables();
  mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
  mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
  mStateVariables = mpModel->GetStdVecStateVariables();

  mCurrentMRMS = mrms(tmp_state_variables, mStateVariables);
  mMRMSBuffer.push_back(log(mCurrentMRMS));

  if(mCurrentMRMS < mCurrentMinimalMRMS || mCurrentMinimalMRMS == DOUBLE_UNSET)
    mCurrentMinimalMRMS = mCurrentMRMS;

  if(mPace % mMRMSBufferSize == 0 && mCurrentMinimalMRMS!=DOUBLE_UNSET){
    mPreviousMinimalMRMSs.push_back(mCurrentMinimalMRMS);
    mPreviousMaximalMRMSs.push_back(mCurrentMinimalMRMS);
    mCurrentMinimalMRMS = DOUBLE_UNSET;
    mCurrentMaximalMRMS = DOUBLE_UNSET;

    if(mPreviousMinimalMRMSs.full()){
      std::vector<double> log_minima, log_maxima;
      log_minima.reserve(mPreviousMinimalMRMSs.size());
      log_maxima.reserve(mPreviousMinimalMRMSs.size());

      for(auto val : mPreviousMinimalMRMSs)
        log_minima.push_back(log(val));

      for(auto val : mPreviousMaximalMRMSs)
        log_maxima.push_back(log(val));

      mPreviousMinimalMRMSsPMCC = CalculatePMCC(log_minima);
      mPreviousMaximalMRMSsPMCC = CalculatePMCC(log_maxima);
    }
  }

  /*  Check stopping criteria */
  if(mPreviousMinimalMRMSsPMCC!=DOUBLE_UNSET && mPreviousMaximalMRMSsPMCC!=DOUBLE_UNSET && mTerminateOnConvergence && mCurrentMRMS < mThreshold && mPreviousMinimalMRMSsPMCC != DOUBLE_UNSET && mPreviousMinimalMRMSsPMCC >= -0.05 && std::abs(mPreviousMaximalMRMSsPMCC) < 0.1 && mMRMSBuffer.full()){

    /*  Perform Dickey-Fuller test */

    /* Calculate auto-difference statistics */
    double sum_x = 0;
    double sum_x2 = 0;
    double n = mPreviousMinimalMRMSsWindowSize-1;
    for(unsigned int i = 0; i < n; i++){
      //Calculate auto-difference
      const double diff = (mMRMSBuffer[i+1] - mMRMSBuffer[i])/mMRMSBuffer[i];

      sum_x += diff;
      sum_x2+= diff*diff;
    }

    const double mean = sum_x/n;
    const double std_dev = sqrt(sum_x2/n - sum_x*sum_x/(n*n));

    const double test_statistic = std::abs(mean/(std_dev/n));

    // std::cout << "test_statistic is " << test_statistic << "\n";
    if(test_statistic<0.1){
      mHasTerminated = true;
      std::cout <<"terminating\n";
    }
    else
      mHasTerminated = false;
  }

  return mHasTerminated;
}

void Simulation::WritePaceToFile(std::string dirname, std::string filename, double sampling_timestep, bool update_vars){
  mpModel->SetForceReset(true);
  OdeSolution solution = mpModel->Compute(0, mPeriod, 1);
  solution.CalculateDerivedQuantitiesAndParameters(mpModel.get());
  solution.WriteToFile(dirname, filename, "ms", 1, false, 20, true);
  if(!update_vars)
    mpModel->SetStateVariables(mStateVariables);
  else
    mStateVariables = mpModel->GetStdVecStateVariables();
}

void Simulation::WriteStatesToFile(boost::filesystem::path dir, std::string filename){
  boost::filesystem::create_directories(dir);
  boost::filesystem::path filepath = (dir / boost::filesystem::path(filename));
  std::cout << filepath.string() << "\n";
  std::ofstream f_out(filepath.string());
  std::vector<double> states = GetStateVariables();
  std::vector<std::string> var_names = mpModel->GetSystemInformation()->rGetStateVariableNames();

  f_out << std::setprecision(20);

  for(std::string var_name : var_names){
    f_out << var_name << " ";
  }
  f_out << "\n";

  for(auto state_var : states){
    f_out << state_var << " ";
  }
  f_out << "\n";
  f_out.close();
}

double Simulation::GetMRMS(bool update){
  if(!mTerminateOnConvergence){
    std::vector<double> last_variables = mpModel->GetStdVecStateVariables();
    mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
    mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
    std::vector<double> new_variables = mpModel->GetStdVecStateVariables();
    mpModel->SetStateVariables(last_variables);
    return mrms(last_variables, new_variables);
  }
  else
    return mCurrentMRMS;
}

bool Simulation::HasTerminated(){
  return mHasTerminated;
}

OdeSolution Simulation::GetPace(double sampling_timestep, bool update_vars){
  mStateVariables = mpModel->GetStdVecStateVariables();
  const OdeSolution solution = mpModel->Compute(0, mPeriod, sampling_timestep);
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
  mStateVariables = states;
}

void Simulation::SetThreshold(double threshold){
  /* The threshold must be non-negative */
  if(threshold<0){
    throw std::exception();
  }
  mThreshold = threshold;
}

double Simulation::GetApd(double percentage, bool update_vars){
  // Need to use fine sampling timestep to be accurate, but has the effect of increasing
  // the solver tolerance and thus the shape of the action potential compared to lower tolerances
  if(!mpModel)
    return DOUBLE_UNSET;
  mpModel->SetMinimalReset(false);
  mStateVariables = mpModel->GetStdVecStateVariables();
  OdeSolution solution = mpModel->Compute(0, mPeriod, 0.01);
  mpStimulus->SetStartTime(0);
  solution.CalculateDerivedQuantitiesAndParameters(mpModel.get());
  const std::vector<double> times = solution.rGetTimes();
  const std::vector<double> voltages = solution.GetAnyVariable("membrane_voltage");
  if(times.size()==0)
    return DOUBLE_UNSET;
  CellProperties cell_props(voltages, times);

  if(update_vars)
    mStateVariables = mpModel->GetStdVecStateVariables();
  else
    mpModel->SetStateVariables(mStateVariables);
  const double return_val = cell_props.GetAllActionPotentialDurations(90).front();
  return return_val;
}

std::vector<double> Simulation::GetVoltageTrace(double sampling_timestep, bool update_vars){
  OdeSolution solution = GetPace();
  solution.CalculateDerivedQuantitiesAndParameters(mpModel.get());
  if(!update_vars)
    mpModel->SetStateVariables(mStateVariables);
  return solution.GetAnyVariable("membrane_voltage");
}

void Simulation::SetIKrBlock(double block){
  if(block > 1 || block < 0){
    EXCEPTION("Tried setting GKrConductance to an invalid amount");
  }

  const std::vector<std::string> parameter_names = mpModel->rGetParameterNames();

  mIKrBlock = block;

  {
    const std::string GKr_parameter_name = "membrane_rapid_delayed_rectifier_potassium_current_conductance";

    if(std::count_if(parameter_names.begin(), parameter_names.end(), [&](std::string name) -> bool {return name==GKr_parameter_name;}) > 0){
      // Parameter exists
      if(mDefaultGKr == DOUBLE_UNSET)
        mDefaultGKr = mpModel->GetParameter(GKr_parameter_name);
      if(block==0)
        mpModel->SetParameter(GKr_parameter_name, mDefaultGKr);
      else
        mpModel->SetParameter(GKr_parameter_name, mDefaultGKr*(1-block));
      return;
    }
  }

  const std::string GKr_scaling_factor_name = "membrane_rapid_delayed_rectifier_potassium_current_conductance_scaling_factor";
  if(std::count_if(parameter_names.begin(), parameter_names.end(), [&](std::string name) -> bool {return name==GKr_scaling_factor_name;}) > 0){
    if(mDefaultGKr == DOUBLE_UNSET)
      mDefaultGKr = mpModel->GetParameter(GKr_scaling_factor_name);
    if(block==0)
      mpModel->SetParameter(GKr_scaling_factor_name, mDefaultGKr);
    else
      mpModel->SetParameter(GKr_scaling_factor_name, mDefaultGKr*(1-block));
    return;
  }
  else{
    EXCEPTION("Couldn't find parameter to adjust GKr");
  }
}

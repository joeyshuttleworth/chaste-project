#include "Simulation.hpp"
#include <iomanip>
#include <algorithm>

// These header files are generated from the cellml files provided at github.com/chaste/cellml
#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_epiCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "decker_2009Cvode.hpp"
#include "ohara_rudy_cipa_v1_2017Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"
#include "iyer_2004Cvode.hpp"
#include "ToRORd_dynCl_epiCvode.hpp"
#include "hund_rudy_2004Cvode.hpp"

// Analytic models
#include "decker_2009_analytic_voltageCvode.hpp"
#include "hund_rudy_2004_analytic_voltageCvode.hpp"
#include "iyer_2004_analytic_voltageCvode.hpp"
#include "ohara_rudy_2011_epi_analytic_voltageCvode.hpp"
#include "ohara_rudy_cipa_2017_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2006_epi_analytic_voltageCvode.hpp"
#include "ten_tusscher_2004_epi_analytic_voltageCvode.hpp"
#include "ToRORd_dyn_chloride_epi_analytic_voltageCvode.hpp"

Simulation::Simulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path, double _tol_abs, double _tol_rel) : mpModel(_p_model), mPeriod(_period), mTolAbs(_tol_abs), mTolRel(_tol_rel){
  mFinished = false;
  mpStimulus = mpModel->UseCellMLDefaultStimulus();
  mpStimulus->SetStartTime(0);
  //There's no need to be near the second stimulus because Solve is called for
  //each pace
  mpStimulus->SetPeriod(2*mPeriod);
  mpModel->SetMaxSteps(1e5);
  mpModel->SetMaxTimestep(1000);
  mNumberOfStateVariables = mpModel->GetSystemInformation()->rGetStateVariableNames().size();
  if(input_path.length()>=1){
    LoadStatesFromFile(mpModel, input_path);
  }
  mStateVariables = mpModel->GetStdVecStateVariables();
  SetTolerances(_tol_abs, _tol_rel);
}

Simulation::~Simulation(){
  mpModel->SetStateVariables(mStateVariables);
  mpStimulus->SetPeriod(mPeriod);
  mpModel->SetForceReset(true);

  // Return IKr block to its original value
  if(mDefaultGKr!=DOUBLE_UNSET)
    SetIKrBlock(mDefaultGKr);
}

bool Simulation::RunPaces(int max_paces){
  RunPace();
  mpModel->SetForceReset(false);
  for(int i = 1; i < max_paces; i++){
    if(RunPace() || mFinished)
      return true;
  }
  return false;
}

bool Simulation::RunPace(){
  mPaces++;
  if(mFinished)
    return false;
  if(mTerminateOnConvergence){
    /*Solve in two parts*/
    std::vector<double> tmp_state_variables = mpModel->GetStdVecStateVariables();
    mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
    mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
    std::vector<double> mStateVariables = mpModel->GetStdVecStateVariables();
    mCurrentMrms = mrms(tmp_state_variables, mStateVariables);
    if(mCurrentMrms < mThreshold){
       mFinished = true;
       mStateVariables = mpModel->GetStdVecStateVariables();
       std::cout << "finished after " << mPaces << " paces \n";
      return true;
    }
  }
  else{
    mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
    mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
  }
  mStateVariables = mpModel->GetStdVecStateVariables();
  return false;
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

double Simulation::GetMrms(bool update){
  if(!mTerminateOnConvergence){
    std::vector<double> last_variables = mpModel->GetStdVecStateVariables();
    mpModel->SolveAndUpdateState(0, mpStimulus->GetDuration());
    mpModel->SolveAndUpdateState(mpStimulus->GetDuration(), mPeriod);
    std::vector<double> new_variables = mpModel->GetStdVecStateVariables();
    mpModel->SetStateVariables(last_variables);
    return mrms(last_variables, new_variables);
  }
  else
    return mCurrentMrms;
}

bool Simulation::IsFinished(){
  return mFinished;
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

  const std::string GKr_parameter_name = "membrane_rapid_delayed_rectifier_potassium_current_conductance";

  const std::vector<std::string> parameter_names = mpModel->rGetParameterNames();

  if(std::count_if(parameter_names.begin(), parameter_names.end(), [&](std::string name) -> bool {return name==GKr_parameter_name;}) > 0){
    // Parameter exists
    if(mDefaultGKr == DOUBLE_UNSET)
      mDefaultGKr = mpModel->GetParameter(GKr_parameter_name);
    mpModel->SetParameter(GKr_parameter_name, mDefaultGKr*(1-block));
    return;
  }

  const std::string GKr_scaling_factor_name = "membrane_rapid_delayed_rectifier_potassium_current_conductance_scaling_factor";
  if(std::count_if(parameter_names.begin(), parameter_names.end(), [&](std::string name) -> bool {return name==GKr_scaling_factor_name;}) > 0){
    if(mDefaultGKr == DOUBLE_UNSET)
      mDefaultGKr = mpModel->GetParameter(GKr_scaling_factor_name);
    mpModel->SetParameter(GKr_scaling_factor_name, mDefaultGKr*(1-block));
    return;
  }
  else{
    EXCEPTION("Couldn't find parameter to adjust GKr");
  }
}



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

#include "beeler_reuter_model_1977Cvode.hpp"
#include "ten_tusscher_model_2004_epiCvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "SimulationTools.hpp"

class Simulation
{
private:
protected:
  bool finished;
  boost::shared_ptr<AbstractCvodeCell> p_model;
  unsigned int number_of_state_variables;
  std::vector<double> state_variables;
  std::ofstream output_file;
  double period = 1000;
  // double TolAbs = 1e-8, TolRel = 1e-8;
  double TolAbs;
  double TolRel;
  double sampling_timestep = 1;
  double current_mrms = NAN;
  const double threshold = 1.8e-07;
  boost::shared_ptr<RegularStimulus> p_stimulus;
public:
  Simulation(){
    return;
  }
  Simulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-7, double _tol_rel=1e-7) : p_model(_p_model), period(_period), TolAbs(_tol_abs), TolRel(_tol_rel){
    finished = false;
    p_stimulus = p_model->UseCellMLDefaultStimulus();
    p_stimulus->SetStartTime(0);
    //There's no need to be near the second stimulus because Solve is called for
    //each pace
    p_stimulus->SetPeriod(2*period);
    p_model->SetMaxSteps(1e5);
    p_model->SetMaxTimestep(1000);
    p_model->SetTolerances(TolAbs, TolRel);
    p_model->SetMinimalReset(false); //Not sure if this is needed
    number_of_state_variables = p_model->GetSystemInformation()->rGetStateVariableNames().size();
    state_variables = p_model->GetStdVecStateVariables();
    if(input_path.length()>=1){
      LoadStatesFromFile(p_model, input_path);
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
    if(finished)
      return false;
    /*Solve in two parts*/
    std::vector<double> tmp_state_variables = p_model->GetStdVecStateVariables();
    p_model->SolveAndUpdateState(0, p_stimulus->GetDuration());
    p_stimulus->SetPeriod(period*2);
    p_model->SolveAndUpdateState(p_stimulus->GetDuration(), period);
    p_stimulus->SetPeriod(period);
    std::vector<double> new_state_variables = p_model->GetStdVecStateVariables();
    current_mrms = mrms(tmp_state_variables, new_state_variables);
    p_model->SetStateVariables(new_state_variables);
    // p_model->SetVoltage(p_model->CalculateAnalyticVoltage());
    if(current_mrms < threshold){
      finished = true;
      return true;
    }
    return false;
  }

  /**Output a pace to file*/
  void WriteToFile(std::string filename){
    OdeSolution solution = p_model->Compute(0, p_stimulus->GetDuration(), 1);
    solution.WriteToFile(filename, filename, "ms");
    solution = p_model->Compute(p_stimulus->GetDuration(), period, 1);
    solution.WriteToFile(filename, filename, "ms");
    return;
  }

  double GetMrms(){
    if(finished)
      return NAN;
    else
      return current_mrms;
  }
  bool is_finished(){
    return finished;
  }
  std::vector<double> GetStateVariables(){
    return p_model->GetStdVecStateVariables();
  }
};

#endif

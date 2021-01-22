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
protected:
  bool mFinished;
  boost::shared_ptr<AbstractCvodeCell> mpModel;
  unsigned int mNumberOfStateVariables;
  std::vector<double> mStateVariables;
  std::ofstream mOutputFile;
  double mPeriod = 1000;
  double mTolAbs;
  double mTolRel;
  double mCurrentMrms = NAN;
  double mThreshold = 1e-08;
  boost::shared_ptr<RegularStimulus> mpStimulus;
  bool mTerminateOnConvergence = true;
public:
  Simulation(){
    return;
  }

  Simulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-8, double _tol_rel=1e-8);

  ~Simulation();

  /* Run paces until max_paces is exceeded or the model reaches a steady state */
  bool RunPaces(int);

  void SetTolerances(double atol, double rtol){
    if(mpModel)
      mpModel->SetTolerances(atol, rtol);
    mTolAbs = atol;
    mTolRel = rtol;
  }

  bool RunPace();

  /**Output a pace to file*/
  void WritePaceToFile(std::string dirname, std::string filename, double sampling_timestep = 1, bool update_variables=false);

  void WriteStatesToFile(boost::filesystem::path dirname, std::string filename);

  OdeSolution GetPace(double sampling_timestep = 1, bool update_vars=false);

  double GetMrms();

  std::vector<double> GetStateVariables();

 // Setters and getters
  void SetThreshold(double threshold);

  void SetStateVariables(std::vector<double> states);

  double GetApd(double percentage = 90, bool update_vars=false);

  std::vector<double> GetVoltageTrace(double sampling_timestep = 1, bool update_vars=false);

  void SetTerminateOnConvergence(bool b){mTerminateOnConvergence=b;}

  bool IsFinished();

  double GetMrms(bool update=false);

  boost::shared_ptr<AbstractCvodeCell> GetModel(){return mpModel;}

};

#endif

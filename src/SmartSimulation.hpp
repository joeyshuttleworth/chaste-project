#ifndef SMART_SIMULATION_HPP
#define SMART_SIMULATION_HPP

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include "Simulation.hpp"


class SmartSimulation : public Simulation{


  SmartSimulation(boost::shared_ptr<AbstractCvodeCell> _p_model, double _period, std::string input_path = "", double _tol_abs=1e-7, double _tol_rel=1e-7, _buffer_size = 200, _extrapolation_constant = 1){
    mBufferSize = _buffer_size;
    mExtrapolationConstant = _extrapolation_constant;
    Simulation(_p_model, _period, input_path, _tol_abs, _tol_rel);
  }


private:
  unsigned int mBufferSize;
  double  mExtrapolationConstant;
  boost::circular_buffer<std::vector<double>>  states_buffer;
  boost::circular_buffer<double> mrms_buffer;
  unsigned int jumps = 0;
  unsigned int max_jumps = 100;
  std::vector<double> mSafeStateVariables;
  unsigned int pace = 0;
  std::ofstream errors;

  bool ExtrapolateState(unsigned int state_index);
  bool ExtrapolateStates();
}

#endif

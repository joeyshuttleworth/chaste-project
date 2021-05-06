#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include "Simulation.hpp"
#include "SmartSimulation.hpp"
#include "SimulationTools.hpp"
#include "CellProperties.hpp"
#include <iomanip>

#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>

class TestOneExtrapolation : public CxxTest::TestSuite
{
private:
public:
  void TestRun(){
    int max_paces = get_max_paces();
    auto extrapolation_constants = get_extrapolation_constants();
    auto buffer_sizes = get_buffer_sizes();
    auto periods = get_periods();
    auto blocks = get_IKr_blocks();
    auto models = get_analytic_models();

    std::cout << std::setprecision(20);

    std::ofstream output("TestOneExtrapolation_output.dat");

    for(auto model : models){
      for(auto e_c : extrapolation_constants){
        for(auto bs : buffer_sizes){
          for(auto period : periods){
            for(auto block : blocks){
              // Just change IKrBlock
              OutputScore("IKrBlock", model, period, block, e_c, bs, output);
            }
          }
        }
      }
    }
  }
};

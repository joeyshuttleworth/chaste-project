#include <cxxtest/TestSuite.h>
#include <sstream>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"
#include "SimulationTools.hpp"
#include <boost/filesystem.hpp>
#include <fstream>

class TestGroundTruthSimulation : public CxxTest::TestSuite
{
public:
  std::string username;
  void TestTolerances(){
    const std::vector<double> tolerances={1e-8, 1e-6, 1e-04, 1e-10, 1e-12};
    const std::vector<double> periods={1000};
    const std::vector<double> IKrBlocks={0};
    const std::string filename_suffix = "test_tolerances";

    auto models = get_models();

    const int default_max_paces = 10000;
    int paces = get_max_paces();

    paces = paces==INT_UNSET?default_max_paces:paces;

    for(auto tolerance : tolerances){
      for(auto IKrBlock : IKrBlocks){
        for(auto period : periods){
           for(auto model : models){
             compare_error_measures(paces, model, period, IKrBlock, tolerance, filename_suffix);
          }
        }
      }
    }
  }
};

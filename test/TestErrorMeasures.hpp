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

class TestErrorMeasures : public CxxTest::TestSuite
{
public:
  std::string username;
  const int default_max_paces = 10000;
  void TestErrors(){
    int paces = get_max_paces();
    paces = paces==INT_UNSET?default_max_paces:paces;

    const std::vector<double> periods= get_periods();
    const std::vector<double> IKrBlocks= get_IKr_blocks();

    std::vector<boost::shared_ptr<AbstractCvodeCell>> models = get_models();

    const std::string filename_suffix = "error_measures";
    const double tolerance = 1e-8;

    boost::shared_ptr<RegularStimulus> p_stimulus;
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

    for(auto model : models){
      for(auto period : periods){
        for(auto IKrBlock : IKrBlocks){
          compare_error_measures(paces, model, period, IKrBlock, tolerance, filename_suffix);
        }
      }
    }
  }
};

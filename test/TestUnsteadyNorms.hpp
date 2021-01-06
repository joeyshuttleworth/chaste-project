/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "TenTusscher2006EpiCvode.hpp"
#include "FakePetscSetup.hpp"
#include <fstream>

/*Make sure the sampling timestep divides the stimulus duration*/

class TestUndsteadyNorms : public CxxTest::TestSuite
{
public:
    void TestSimulation()
    {
#ifdef CHASTE_CVODE
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellTenTusscher2006EpiFromCellMLCvode(p_solver, p_stimulus));
        boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
	const double period = 500.0;
	const unsigned int number_of_variables = p_model->GetNumberOfStateVariables();
	  
        p_regular_stim->SetPeriod(period);
     	p_model->SetTolerances(1e-9,1e-9);

	double max_timestep = p_regular_stim->GetDuration()/2;

        p_model->SetMaxTimestep(max_timestep);

        double steps = 100;
	
	OdeSolution *current_solution = NULL, *last_solution = NULL;
	std::ofstream changes_file;
	std::string username = std::string(getenv("USER"));
	
	changes_file.open("/tmp/"+username+"/changes.ssv");
	
	for(int i=0; i < steps; i++){
	  double start_time = i*period;
	  double end_time   = start_time + period;
	  
	  if(last_solution)
	    delete(last_solution);
	  if(current_solution)
	    last_solution = current_solution;

	  current_solution  = new OdeSolution;
	  *current_solution = p_model->Compute(start_time, end_time, max_timestep); 

	  /*Set the initial values to be the terminal values of the last solution*/
	  p_model->SetStateVariables(current_solution->rGetSolutions().back());

	  if(i>1){
	    /*Calculate the change between the state variables over the last two beats using the Euclidean norm*/
	    double sum_2 = 0;
	    double sum_mrms = 0;
	    std::vector<double> last_times = last_solution->rGetTimes();
	    std::vector<double> current_times = current_solution->rGetTimes();
	    for(unsigned int j = 0; j < number_of_variables; j++){
	      for(unsigned int k = 0; k < current_solution->GetNumberOfTimeSteps(); k++){
          sum_2 = sum_2 + pow(current_solution->rGetSolutions()[k][j] - last_solution->rGetSolutions()[k][j], 2);
	      }
	    }
	    changes_file << sqrt(sum_2) << " " << sqrt(sum_mrms / number_of_variables) << "\n";
	  }
	}
#else
        std::cout << "Cvode is not enabled.\n";
#endif
	}
};

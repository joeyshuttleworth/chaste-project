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

#include "CellProperties.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/circular_buffer.hpp>
#include <string>
#include <sstream>
#include <iostream>

const std::vector<std::string> model_names = {"beeler_reuter_model_1997", "ten_tusscher_model_2004_epi", "ohara_rudy_2011_endo", "shannon_wang_puglisi_weber_bers_2004"};

void RunSimulation(boost::shared_ptr<AbstractCvodeCell>, unsigned int paces, double tolerances);

int LoadStatesFromFile(boost::shared_ptr<AbstractCvodeCell>, std::string file_path);

void OutputVariablesToFile(boost::shared_ptr<AbstractCvodeCell>, std::string file_path);

std::vector<double> GetNthVariable(std::vector<std::vector<double>>, unsigned int);

std::vector<double> cGetNthVariable(boost::circular_buffer<std::vector<double>>, unsigned int);

double mrms(std::vector<double>, std::vector<double>);

double TwoNorm(std::vector<double>, std::vector<double>);

double mrmsTrace(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

double TwoNormTrace(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

double CalculateAPD(boost::shared_ptr<AbstractCvodeCell>, double, double, double);

std::vector<double> FitExponential(std::vector<double> x_vals, std::vector<double> y_vals);

double CalculatePace2Norm(boost::shared_ptr<AbstractCvodeCell> p_model, std::vector<double> first_states, std::vector<double> second_states, double period, double duration);

double CalculatePaceMrms(boost::shared_ptr<AbstractCvodeCell> p_model, std::vector<double> first_states, std::vector<double> second_states, double period, double duration);

void WriteStatesToFile(std::vector<double> states, std::ofstream &f_out);

std::vector<std::vector<double>> GetPace(std::vector<double> initial_conditions, boost::shared_ptr<AbstractCvodeCell> p_model, double period, double duration); 

double CalculatePMCC(std::vector<double>, std::vector<double>);

template<typename Container>
double CalculatePMCC(Container values){
  const unsigned int N = values.size();
  // const double sum_x = N*(N-1)/2;
  //const double sum_x2 = (N-1)*N*(2*N-1)/6;

  if(values.size()<=2 || values.size() <= 2){
    return -NAN;
  }
  double sum_x = 0, sum_x2 = 0, sum_y = 0, sum_y2 = 0, sum_xy = 0;

  for(unsigned int i = 0; i < N; i++){
    sum_x += i;
    sum_x2+= i*i;
    sum_y += values[i];
    sum_y2+= values[i]*values[i];
    sum_xy+= i*values[i];
  }

  double pmcc = (N*sum_xy - sum_x*sum_y)/sqrt((N*sum_x2 - sum_x*sum_x)*(N*sum_y2 - sum_y*sum_y));
  return pmcc;
}



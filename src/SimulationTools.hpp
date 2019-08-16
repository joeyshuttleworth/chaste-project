#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <string>
#include <sstream>
#include <iostream>

const std::vector<std::string> model_names = {"beeler_reuter_model_1997", "ten_tusscher_model_2004_epi", "ohara_rudy_2011_endo", "shannon_wang_puglisi_weber_bers_2004"};

void RunSimulation(boost::shared_ptr<AbstractCvodeCell>, unsigned int paces, double tolerances);

int LoadStatesFromFile(boost::shared_ptr<AbstractCvodeCell>, std::string file_path);

void OutputVariablesToFile(boost::shared_ptr<AbstractCvodeCell>, std::string file_path);

std::vector<double> GetNthVariable(std::vector<std::vector<double>>, unsigned int);

double mrms(std::vector<double>, std::vector<double>);

double TwoNorm(std::vector<double>, std::vector<double>);

double mrmsTrace(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

double TwoNormTrace(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

double CalculateAPD(boost::shared_ptr<AbstractCvodeCell>, boost::shared_ptr<RegularStimulus>, double);

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

void RunSimulation(boost::shared_ptr<AbstractCvodeCell>, unsigned int paces, double tolerances);

void LoadStatesFromFile(boost::shared_ptr<AbstractCvodeCell>, std::string file_path);

void OutputVariablesToFile(boost::shared_ptr<AbstractCvodeCell>, std::string file_path);

std::vector<double> GetNthVariable(std::vector<std::vector<double>>*, unsigned int);

double mrms(std::vector<double>, std::vector<double>);

double TwoNorm(std::vector<double>, std::vector<double>);

double mrmsTrace(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

double TwoNormTrace(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

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

class Simulation
{
private:
  std::vector<std::vector<double>> results;
  boost::shared_ptr<AbstractCvodeCell> p_model;
  double tol_abs = 1e-12, tol_rel = 1e-12;
  std::ofstream output_file;
  double period;
  double sampling_timestep = 1;
public:
  Simulation(double period, std::string file_path = "", std::string model_name, double tol_abs = 1e-12, double tol_rel = 1e-12), double sampling_timestep = 1){
    boost::shared_ptr<AbstractCvodeCell> p_stimulus;
    if(model_name == "ohara_rudy_2011_endo"){
      p_model.reset(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, p_stimulus));
    }
    p_stimulus = p_model->UseCellMLDefaultStimulus();
    p_stimulus->SetStartTime(0);
    p_stimulus->SetPeriod(period);
    p_model->SetMaxSteps(1e5);
    p_model->SetMaxTimestep(1000);
    TS_ASSERT(period < p_stimulus->GetDuration());
    if(file_path != ""){
      output_file.open(file_path);
      boost::filesystem::create_director(file_path);
      output_file.open(file_path);
      if(!output_file.is_open()){
	std::cout << "Failed to open file!\n";
	return;
      }
      else{
	output_file << p_model->GetSystemInformation()->GetSystemName() << " ";
	output_file << tol_abs << " " << tol_rel << " " << period;
      }
    }

  void Run(unsigned int paces){
    for(unsigned int i = 0; i < paces; i++){
      /*Solve in two parts*/
      p_model->Solve(0, p_stimulus->GetDuration());
      p_model->Solve(p_stimulus->GetDuration(), period);
    }
    return;
  }
  
  void RunAndSave(unsigned int paces){
    for(unsigned int i = 0; i < paces; i++){
      OdeSolution solution = p_model->compute(0, p_stimulus->GetDuration(), sampling_timestep);
      std::vector<std::vector<double>> states = solution->rGetSolutions();
      results.insert(results.end(), states.begin(), states.end());
      OdeSolution solution = p_model->compute(0, p_stimulus->GetDuration(), sampling_timestep);
      std::vector<std::vector<double>> states = solution->rGetSolutions();      
      results.insert(results.end(), states.begin(), states.end());
    }
    return;
  }

  void WriteToFile(std::string){
    std::vector<std::string> state_names = p_model->rGetStateVariableNames();
    for(unsigned int i = 0; i < state_names.size(); i++){
      output_file << state_names[i];
    }
    output_file << "\n";
    for(unsigned int i = 0; i < results.size(); i++){
      for(unsigned int j = 0; j < results[0].size(); j++){
	output_file<<results[i][j]<<",";
      }
      output_file << "\n";
    }
  }
  
  void SetResults(std::vector<std::vector<double>> r){
    results = r;
  }
  
  std::vector<std::vector<double>> GetResults(){
    return results;
  }
  
  void ClearResults(){
    results.clear();
  }
}

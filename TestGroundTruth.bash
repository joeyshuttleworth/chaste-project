#!/bin/bash

# Run TestGroundTruthSimulation in multiple processes using the command line arguments
declare -a models=("shannon_wang_puglisi_weber_bers_2004_model_updated" "beeler_reuter_model_1977" "Tomek2020epi" "IyerMazhariWinslow2004" "HundRudy2004_units" "decker_2009" "tentusscher_model_2004_epi" "tentusscher_model_2006_epi" "ohara_rudy_2011_epi" "ohara_rudy_cipa_v1_2017" "Tomek2020epi_analytic_voltage" "IyerMazhariWinslow2004_analytic_voltage" "HundRudy2004_analytic_voltage_units" "decker_2009_analytic_voltage" "tentusscher_model_2004_epi_analytic_voltage" "tentusscher_model_2006_epi_analytic_voltage" "ohara_rudy_2011_epi_analytic_voltage" "ohara_rudy_cipa_2017_analytic_voltage")

periods=(1000 500 750 1250)
IKrBlocks=(0 0.25 0.5);

TestGroundTruth_bin=$1

mkdir log
mkdir log/TestGroundTruth

for IKrBlock in ${IKrBlocks[@]}; do
    for period in ${periods[@]}; do
        for model in ${models[@]}; do
            ${TestGroundTruth_bin}/projects/chaste-project/test/TestGroundTruthSimulation --models $model --periods $period --IKrBlocks $IKrBlocks & >> log/TestGroundTruth/${model}_${period}_${IKrBlock}.log;
        done;
    done;
done;

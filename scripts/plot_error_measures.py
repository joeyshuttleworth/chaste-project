#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

test_dir =  os.path.expanduser("~/Chaste-from-server/testoutput")

def plot_error_measure(measure, period=1000, IKrBlock=0):
    initial_period = 500 if period == 1000 else 500
    initial_block = 0.5 if IKrBlock == 0 else 0
    for dir in os.listdir(test_dir):
        if re.search("_{}ms_{}_percent".format(period, int(100*IKrBlock)), dir):
            model_name = re.search("([A-Z|a-z|0-9|_]*)_{}ms_{}_percent".format(period, int(IKrBlock*100)), dir).group(1)
            print(dir)
            file_path = os.path.join(os.path.join(test_dir, dir, "error_measures_1e-08.dat"))
            # apd_file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace = True)
                y_vals = np.log10(df[measure])  # df[['MRMS', '2-Norm', "Trace-2-Norm", "Trace-MRMS"]].values)
                apds = df["APD"].values

                apd_file_path = os.path.join(os.path.join(os.path.expanduser("~/Chaste-from-server"), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
                true_apd = pd.read_csv(apd_file_path).values[-1,0]

                apd_errors = np.log10(abs(apds - true_apd))

                data_to_plot = np.array([[apd_error, y] for [apd_error, y] in zip(apd_errors, y_vals) if math.isfinite(apd_error)])
                plt.plot(data_to_plot[:,0], data_to_plot[:,1], label=model_name, markevery=1000)
            else:
                print("{} Does not exist".format(file_path))
    plt.xlabel("log10 absolute APD90 error")
    plt.ylabel("{} error between sucessive paces".format(measure))
    fig = plt.gcf()
    fig.set_size_inches(12, 10)
    plt.savefig("compare_measures_{}.pdf".format(measure))


def plot_half_log_error_measure(measure, model_name):
    print(model_name)
    for dir in os.listdir(test_dir):
        if re.search("{}_1000ms_0_percent_block".format(model_name), dir):
            print(dir)
            file_path = os.path.join(os.path.join(test_dir, dir, "error_measures_1e-08.dat"))
            # apd_file_path = os.path.join(os.path.join(test_dir, dir, "apds_using_groundtruth_1e-12.dat"))
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace = True)
                y_vals = np.log10(df[measure])  # df[['MRMS', '2-Norm', "Trace-2-Norm", "Trace-MRMS"]].values)
                apds = df["APD"].values
                apd_errors = [np.log10(abs(apd - apds[-1])) for apd in apds]
                # true_apd = pd.read_csv(apd_file_path).values[-1,0]
                plt.plot(y_vals[0:3000], markevery=1000)
            else:
                print("{} Does not exist".format(file_path))

    plt.xlabel("pace")
    plt.ylabel("log_10 pace-to-pace MRMSE".format(measure))
    # plt.title("{} various scenarios".format(model_name))
    plt.show()
    # plt.savefig("compare_scenarios_{}_{}".format(measure,model_name))


def main():
    measures = ["MRMS", "2-Norm", "Trace-MRMS", "Trace-2-Norm"]
    # plot_error_measure("MRMS")
    # plot_error_measure("2-Norm")
    # plot_error_measure("Trace-MRMS")
    # plot_error_measure("Trace-2-Norm")

    # compare_scenarios("MRMS", "Tomek2020epi_analytic_voltage")
    # plot_half_log_error_measure("MRMS", "tentusscher_model_2004_epi_analytic_voltage")

    for measure in measures:
        print("using measure {}".format(measure))
        # compare_scenarios(measure, "HundRudy2004_analytic_voltage_units")
        # compare_scenarios(measure, "IyerMazhariWinslow2004")
        plot_error_measure(measure)

if __name__ == "__main__":
    main()

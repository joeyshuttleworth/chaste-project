#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_error_measure(measure):
    for dir in os.listdir("testoutput/"):
        if re.search("_1000ms_0_percent", dir):
            model_name = re.search("([a-z|0-9|_]*)_1000ms_0_percent", dir).group(1)
            print(dir)
            file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "error_measures_1e-10.dat"))
            apd_file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace = True)
                y_vals = np.log10(df[measure])  # df[['MRMS', '2-Norm', "Trace-2-Norm", "Trace-MRMS"]].values)
                apds = df["APD"].values
                apd_errors = [np.log10(abs(apd - apds[-1])) for apd in apds]
                true_apd = pd.read_csv(apd_file_path).values[-1,0]
                plt.plot(apd_errors, y_vals, label=model_name)
            else:
                print("{} Does not exist".format(file_path))
    plt.legend()
    plt.xlabel("log10 error in APD90")
    plt.ylabel("{} error between sucessive paces".format(measure))
    # plt.savefig("compare_measures_{}".format(measure))
    plt.show()


def compare_scenarios(measure, model_name):
    for dir in os.listdir("testoutput/"):
        print(dir)
        if re.search(model_name, dir):
            print(dir)
            file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "error_measures_1e-10.dat"))
            apd_file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace = True)
                y_vals = np.log10(df[measure])  # df[['MRMS', '2-Norm', "Trace-2-Norm", "Trace-MRMS"]].values)
                apds = df["APD"].values
                apd_errors = [np.log10(abs(apd - apds[-1])) for apd in apds]
                true_apd = pd.read_csv(apd_file_path).values[-1,0]
                plt.plot(apd_errors, y_vals)
            else:
                print("{} Does not exist".format(file_path))

    plt.xlabel("log10 error in APD90")
    plt.ylabel("{} error between sucessive paces".format(measure))
    plt.title("{} various scenarios".format(model_name))
    plt.show()
    # plt.savefig("compare_scenarios_{}_{}".format(measure,model_name))


def main():
    measures = ["MRMS", "2-Norm", "Trace-MRMS", "Trace-2-Norm"]
    # plot_error_measure("MRMS")
    # plot_error_measure("2-Norm")
    # plot_error_measure("Trace-MRMS")
    # plot_error_measure("Trace-2-Norm")

    for measure in measures:
        print("using measure {}".format(measure))
        compare_scenarios(measure, "ohara_rudy_2011")
        compare_scenarios(measure, "decker_2009")
        plot_error_measure(measure)

if __name__ == "__main__":
    main()

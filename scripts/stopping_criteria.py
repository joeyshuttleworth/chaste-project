#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

def plot_moving_average(measure, model):
    window_size = 20
    for dir in os.listdir("testoutput/"):
        if re.search("{}_1000ms_0_percent".format(model), dir):
            model_name = re.search("([a-z|0-9|_]*)_1000ms_0_percent", dir).group(1)
            plot = False
            tol = "1e-10"
            file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "test_tolerances_{}.dat".format(tol)))
            apd_file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace = True)
                y_vals = np.log10(df[measure])
                print("y_vals length is {}".format(len(y_vals)))
                apds = df["APD"].values
                true_apd = pd.read_csv(apd_file_path).values[-1,0]
                apd_errors = [np.log10(abs(apd - true_apd)) for apd in apds]
                average_label = "rolling_average_size_{}".format(window_size)
                print(average_label)
                df["average"] = np.log10(df[measure].rolling(window=window_size).mean())
                df["pmcc"] = np.log10(df[measure].rolling(window=window_size).mean())
                plt.plot(y_vals, label=tol)
                df["average"].plot(label=average_label)
                print("plotted")
                plot=True
            if plot:
                plt.legend()
                plt.title(model_name)
                plt.xlabel("log10 error in APD90")
                plt.ylabel("{} error between sucessive paces".format(measure))
                # plt.savefig("compare_measures_{}".format(measure))
                plt.show()

def main():
    measures = ["MRMS"]
    models   = ["tentusscher_model_2006_epi_analytic", "tentusscher_model_2004_epi", "ohara_rudy_2011"]
    for model in models:
        for measure in measures:
            plot_moving_average(measure, model)

if __name__ == "__main__":
    main()










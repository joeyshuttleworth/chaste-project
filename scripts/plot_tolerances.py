#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_tolerances(measure, model):
    for dir in os.listdir("testoutput/"):
        if re.search("{}_1000ms_0_percent".format(model), dir):
            model_name = re.search("([A-Z|a-z|0-9|_]*)_1000ms_0_percent", dir).group(1)
            plot = False
            for tol in ["0.0001", "1e-12", "1e-06", "1e-08", "1e-10"]:
                file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "test_tolerances_{}.dat".format(tol)))
                apd_file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
                if os.path.exists(file_path):
                    df = pd.read_csv(file_path, delim_whitespace = True)
                    y_vals = np.log10(df[measure])
                    print("y_vals length is {}".format(len(y_vals)))
                    apds = df["APD"].values
                    true_apd = pd.read_csv(apd_file_path).values[-1,0]
                    apd_errors = [np.log10(abs(apd - true_apd)) for apd in apds]
                    plt.plot(apd_errors, y_vals, label=tol)
                    print("plotted")
                    plot=True
            if plot:
                plt.legend()
                plt.title(model_name)
                plt.xlabel("log10 error in APD90")
                plt.ylabel("{} error between sucessive paces".format(measure))
                plt.show()

def main():
    measures = ["MRMS", "2-Norm", "Trace-MRMS", "Trace-2-Norm"]
    models = ["Tomek2020epi_analytic_voltage", "IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]

    for model in models:
        for measure in measures:
            plot_tolerances(measure, model)

if __name__ == "__main__":
    main()










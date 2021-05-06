#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

def plot_tolerances(measure, model):
    for dir in os.listdir(os.path.expanduser("~/Chaste-from-server/testoutput/")):
        if re.search("{}_1000ms_0_percent".format(model), dir):
            print(dir)
            model_name = re.search("([A-Z|a-z|0-9|_]*)_1000ms_0_percent_block", dir).group(1)
            print(model_name)
            plot = False
            for tol in ["0.0001", "1e-06", "1e-08", "1e-10", "1e-12"]:
                file_path = os.path.join(os.path.join(os.path.expanduser("~/Chaste-from-server"), "testoutput", dir, "test_tolerances_{}.dat".format(tol)))
                print(file_path)
                apd_file_path = os.path.join(os.path.join(os.path.expanduser("~/Chaste-from-server"), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
                if os.path.exists(file_path):
                    df = pd.read_csv(file_path, delim_whitespace = True)
                    y_vals = np.log10(df[measure])
                    print("y_vals length is {}".format(len(y_vals)))
                    apds = df["APD"].values
                    true_apd = pd.read_csv(apd_file_path).values[-1,0]
                    apd_errors = [np.log10(abs(apd - true_apd)) for apd in apds]
                    # Join values together
                    data_to_plot = np.stack([apd_errors,y_vals], axis=1)
                    data_to_plot = np.array([val for val in data_to_plot if math.isfinite(val[0])])
                    plt.plot(data_to_plot[:,0], data_to_plot[:,1], label="abstol=reltol={}".format(tol))
                    print("plotted")
                    plot=True
            if plot:
                plt.legend()
                plt.title(model_name)
                plt.xlabel("log10 error in APD90")
                plt.ylabel("pace-to-pace MRMSE".format(measure))
                plt.legend()
                fig = plt.gcf()
                fig.set_size_inches(12, 10)
                plt.savefig("tolerance_plot_{}.pdf".format(model))
                plt.close()

def main():
    measures = ["MRMS"] #"2-Norm", "Trace-MRMS", "Trace-2-Norm"]
    models = ["Tomek2020epi_analytic_voltage", "IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]

    for model in models:
        for measure in measures:
            plot_tolerances(measure, model)

if __name__ == "__main__":
    main()

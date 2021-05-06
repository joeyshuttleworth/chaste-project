#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import adfuller

def plot_moving_average(measure, model_name="", period=1000, block=0, to_plot=True):
    final_apd_error = None
    stopping_pace=None
    minima_window_size=200
    auto_difference_window_size = 50
    window_size = 50
    test_dir = os.path.expanduser("~/Chaste-from-server/testoutput/")
    for dir in os.listdir(os.path.expanduser("~/Chaste-from-server/testoutput/")):
        search_string = "{}_{}ms_{}_percent_block".format(model_name, period, block)
        if re.search(search_string, dir):
            model_name = re.search("([a-z|A-Z|0-9|_]*)_{}ms_{}_percent".format(period,block), dir).group(1)
            print(model_name)
            plot = False
            tol="1e-08"
            file_path = os.path.join(os.path.join(test_dir, dir, "error_measures_{}.dat".format(tol)))
            apd_file_path = os.path.join(os.path.join(test_dir, dir, "apds_using_groundtruth_1e-12.dat"))
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace = True)
                y_vals = np.log10(df[measure])
                apds = df["APD"].values
                # apd_differences=np.log(np.diff(np.abs(apds)))
                true_apd = pd.read_csv(apd_file_path).values[-1,0]
                apd_errors = [np.log10(abs(apd - true_apd)) for apd in apds]
                average_label = "rolling_average_size_{}".format(window_size)
                df["logged_measure"] = np.log10(df[measure])
                df["average"] = df["logged_measure"].rolling(window=window_size).mean()
                df["variance"] = df["logged_measure"].rolling(window=window_size).std()

                minima_windows = np.array_split(df["logged_measure"].values, int(len(df)/minima_window_size))
                minima=[]

                for window in minima_windows:
                    minima.append(min(window))

                windows = df["average"].rolling(window=200)
                adfs  = []

                fig, ax1 = plt.subplots()
                # ax1.plot(y_vals, label=tol)
                ax1.plot(y_vals, label=tol)

                ax2 = ax1.twinx()

                # auto_differences = pd.DataFrame(np.diff(df["average"]), columns=["Auto Differences"])
                auto_differences = pd.DataFrame(np.diff(df["logged_measure"]), columns=["Auto Differences"])
                auto_differences = auto_differences.rolling(window=auto_difference_window_size)

                adfs = [math.nan if abs(val)>5 else val for val in adfs]
                # ax2.plot(adfs, label="ADFuller test", color="red", alpha=0.75)

                test_statistic   = auto_differences.mean()/(auto_differences.std()/np.sqrt(auto_difference_window_size))
                test_statistic   = [math.nan if abs(val)>5 else val for val in test_statistic.values]
                ax2.plot(test_statistic, label="auto-difference test", color="black", alpha=0.75)
                print(len(test_statistic))


                minima_to_consider = 20
                for i in range(minima_window_size*minima_to_consider, len(test_statistic)):
                    # if df["logged_measure"].values[i] > -5:
                    #     continue
                    finish = True
                    if abs(test_statistic[i]) > 0.2:
                            finish = False
                            continue
                    if finish:
                        minima_window_index = int(i / minima_window_size)
                        if minima_window_index >= len(minima):
                            continue
                        vals = [minima[k] for k in range(minima_window_index - minima_to_consider, minima_window_index+1) if math.isfinite(minima[k])]
                        pmcc = pearsonr(range(0,len(vals)), vals)[0]
                        if pmcc < -0:
                            finish = False
                        else:
                            stopping_pace = i
                            ax1.axvline(stopping_pace, label="terminated at")
                            break

                print("plotted")
                plot=True
            if plot and to_plot:
                ax1.legend()
                ax2.legend()
                plt.title("compare_measures_{}_{}_{}ms_{}_percent_block".format(measure,model_name, period, block))
                plt.xlabel("log10 error in APD90")
                ax1.set_ylabel("{} error between sucessive paces".format(measure))
                fig.set_size_inches(16, 10)
                # plt.savefig("compare_measures_{}_{}_{}ms_{}_percent_block".format(measure,model_name, period, block))
                plt.show()
                final_apd_error = apd_errors[stopping_pace]
    return [stopping_pace, final_apd_error]



def main():
    measures = ["MRMS"]
    # models = ["IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]
    periods = [1000, 1250, 750, 500]
    blocks  = [0, 25, 50]
    models = ["Tomek2020epi_analytic_voltage", "IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]
    outputfile = open("stopping_criteria.out", "w")

    plot_moving_average("MRMS", "Tomek2020epi_analytic_voltage", to_plot=True)

    outputfile.write("model period IKrBlock initial_conditions_period initial_conditions_IKrBlock terminal_pace APD90_error_at_terminal_pace\n")
    for period in periods:
        for block in blocks:
                for model in models:
                    # plot_moving_average(measure, model)
                    icperiod = 500 if period == 1000 else 500
                    icblock  = 0.5 if period == 0 else 0.5
                    [stopping_pace, final_apd_error] = plot_moving_average("MRMS", model, period, block, True)
                    print("Stopped at {} with APD90 error {}".format(stopping_pace, final_apd_error))
                    outputfile.write("{} {} {} {} {} {} {}\n".format(model, period, block, icperiod, icblock, stopping_pace, final_apd_error))
    outputfile.close()

if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import adfuller

def plot_moving_average(measure, model_name=""):

    auto_difference_window_size = 200
    window_size = 50
    for dir in os.listdir("testoutput/"):
        if re.search("{}_1000ms_0_percent".format(model_name), dir):
            model_name = re.search("([a-z|A-Z|0-9|_]*)_1000ms_0_percent", dir).group(1)
            plot = False
            tol="1e-08"
            file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "error_measures_{}.dat".format(tol)))
            apd_file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace = True)
                y_vals = np.log10(df[measure])
                print("y_vals length is {}".format(len(y_vals)))
                apds = df["APD"].values
                # apd_differences=np.log(np.diff(np.abs(apds)))
                true_apd = pd.read_csv(apd_file_path).values[-1,0]
                apd_errors = [np.log10(abs(apd - true_apd)) for apd in apds]
                average_label = "rolling_average_size_{}".format(window_size)
                print(average_label)
                df["average"] = np.log10(df[measure].rolling(window=window_size).mean())

                windows = df["average"].rolling(window=200)
                pmccs = []
                adfs  = []
                for window in windows:
                    is_finite = False if sum([1 if not math.isfinite(val) else 0 for val in window]) > 0 else True
                    if len(window)>2 and is_finite:
                        pmcc = pearsonr(range(0, len(window)), window)[0]
                        # adf = adfuller(window, maxlag=10)[0]
                    else:
                        # adf  = math.nan
                        pmcc = math.nan
                    pmccs.append(pmcc)
                    # adfs.append(adf)

                fig, ax1 = plt.subplots()
                # ax1.plot(y_vals, label=tol)
                ax1.plot(y_vals, label=tol)
                ax1.plot(apd_errors,label="apd_errors")
                # ax1.plot(apd_differences, label="apd differences")
                df["average"].plot(label=average_label, ax=ax1)

                ax2 = ax1.twinx()

                # auto_differences = pd.DataFrame(np.diff(df["average"]), columns=["Auto Differences"])
                auto_differences = pd.DataFrame(np.diff(df["average"]), columns=["Auto Differences"])
                auto_differences = auto_differences.rolling(window=auto_difference_window_size)

                test_statistic = pmccs
                ax2.plot(test_statistic, label="rolling pmcc", color="pink", alpha=0.75)

                test_statistic = adfs
                ax2.plot(test_statistic, label="ADFuller test", color="red", alpha=0.75)

                test_statistic   = auto_differences.mean()/(auto_differences.std()/np.sqrt(auto_difference_window_size))
                test_statistic   = [math.nan if abs(val)>5 else val for val in test_statistic.values]
                ax2.plot(test_statistic, label="autodifference test", color="black", alpha=0.75)

                auto_differences = pd.DataFrame(np.diff(df["average"]), columns=["Auto Differences"])
                auto_differences = auto_differences.rolling(window=int(auto_difference_window_size/2))
                test_statistic2   = auto_differences.mean()/(auto_differences.std()/np.sqrt(int(auto_difference_window_size/2)))
                test_statistic2   = [math.nan if abs(val)>5 else val for val in test_statistic2.values]

                # Now find where our criteria would have stopped
                pmccs_to_check = 5
                for i in range(0, len(pmccs)):
                    finish = True
                    if pmccs[i] > -0.2 and i + 5 < df.size:
                        for j in range(i+1, i + pmccs_to_check):
                            if pmccs[j] < -0.75:
                                finish=False
                                break
                        if not finish:
                            break
                        else:
                            # Check the test_statistic
                            if test_statistic[i] < -0.8:
                                finish = False
                            elif test_statistic2[i] < -0.8:
                                finish = False
                                # Check auto difference test with smaller window size
                    else:
                        finish=False

                    if finish:
                        print("finished at pace {}".format(i))
                        ax1.axvline(i, label="terminated at")
                        break

                print("plotted")
                plot=True
            if plot:
                ax1.legend()
                ax2.legend()
                plt.title(model_name)
                plt.xlabel("log10 error in APD90")
                ax1.set_ylabel("{} error between sucessive paces".format(measure))
                ax2.set_ylabel("test statistic with window_size {}".format(auto_difference_window_size))
                fig.set_size_inches(16, 10)
                # plt.savefig("compare_measures_{}_{}".format(measure,model_name))
                plt.show()


def main():
    measures = ["MRMS"]
    models = ["Tomek2020epi_analytic_voltage", "IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]
    # models =  ["IyerMazhariWinslow2004_analytic_voltage"]
    for measure in measures:
        for model in models:
            plot_moving_average(measure, model)

if __name__ == "__main__":
    main()

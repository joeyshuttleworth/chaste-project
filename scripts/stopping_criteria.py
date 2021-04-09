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
    stopping_pace=0
    auto_difference_window_size = 200
    window_size = 50
    for dir in os.listdir("testoutput/"):
        if re.search("{}_{}ms_{}_percent".format(model_name, period, block), dir):
            model_name = re.search("([a-z|A-Z|0-9|_]*)_{}ms_{}_percent".format(period,block), dir).group(1)
            print(model_name)
            plot = False
            tol="1e-08"
            file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "error_measures_{}.dat".format(tol)))
            apd_file_path = os.path.join(os.path.join(os.getcwd(), "testoutput", dir, "apds_using_groundtruth_1e-12.dat"))
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

                minimum_window_size = 50
                minimums = df["logged_measure"].rolling(window=minimum_window_size).min()

                minima_windows = np.split(df["logged_measure"].values, int(len(df)/minimum_window_size))
                minima=[]

                for window in minima_windows:
                    minima.append(min(window))

                windows = df["average"].rolling(window=200)
                adfs  = []
                mean_pmccs = []
                for window in windows:
                    is_finite = False if sum([1 if not math.isfinite(val) else 0 for val in window]) > 0 else True
                    if len(window)>2 and is_finite:
                        pmcc = pearsonr(range(0, len(window)), window)[0]
                        # adf = -adfuller(window, regression="nc")[0]
                    else:
                        # adf  = math.nan
                        pmcc = math.nan
                    mean_pmccs.append(pmcc)
                    # adfs.append(adf)

                var_pmccs = []
                windows = df["variance"].rolling(window=100)
                for window in windows:
                    is_finite = False if sum([1 if not math.isfinite(val) else 0 for val in window]) > 0 else True
                    if len(window)>2 and is_finite:
                        pmcc = pearsonr(range(0, len(window)), window)[0]
                    else:
                        pmcc = math.nan
                    var_pmccs.append(pmcc)


                fig, ax1 = plt.subplots()
                # ax1.plot(y_vals, label=tol)
                ax1.plot(y_vals, label=tol)
                ax1.plot(apd_errors,label="apd_errors")
                # ax1.plot(apd_differences, label="apd differences")
                ax1.plot(minimums, label="rolling minimums")
                df["average"].plot(label=average_label, ax=ax1)
                # df["variance"].plot(label="rolling standard deviation", ax=ax1)

                ax2 = ax1.twinx()

                # auto_differences = pd.DataFrame(np.diff(df["average"]), columns=["Auto Differences"])
                auto_differences = pd.DataFrame(np.diff(df["average"]), columns=["Auto Differences"])
                auto_differences = auto_differences.rolling(window=auto_difference_window_size)

                ax2.plot(mean_pmccs, label="pmcc of rolling average", color="pink", alpha=0.75)
                # ax2.plot(var_pmccs, label="pmcc of rolling standard deviation", color="purple", alpha=0.75)


                adfs = [math.nan if abs(val)>5 else val for val in adfs]
                # ax2.plot(adfs, label="ADFuller test", color="red", alpha=0.75)

                test_statistic   = auto_differences.mean()/(auto_differences.std()/np.sqrt(auto_difference_window_size))
                test_statistic   = [math.nan if abs(val)>5 else val for val in test_statistic.values]
                ax2.plot(test_statistic, label="Dickey Fuller Test Statistic", color="black", alpha=0.75)

                auto_differences = pd.DataFrame(np.diff(df["average"]), columns=["Auto Differences"])
                auto_differences = auto_differences.rolling(window=int(auto_difference_window_size/2))
                test_statistic2   = auto_differences.mean()/(auto_differences.std()/np.sqrt(int(auto_difference_window_size/2)))
                test_statistic2   = [math.nan if abs(val)>5 else val for val in test_statistic2.values]

                # Now find where our criteria would have stopped
                pmccs_to_check = 5
                dfs_to_check = 5
                for i in range(minimum_window_size, len(mean_pmccs)):
                    finish = True
                    if df["logged_measure"].values[i] > -5:
                        finish=False
                    elif mean_pmccs[i] > -0.2 and i + pmccs_to_check < len(mean_pmccs):
                        for j in range(i+1, i + pmccs_to_check):
                            if mean_pmccs[j] < -0.5:
                                finish=False
                                break
                        if not finish:
                            break
                        else:
                            # Check the test_statistic
                            for j in range(0, dfs_to_check):
                                if j + dfs_to_check > len(mean_pmccs):
                                    break
                                if test_statistic[i] < -0.2:
                                    finish = False
                                    break
                                # Check auto difference test with smaller window size
                                elif test_statistic2[i] < -0.2:
                                    break
                                    finish = False
                    else:
                        finish=False

                    if finish:
                        minima_window_index = int(i / minimum_window_size)
                        # finally check last five minima
                        if minima_window_index - 10 < 0 or minima_window_index >= len(minima):
                            continue
                        vals = [minima[k] for k in range(minima_window_index-10, minima_window_index+1)]
                        if pearsonr(range(0,len(vals)), vals)[0] < -0.1:
                            finish = False
                        else:
                            stopping_pace = i+max(pmccs_to_check, dfs_to_check)
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
                ax2.set_ylabel("test statistic with window_size {}".format(auto_difference_window_size))
                fig.set_size_inches(16, 10)
                plt.savefig("compare_measures_{}_{}_{}ms_{}_percent_block".format(measure,model_name, period, block))
                # plt.show()

            final_apd_error = apd_errors[stopping_pace]
            return [stopping_pace, final_apd_error]


def main():
    measures = ["MRMS"]
    models = ["IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]
    periods = [1000, 1250, 750, 500]
    blocks  = [0, 25, 50]
    # models = ["Tomek2020epi_analytic_voltage", "IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]
    outputfile = open("stopping_criteria.out", "w")

    plot_moving_average("MRMS", "decker_2009_analytic_voltage", to_plot=True)

    outputfile.write("model period IKrBlock initial_conditions_period initial_conditions_IKrBlock terminal_pace APD90_error_at_terminal_pace\n")
    for period in periods:
        for block in blocks:
                for model in models:
                    # plot_moving_average(measure, model)
                    icperiod = 500 if period == 1000 else 500
                    icblock  = 0.5 if period == 0 else 0.5
                    [stopping_pace, final_apd_error] = plot_moving_average("MRMS", model, period, block, True)
                    print("Stopped at {} with APD90 error {}".format(stopping_pace, 10**final_apd_error))
                    outputfile.write("{} {} {} {} {} {} {}\n".format(model, period, block, icperiod, icblock, stopping_pace, 10**final_apd_error))
    outputfile.close()

if __name__ == "__main__":
    main()

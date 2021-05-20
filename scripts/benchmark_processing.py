#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import seaborn as sns

def summarise_model(model_name):
    print(model_name)
    file_path = os.path.join(os.path.join(os.path.expanduser("~"), "Chaste-from-server", "testoutput", "TestBenchmark", "{}_results.dat".format(model_name)))

    df_results = pd.read_csv(file_path, delim_whitespace=True)
    df_results = df_results[df_results.model_name==model_name]

    extrapolation_constants = np.unique(df_results["extrapolation_constant"].values)
    extrapolation_constants = [float(val) for val in extrapolation_constants]

    buffer_sizes = np.unique(df_results["buffer_size"].values)
    buffer_sizes = [float(val) for val in buffer_sizes]

    periods = np.unique(df_results["period"].values)
    periods = [float(val) for val in periods]

    blocks = np.unique(df_results["IKrBlock"].values)
    blocks = [float(val) for val in blocks]

    total_scores = []

    # Get reference rows
    reference_df=df_results[df_results.extrapolation_constant==0]
    reference_df=reference_df[reference_df.buffer_size==100]
    reference_df=reference_df[["score", "ic_period", "ic_block", "period", "IKrBlock"]]

    # Add fraction paces saved column
    reference_df=reference_df.set_index(["ic_period", "ic_block", "period", "IKrBlock"])
    normalised_paces_saved_list=[]
    # print scores for each set of parameters
    for e_c in extrapolation_constants:
        for buffer_size in buffer_sizes:
            rows = df_results[df_results.extrapolation_constant == e_c]
            rows = rows[rows.buffer_size == buffer_size]

            rows = rows.set_index(["ic_period", "ic_block", "period", "IKrBlock"])

            scores = rows["score"]

            paces_saved = (reference_df- rows)["score"]
            normalised_paces_saved = ((reference_df- rows)/reference_df)["score"]
            normalised_paces_saved_list.append([buffer_size, e_c, normalised_paces_saved])

            number_of_scenarios = len(scores)
            total_score = sum(scores)
            jumps_used = sum(rows["jumps_used"])
            total_scores.append([e_c, int(buffer_size), jumps_used, total_score, normalised_paces_saved.quantile(0), normalised_paces_saved.quantile(0.25), normalised_paces_saved.quantile(0.5), normalised_paces_saved.quantile(0.75), normalised_paces_saved.quantile(1), paces_saved.mean()])

    total_scores = pd.DataFrame(total_scores, columns=["lmbda", "n", "jumps used", "paces_saved", "min", "lq", "median", "uq", "max", "mean"])

    total_scores = total_scores[total_scores.paces_saved != 0]

    reference_score = total_scores[total_scores.lmbda==0]
    reference_score = reference_score[reference_score.n==100]["paces_saved"].values[0]
    print("reference_score {}".format(reference_score))

    total_scores["percentage_paces_saved"] = (100*(reference_score - total_scores["paces_saved"])/reference_score)

    total_scores["paces_saved"] = reference_score - total_scores["paces_saved"]

    total_scores.sort_values(by=['paces_saved'], inplace=True, ascending=False)
    # plot boxplots for the top 3 scores
    scores_to_plot=total_scores.iloc[0:5]

    df_to_plot=[]
    column_names=[]
    for [index, settings] in scores_to_plot.iterrows():
        lmbda = settings["lmbda"]
        n    = settings["n"]
        # setting_results = df_results[df_results.buffer_size==n]
        # setting_results = setting_results[setting_results.extrapolation_constant==lmda]["normalised_paces_saved"].values

        setting_results=next(filter(lambda x: x[0]==n and x[1] == lmbda, normalised_paces_saved_list))[2]

        column_names.append("({}, {})".format(int(n),lmbda))
        df_to_plot.append(setting_results)

    df_to_plot=np.stack(df_to_plot).transpose()
    df_to_plot = pd.DataFrame(df_to_plot, columns=column_names)



    # make boxplot with Seaborn
    bplot=sns.boxplot(
                 data=df_to_plot,
                 width=0.5,
                 palette="colorblind")

    # add stripplot to boxplot with Seaborn
    bplot=sns.stripplot(data=df_to_plot,
                   jitter=True,
                   marker='o',
                   alpha=0.5,
                   color='black')
    plt.ylabel("proportion of paces saved")
    plt.xlabel("(n, λ)")
    plt.ylim(0,1)
    plt.savefig("{}_boxplots.pdf".format(model_name))

    max_final_mrms=0
    max_apd_range=0
    max_reference_mrms=0
    # Check all of the APDs are close
    for block in blocks:
        for period in periods:
            # max_final_mrms=0
            # max_apd_range=0
            # max_reference_mrms=0

            rows = df_results[df_results.period==period]
            rows = rows[rows.IKrBlock==block]
            rows = rows[rows.APD90.notnull()]
            apd_vals = rows["APD90"].values
            apd_min = apd_vals.min()
            apd_max = apd_vals.max()

            print("block={}, period={}, APD range was {}, min paces taken = {}".format(block, period, apd_max-apd_min, rows["score"].min()))

            new_max_apd_range=abs(apd_max - apd_min)

            if new_max_apd_range > max_apd_range:
                min_row = rows[rows.APD90==apd_min]
                max_row = rows[rows.APD90==apd_max]
                max_apd_range=new_max_apd_range

            mrms_errors = rows["last_mrms"].values
            max_final_mrms = max([max(mrms_errors), max_final_mrms])
            reference_errors = rows["reference_mrms"].values
            max_reference_mrms = max([max(reference_errors), max_reference_mrms])

    print("biggest final mrms was {}".format(max_final_mrms))
    print("biggest reference mrms was {}".format(max_reference_mrms))
    print("max apd range was {}".format(max_apd_range))
    if(max_apd_range > 0.5):
        print("Bifurcation detected:")
        print(min_row[["ic_period", "ic_block", "period", "IKrBlock", "extrapolation_constant", "buffer_size", "APD90"]])
        print(max_row[["ic_period", "ic_block", "period", "IKrBlock", "extrapolation_constant", "buffer_size", "APD90"]])


    # return the total_scores data_frame

    # Print a latex table to file
    with open("{}_scores_table.tex".format(model_name), "w") as f:
        table = total_scores[["n", "lmbda", "percentage_paces_saved", "jumps used", "min", "lq", "median", "uq", "max"]].iloc[0:10].to_latex(index=False, float_format="%.2f")
        f.write(table)
        f.write("\\biggest final mrms was {}\\".format(max_final_mrms))
        f.write("max apd range was {}\\".format(max_apd_range))

    return([df_to_plot.iloc[:,0], total_scores.iloc[0]])

if __name__=="__main__":
    # models = ["Tomek2020epi_analytic_voltage", "IyerMazhariWinslow2004_analytic_voltage", "HundRudy2004_analytic_voltage_units", "decker_2009_analytic_voltage", "tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]
    models = ["tentusscher_model_2004_epi_analytic_voltage", "tentusscher_model_2006_epi_analytic_voltage", "ohara_rudy_2011_epi_analytic_voltage", "ohara_rudy_cipa_2017_analytic_voltage"]
    # models = ["ohara_rudy_2011_epi_analytic_voltage"]
    # models = ["decker_2009_analytic_voltage"]

    total_scores = []
    boxplot_data = []
    for model in models:
        normalised_data_points, total_score = summarise_model(model)
        # Get best score
        total_scores.append(total_score)
        boxplot_data.append(normalised_data_points)
    # Construct new data frame with best results from each model to compare

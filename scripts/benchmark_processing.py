#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def summarise_model(model_name):
    file_path = os.path.join(os.path.join(os.path.expanduser("~"), "Chaste-from-server", "testoutput", "TestBenchmark", "{}_results.dat".format(model_name)))

    print("filepath is ".format(file_path))

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
    print(df_results.columns)
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

            scores = rows["score"].values

            paces_saved = (reference_df- rows)["score"]
            normalised_paces_saved = ((reference_df- rows)/reference_df)["score"]
            normalised_paces_saved_list.append([buffer_size, e_c, normalised_paces_saved])

            number_of_scenarios = len(scores)
            total_score = sum(scores)
            jumps_used = sum(rows["jumps_used"])
            total_scores.append([e_c, int(buffer_size), jumps_used, total_score, total_score/sum(reference_df["score"]), normalised_paces_saved.quantile(0), normalised_paces_saved.quantile(0.25), normalised_paces_saved.quantile(0.5), normalised_paces_saved.quantile(0.75), normalised_paces_saved.quantile(1), paces_saved.mean()])

    total_scores = pd.DataFrame(total_scores, columns=["lmbda", "n", "jumps used", "paces_saved", "normalised_paces_saved", "min", "lq", "median", "uq", "max", "mean"])

    total_scores = total_scores[total_scores.paces_saved != 0]

    reference_score = total_scores[total_scores.lmbda==0]
    reference_score = reference_score[reference_score.n==100]["paces_saved"].values[0]
    print("reference_score {}".format(reference_score))

    total_scores["percentage_paces_saved"] = (100*(reference_score - total_scores["paces_saved"])/reference_score).astype(int)

    total_scores["paces_saved"] = reference_score - total_scores["paces_saved"]

    total_scores.sort_values(by=['paces_saved'], inplace=True, ascending=False)
    # plot boxplots for the top 3 scores
    scores_to_plot=total_scores.iloc[0:8]

    df_to_plot=[]
    column_names=[]
    for [index, settings] in scores_to_plot.iterrows():
        lmbda = settings["lmbda"]
        n    = settings["n"]
        # setting_results = df_results[df_results.buffer_size==n]
        # setting_results = setting_results[setting_results.extrapolation_constant==lmda]["normalised_paces_saved"].values

        setting_results=next(filter(lambda x: x[0]==n and x[1] == lmbda, normalised_paces_saved_list))[2]

        if(len(setting_results)!=36):
            print(len(setting_results), n, lmbda)
        column_names.append("n={}, lamda={}".format(n,lmbda))
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

    plt.show()

    max_final_mrms=0
    max_apd_range=0
    # Check all of the APDs are close
    for block in blocks:
        for period in periods:
            rows = df_results[df_results.period==period]
            rows = rows[rows.IKrBlock==block]
            apd_vals = rows["APD90"].values
            apd_min = min([min(apd_vals)])
            apd_max = max([max(apd_vals)])

            max_apd_range=max([max_apd_range, abs(apd_max - apd_min)])

            mrms_errors = rows["last_mrms"].values
            max_final_mrms = max([max(mrms_errors), max_final_mrms])

    print("biggest final mrms was {}".format(max_final_mrms))
    print("max apd range was {}".format(max_apd_range))

if __name__=="__main__":
    summarise_model("tentusscher_model_2006_epi_analytic_voltage")
    # summarise_model("tentusscher_model_2006_epi_analytic_voltage")
    # summarise_model("tentusscher_model_2006_epi_analytic_voltage")

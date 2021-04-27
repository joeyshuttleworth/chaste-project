import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

    # print scores for each set of parameters
    for e_c in extrapolation_constants:
        for buffer_size in buffer_sizes:
            rows = df_results[df_results.extrapolation_constant == e_c]
            rows = rows[rows.buffer_size == buffer_size]
            scores = rows["score"].values
            number_of_scenarios = len(scores)
            total_score = sum(scores)
            jumps_used = sum(rows["jumps_used"])
            total_scores.append([e_c, int(buffer_size), jumps_used, total_score])

    total_scores = pd.DataFrame(total_scores, columns=["lmbda", "n", "jumps used", "paces_saved"])

    total_scores = total_scores[total_scores.paces_saved != 0]
    # Find reference result - this was given by n=100, lambda=0
    row=total_scores[total_scores.lmbda==0]
    row=row[total_scores.n==100]
    reference_score = row["paces_saved"].values[0]
    print("reference_score {}".format(reference_score))

    total_scores["paces_saved"] = reference_score - total_scores["paces_saved"]

    total_scores.sort_values(by=['paces_saved'], inplace=True)
    print(total_scores)

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

    print("biggest final mrms was ".format(max_final_mrms))
    print("max apd range was {}".format(max_apd_range))

if __name__=="__main__":
    summarise_model("tentusscher_model_2006_epi_analytic_voltage")
    # summarise_model("tentusscher_model_2006_epi_analytic_voltage")
    # summarise_model("tentusscher_model_2006_epi_analytic_voltage")

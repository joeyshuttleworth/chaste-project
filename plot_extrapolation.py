#!/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style

matplotlib.style.use('classic')

def make_plot(model_name):
    numeric_dir = "~/extrapolation_data_numeric/"
    dir = "~/extrapolation_data/"


    dir = dir + model_name + "/TestExtrapolation/"
    numeric_dir = numeric_dir + model_name + "/TestExtrapolation/"

    smart_data = pd.read_csv(dir + "smart.dat", delim_whitespace=True)
    brute_data = pd.read_csv(dir + "bruteforce.dat", delim_whitespace=True)
    numeric_brute_data = pd.read_csv(numeric_dir + "bruteforce.dat", delim_whitespace=True)

    brute_paces = brute_data['pace']
    smart_paces = smart_data['pace']

    for state_var in smart_data:
        smart_var = smart_data[state_var]
        brute_var = brute_data[state_var]

        plt.plot(smart_paces, smart_var)
        plt.plot(brute_paces, brute_var)
        # plt.axhline(smart_var.values[-1], c="burlywood", linestyle='--')
        # plt.axhline(brute_var.values[-1], c="burlywood", linestyle='--')
        plt.axhline(numeric_brute_data[state_var].values[-1], linestyle ='--')
        # plt.legend(["Extrapolation Method", "Brute Force Method", "Terminal value from extrapolation method", "Terminal value from brute force method"])
        plt.title(model_name + " " + state_var)
        plt.show()

if __name__=="__main__":
    make_plot("ohara_rudy_cipa_v1_2017")

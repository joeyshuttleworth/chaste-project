#!/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style
import csv
import os.path
import numpy as np

matplotlib.style.use('classic')

def make_plot(model_name):
    numeric_dir = "~/Chaste/data"
    dir = "~/Chaste/data/"

    with open(os.path.expanduser('~/Chaste/data/') + model_name + "_analytic/500JumpParameters.dat") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        lines = [row for row in csv_reader]

    jump_pace = int(lines[0][0])
    buffer_size = int(lines[0][1])
    extrapolation_coefficient = float(lines[0][2])

    extrap_parameters = [[float(cell) for cell in line] for line in lines[1:]]

    dir = dir + model_name + "/TestExtrapolation/"
    numeric_dir = numeric_dir + model_name + "/TestExtrapolation/"

    smart_data = pd.read_csv(dir + "smart.dat", delim_whitespace=True)
    brute_data = pd.read_csv(dir + "bruteforce.dat", delim_whitespace=True)

    # brute_paces = brute_data['pace']
    # smart_paces = smart_data['pace']

    line_no = 0
    for state_var in smart_data.drop(['pace'], axis=1):
        print("plotting extrapolation of {}".format(state_var))
        if line_no == len(extrap_parameters):
            break
        smart_var = smart_data[state_var]
        brute_var = brute_data[state_var]

        state_no = int(extrap_parameters[line_no][0])
        tau = -1/extrap_parameters[line_no][1]
        beta = -1/tau
        alpha = extrap_parameters[line_no][2]

        V_diff = np.exp(alpha  - buffer_size/tau + 1/tau)/(np.exp(1/tau) - 1);

        if smart_var.values[jump_pace-1] - smart_var.values[jump_pace - buffer_size-1] < 0:
            V_diff = -V_diff;

        V = brute_var[jump_pace-1] + V_diff*(1 - np.exp(-(np.array(range(-jump_pace+1, len(brute_var) - jump_pace + 1))/tau)))

        plt.plot(smart_var)
        plt.plot(brute_var)
        plt.plot(range(0 ,len(brute_data)), V, linestyle = '--', color="orange")
        plt.axhline(brute_var.values[-1], linestyle ='--')
        plt.axhline(smart_var[jump_pace-2] + V_diff, color="grey", linestyle = "-.")
        plt.plot(jump_pace, smart_var[jump_pace-1] + V_diff, label="point", markersize=5, color='red', marker='x')

        bbox = []
        loc  = []
        if V_diff > 0:
            bbox=(1,0)
            loc="lower right"
        else:
            bbox=(1,1)
            loc="upper right"


        plt.legend(["extrapolation method solution", "brute force solution", "fitted curve", "terminal value from brute force method on numeric model", "predicted terminal value from fitted curve", "point extrapolated to"], prop={'size': 5}, bbox_to_anchor=bbox, loc=loc)
        plt.title(model_name + " " + state_var)
        if not os.path.exists("data"):
            os.mkdir("data")
        plt.savefig("data/"+model_name+"_"+state_var+".pdf", format="pdf")
        plt.clf()
        line_no += 1

if __name__=="__main__":
    make_plot("ohara_rudy_cipa_v1_2017")
    make_plot("tentusscher_model_2006_epi")

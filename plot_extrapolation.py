#!/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style
import csv
import os.path
import numpy as np

# matplotlib.style.use('classic')

def make_plot(model_name, jump_number=0):
    print("Plotting extrapolation for {}".format(model_name))
    dir = "~/Chaste/testoutput/"

    with open(os.path.expanduser(dir) + model_name + "/TestExtrapolationMethod/750JumpParameters{}.dat".format(jump_number)) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        lines = [row for row in csv_reader]

    jump_pace = int(lines[0][0])
    buffer_size = int(lines[0][1])
    extrapolation_coefficient = float(lines[0][2])

    extrap_parameters = [[float(cell) for cell in line[2:]] for line in lines[1:]]

    dir = dir + model_name + "/TestExtrapolationMethod/"

    smart_data = pd.read_csv(dir + "smart.dat", delim_whitespace=True)
    brute_data = pd.read_csv(dir + "bruteforce.dat", delim_whitespace=True)

    for i in range(1, len(lines)):
        var_index = int(lines[i][0])
        var_name  = lines[i][1]
        smart_var = smart_data[var_name]
        brute_var = brute_data[var_name]

        state_no = int(lines[i][0])
        tau = -1/extrap_parameters[i-1][0]
        beta = -1/tau
        alpha = extrap_parameters[i-1][1]
        vinf = extrap_parameters[i-1][2]

        print("plotting extrapolation of {}".format(var_name))
        print(alpha, tau, vinf)

        V_diff = np.exp(alpha  - buffer_size/tau + 1/tau)/(np.exp(1/tau) - 1);

        if smart_var.values[jump_pace-1] - smart_var.values[jump_pace - buffer_size-1] < 0:
            V_diff = -V_diff;

        V = brute_var[jump_pace-1] + V_diff*(1 - np.exp(-(np.array(range(-jump_pace+1, len(brute_var) - jump_pace + 1))/tau)))

        # Predict the value at the start of our buffer and see how far off this value is
        v0 = V[jump_pace-buffer_size-1]

        V_total_difference = np.exp(alpha)/(1 - np.exp(-1/tau))

        check_val = np.abs((brute_var[jump_pace-buffer_size-1] - v0)/ V_total_difference)
        predicted_val = 1 - np.exp(-buffer_size/tau)
        print("check_val, predicted_val = {}, {}".format(check_val, predicted_val))

        plt.plot(smart_var)
        plt.plot(brute_var)
        plt.plot(range(0 ,len(brute_data)), V, linestyle = '--', color="orange")
        plt.axhline(brute_var.values[-1], linestyle ='--')
        plt.axhline(vinf, color="grey", linestyle = "-.")
        plt.plot(jump_pace, vinf, label="point", markersize=5, color='red', marker='x')
        plt.plot(jump_pace - buffer_size, v0, markersize=5, color="purple", marker="x")

        buffer_range = range(jump_pace-buffer_size-1, jump_pace-1)
        plt.plot(buffer_range, smart_var[buffer_range] ,marker="x", color="green", linestyle=None)

        bbox = []
        loc  = []
        if V_diff > 0:
            bbox=(1,0)
            loc="lower right"
        else:
            bbox=(1,1)
            loc="upper right"


        plt.legend(["extrapolation method solution", "brute force solution", "fitted curve", "terminal value from brute force method on numeric model", "predicted terminal value from fitted curve", "point extrapolated to", "predicted start of buffer", "Points used for extrapolation"], prop={'size': 10}, bbox_to_anchor=bbox, loc=loc)
        plt.title(model_name + " " + var_name)
        if not os.path.exists("plots"):
            os.mkdir("plots")
        plt.show()
        # plt.savefig("plots/"+model_name+"_"+var_name+".pdf", format="pdf")

        plt.clf()

def compare_p2p_errors(model_name):
    dir = "~/Chaste/testoutput/"
    dir = dir + model_name + "/TestExtrapolationMethod/"
    smart_df = pd.read_csv(dir + "smart.dat", delim_whitespace=True)
    smart_df = smart_df.drop("pace", 1)

    brute_df = pd.read_csv(dir + "bruteforce.dat", delim_whitespace=True)
    brute_df = brute_df.drop("pace", 1)
    np.log10(smart_df['mrms']).plot(label="extrapolation method", title = "{} pace to pace mixed root mean square error comparison".format(model_name))
    np.log10(brute_df['mrms']).plot(label="brute force method")
    plt.show()

    # for column in df:
    #     df[column].plot(title=column)
    #     plt.show(


if __name__=="__main__":
    # compare_p2p_errors("IyerMazhariWinslow2004")
    # make_plot("IyerMazhariWinslow2004")
    compare_p2p_errors("HundRudy2004_units")
    make_plot("HundRudy2004_units", 2)
    # compare_p2p_errors("tentusscher_model_2004_epi")
    # make_plot("tentusscher_model_2004_epi")
    # make_plot("ohara_rudy_2011_epi")
    # compare_p2p_errors("Tomek2020epi")
    # make_plot("Tomek2020epi", 0)
    # compare_p2p_errors("decker_2009")
    # make_plot("decker_2009", 0)

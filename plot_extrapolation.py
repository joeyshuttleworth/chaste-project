#!/bin/env python3
"""Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style
import csv
import os.path
import numpy as np

matplotlib.style.use('classic')

def make_plot(model_name):
    numeric_dir = "~/extrapolation_data_numeric/"
    dir = "~/extrapolation_data/"

    with open(os.path.expanduser('~/extrapolation_data/') + model_name + "/500JumpParameters.dat") as csv_file:
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
    numeric_brute_data = pd.read_csv(numeric_dir + "bruteforce.dat", delim_whitespace=True)

    brute_paces = brute_data['pace']
    smart_paces = smart_data['pace']

    line_no = 0
    for state_var in smart_data.drop(['pace'], axis=1):
        if line_no == len(extrap_parameters):
            break
        smart_var = smart_data[state_var]
        brute_var = brute_data[state_var]

        state_no = int(extrap_parameters[line_no][0])
        tau = -1/extrap_parameters[line_no][1]
        beta = -1/tau
        alpha = extrap_parameters[line_no][2]

        V_diff = np.exp(alpha  - buffer_size/tau)/(1 - np.exp(-1/tau));


        if smart_var.values[jump_pace-1] - smart_var.values[jump_pace - buffer_size-1] < 0:
            V_diff = -V_diff;

        V = brute_var[jump_pace-1] + V_diff*(1 - np.exp(-(np.array(range(-jump_pace+1, len(brute_var) - jump_pace + 1))/tau)))

        plt.plot(smart_paces, smart_var)
        plt.plot(brute_paces, brute_var)
        plt.plot(range(0 ,len(brute_data)), V, linestyle = '--', color="orange")
        plt.axhline(numeric_brute_data[state_var].values[-1], linestyle ='--')
        plt.axhline(smart_var[jump_pace-2] + V_diff, color="grey", linestyle = "-.")

        plt.plot(jump_pace, smart_var[jump_pace-1] + V_diff, label="point", markersize=5, color='red', marker='x')
        plt.legend(["extrapolation method solution", "brute force solution", "fitted curve", "terminal value from brute force method on numeric model", "predicted terminal value from fitted curve", "point extrapolated to"], prop={'size': 5})
        plt.title(model_name + " " + state_var)
        if not os.path.exists("data"):
            os.mkdir("data")
        plt.savefig("data/"+model_name+"_"+state_var+".eps", format="eps")
        plt.clf()
        line_no += 1

if __name__=="__main__":
    make_plot("ohara_rudy_cipa_v1_2017")
    make_plot("tentusscher_model_2006_epi")

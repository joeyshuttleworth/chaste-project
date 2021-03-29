#/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import os
import re

def main():
    os.chdir("TestAlgebraicVoltage")
    model_list = []
    for filename in os.listdir():
        split_filename = os.path.splitext(filename)
        if split_filename[1]==".dat" and re.search("block_final",split_filename[0]):
            print(filename)
            model_name = re.search("([a-z|A-Z|0-9|_]*)_final", split_filename[0]).group(1)
            print(filename)
            model_list.append(filename[0])
            # plot final traces
            df = pd.read_csv(filename, delim_whitespace=True)
            print(df)
            if "membrane_voltage(mV)" in df.columns:
                df["membrane_voltage(mV)"].plot(label="original_model")
            else:
                df["membrane_voltage(millivolt)"].plot(label="original_model")

            df = pd.read_csv(model_name+"_analytic_final.dat", delim_whitespace=True)
            if "membrane_voltage(mV)" in df.columns:
                df["membrane_voltage(mV)"].plot(label="analytic_model")
            else:
                df["membrane_voltage(millivolt)"].plot(label="analytic_model")

            plt.title("{} final pace".format(model_name))
            plt.show()
            plt.legend()
    # for model in model_list:
if __name__ == "__main__":
    main()

import pandas as pd
import numpy as np

# Calculate the mrms between pandas.DataFrames x and y (x is in the denominator)
def mrms(X, Y):
    N = len(X.index)
    if N != len(Y.index):
        return NaN
    else:
        list_of_summands = [int((x - y)**2 / (1 + np.abs(x))**2) for x,y in zip(X.values,Y.values)]
        return(np.sqrt(sum(list_of_summands)/N))

# Test this works
if __name__ == "__main__":
    x = pd.DataFrame(range(0,100))
    y = pd.DataFrame(range(100,200))
    print(mrms(x,y))



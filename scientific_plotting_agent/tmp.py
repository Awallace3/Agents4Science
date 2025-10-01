import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pprint import pprint as pp
import qcelemental as qcel

def isolate_curve():
    df = pd.read_pickle("./combined_df_4569.pkl")
    pp(df.columns.to_list())
    return


def main():
    isolate_curve()
    return


if __name__ == "__main__":
    main()

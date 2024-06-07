from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import scale
from scipy.stats import bootstrap
import numpy as np


def normalize(x_list):
    x_scaled = scale(x_list)
    return x_scaled


def bootstrapping(x_list):
    res = bootstrap((x_list,), np.median)
    standard_error = res.standard_error
    median = np.median(res.bootstrap_distribution)
    return [median, standard_error]

"""Performing Bias calculations."""

from typing import List
import numpy as np
import pylab as plt
import ccdproc
from astropy.nddata import CCDData
import funcs
from exposure_type import ExpType


def create_master_bias(biases: List[CCDData]) -> CCDData:
    """Creates the master bias."""
    #funcs.oscan_and_trim(biases)
    biases = ccdproc.Combiner(biases)
    master_bias = biases.average_combine()
    return master_bias

def do_biases(bias_directory: str):
    """Performs all the steps to create a master bias."""
    bias_list = funcs.read_in_file_type(bias_directory, ExpType.bias)
    return create_master_bias(bias_list)

if __name__ == '__main__':
    TEST_DIRECTORY = '/home/tlambert/Downloads/g_band/BIAS'
    bias = do_biases(TEST_DIRECTORY)
    imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
    bias_min, bias_max, bias_mean, bias_std = imstats(np.asarray(bias))
    plt.figure(figsize=(15, 15))
    plt.imshow(bias, vmax=bias_mean + 4*bias_std, vmin=bias_mean - 4*bias_std)
    plt.show()

"""Performing Bias calculations."""

import numpy as np
import pylab as plt
import ccdproc
from ccdproc import ImageFileCollection
from astropy.nddata import CCDData
from typing import List
import funcs


def read_in_biases(bias_directory: str) -> List[CCDData]:
    """Returns all the biases in a given directory."""
    images = ImageFileCollection(bias_directory, keywords='*')
    bias_list = []
    for hdu, fname in images.hdus(exptype='Bias', return_fname=True):
        meta = hdu.header
        meta['filename'] = fname
        bias_list.append(ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu"))
    return bias_list

def create_master_bias(biases: List[CCDData]) -> CCDData:
    """Creates the master bias."""
    funcs.oscan_and_trim(biases)
    biases = ccdproc.Combiner(biases)
    master_bias = biases.average_combine()
    return master_bias

def do_biases(bias_directory: str):
    """Performs all the steps to create a master bias."""
    bias_list = read_in_biases(bias_directory)
    return create_master_bias(bias_list)


if __name__ == '__main__':
    TEST_DIRECTORY = '/home/tlambert/Downloads/g_band/BIAS'
    master_bias = do_biases(TEST_DIRECTORY)
    imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
    bias_min, bias_max, bias_mean, bias_std = imstats(np.asarray(master_bias))
    plt.figure(figsize=(15, 15))
    plt.imshow(master_bias, vmax=bias_mean + 4*bias_std, vmin=bias_mean - 4*bias_std)
    plt.show()

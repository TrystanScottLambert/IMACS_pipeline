"""Flat reduction."""

from typing import List, Tuple
import pylab as plt
import numpy as np
import ccdproc
import astropy.units as u
from astropy.nddata import CCDData
import funcs
from exposure_type import ExpType
from bias_reduction import do_biases

GAIN = 0.83 * u.electron / u.adu

def make_master_flat(flat_list: List[CCDData], master_bias: CCDData) -> Tuple[CCDData, CCDData]:
    """Goes through list of flats and corrects them for bias."""
    for flat in flat_list:
        flat = ccdproc.subtract_bias(flat, master_bias, add_keyword={'calib': 'subtracted bias'})
    #funcs.oscan_and_trim(flats)
    flat_combiner = ccdproc.Combiner(flat_list)
    scaling_func = lambda arr: 1/np.ma.average(arr)
    flat_combiner.scaling = scaling_func
    master_flat = flat_combiner.median_combine(median_func=funcs.bn_median)
    master_flat.header = flat_list[0].meta
    master_flat_electron = ccdproc.gain_correct(master_flat, gain=GAIN)
    return master_flat, master_flat_electron


def do_flats(flat_directory: str, master_bias: CCDData) ->  Tuple[CCDData, CCDData]:
    """performs the flat reduction."""
    flats = funcs.read_in_file_type(flat_directory, ExpType.flat)
    return make_master_flat(flats, master_bias)

if __name__ == '__main__':
    imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
    FLAT_DIRECTORY = '/home/tlambert/Downloads/g_band/FLAT'
    BIAS_DIRECTORY = '/home/tlambert/Downloads/g_band/BIAS'
    master_bias = do_biases(BIAS_DIRECTORY)
    master_flat, master_flat_electron = do_flats(FLAT_DIRECTORY, master_bias)

    f_min, f_max, f_mean, f_std = imstats(np.asarray(master_flat))
    plt.figure(figsize=(15, 15))
    plt.imshow(master_flat, vmin=f_mean-5*f_std, vmax=f_mean+5*f_std)
    plt.show()

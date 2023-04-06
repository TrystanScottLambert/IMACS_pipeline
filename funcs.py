"""Convience functions for reductions."""

import bottleneck as bn
import numpy as np
from typing import List
from astropy.modeling import models
from astropy.nddata import CCDData
import ccdproc
from ccdproc import ImageFileCollection

def oscan_and_trim(image_list):
    """
    Remove overscan and trim a list of images. The original list is replaced by a list of images
    with the changes applied.
    """
    for idx, img in enumerate(image_list):
        oscan = ccdproc.subtract_overscan(img, img[:, 1025:1088], add_keyword={'oscan_sub': True, 'calstat': 'O'}, model=models.Polynomial1D(1))
        image_list[idx] = ccdproc.trim_image(oscan[:, :1088], add_keyword={'trimmed': True, 'calstat': 'OT'})


def med_over_images(masked_arr, axis=0):
    """
    Calculate median pixel value along specified axis
    
    Uses bottleneck.nanmedian for speed
    """

    dat = masked_arr.data.copy()
    dat[masked_arr.mask] = np.NaN
    return bn.nanmedian(dat, axis=axis)

def bn_median(masked_array, axis=None):
    """
    Perform fast median on masked array
    
    Parameters
    ----------
    
    masked_array : `numpy.ma.masked_array`
        Array of which to find the median.
    
    axis : int, optional
        Axis along which to perform the median. Default is to find the median of
        the flattened array.
    """
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))

def read_in_file_type(type_directory: str, exposure_type: str) -> List[CCDData]:
    """Returns all the biases in a given directory."""
    images = ImageFileCollection(type_directory, keywords='*')
    type_list = []
    for hdu, fname in images.hdus(exptype=exposure_type, return_fname=True):
        meta = hdu.header
        meta['filename'] = fname
        type_list.append(ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu"))
    return type_list

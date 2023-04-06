"""Science reduction."""

import glob
from typing import List
from tqdm import tqdm
from astropy.io import fits
from astropy.nddata import CCDData
import ccdproc
import astropy.units as u
import funcs
from exposure_type import ExpType
from flat_reduction import do_flats, GAIN
from bias_reduction import do_biases

READ_NOISE = 4.90 * u.electron

def calibrate_science(science_directory: str, master_flat_electron:CCDData, master_bias:CCDData) -> List[CCDData]:
    """Performs the Science reduction."""
    science_images = funcs.read_in_file_type(science_directory, ExpType.science)
    #funcs.oscan_and_trim(science_images)
    science_calibrated = []
    for science_image in tqdm(science_images):
        science_bias = ccdproc.subtract_bias(science_image, master_bias, add_keyword={'calib': 'subtracted bias'})
        science_gain = ccdproc.gain_correct(science_bias, gain=GAIN)
        science_flat = ccdproc.flat_correct(science_gain, master_flat_electron)
        science_calibrated.append(science_flat)
    return science_images

def write_list_to_fits(calibrated_list: List[CCDData], outfiles: List[str]) -> None:
    """Writes the calibrated science images as fits files."""
    for i, image in enumerate(calibrated_list):
        image_hdu_list = image.to_hdu()
        image_hdu_list.writeto(outfiles[i])

def do_science(science_directory: str, master_flat_electron: CCDData, master_bias: CCDData) -> None:
    """Perform the full science rediuction"""
    infiles = glob.glob(f'{science_directory}/*.fits')
    infiles.sort()
    outfiles = [file.split('.fits')[0] + '.cal.fits' for file in infiles]
    science_images = calibrate_science(science_directory, master_flat_electron, master_bias)
    write_list_to_fits(science_images, outfiles)


if __name__ == '__main__':
    SCIENCE_DIRECTORY = '/home/tlambert/Downloads/g_band/SCIENCE'
    FLAT_DIRECTORY = '/home/tlambert/Downloads/g_band/FLAT'
    BIAS_DIRECTORY = '/home/tlambert/Downloads/g_band/BIAS'
    master_bias = do_biases(BIAS_DIRECTORY)
    master_flat, master_flat_electron = do_flats(FLAT_DIRECTORY, master_bias)
    do_science(SCIENCE_DIRECTORY, master_flat_electron, master_bias)

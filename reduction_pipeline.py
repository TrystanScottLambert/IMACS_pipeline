"""
IMACS Pipeline
"""

from pathlib import Path
import ccdproc as ccdp
from astropy.nddata import CCDData
from astropy.stats import mad_std
import numpy as np


def inv_median(val):
    """Inverted median for scaling flats."""
    return 1 / np.median(val)

def get_trim_index(ccd_object: CCDData) -> int:
    """Determines the index needed for trimming by reading the biassec in the header."""
    idx =  int(ccd_object.header['biassec'].split('[')[-1].split(':')[0]) -1
    return idx


class ExpType:
    """Exposure types."""
    def __init__(self, input_directory: str, output_directory: str):
        """Directory of the exposure type. Use 'SCIENCE' and not 'SCIENCE/'."""
        self.path = Path(input_directory)
        self.files = ccdp.ImageFileCollection(self.path)
        self.calibrated_data = Path(output_directory)
        self.calibrated_data.mkdir(exist_ok=True)

    def _subtract_overscan(self, ccd_object: CCDData, outfile: str) -> CCDData:
        """Subtracts and trims the overscan."""
        ccd = ccdp.subtract_overscan(
            ccd_object, overscan=ccd_object[:, get_trim_index(ccd_object):])
        ccd = ccdp.trim_image(ccd[:, :get_trim_index(ccd_object)])
        ccd.write(self.calibrated_data / outfile)


class Biases(ExpType):
    """Bias Class"""

    def make_master_bias(self) -> CCDData:
        """Makes the master bias."""
        for ccd, file_name in self.files.ccds(
            EXPTYPE='Bias', ccd_kwargs={'unit': 'adu'}, return_fname=True):
            self._subtract_overscan(ccd, 'bias-' + file_name)

        reduced_biases = ccdp.ImageFileCollection(self.calibrated_data)
        calibrated_biases = reduced_biases.files_filtered(EXPTYPE='Bias', include_path=True)
        combined_bias = ccdp.combine(calibrated_biases,
                             method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e6
                            )
        combined_bias.meta['combined'] = True
        combined_bias.write(self.calibrated_data / 'master_bias.fits')
        return combined_bias


class Flats(ExpType):
    """Flat class"""

    def make_master_flat(self) -> CCDData:
        """Creates the master flat."""
        for ccd, file_name in self.files.ccds(
            EXPTYPE='Flat', ccd_kwargs={'unit': 'adu'}, return_fname=True):
            self._subtract_overscan(ccd, 'flat-' + file_name)

        reduced_images = ccdp.ImageFileCollection(self.calibrated_data)
        calibrated_flats = reduced_images.files_filtered(EXPTYPE='Flat', include_path=True)
        combined_flats = ccdp.combine(calibrated_flats,
                             method='average', scale=inv_median,
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e6
                            )

        combined_flats.meta['combined'] = True
        combined_flats.write(self.calibrated_data / 'master_flat.fits')
        return combined_flats


class Objects(ExpType):
    """Science expsoures."""
    def reduce_science_images(self, master_bias:CCDData, master_flat: CCDData) -> None:
        """Reduces each individual science image in the directory."""
        for ccd, file_name in self.files.ccds(
            ExpType='Object', return_fname=True, ccd_kwargs=dict(unit='adu')):

            reduced = ccdp.subtract_overscan(
                ccd, overscan=ccd[:, get_trim_index(ccd):], median=True)
            reduced = ccdp.trim_image(reduced[:, :get_trim_index(ccd)])
            reduced = ccdp.subtract_bias(reduced, master_bias)
            reduced = ccdp.flat_correct(reduced, master_flat)
            reduced.write(self.calibrated_data / ('science-'+file_name))

def reduce_images(bias_directory: str, flat_directory: str, science_directory: str) -> None:
    """Performs basic reduction on all the science images."""
    biases = Biases(bias_directory, 'reduced')
    flats = Flats(flat_directory, 'reduced')
    objects = Objects(science_directory, 'reduced')

    master_bias = biases.make_master_bias()
    master_flat = flats.make_master_flat()
    objects.reduce_science_images(master_bias, master_flat)


if __name__ == '__main__':
    BIAS_DIRECTORY = 'BIAS'
    FLAT_DIRECTORY = 'FLAT'
    SCIENCE_DIRECTORY = 'SCIENCE'

    reduce_images(BIAS_DIRECTORY, FLAT_DIRECTORY, SCIENCE_DIRECTORY)

""" Script to combine all images in the SCIENCE FOLDER """

import glob
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd


def combine_images(file_pattern: str, output_folder='') -> tuple[np.ndarray, WCS]:
	"""
	Combines all images for the given string pattern.
	Eg. a string pattern of *.fits will combine all fits files
	whilst a string patter of *wcs.fits will combine all files
	ending in 'wcs.fits'. 

	Optinal output folder. If not specified then file will be written in 
	same directory as where the code is run. Alternatively an output folder 
	could be given, i.e. REDUCED/ or /home/user/Desktop/
	"""

	if os.path.isfile(f'{output_folder}coadded_image.fits'):
		raise FileExistsError('coadded_image.fits already exists. Remove, move, or rename before running.')

	files = glob.glob(file_pattern)
	if len(files) == 0:
		raise ValueError('No files have been selected.')
	hdus = [fits.open(file)[0] for file in files]
	wcs_out, shape_out = find_optimal_celestial_wcs(hdus)
	array, _ = reproject_and_coadd(hdus,wcs_out, shape_out=shape_out,reproject_function=reproject_interp)


	header = hdus[0].header
	header = header.update(wcs_out.to_header())


	fits.writeto(f'{output_folder}coadded_image.fits', array, header)
	return wcs_out

if __name__ == '__main__':
	FOLDER = '/home/tlambert/Downloads/g_band/REDUCED/aligned/*science*wcs*.fits'
	wcs_test = combine_images(FOLDER)

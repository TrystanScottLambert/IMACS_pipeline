"""
IMACS Pipeline
"""

from pathlib import Path
import ccdproc as ccdp
from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy import units as u
import numpy as np

BIAS_DIRECTORY = 'BIAS'
FLAT_DIRECTORY = 'FLAT'
SCIENCE_DIRECTORY = 'SCIENCE'


# BIAS #
calibrated_data = Path('.', 'example-reduced')
calibrated_data.mkdir(exist_ok=True)

bias_path = Path(BIAS_DIRECTORY)
bias_files = ccdp.ImageFileCollection(bias_path)

# remove the overscan region and trim.
for ccd, file_name in bias_files.ccds(EXPTYPE='Bias',            # Just get the bias frames
                                 ccd_kwargs={'unit': 'adu'}, # CCDData requires a unit for the image if 
                                                             # it is not in the header
                                 return_fname=True           # Provide the file name too.
                                ):
        # Subtract the overscan
    ccd = ccdp.subtract_overscan(ccd, overscan=ccd[:, 1034:], median=True)
    
    # Trim the overscan
    ccd = ccdp.trim_image(ccd[:, :1024])
    
    # Save the result
    ccd.write(calibrated_data / ('bias-' + file_name))

# create the master bias.

reduced_biases = ccdp.ImageFileCollection(calibrated_data)

calibrated_biases = reduced_biases.files_filtered(EXPTYPE='Bias', include_path=True)

combined_bias = ccdp.combine(calibrated_biases,
                             method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e6
                            )

combined_bias.meta['combined'] = True

combined_bias.write(calibrated_data / 'combined_bias.fits')

############################################################

# Flat 

# Remove the overscan and trim

flat_raw = Path(FLAT_DIRECTORY)
ifc_raw = ccdp.ImageFileCollection(flat_raw)

for ccd, file_name in ifc_raw.ccds(EXPTYPE='Flat',            # Just get the bias frames
                                         ccd_kwargs={'unit': 'adu'}, # CCDData requires a unit for the image if 
                                                                     # it is not in the header
                                         return_fname=True           # Provide the file name too.
                                        ):    
    # Subtract the overscan
    ccd = ccdp.subtract_overscan(ccd, overscan=ccd[:, 1034:], median=True)
    
    # Trim the overscan
    ccd = ccdp.trim_image(ccd[:, :1024])

    # Save the result; there are some duplicate file names so pre-pend "flat"
    ccd.write(calibrated_data / ('flat-' + file_name))

# create the master Flat.

def inv_median(a):
    return 1 / np.median(a)
reduced_images = ccdp.ImageFileCollection(calibrated_data)

calibrated_flats = reduced_images.files_filtered(EXPTYPE='Flat', include_path=True)

combined_flats = ccdp.combine(calibrated_flats,
                             method='average', scale=inv_median,
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e6
                            )

combined_flats.meta['combined'] = True

combined_flats.write(calibrated_data / 'combined_flat.fits')

####################################################################

# Science image reduction

science_raw = ccdp.ImageFileCollection(SCIENCE_DIRECTORY)

# These two lists are created so that we have copies of the raw and calibrated images
# to later in the notebook. They are not ordinarily required.
all_reds = []
science_ccds = []
for light, file_name in science_raw.ccds(EXPTYPE='Object', return_fname=True, ccd_kwargs=dict(unit='adu')):
    science_ccds.append(light)
    
    reduced = ccdp.subtract_overscan(light, overscan=light[:, 1034:], median=True)
    reduced = ccdp.trim_image(reduced[:, :1024])

    reduced = ccdp.subtract_bias(reduced, combined_bias)

    reduced = ccdp.flat_correct(reduced, combined_flats)
    all_reds.append(reduced)

    reduced.write(calibrated_data / ('science-' + file_name))


import os
import montage_wrapper as montage
from glob import glob
from astropy.io import fits

# FITS to combine should be added into this directory
path = "/media/eric/Data_3/NGC1313/HST/"
raw_path = os.path.join(path, "raw_mosaic_fits")
outpath = "/media/eric/Data_3/NGC1313/HST/final_mosaic"

if not os.path.exists(raw_path):
    os.mkdir(raw_path)

# Grab all of the total frames, extract the SCI extension
# into the raw_mosaic_fits folder
totals = glob(os.path.join(path, "*total"))

for fold in totals:
    filename = os.path.join(fold, fold.split("/")[-1]+"_drz.fits")
    hdu = fits.open(filename)

    sci_hdu = hdu["SCI"]
    out_name = os.path.join(raw_path, fold.split("/")[-1]+"_sci.fits")
    sci_hdu.writeto(out_name)

# Create the mosaic
montage.mosaic(raw_path, outpath, background_match=False)

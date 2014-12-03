
'''
Combine all HGBS filament results into a single FITS file.
'''

from astropy.io import fits

import sys
import os

folder = sys.argv[1]

files = [f for f in os.listdir(folder) if os.path.isfile(f)
         and f[-4:] == "fits"
         and f[:3] != "deg"]

header = fits.Header()
header["COMMENT"] = "Contains data on all HGPS filaments at 250 microns."

table_list = [fits.PrimaryHDU(header=header)]
for f in files:
    # Read table in
    hdu = fits.open(f)
    hdu[1].name = f
    table_list.append(hdu[1])

hdulist = fits.HDUList(table_list)
hdulist.writeto("hgps_250_all.fits")

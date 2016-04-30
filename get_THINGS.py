
from bs4 import BeautifulSoup
import urllib2
import re
import os
import shutil

from astroquery.utils import download_list_of_fitsfiles

'''
When you want the THINGS data, but don't want to click on every link
(and wait...)
'''

# base = ("https://bulk.cv.nrao.edu/littlethings/")
base = "http://www.mpia.de/THINGS/"

# Open the page; convert to string
page = urllib2.urlopen(base + "Data.html")
html_string = page.read().decode('utf-8')

fits_reg = '.\FITS'

# Now soup
soup = BeautifulSoup(html_string, 'lxml')

files = soup.findAll(href=re.compile(fits_reg))

full_files = [os.path.join(base, f['href']) for f in files]

output_dir = "/media/eric/Data_3/VLA/THINGS/"

while True:
    try:
        download_list_of_fitsfiles(full_files, output_directory=output_dir,
                                   verbose=True, save=True, overwrite=False)
        break
    except urllib2.URLError:
        pass

# Now let's go through them and sort into a reasonable folder structure
saved_files = [os.path.join(output_dir, f['href'].split("/")[1]) for f in
               files]

for f in saved_files:
    if "_NA_" in f:
        name = f.split("_NA_")[0]
    elif "_RO_" in f:
        name = f.split("_NA_")[0]
    else:
        print(f + "has no NA or RO?? Skipping.")

    folder_name = os.path.join(output_dir, name)

    # If we haven't gotten to this galaxy yet, make a new folder
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)

    # There's an extra _ in front for some reason??
    shutil.move("_" + f, f)

    # Now move that file
    shutil.move(f, folder_name)

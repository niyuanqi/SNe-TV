import glob
import astroalign as aa
from astropy.io import fits

filenames = glob.glob("crop/*.fits")
filenames.sort()

i = 0 #reference image number
hdulist = fits.open(filenames[i])
target_image = hdulist[0].data.byteswap().newbyteorder()

for i in range(len(filenames)):
    hdulist = fits.open(filenames[i])
    source_image = hdulist[0].data.byteswap().newbyteorder()

    try:
        aligned_image, footprint = aa.register(source_image, target_image, detection_sigma=3.0)
        hdulist[0].data = aligned_image.byteswap().newbyteorder()
        hdulist.writeto('align/'+filenames[i].split('/')[-1])
        hdulist.close()

        print(filenames[i], "SUCCESS")
    except:
        hdulist.close()
        print(filenames[i], "FAILED")
        

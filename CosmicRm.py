from astropy.io import fits
import glob
import matplotlib.pyplot as plt
from PIL import Image
import astroscrappy

Bs, Vs, Is = np.loadtxt("rgblist.txt", unpack=True, usecols=(0, 1, 2), dtype=str)

for i in range(len(bs)):
    hdulist = fits.open(Bs[i])
    image = hdulist[0].data
    mask, image = astroscrappy.detect_cosmics(image, verbose=True)
    hdulist[0].data = image
    hdulist.writeto('cosrm/'+filename.split('/')[-1])
    hdulist.close()

    hdulist = fits.open(Vs[i])
    image = hdulist[0].data
    mask, image = astroscrappy.detect_cosmics(image, verbose=True)
    hdulist[0].data = image
    hdulist.writeto('cosrm/'+filename.split('/')[-1])
    hdulist.close()

    hdulist = fits.open(Is[i])
    image = hdulist[0].data
    mask, image = astroscrappy.detect_cosmics(image, verbose=True)
    hdulist[0].data = image
    hdulist.writeto('cosrm/'+filename.split('/')[-1])
    hdulist.close()

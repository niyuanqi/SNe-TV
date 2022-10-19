from astropy.io import fits
import glob
import matplotlib.pyplot as plt
from PIL import Image
import sys

from ObjData import*
from SNAP.Analysis.LCRoutines import*
from SNAP.MagCalc import*

#function: load image
def load_im(filename):
    from astropy.io import fits

    hdulist = fits.open(filename)
    image = hdulist[0].data
    hdulist.close()
    return image

#function: crop image to size
def crop_im(image, rad, fac=1.0):
    sh = image.shape
    image = image[int(sh[0]/2-rad*fac):int(sh[0]/2+rad*fac),
                  int(sh[1]/2-rad):int(sh[1]/2+rad)]
    return image

#function: normalize an image
def norm_im(image, vmin=0.1, vmax=1):
    norm_im = np.power((image - vmin)/(vmax-vmin), 0.5)
    norm_im[norm_im>1] = 1.0
    norm_im[norm_im<0] = 0.0
    return norm_im

#click listener for reference stars
def ref_click(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(ix, iy)
    plt.Circle((ix, iy), 0.2, color='r')

    global ref_coords
    ref_coords.append((ix, iy))

    if len(ref_coords) == 5:
        fig.canvas.mpl_disconnect(cid)
    return ref_coords

#function: distance metric on images
def dist(x1, y1, x2, y2):
    #Euclidean distance
    return np.sqrt(np.square(x1-x2)+np.square(y1-y2))

#function: photometric aperture at (x0,y0) from r1 to r2
def ap_get(image, x0, y0, r1, r2):
    xaxis = np.arange(max([0,x0-r2]), min(image.shape[1],x0+r2+1), dtype=int)
    yaxis = np.arange(max([0,y0-r2]), min(image.shape[0],y0+r2+1), dtype=int)
    api = np.array([image[y][x] for x in xaxis for y in yaxis if (dist(x0,y0,x,y)<=r2 and dist(x0,y0,x,y)>=r1)])
    apx = np.array([x for x in xaxis for y in yaxis if (dist(x0,y0,x,y)<=r2 and dist(x0,y0,x,y)>=r1)])
    apy = np.array([y for x in xaxis for y in yaxis if (dist(x0,y0,x,y)<=r2 and dist(x0,y0,x,y)>=r1)])
    return api, apx, apy

#function: Gaussian
def S2gaus((x, y), A, sig, x0, y0):
    from scipy.signal import gaussian

    r = np.sqrt(np.square(x-x0)+np.square(y-y0))
    g = A*np.exp(-np.power(r, 2.) / (2 * np.power(sig, 2.)))
    return g.ravel()

#function: measure reference star psf
def meas_psf(im, x0, y0):
    from scipy.optimize import curve_fit

    #get centered fit box
    intens, x, y = ap_get(im, x0, y0, 0, 15)
    mask = np.logical_not(np.isnan(intens))
    #get an approximate fix on position
    x0 = np.nansum(intens[mask]*x[mask])/np.nansum(intens[mask])
    y0 = np.nansum(intens[mask]*y[mask])/np.nansum(intens[mask])
    #get centered fit box
    intens, x, y = ap_get(im, x0, y0, 0, 15)
    mask = np.logical_not(np.isnan(intens))

    est = [max(intens), 5, x0, y0]
    popt, pcov = curve_fit(S2gaus, (x[mask], y[mask]), intens[mask],
                           sigma=np.sqrt(np.absolute(intens[mask])), p0=est)
    return popt[1], popt[0]*popt[1]**2


#function: main
if __name__ == '__main__':
    import argparse
    import glob
    
    #load rgb list
    Bs, Vs, Is = np.loadtxt("rgblist.txt", unpack=True, usecols=(0, 1, 2), dtype=str)

    #command line arguments
    parser = argparse.ArgumentParser(description="Compile list of PSFs from rgb list.")
    parser.add_argument("-w", "--width", type=int, default=1100, help="width of composed image [pixels]")
    parser.add_argument("-f", "--factor", type=int, default=1.0, help="ratio of height to width")
    parser.add_argument("-r", "--reference", type=int, default=0, help="reference image number")
    args = parser.parse_args()

    #image composition
    radius = args.width
    fac = args.factor
    i = args.reference

    #load reference image
    b = load_im(glob.glob("cosrm/*"+Bs[i]+"*")[0])
    b = crop_im(b, radius, fac=fac)
    vmin = np.median(b) - 0.03*np.std(b)
    vmax = np.median(b) + 3*np.std(b)
    b = norm_im(b, vmin, vmax)
    g = load_im(glob.glob("cosrm/*"+Vs[i]+"*")[0])
    g = crop_im(g, radius, fac=fac)
    vmin = np.median(g) - 0.03*np.std(g)
    vmax = np.median(g) + 3*np.std(g)
    g = norm_im(g, vmin, vmax)
    r = load_im(glob.glob("cosrm/*"+Is[i]+"*")[0])
    r = crop_im(r, radius, fac=fac)
    vmin = np.median(r) - 0.03*np.std(r)
    vmax = np.median(r) + 3*np.std(r)
    r = norm_im(r, vmin, vmax)
    #create rgb image
    rgb = (np.dstack((r,g,b)) * 254.999).astype(np.uint8)
    img = Image.fromarray(rgb)
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)
    ax.imshow(img, origin='lower', interpolation='bilinear')

    #get reference stars
    ref_coords = []   
    #find saturation stars
    print("Click 5 reference stars for PSF matching")
    cid = fig.canvas.mpl_connect('button_press_event', ref_click)
    plt.show()

    #function: all of the above
    def do_all(filename):
        im = load_im(filename)
        im = crop_im(im, radius, fac=fac)
        #measure psf of source image
        im_sig = 0
        im_vol = 0
        for i in range(len(ref_coords)):
            sig, vol = meas_psf(im, *ref_coords[i])
            im_sig += sig/len(ref_coords)
            im_vol += vol/len(ref_coords)
        return im, im_sig, im_vol

    #measure fwhm for each image
    sigs = []
    vols = []
    for i in range(len(Bs)):
        print "Doing", Bs[i]
        b, b_sig, b_vol = do_all(glob.glob("cosrm/*"+Bs[i]+"*")[0])
        g, g_sig, g_vol = do_all(glob.glob("cosrm/*"+Vs[i]+"*")[0])
        r, r_sig, r_vol = do_all(glob.glob("cosrm/*"+Is[i]+"*")[0])
        sigs.append([b_sig, g_sig, r_sig])
        vols.append([b_vol, g_vol, r_vol])
    sigs = np.array(sigs)
    vols = np.array(vols)

    #display psf histogram
    figh = plt.figure(figsize=(5,5))
    ax2 = figh.add_subplot(111)
    ax2.hist(sigs.T[0], color='b', density=True, bins=100, alpha=0.2)
    ax2.hist(sigs.T[1], color='g', density=True, bins=100, alpha=0.2)
    ax2.hist(sigs.T[2], color='r', density=True, bins=100, alpha=0.2)
    plt.show()

    #display brightness histogram
    figh = plt.figure(figsize=(5,5))
    ax2 = figh.add_subplot(111)
    ax2.hist(vols.T[0], color='b', density=True, bins=100, alpha=0.2)
    ax2.hist(vols.T[1], color='g', density=True, bins=100, alpha=0.2)
    ax2.hist(vols.T[2], color='r', density=True, bins=100, alpha=0.2)
    plt.show()

    #output PSFs list
    np.savetxt("siglist.txt", sigs)
    #output star brightness list (optional)
    np.savetxt("vollist.txt", vols)

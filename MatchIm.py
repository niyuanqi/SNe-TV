
#function: load image
def load_im(filename):
    from astropy.io import fits

    hdulist = fits.open(filename)
    image = hdulist[0].data
    hdulist.close()
    return image

#function: mask saturation artefact
def mask_sat(image, x, y):
    image[y+20:y+500, x-12:x+13] = 0.0
    image[y-500:y-20, x-12:x+13] = 0.0
    image[y-20:y+20, x-5:x+6] = 0.0
    return image

#function: crop image to size
def crop_im(image, rad, fac=1.0):
    sh = image.shape
    image = image[int(sh[0]/2-rad*fac):int(sh[0]/2+rad*fac),
                  int(sh[1]/2-rad):int(sh[1]/2+rad)]
    return image

#function: gaussian blur an image
def blur_im(image, sig, sig0):
    from scipy.ndimage import gaussian_filter
    
    dsig = np.sqrt(sig0**2 - sig**2)
    return gaussian_filter(image, dsig)

#function: normalize an image
def norm_im(image, vmin=0.1, vmax=1, power):
    norm_im = np.power((image - vmin)/(vmax-vmin), power)
    norm_im[norm_im>1] = 1.0
    norm_im[norm_im<0] = 0.0
    return norm_im

#function: measure sky information
def meas_hist(image, x, y, r):
    from SNAP.Photometry import ap_get

    skyi, skyx, skyy = ap_get(image, x, y, 0, r)
    iqr = np.percentile(skyi, 75) - np.percentile(skyi, 50)
    return np.median(skyi), iqr


if __name__ == '__main__':
    import glob
    import matplotlib.pyplot as plt
    from PIL import Image
    import sys
    
    #option selector argument
    sysarg = sys.argv[1]
    # 0 -> do all images in rgb list
    # 1 -> just plot one image, don't save anything
    # 2 -> same as 1, but shows histograms

    #Adjust the following parameters:
    ########################################
    #image composition
    radius = 1100
    fac = 1.

    #select images by PSF
    clip = 80 #percentile

    #histogram sampling patch for the object of interest (usually host galaxy)
    objx = 1025
    objy = 602
    objr = 100 #radius of patch
    #histogram sampling patch for the sky background (exclude stars)
    skyx = 979
    skyy = 1669
    skyr = 70 #radius of patch
    #normalization factors
    smax = 1.0 #max sigma of sky background, sets black-point of image
    amin = 1.3 #minimum sigma of object that you want to show
    amax = 1.8 #max sigma of object to set the image white-point
    power = 0.5 #power of intensity scale (lower -> favor midtones/highlights)

    #load rgb list
    Bs, Vs, Is = np.loadtxt("rgblist.txt", unpack=True, usecols=(0, 1, 2), dtype=str)
    ts = np.loadtxt("rgblist.txt", usecols=(3,), dtype=float)
    #load psf list and mask best images
    sigs = np.loadtxt("siglist.txt")
    clip = np.percentile(sigs.flatten(), 80)
    smask = sigs <= clip
    #load brightness list (optional: can be used to scale image below)
    vols = np.loadtxt("vollist.txt")

    #function: all of the above
    def do_all(filename, sig, vol, wb=1.0):
        im = load_im(filename)
        im = crop_im(im, radius, fac=fac)
        if sig < clip:
            #blur to the PSF clip value
            im = blur_im(im, sig, clip)
        im = mask_sat(im, 1629, 1629)
        obj, std = meas_hist(im, objx, objy, objr)
        sky, noise = meas_hist(im, skyx, skyy, skyr)

        #scale histogram using reference stars
        #print "vmin", sky+amin*noise, vol*0.0007
        #vmin = max(sky+amin*noise, vol*0.0007*wb)
        #im = norm_im(im, vmin, amax*vol*wb, power)

        #scale histogram using object
        print "vmin", sky+smax*noise, (obj-amin*std)*wb
        vmin = max(sky+smax*noise, (obj-amin*std)*wb)
        im = norm_im(im, vmin, (obj+amax*std)*wb, power)
        return im
    
    #color composite each set of images
    for i in range(len(Bs)):
        #only do images with good enough PSF
        if smask[i].all():
            print "Doing", i, Bs[i]
        
            #process rgb images
            b = do_all(glob.glob("cosrm/*"+Bs[i]+"*")[0], sigs[i][0], vols[i][0], 1.3)
            g = do_all(glob.glob("cosrm/*"+Vs[i]+"*")[0], sigs[i][1], vols[i][1], 1.15)
            r = do_all(glob.glob("cosrm/*"+Is[i]+"*")[0], sigs[i][2], vols[i][2], 1.0)
            #create rgb image
            rgb = (np.dstack((r,g,b)) * 254.999).astype(np.uint8)
            img = Image.fromarray(rgb)
            #plot image
            fig = plt.figure(figsize=(7,7))
            fig.patch.set_visible(False)
            ax = fig.add_subplot(111)
            ax.imshow(img, origin='lower', interpolation='bilinear')
            
            #plot figure labels
            ##################################################
            # Edit parameters below based on your composition

            #arrow indicating SN position
            ax.arrow(b.shape[1]/2 + radius/6, b.shape[0]/2 - radius/6,
                     -radius/11, radius/11, head_width=20, color='w')
            #epoch: date and time
            ax.text(80, b.shape[0]-120, Bs[i][2:8]+' - '+Bs[i][9:11]+':'+Bs[i][11:13], color='w')
            #scale indicator
            ax.plot([70, 220], [70, 70], linewidth=2, color='w')
            ax.text(120, 100, "1'", color='w')
            #compass
            ax.plot([b.shape[1]-300, b.shape[1]-150],
                    [b.shape[0]-300, b.shape[0]-300], linewidth=2, color='w')
            ax.text(b.shape[1]-320, b.shape[0]-120, "N", color='w')
            ax.plot([b.shape[1]-300, b.shape[1]-300],
                    [b.shape[0]-300, b.shape[0]-150], linewidth=2, color='w')
            ax.text(b.shape[1]-120, b.shape[0]-320, "W", color='w')
            ##################################################
            plt.axis("off")
            extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        
    
            #output options
            if sysarg == "0":
                plt.savefig('rgbout/'+str(i)+'.png',  bbox_inches=extent)
            elif sysarg == "1":
                plt.show()
            elif sysarg == "2":
                figh = plt.figure(figsize=(5,5))
                ax2 = figh.add_subplot(111)
                rhist = r.flatten()[np.logical_not(np.isnan(r.flatten()))]
                ghist = g.flatten()[np.logical_not(np.isnan(g.flatten()))]
                bhist = b.flatten()[np.logical_not(np.isnan(b.flatten()))]
                ax2.hist(rhist, color='r', density=True, bins=100, alpha=0.3)
                ax2.hist(ghist, color='g', density=True, bins=100, alpha=0.3)
                ax2.hist(bhist, color='b', density=True, bins=100, alpha=0.3)
                plt.show()
            else:
                print "invalid output option"
    

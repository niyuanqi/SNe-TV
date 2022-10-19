
#function: compile rgb list of files
def rgblist(ss, ts, dt):
    """
    s: [filename ids lists for r, g, b]
    t: [time lists in days for r, g, b]
    dt: time separation upper limit for r g b channels [days]
    """

    Bs = [] #b identifier
    Bts = [] #b times
    Vs = [] #g identifier
    Vts = [] #g times
    Is = [] #r identifier
    Its = [] #r times
    #curate list based on the B-band
    for i in range(len(tl[0])):
        #evaluate time differences within 5 mins
        Vmask = np.absolute(tl[1] - tl[0][i]) < dt
        Imask = np.absolute(tl[2] - tl[0][i]) < dt
        #only add to list if all bands are observed
        if np.sum(Vmask) > 0 and np.sum(Imask) > 0:
            #only add if none of these images have been added before
            Vclose = ""
            Vcloses = ssl[1][Vmask]
            Vtcloses = tl[1][Vmask]
            for j in range(len(Vcloses)):
                if Vcloses[j] not in Vs:
                    Vclose = Vcloses[j]
                    Vtclose = Vtcloses[j]
                    break
            Iclose = ""
            Icloses = ssl[2][Imask]
            Itcloses = tl[2][Imask]
            for j in range(len(Icloses)):
                if Icloses[j] not in Is:
                    Iclose = Icloses[j]
                    Itclose = Itcloses[j]
                    break
            #check if both unadded V and I are found
            if Vclose != "" and Iclose != "":
                Bs.append(ssl[0][i])
                Bts.append(tl[0][i])
                Vs.append(Vclose)
                Vts.append(Vtclose)
                Is.append(Iclose)
                Its.append(Itclose)
            
            #print Bs[-1], Vs[-1], Is[-1]

    #save filelist
    np.savetxt("rgblist.txt", np.array([Bs, Vs, Is, Bts]).T)

if __name__ == '__main__':
    import glob
    import argparse
    from astropy.io import fits
    from astropy.time import Time

    #command line arguments
    parser = argparse.ArgumentParser(description="Compile list of rgb image triplets from random folder of BVi-band images.")
    parser.add_argument("filenames", nargs='*', help="files to include")
    parser.add_argument("-y", "--year", type=int, default=2020, help="SN discovery year")
    parser.add_argument("-dt", "--dt", type=float, default=0.004, help="rgb max time difference")
    args = parser.parse_args()

    filenames = args.filenames.sort()
    filenames.sort()
    year = args.year
    dt = args.dt

    ts = [[], [], []]
    ss = [[], [], []]
    #set reference time
    t_ref = str(year)+"-01-01T00:00:00.000"
    t_ref = Time(t_ref, format='isot', scale='utc')
    #load each file
    for filename in range(len(filenames)):
        hdulist = fits.open(filename)
        image = hdulist[0].data
        header = hdulist[0].header
        time = Time(header['DATE-OBS'], format='isot', scale='utc')
        time = float((time - t_ref).value) #time in days of year
        
        if 'B' in filename:
            ss[0].append(filename.split('/')[-1])
            ts[0].append(time)
        elif 'V' in filename:
            ss[1].append(filename.split('/')[-1])
            ts[1].append(time)
        elif 'I' in filename:
            ss[2].append(filename.split('/')[-1])
            ts[2].append(time)
        else:
            print("Unknown band", filename)

    #make rgb list in txt
    rgblist(ss, ts, dt)

#################################################################
# Name:     ObjData.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Important configuration data.                       #
#################################################################

#essential modules
import numpy as np

#window of interest
t1 = 87
t2 = 200
tdet = 87
t1_early = 85
t2_early = 94.5
limcuts = [20,20,20]

#supernova basic data
name = "KSP-N3923-2_2018ku"
RA = 177.757515
DEC = -28.744022

#kasfile = "kasen.CN_190822.txt"
kasfile = "kasen-mc-200915/kasen.exp.CN_200930.txt"

#background star
star_mag = 22.39
star_err = 0.08

#light curve files
Bfile = "../KSP-N3923-2-2018ku.B.lc.CN_190103.txt"
Vfile = "../KSP-N3923-2-2018ku.V.lc.CN_190103.txt"
Ifile = "../KSP-N3923-2-2018ku.I.lc.CN_190103.txt"
Ifile = "../KSP-N3923-2-2018ku.I.lcapmulti.CN_200915.txt"
files = [Bfile, Vfile, Ifile]
#S correction files
Bscorrfile = "../Scorrect/KSP-N3923-2-2018ku.B.scorr.CN_190207.txt"
Iscorrfile = "../Scorrect/KSP-N3923-2-2018ku.I.scorr.CN_190207.txt"
firstspec = 91.0
#observed bands
band = ['B','V','i']
#band labels
Band = ['B','V','I']
#snoopy bands
snband_all = ['B', 'V', 'g', 'r', 'i']
#snband = ['B','V','i']
#parameters for template fitting
#SNooPy format light curve
sn_file_all = "KSP-N3923-2_2018ku_final.txt"
#sn_file_all = "KSP-N3923-2_2018ku_all.txt"
#sn_file = "KSP-N3923-2_2018ku.txt"
#Phillips relation dataset (including sBV) from Burns
ph_file = "Phillips_dm15+st_results.dat"
#redshift
#z = 0.005767
#zerr = 0.000667
z = 0.005801
zerr = 0.000030
#Extinction coefficient (galactic) in each band S & F (2011)
#177h45m27.02s -28d44m38.48s Equ J2000
#EBVgal = 0.0934 is from S&F 2011, but their coefs are calibrated in SFD1998
EBVgal = 0.1086 #SFD1998

#Chris Ni:
#220122 the above ra and dec query may be mistaken?
#11h51m01.80s -28d44m38.48s Equ J2000
#EBVgal = 0.0839 #SFD1998

#CTIO B, CTIO V, SDSS i
Coefs = np.array([3.641, 2.682, 1.698])
#CTIO U, CTIO B, CTIO V, SDSS g, SDSS r, SDSS i
#snCoefs = np.array([4.107, 3.641, 2.682, 3.303, 2.285, 1.698])

#parameters from template fitting (for early light curve fitting)
#maximum epoch
Tmax = 103.272
Tmaxerr = 0.021
#measured dm15
dm = 1.117
sbv = 0.797
#distance modulus
DM = 31.75
DMerr = 0.08
#luminosity distance 
Dl = np.power(10, 0.2*(DM+5))
Dlerr = Dl*np.log(10)*0.2*DMerr
#max observed magnitudes
max_obs = np.array([12.808, 12.798, 13.259])
#max absolute magnitudes
max_abs = np.array([-19.319, -19.226, -18.614])
#max times in rest frame
max_times = np.array([-0.858736, -0.965898, -3.930457])

#early binned time series data files
binBfile = "KSP-N3923-2-2018ku.B.color.CN_190103.txt"
binVfile = "KSP-N3923-2-2018ku.V.color.CN_190103.txt"
binIfile = "KSP-N3923-2-2018ku.I.color.CN_190103.txt"
binIfile = "KSP-N3923-2-2018ku.I.colorapmulti.CN_200915.txt"
binfiles = [binBfile, binVfile, binIfile]

#parameters from early light curve fitting (for Arnett fitting)
#t0 = -16.231
#t0err = 0.015
#t0obs = -16.326
#t0obserr = 0.015
t0 = -16.180
t0err = 0.020
t0obs = -16.274
t0obserr = 0.020

#parameters from spectra
vej = 11.43
vej_err = 0.12

#bolometric light curve
#bolofile = "KSP-N3923-2_2018ku.uvoirbolo.CN_190822.txt"
bolofile = "KSP-N3923-2_2018ku.uvoirbolo.CN_210204.txt"
#bololbfile = "KSP-N3923-2_2018ku.ukbolo.txt"
#bolometric light curve from template warping
#bolotempfile = "KSP-N3923-2_2018ku.bolotemp.CN_190822.txt"
bolotempfile = "KSP-N3923-2_2018ku.bolotemp.CN_210204.txt"

#parameters from Arnett fitting (for Kasen modelling)
#explosion parameters
#m_ni = 0.785 #solar mass
#e_51 = 0.616 #x10^51 ergs
#m_ni_err = 0.044
#e_51_err = 0.036

#m_ni = 0.819 #solar mass
#e_51 = 0.642 #x10^51 ergs
#m_ni_err = 0.048
#e_51_err = 0.040

#2021 uv issue
m_ni = 0.799 #solar mass
e_51 = 0.626 #x10^51 ergs
m_ni_err = 0.044
e_51_err = 0.037

#parameters from Polin model
m_he = 1.05
e_he = (3./10.)*(m_he*2.e33)*(vej*1.e8)**2/1.e51

limBfile = "kasen-mc/KSP-N3923-2-2018ku.B.lcbin.CN_190103.S3.txt"
limVfile = "kasen-mc/KSP-N3923-2-2018ku.V.lcbin.CN_190103.S3.txt"
limIfile = "kasen-mc/KSP-N3923-2-2018ku.I.lcbin.CN_190103.S3.txt"
limfiles = [limBfile, limVfile, limIfile]

#function: load final light curve for this source
def load_final_lc(tlim1=tdet, tlim2=600, SNthres=3.0, binned=False, retSN=False, retlims=False, Bplat=True):
    #essential imports
    from SNAP.Analysis.LCRoutines import*
    from SNAP.Analysis.FilterCorrect import*
    from SNAP.Analysis.Cosmology import*
    
    #load light curve file
    if binned:
        #get light curve for SN>SNthres
        t, M, M_err, F, SN, Mlim, ra, dec, ss = LCload(binfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=SNthres, racols=2, deccols=3, scols=9, flags=['-99.99999'], mode='multi')
    else:
        #get light curve for SN>SNthres
        t, M, M_err, F, SN, Mlim, ra, dec, ss = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=SNthres, racols=2, deccols=3, scols=9, flags=['-99.99999'], mode='multi')
    #crop window in data
    for i in range(len(M)):
        t[i], M[i], M_err[i], F[i], SN[i], Mlim[i] = LCcrop(t[i], tlim1, tlim2, M[i], M_err[i], F[i], SN[i], Mlim[i])
    #flux corrections for lcs made before 190215 for flux_0
    fluxcorr = [4063.0/4260.0, 3636.0/3640.0, 3631.0/2550.0]
    for i in [0,1]:
        F[i] = F[i]*fluxcorr[i]
    
    #subtract background star from i-band
    flux_sub, err_sub = Flux_subMag(F[2], F[2]/SN[2], star_mag, star_err)
    F[2] = flux_sub
    SN[2] = flux_sub/err_sub
    M[2], M_err[2] = Mag_subMag(M[2], M_err[2], star_mag, star_err)
    Mlim[2] = Mag_addMag(Mlim[2], star_mag)
    imask = SN[2] > SNthres
    t[2], M[2], M_err[2], F[2], SN[2], Mlim[2] = t[2][imask], M[2][imask], M_err[2][imask], F[2][imask], SN[2][imask], Mlim[2][imask]
    #saturated I-band replacement
    t1, M1, M1_err, F1, SN1, Mlim1 = LCcrop(t[2], tlim1, 89, M[2], M_err[2], F[2], SN[2], Mlim[2])
    Isatfile = ["KSP-N3923-2-2018ku.I.lc.CN_190103.txt"]
    ti, Mi, Mi_err, Fi, SNi, Mlimi, rai, deci, ssi = LCload(Isatfile, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=SNthres, racols=2, deccols=3, scols=9, flags=['-99.99999'], mode='multi')
    Fi[0] = Fi[0]*fluxcorr[2] #correct the flux for lcs made before 190215 for flux_0
    t2, M2, M2_err, F2, SN2, Mlim2 = LCcrop(ti[0], max(tlim1, 89), min(tlim2, 180), Mi[0], Mi_err[0], Fi[0], SNi[0], Mlimi[0])
    t3, M3, M3_err, F3, SN3, Mlim3 = LCcrop(t[2], 180, tlim2, M[2], M_err[2], F[2], SN[2], Mlim[2])
    t[2] = np.concatenate([t1, t2, t3])
    M[2] = np.concatenate([M1, M2, M3])
    M_err[2] = np.concatenate([M1_err, M2_err, M3_err])
    F[2] = np.concatenate([F1, F2, F3])
    SN[2] = np.concatenate([SN1, SN2, SN3])
    Mlim[2] = np.concatenate([Mlim1, Mlim2, Mlim3])
    #get noise in flux
    F_err = [F[i]/SN[i] for i in range(len(t))]
    
    #filter correct B-band limiting magnitudes
    t, Mlim = BVcorrectLim(t, Mlim, mBVr=0.639)
    #perform filter corrections for detections after the first detection
    tc, Mc, Mc_err = [],[],[]
    for i in range(len(M)):
        SNci, tci, Mci, Mci_err = LCcrop(SN[i], 3, max(SN[i])+1, t[i], M[i], M_err[i])
        #don't correct nebular light curve
        tci, Mci, Mci_err = LCcrop(tci, tdet, min(tlim2, 180), Mci, Mci_err)
        tc.append(tci)
        Mc.append(Mci)
        Mc_err.append(Mci_err)
    #S correct B band
    tcorr, scorr, scorr_err = np.loadtxt(Bscorrfile, unpack=True, comments=';')
    tc, Mc, Mc_err = SBcorrectMag(tc,Mc,Mc_err,tcorr, scorr, tdiv=firstspec,
                                  SBVega=-0.00238,mBVr=0.639,mBVrerr=0.040,
                                  interp='lin', Sinterp='GP', Scorr_err=scorr_err)
    #S correct i band
    tcorr, scorr, scorr_err = np.loadtxt(Iscorrfile, unpack=True, comments=';')
    tc, Mc, Mc_err = SIcorrectMag(tc,Mc,Mc_err,tcorr, scorr, tdiv=firstspec,
                                  SIVega=0.09287,
                                  interp='lin', Sinterp='GP', Scorr_err=scorr_err)
    #slip corrected data into light curve
    for i in range(len(M)):
        #after first detection
        mask = np.logical_and(t[i] > tdet, t[i] < min(tlim2, 180))
        mask = np.logical_and(mask, SN[i] > 3)
        #replace with corrected data
        M[i][mask] = Mc[i]
        M_err[i][mask] = Mc_err[i]
        fc, fc_err = Mag_toFlux(Band[i], Mc[i], Mc_err[i])
        F[i][mask] = fc*1e6
        F_err[i][mask] = fc_err*1e6
    #check if omit B-band plateau
    if not Bplat:
        i = 0
        mask = np.logical_or(t[i]<tdet, t[i]>tdet+1)
        t[i], M[i], M_err[i], F[i], F_err[i], SN[i], Mlim[i] = t[i][mask], M[i][mask], M_err[i][mask], F[i][mask], F_err[i][mask], SN[i][mask], Mlim[i][mask]
    
    #return finalized light curve
    retlist = [t, M, M_err, F, F_err, Mlim]
    if retSN:
        retlist += [SN]
    if retlims:
        #get limiting magnitudes
        tl, Ml, Ml_err, Fl, SNl, Mliml, ral, decl, ssl = LCload(limfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-89.0, racols=2, deccols=3, scols=9, flags=['-99.99999'], mode='multi')
        #get corresponding magnitudes for i band SN calculation
        if binned:
            td, Md, Md_err, Fd, SNd, Mlimd, rad, decd, ssd = LCload(binfiles, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-89.0, racols=2, deccols=3, scols=9, flags=['-99.99999'], mode='multi')
        else:
            td, Md, Md_err, Fd, SNd, Mlimd, rad, decd, ssd = LCload(files, tcol=0, magcols=6, errcols=7, fluxcols=4, SNcols=5, limcols=8, SNthres=-89.0, racols=2, deccols=3, scols=9, flags=['-99.99999'], mode='multi')
        flux_sub, err_sub = Flux_subMag(Fd[2], Fd[2]/SNd[2], star_mag, star_err)
        SNl[2] = flux_sub/err_sub
        SNl[2] = SNl[2][np.logical_and(td[2]>=tl[2][0], td[2]<=tl[2][-1])]
        for i in range(len(Ml)):
            mask = np.logical_and(tl[i]>tlim1, tl[i]<tlim2)
            mask = np.logical_and(mask, SNl[i] < SNthres)
            tl[i], Mliml[i] = tl[i][mask], Mliml[i][mask]
        #filter correct B-band limiting magnitudes
        tl, Mliml = BVcorrectLim(tl, Mliml, mBVr=0.639)
        Mliml[2] = Mag_addMag(Mliml[2], star_mag)
        retlist += [tl, Mliml]
    #returned list
    return retlist

#function: source values at binned epoch
def return_bin_values():
    #essential imports
    from SNAP.Analysis.LCRoutines import*
    from SNAP.Analysis.FilterCorrect import*
    from SNAP.Analysis.Cosmology import*
    
    #Binned image data
    bint = np.array([[87.517990], [87.509867], [87.511320]])
    binterr = np.array([[0.044], [0.035], [0.035]])
    binF = np.array([[10.984529], [35.863009], [30.871249]])
    binSN = np.array([[7.446], [27.30], [16.31]])
    binM = np.array([[21.472], [20.016], [20.176]])
    binMerr = np.array([[0.149], [0.040], [0.067]])
    binMlim = np.array([[22.370], [22.426], [22.639]])
    bint, binMlim = BVcorrectLim(bint, binMlim, mBVr=0.639)
    #subtract background star from i-band
    flux_sub, err_sub = Flux_subMag(binF[2], binF[2]/binSN[2], star_mag, star_err)
    binF[2] = flux_sub
    binSN[2] = flux_sub/err_sub
    binM[2], binMerr[2] = Mag_subMag(binM[2], binMerr[2], star_mag, star_err)
    #S correct B band
    tcorr, scorr, scorr_err = np.loadtxt(Bscorrfile, unpack=True, comments=';')
    bint, binM, binMerr = SBcorrectMag(bint,binM,binMerr,tcorr, scorr,
                                       tdiv=firstspec, SBVega=-0.00238,
                                       mBVr=0.639,mBVrerr=0.040,interp='lin', 
                                       Sinterp='GP', Scorr_err=scorr_err)
    #S correct i band
    tcorr, scorr, scorr_err = np.loadtxt(Iscorrfile, unpack=True, comments=';')
    bint, binM, binMerr = SIcorrectMag(bint,binM,binMerr,tcorr, scorr,
                                       tdiv=firstspec,SIVega=0.09287, interp='lin',
                                       Sinterp='GP', Scorr_err=scorr_err)
    return bint, binterr, binM, binMerr

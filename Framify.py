import matplotlib.pyplot as plt
import glob
import numpy as np

#get sorted frame list
filenames = glob.glob("rgbout/*")
nums = [int((fn.split('.')[0]).split('/')[-1]) for fn in filenames]
filenames = np.array([x for _, x in sorted(zip(nums, filenames))])
ixs = sorted(nums) #indices
#load times list
t_ims = np.loadtxt("rgblist.txt", usecols=(3,), dtype=float)
#times of processed frames
t_ims = np.array([t_ims[ix] for ix in ixs])

#make time axis
dt = np.min(t_ims[1:] - t_ims[0:-1])
#fill in the gaps
################################################
tol = 45 #frame skip tolerance
################################################

#output filename padding
def zpad(num, l):
    out = str(num)
    for i in range(l - len(out)):
        out = "0" + out
    return out

#save first frame
num = 0
######################################################
#load image and make final cut to get even dimensions
img = plt.imread(filenames[0])
img_cropped = img[:, 2:-1]
######################################################
#tick the frame counter
outnum = zpad(num, 4)
num += 1
#save the frame
plt.imsave('frames/'+outnum+".png", img_cropped)

#elapse the timelist
for i in range(len(t_ims)-1):
    print num
    #check if the last bin was too long
    if t_ims[i+1] - t_ims[i] > tol*dt:
        #add in pause frames
        npause = int(np.floor((t_ims[i+1] - t_ims[i])/(tol*dt)))
        for n in range(npause):
            #tick the frame counter
            outnum = zpad(num, 4)
            num += 1
            #save the frame
            plt.imsave('frames/'+outnum+".png", img_cropped)
    ######################################################
    #load image and make final cut to get even dimensions
    img = plt.imread(filenames[i+1])
    img_cropped = img[:, 2:-1]
    ######################################################
    #tick the frame counter
    outnum = zpad(num, 4)
    num += 1
    #save the frame
    plt.imsave('frames/'+outnum+".png", img_cropped)

#time step frame speed
def frame_speed(dt, fps):
    dtsec = dt*24*60*60 #time step [s]
    ratio = dtsec*fps #time lapse [s] per s
    return ratio #speed X
def vid_len(num, fps):
    vidlen = (num+1.)/fps #s
    return vidlen #video length

print "Number of frames:", num
print "Video speed [X] at 45 fps:", frame_speed(dt, 45)
print "Video length [s] at 45 fps:", vid_len(num, 45)

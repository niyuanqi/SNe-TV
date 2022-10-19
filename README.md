# SNeTV
SuperNova(e) and Transients Visualizer

I wrote this code pretty fast in order to compose 3-channel supernova images into time-evolving videos that I can show off when giving talks. Given that, you should be able to get good results with the suite of programs that I've included in  this package if you follow these steps with your own tweaks as needed. There doesn't seem to exist any other program for doing this kind of task that prioritized the aesthetics and artistic composition of the output video.


1. Make a directory for your project, and do everything in it.

2. Setting up your raw images: Crop the images that you'd like to include in your video and move them to a subdirectory called "crops". You don't need to decide your final video composition yet, just crop your image for processing efficiency to a region that contains just about everything you might want to show, including the supernova/transient, its host galaxy, and make sure to include numerous foreground field stars.

2. Aligning the stars in the images: Make a subdirectory called "align", and run AlignIm.py in order to register your images and align their coordinates based on the foreground field stars. It requires the Astroalign (https://astroalign.quatrope.org/en/latest/) package, which requires Python3 (The rest of the steps are in Python2). You should select a not terrible image for the target in the code. A handful of images may fail if the image quality is too poor, and that's fine.

3. Compile list of images that will become RGB: Use RGBSetup.py to compile the RGB list, rgblist.txt, whose four columns are bname, gname, rname, time. These are the files that will go into the final video.

4. Automatically remove cosmic rays: Make a subdirectory called "cosrm", and run Cosrm.py in order to remove cosmic rays. It requires the Astro-SCRAPPY package (https://astroscrappy.readthedocs.io/en/latest/).

5. Measure image PSFs: Use PSFMsr.py to compose your final image region and select 5 reference stars by clicking close to their centers within the region. This will create a psflist.txt that evaluates the seeing of each image for final matching.

6. Make matched RGB images: Make a subdirectory called "rgbout" and use MatchIm.py to perform the image matching and final histogram scaling. Note, you need to edit the parameters in this code to match the histograms appropriately according to what you have composed.

7. Make video frames: Make a subdirectory called "frames" and use Framify.py to make equally spaced video frames. Note you need to choose how many frames you allow the code to skip in the code when there is a gap in the time series. Also, check what the dimensions of the rgbout/ images are, because the video frames need to be even dimension on both sides. Do some last minute cropping in Framify.py otherwise.

Now you can make your video. I prefer to use ffmpeg. Run this, and change the resolution to the pixel scale of your frames.

ffmpeg -r 45 -f image2 -s 546x546 -i frames/%04d.png -vcodec libx264 -crf 10 -pix_fmt yuv420p out45.mp4

If there is a bit of white border from the last step, do a bit of cropping to get rid of it.

ffmpeg -i out45.mp4 -filter:v "crop=536:536:5:5" out45-crop.mp4

If you want to add any colored border for effect, use the following.

ffmpeg -i out45-crop.mp4 -filter_complex "[0]pad=w=iw+20:h=ih+20:x=10:y=10:color=orange" out45-border.mp4
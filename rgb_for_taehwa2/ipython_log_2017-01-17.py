########################################################
# Started Logging At: 2017-01-17 14:06:08
########################################################

import reproject
get_ipython().magic('paste')
import numpy as np
import copy
import matplotlib
from astropy.io import fits
import reproject

import os

b1 = fits.open('AG-Laboca-Planck.49.5.fits')
b2 = fits.open('W51_90cm_CDBB.fits')
b3 = fits.open('grs_13CO_max_30to90kms.fits')

hdr = fits.Header.fromtextfile('hdr4096.hdr')

b1d,_ = reproject.reproject_interp(b1, hdr)
b2d,_ = reproject.reproject_interp(b2, hdr)
b3d,_ = reproject.reproject_interp(b3, hdr)
b2
b2[0]
b2[1].data
b2[0].data
b2 = fits.open('W51_90cm_CDBB.fits')[0]
b2.data = b2.data.squeeze()
b2
b2d,_ = reproject.reproject_interp(b2, hdr)
b2.data.shape
b2d,_ = reproject.reproject_interp(b2, hdr)
hdr
get_ipython().magic('debug')
get_ipython().magic('pinfo reproject.reproject_interp')
get_ipython().magic('debug')
get_ipython().magic('paste')
import numpy as np
import copy
import matplotlib
from astropy.io import fits
from astropy import wcs
import reproject

import os

b1 = fits.open('AG-Laboca-Planck.49.5.fits')
b2 = fits.open('W51_90cm_CDBB.fits')[0]
b3 = fits.open('grs_13CO_max_30to90kms.fits')

hdr = fits.Header.fromtextfile('hdr4096.hdr')

b1d,_ = reproject.reproject_interp(b1, hdr)
b2d,_ = reproject.reproject_interp((b2.data, wcs.WCS(b2.header).celestial), hdr)
b3d,_ = reproject.reproject_interp(b3, hdr)

get_ipython().magic('run rgb_overview_big.py')
import pylab as pl
pl.imshow(rgb)
rgb.shape
rgb = np.array([red,green,blue]).T.swapaxes(0,1)
pl.imshow(rgb)
rgb
linearize(red)
linearize(red).max()
get_ipython().magic('paste')
#red = logscale(linearize(b2d, -0.05, 0.5))
red = linearize(b2d)
green = (linearize(b3d, 0, 7))
# blue = logscale(linearize(b1d, 0.05, 4))
blue = linearize(b1d)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
get_ipython().magic('paste')
red = linearize(b2d, -0.05, 0.5)
green = (linearize(b3d, 0, 7))
# blue = logscale(linearize(b1d, 0.05, 4))
blue = linearize(b1d, 0.05, 4)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
logscale(b2d)
get_ipython().magic('paste')
red = logscale(b2d, xmin=-0.05, xmax=0.5)
green = (linearize(b3d, 0, 7))
blue = logscale(b1d, xmin=0.05, xmax=4)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
get_ipython().magic('paste')
red = logscale(b2d, xmin=0.05, xmax=0.5)
green = (linearize(b3d, 0, 7))
blue = logscale(b1d, xmin=0.15, xmax=4)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
pl.imshow(logscale(b2d, xmin=-0.05, xmax=0.5))
pl.imshow(logscale(b2d, xmin=-0.00, xmax=0.5))
pl.imshow(logscale(b2d, xmin=-0.01, xmax=0.5))
pl.imshow(logscale(b1d, xmin=0.15, xmax=4))
pl.imshow(logscale(b1d, xmin=0.15, xmax=5))
pl.imshow(logscale(b1d, xmin=0.15, xmax=6))
pl.imshow(logscale(b1d, xmin=0.25, xmax=6))
pl.imshow(logscale(b1d, xmin=0.35, xmax=6))
pl.imshow(logscale(b1d, xmin=0.35, xmax=7))
pl.imshow(logscale(b1d, xmin=0.4, xmax=7))
pl.imshow(logscale(b1d, xmin=0.4, xmax=9))
pl.imshow(logscale(b1d, xmin=0.4, xmax=10))
pl.imshow(linearize(b1d, xmin=0.4, xmax=10))
pl.imshow(linearize(b1d,0,10) * 0.8 + logscale(b1d, xmin=10, xmax=50)*0.2)
(linearize(b1d,0,10) * 0.8 + logscale(b1d, xmin=10, xmax=50)*0.2).max()
np.nanmax((linearize(b1d,0,10) * 0.8 + logscale(b1d, xmin=10, xmax=50)*0.2))
logscale(b1d, xmin=10, xmax=50)
logscale(b1d, xmin=10, xmax=50).max()
get_ipython().magic('paste')
red = logscale(b2d, xmin=0.05, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = logscale(b1d, xmin=0.4, xmax=10, toint=False)
blue = linearize(b1d,0,10) * 0.8 + logscale(b1d, xmin=10, xmax=50, toint=False)*0.2

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
pl.imshow(blue)
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
np.nanmax(blue)
np.nanmax(logscale(b1d, xmin=5, xmax=50, toint=False))
get_ipython().magic('paste')
def linearize(x, xmin=None, xmax=None, truncate=True):
    if np.isscalar(x):
        return x
    else:
        if xmin is None:
            xmin = np.nanmin(x)
        if xmax is None:
            xmax = np.nanmax(x)
        if truncate:
            x = np.copy(x)
            x[x<xmin] = xmin
            x[x>xmax] = xmax
        return ((x-xmin)/(xmax-xmin))

def logscale(arr, logexp=3.0, toint=True, relinearize=True, **kwargs):
    linarr = linearize(arr, **kwargs)
    if logexp is None:
        logarr = linarr
    else:
        logarr = np.log10(linarr * 10**logexp + 1)
    if relinearize:
        return linearize(logarr)
    elif toint:
        lla = linearize(logarr)*255
        return lla.astype('uint8')
    else:
        return logarr

def expscale(arr, exp=2, toint=True, **kwargs):
    linarr = linearize(arr, **kwargs)
    if toint:
        lla = linearize(linarr**exp)*255
        return lla.astype('uint8')
    else:
        return linarr**exp


red = logscale(b2d, xmin=0.05, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = logscale(b1d, xmin=0.4, xmax=10, toint=False)
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
get_ipython().magic('paste')
red = logscale(b2d, xmin=0.05, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
blue = logscale(b1d, xmin=0.4, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
get_ipython().magic('paste')
red = logscale(b2d, xmin=-0.05, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
blue = logscale(b1d, xmin=0.04, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
get_ipython().magic('paste')
red = logscale(b2d, xmin=0.05, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
blue = logscale(b1d, xmin=0.04, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
get_ipython().magic('paste')
red = logscale(b2d, xmin=0.05, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb)
get_ipython().magic('paste')
red = logscale(b2d, xmin=0.05, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[400:1300,400:1600,:])
pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = logscale(b2d, xmin=0.025, xmax=0.5, toint=False)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
blue = linearize(b1d,0,5) * 0.7 + logscale(b1d, xmin=5, xmax=50, toint=False)*0.3
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
blue = linearize(b1d,0,1) * 0.7 + logscale(b1d, xmin=1, xmax=50, toint=False)*0.3
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
blue = linearize(b1d,0,0.5) * 0.7 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.3
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d,-0.05,0.025) * 0.25 + logscale(b2d, xmin=0.025, xmax=0.5, toint=False)*0.75
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d,-0.05,0.025) * 0.25 + logscale(b2d, xmin=0.025, xmax=0.5, toint=False)*0.75
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
pl.imshow(red)
pl.imshow(red)
get_ipython().magic('paste')
red = linearize(b2d,-0.05,0.01) * 0.25 + logscale(b2d, xmin=0.01, xmax=0.5, toint=False)*0.75
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d,-0.05,0.01) * 0.25 + logscale(b2d, xmin=0.01, xmax=1.5, toint=False)*0.75
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
pl.imshow(red)
pl.imshow(linearize(b2d, 0.01, 1.5))
get_ipython().magic('paste')
red = linearize(b2d, xmin=0.01, xmax=1.5)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d, xmin=0.01, xmax=0.5)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d, xmin=0.01, xmax=0.2)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d, xmin=0.01, xmax=0.3)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d, xmin=0.01, xmax=0.25)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d,-0.05,0.2) * 0.75 + logscale(b2d, xmin=0.2, xmax=1.5, toint=False)*0.25
#red = linearize(b2d, xmin=0.01, xmax=0.2)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
get_ipython().magic('paste')
red = linearize(b2d,-0.05,0.2) * 0.75 + logscale(b2d, xmin=0.2, xmax=0.7, toint=False)*0.25
#red = linearize(b2d, xmin=0.01, xmax=0.2)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
pl.imshow(red)
get_ipython().magic('paste')
red = linearize(b2d, xmin=-0.05, xmax=0.2) * 0.75 + linearize(b2d, xmin=0.2, xmax=0.7)*0.25
#red = linearize(b2d, xmin=0.01, xmax=0.2)
green = (linearize(b3d, 0, 7))
blue = linearize(b1d,0,0.5) * 0.5 + logscale(b1d, xmin=0.5, xmax=50, toint=False)*0.5
#blue = logscale(b1d, xmin=0.10, xmax=10, toint=False)

rgb = np.array([red,green,blue]).T.swapaxes(0,1)

pl.imshow(rgb[300:1300,400:1600,:])
slices = slice(300,1300), slice(400,1600) #[300:1300,400:1600]
rgb_crop = rgb[slices,:]
rgb_crop = rgb[*slices,:]
slices = slice(300,1300), slice(400,1600), slice(None) #[300:1300,400:1600]
rgb_crop = rgb[slices]
im = PIL.Image.fromarray((rgb_crop*255).astype('uint8')[::-1,:])
im.save("RGB_90cm_CO_ATLASGAL.png")
get_ipython().magic('run rgb_overview_big.py')
get_ipython().magic('ls -rt')
get_ipython().system('open *png')
FF = aplpy.FITSFigure(outfn)
FF.show_rgb(outfn)
import aplpy
FF = aplpy.FITSFigure(outfn)
FF.show_rgb(outfn)
FF.set_tick_labels_format('d.d','d.d')
FF.set_tick_labels_format('dd.d','dd.d')
get_ipython().system('locate pubfiguresrc')
get_ipython().system('ln /Users/adam/.matplotlib/pubfiguresrc .')
pl.matplotlib.rc_file('pubfiguresrc')
get_ipython().magic('run rgb_overview_big.py')
get_ipython().system('open *png')
get_ipython().magic('ls *png')
get_ipython().system('git add rgb_overview_big.py')
get_ipython().system('git commit -av')

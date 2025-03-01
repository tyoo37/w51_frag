import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from radio_beam import Beam
from astropy.wcs import WCS
from astropy import coordinates
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from dendrocat.aperture import Ellipse
from astropy import units as u

def cen_freq(center, delta ,alpha=2):
    start = center-delta/2 ; end = center+delta/2
    freqarr = np.logspace(np.log10(start),np.log10(end),20)
    dfreq = freqarr[1:]-freqarr[:-1]
    dfreq = np.append(dfreq,[dfreq[-1]])
    
    integral_up = np.sum(freqarr**(alpha+1)*dfreq)
    integral_down = np.sum(freqarr**(alpha)*dfreq)
    return integral_up/integral_down  

def get_flux(image,peakxy, beam1, wcsNB, pixel_scale, major,minor,pa,issky=False,savedir=None):
    fluxarr = []

    if issky:
        cen_world = peakxy
    else:
        cen_world = wcsNB.wcs_pix2world(peakxy,0)
    num_source = len(cen_world)
   
    
    for i in range(num_source):
        
        #x_cen = peakxy[i][0]
        #y_cen = peakxy[i][1]
        #cen = (x_cen, y_cen)
        """ 
        major = beam1.major
        minor = beam1.minor
        pa = beam1.pa
        area = beam1.sr.value
        
        beamarea = area
        if major_ai is not None:
            major = major_ai[i] / 3600 *u.deg
            minor = minor_ai[i] / 3600 *u.deg
            pa = pa_ai[i]*u.deg
            area = np.pi *(major.value/2*np.pi/180)*(minor.value/2*np.pi/180)
        """    
        beamarea = beam1.sr.value          
        
        positions = coordinates.SkyCoord(cen_world[i,0],cen_world[i,1], frame=wcs.utils.wcs_to_celestial_frame(wcsNB).name,unit=(u.deg,u.deg))
        
        if isinstance(major, (int, float)):
            cutout = Cutout2D(image, positions, 4.0*major, wcs=wcsNB, mode='partial')
            frame = wcs.utils.wcs_to_celestial_frame(cutout.wcs).name
            aperture = Ellipse(positions, major*pixel_scale, minor*pixel_scale, -1*pa, unit=u.deg, frame=frame)
        else:
            cutout = Cutout2D(image, positions, 4.0*major[i], wcs=wcsNB, mode='partial')
            frame = wcs.utils.wcs_to_celestial_frame(cutout.wcs).name
            aperture = Ellipse(positions, major[i]*pixel_scale, minor[i]*pixel_scale, -1*pa[i], unit=u.deg, frame=frame) # pa in degree with anti-clockwise direction
        this_mask = aperture.place(cutout.data, wcs=cutout.wcs)
        pixel_scale_sr = (pixel_scale.value * np.pi/180)**2 # pixel scale in deg^2 -> sr
        fluxarr.append(np.sum(cutout.data[this_mask]/beamarea*pixel_scale_sr)) # Jy/beam / (sr/beam) * (sr/pixel) = Jy/ pixel
        if savedir is not None:
            if not os.path.isdir(savedir):
                os.mkdir(savedir)
            else:
                fig = plt.figure(figsize=(10,10))
                ax1=fig.add_axes([0.2,0.2,0.8,0.8],projection=wcsNB)
                ax1.imshow(cutout.data, origin='lower')
                ax1.imshow(this_mask, origin='lower',alpha=0.1,cmap='gray')
                plt.savefig(savedir+'aper_%04d.png'%i)
                plt.close()
    return fluxarr
def plot_flux_comp(ax, dat_hr_b3b6, dat_lr_b3b6, ind_b3b6,
                  catdat_w51,use_rec_ind=False,
                   aperture='fwhm+beam',isfluxcorr=True,label='w51n',colors=['r','k'],marker='o',size=50):

    catdata_w51 = ascii.read(catdat_w51,data_start=0,format='commented_header', header_start=120,  comment="!")
    sky_ra = catdata_w51['WCS_ACOOR']
    sky_dec = catdata_w51['WCS_DCOOR']
    pix_x = catdata_w51['XCO_P']
    pix_y = catdata_w51['YCO_P']
    no = catdata_w51['NO']

    for i, (band, dat_hr, dat_lr, index) in enumerate(zip(['b3','b6'],dat_hr_b3b6,dat_lr_b3b6, ind_b3b6)):
        if band is 'b3':
            bandindex = 2
        elif band is 'b6':
            bandindex = 3

        fooa = catdata_w51['FOOA0%s'%bandindex]
        foob = catdata_w51['FOOB0%s'%bandindex]
        theta = catdata_w51['THETA0%s'%bandindex]
        afwhm = catdata_w51['AFWHM0%s'%bandindex]
        bfwhm = catdata_w51['BFWHM0%s'%bandindex]
        flux = catdata_w51['FXT_BST0%s'%bandindex]
        flux_g = catdata_w51['FXT_ALT0%s'%bandindex]


        rec_ind_b3 = np.where((np.abs(catdata_w51['GOODM03'])>1)&
                           (np.abs(catdata_w51['SIGNM03'])>1)&
                           (catdata_w51['FXP_BST03']/catdata_w51['FXP_ERR03']>2)&
                           (catdata_w51['FXT_BST03']/catdata_w51['FXT_ERR03']>2)&
                           (catdata_w51['AFWHM03']/catdata_w51['BFWHM03']<2)&
                           (catdata_w51['FOOA03']/catdata_w51['AFWHM03']>1.15))[0]

        rec_ind_b6 = np.where((np.abs(catdata_w51['GOODM02'])>1)&
                           (np.abs(catdata_w51['SIGNM02'])>1)&
                           (catdata_w51['FXP_BST02']/catdata_w51['FXP_ERR02']>2)&
                           (catdata_w51['FXT_BST02']/catdata_w51['FXT_ERR02']>2)&
                           (catdata_w51['AFWHM02']/catdata_w51['BFWHM02']<2)&
                           (catdata_w51['FOOA02']/catdata_w51['AFWHM02']>1.15))[0]

        rec_ind = np.where((np.abs(catdata_w51['GOODM03'])>1)&
                           (np.abs(catdata_w51['SIGNM03'])>1)&
                           (catdata_w51['FXP_BST03']/catdata_w51['FXP_ERR03']>2)&
                           (catdata_w51['FXT_BST03']/catdata_w51['FXT_ERR03']>2)&
                           (catdata_w51['AFWHM03']/catdata_w51['BFWHM03']<2)&
                           (catdata_w51['FOOA03']/catdata_w51['AFWHM03']>1.15)&
                           (np.abs(catdata_w51['GOODM02'])>1)&
                           (np.abs(catdata_w51['SIGNM02'])>1)&
                           (catdata_w51['FXP_BST02']/catdata_w51['FXP_ERR02']>2)&
                           (catdata_w51['FXT_BST02']/catdata_w51['FXT_ERR02']>2)&
                           (catdata_w51['AFWHM02']/catdata_w51['BFWHM02']<2)&
                           (catdata_w51['FOOA02']/catdata_w51['AFWHM02']>1.15))[0]

        skypos = np.vstack((sky_ra,sky_dec)).T

        fitsdata_lr = fits.open(dat_lr)
        hdr_lr = fits.getheader(dat_lr)
        wcs_lr = WCS(hdr_lr,naxis=2)
        scale_lr = wcs_lr.proj_plane_pixel_scales()[0]
        beam_lr = Beam.from_fits_header(hdr_lr)

        if aperture is 'footprint':
            major = fooa/3600
            minor = foob/3600
            pa = 180-theta
        elif aperture is 'fwhm+beam':
            meanbeamsize = (beam_lr.major.value+beam_lr.minor.value)/4
            major = (afwhm/3600+meanbeamsize)
            minor = (bfwhm/3600+meanbeamsize)
            pa = 180-theta
        elif aperture is 'beam':
            major = beam_lr.major.value
            minor = beam_lr.minor.value
            pa = 180-beam_lr.pa.value




        if use_rec_ind:
            skypos = skypos[rec_ind]
            major = major[rec_ind]
            minor = minor[rec_ind]
            pa = pa[rec_ind]
        else:
            skypos = skypos[index]
            major = major[index]
            minor = minor[index]
            pa = pa[index]
        print(major)

        image_lr = fitsdata_lr[0].data
        print(image_lr.shape)
        peakxy_lr = wcs_lr.wcs_world2pix(skypos,0) 
        if len(image_lr.shape)!=2:
            image_lr = fitsdata_lr[0].data[0][0]
        flux_lr = get_flux(image_lr, peakxy_lr, beam_lr, wcs_lr, scale_lr, major/scale_lr.value,minor/scale_lr.value,pa)
        print(flux_lr)




        fitsdata_hr = fits.open(dat_hr)
        hdr_hr = fits.getheader(dat_hr)
        wcs_hr = WCS(hdr_hr,naxis=2)
        scale_hr = wcs_hr.proj_plane_pixel_scales()[0]
        beam_hr = Beam.from_fits_header(hdr_hr)
        image_hr = fitsdata_hr[0].data
        peakxy_hr = wcs_hr.wcs_world2pix(skypos,0) 
        if len(image_hr.shape)!=2:
            image_hr = fitsdata_hr[0].data[0][0]
        print(image_hr.shape)

        flux_hr = get_flux(image_hr, peakxy_hr, beam_hr, wcs_hr, scale_hr, major/scale_hr.value,minor/scale_hr.value,pa)
        print(flux_lr)
        if isfluxcorr:
            for k in range(4):
                if hdr_lr['CTYPE%d'%(k+1)]=='FREQ':
                    hdrind_lr = k+1
                if hdr_hr['CTYPE%d'%(k+1)]=='FREQ':
                    hdrind_hr = k+1
            freq_lr =cen_freq(hdr_lr['CRVAL%d'%hdrind_lr],hdr_lr['CDELT%d'%hdrind_lr])
            freq_hr=cen_freq(hdr_hr['CRVAL%d'%hdrind_hr],hdr_hr['CDELT%d'%hdrind_hr])
            print(freq_lr, freq_hr,hdr_hr['CRVAL4'],hdr_hr['CDELT4'])
            dflux = np.log10(freq_hr/freq_lr)*2.
            flux_corr = np.log10(flux_lr) + dflux
            flux_lr = 10**flux_corr
            print(flux_lr)
        ax.scatter(flux_lr, flux_hr, color=colors[i], label='%s_%s'%(label,band), s=size,marker=marker)
        
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys
from astropy import units as u
from astropy.table import Table
from astropy import stats
import astropy.constants as c
import matplotlib as mpl
import emcee
import corner
import scipy.integrate as integrate
from scipy.optimize import minimize


W51 = '/orange/adamginsburg/w51/TaehwaYoo/'
W51b6 = '/orange/adamginsburg/w51/TaehwaYoo/2015.1.01596.S_W51_B6_LB/continuum_images/'
W51cont='/orange/adamginsburg/w51/TaehwaYoo/b6contfits/'
w51e2_b6_briggs=W51cont+'W51e2_cont_bigbriggs.image.fits'
w51e2_b6_robust0=W51cont+'W51e2_cont_big_robust0.image.fits'
w51e2_b6_uniform=W51cont+'W51e2_cont_biguniform.image.fits'
w51e2_b6_superuniform=W51cont+'W51e2_cont_bigsuperuniform.image.fits'

w51n_b6_briggs = W51cont+'W51n_cont_bigbriggs.image.fits'
w51n_b6_robust0 = W51cont+'w51n_cont_big_robust0.image.fits'
w51n_b6_uniform = W51cont+'W51n_cont_biguniform.image.fits'
w51n_b6_superuniform = W51cont+'W51n_cont_bigsuperuniform.image.fits'
w51n_b6_natural = W51cont+'W51n_cont_bignatural.image.fits'

W51b3 = '/orange/adamginsburg/w51/TaehwaYoo/2017.1.00293.S_W51_B3_LB/may2021_successful_imaging/'

w51n_b3_tt0 = W51b3+'w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
w51n_b3_tt1 = W51+'w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt1.pbcor.fits'
w51n_b3_alpha = W51+'w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.alpha.pbcor.fits'

w51conv = '/orange/adamginsburg/w51/TaehwaYoo/convolved_new/'
w51n_b6_conv = w51conv + 'w51n_cont_bigbriggs.image.convB3_briggs.fits'

w51e_b3_tt0 = W51b3+'w51e2.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
w51e2_b3_tt1 = W51+'w51e2.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt1.pbcor.fits'
w51e2_b3_alpha = W51+'w51e2.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.alpha.pbcor.fits'

w51e2_b6_conv = w51conv + 'w51e2_cont_bigbriggs.image.convB3_briggs.fits'

w51e_b6_almaimf_conv = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/w51e_B6_conv.fits'
w51e_b3_almaimf_conv = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/w51e_B3_conv.fits'


w51n_b3_almaimf_conv = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/w51n_B3_conv.fits'

w51n_b6_almaimf = '/orange/adamginsburg/w51/TaehwaYoo/w51_alma_imf/W51-IRS2_B6_uid___A001_X1296_X187_continuum_merged_12M_robust0_selfcal9_finaliter.image.tt0.pbcor.fits'

w51n_b6_conv = w51conv + 'w51n_new_nocorr_in_area_B6_conv.fits'
#w51n_b3_conv = w51conv + 'w51n_B3_conv.fits'
#w51e_b3_conv = w51conv + 'w51e_B3_conv.fits'
w51e_b6_conv = w51conv + 'w51e_new_nocorr_in_area_B6_conv.fits'


w51e_matched_catalog = '/home/t.yoo/w51/catalogue/dendrogram/dendro_w51e_matched.fits'
w51n_matched_catalog = '/home/t.yoo/w51/catalogue/dendrogram/dendro_w51n_matched.fits'

w51e_b3_noiseregion = '/orange/adamginsburg/w51/TaehwaYoo/w51e_b3_std_sky.reg'
w51e_b6_noiseregion = '/orange/adamginsburg/w51/TaehwaYoo/w51e_b6_std_sky.reg'
w51n_b3_noiseregion = '/orange/adamginsburg/w51/TaehwaYoo/w51n_b3_std_sky.reg'
w51n_b6_noiseregion = '/orange/adamginsburg/w51/TaehwaYoo/w51n_b6_std_sky.reg'
w51e_b6_calibrated_pbcor = '/orange/adamginsburg/w51/TaehwaYoo/w51e2.spw0thru19.14500.robust0.thr0.15mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
w51n_b6_calibrated_pbcor = '/orange/adamginsburg/w51/TaehwaYoo/w51n.spw0thru19.14500.robust0.thr0.1mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'


photometrydir = '/home/t.yoo/w51/catalogue/photometry/'
w51e_b3_flux = photometrydir+'w51e_b3_flux_size.fits'
w51n_b3_flux = photometrydir+'w51n_b3_flux_size.fits'
w51e_b6_flux = photometrydir+'w51e_b6_flux_size.fits'
w51n_b6_flux = photometrydir+'w51n_b6_flux_size.fits'
w51e_b6_conv_flux = photometrydir+'w51e_b6_conv_flux_size.fits'
w51n_b6_conv_flux = photometrydir+'w51n_b6_conv_flux_size.fits'

def kappa(nu, nu0=271.1*u.GHz, kappa0=0.0114*u.cm**2*u.g**-1, beta=1.75):
    """
    Compute the opacity $\kappa$ given a reference frequency (or wavelength)
    and a power law governing the opacity as a fuction of frequency:
    $$ \kappa = \kappa_0 \left(\\frac{\\nu}{\\nu_0}\\right)^{\\beta} $$
    The default kappa=0.0114 at 271.1 GHz comes from extrapolating the
    Ossenkopf & Henning 1994 opacities for the thin-ice-mantle, 10^6 year model
    anchored at 1.0 mm with an assumed beta of 1.75.
    Parameters
    ----------
    nu: astropy.Quantity [u.spectral() equivalent]
        The frequency at which to evaluate kappa
    nu0: astropy.Quantity [u.spectral() equivalent]
        The reference frequency at which $\kappa$ is defined
    kappa0: astropy.Quantity [cm^2/g]
        The dust opacity per gram of H2 along the line of sight.  Because of
        the H2 conversion, this factor implicitly includes a dust to gas ratio
        (usually assumed 100)
    beta: float
        The power-law index governing kappa as a function of nu
    """
    return (kappa0*(nu.to(u.GHz,u.spectral())/nu0.to(u.GHz,u.spectral()))**(beta)).to(u.cm**2/u.g)


def BB(freq, temp):
    B_nu = (2 * freq**3 *c.h / (c.c**2) * 1 / (np.e**(c.h*freq/(c.k_B*temp))-1))
    return B_nu
#def find_opt_thick_radius(rarr, kappa, rho):
#    
#    kappa*rho*

def get_flux_aperture(theta, rarr, freqb3, freqb6, dist=5.41*u.kpc, r_0 = 1*u.au, verbose=False ):
    rmax = rarr[-1]
    dr = rarr[1:] - rarr[:-1]
    dr = np.append(dr, rarr[-1])
    T1, T2, logrho_0, alpha = theta
    rho_0 = 10**logrho_0 *u.g / u.cm**3
    #print(rho_0)
    
    rho_r =  rho_0 * (rarr/r_0)**(-alpha)
    
    kappa_b3 = kappa(freqb3)
    kappa_b6 = kappa(freqb6)
    #print('ho',(rho_0*kappa_b3*r_0).to(u.cm/u.cm))

    au_to_cm = (1*u.au).to(u.cm)
    #print((kappa_b3*rho_0*r_0).to(u.au/u.au), alpha)
    r_b3thick = r_0 * (1+(1-alpha)/kappa_b3/rho_0/r_0)**(1/(1-alpha))
    r_b6thick = r_0 * (1+(1-alpha)/kappa_b6/rho_0/r_0)**(1/(1-alpha))

    #r_b3thick = (1-alpha + kappa_b3.value * rho_0.value * r_0.value * au_to_cm.value) / (kappa_b3.value * rho_0.value * (r_0.value * au_to_cm.value)**alpha)**(1/(1-alpha)) / au_to_cm.value
    #r_b6thick = (1-alpha + kappa_b6.value * rho_0.value * r_0.value * au_to_cm.value) / (kappa_b6.value * rho_0.value * (r_0.value * au_to_cm.value)**alpha)**(1/(1-alpha)) / au_to_cm.value
   # print('r_b3thick', r_b3thick)
   # print('r_b6thick', r_b6thick)
    
    BBthick_b3 = BB(freqb3, T1*u.K)
    BBthick_b6 = BB(freqb6, T1*u.K)
    
    BBthin_b3 = BB(freqb3, T2*u.K)
    BBthin_b6 = BB(freqb6, T2*u.K)
    
    flux_b3_arr = []
    flux_b6_arr =  []
    flux_b3 = flux_b6 = 0*u.Jy
    for i,r in enumerate(rarr):
        if verbose:
            print('verbose turned on')
            print(i,r)
        int_uplim = np.sqrt(rmax**2-r**2)
        #int_lolim_b3 = np.max((np.sqrt(r_b3thick.value**2-r**2),np.sqrt(1-r**2)))
        #int_lolim_b6 = np.max((np.sqrt(r_b6thick.value**2-r**2),np.sqrt(1-r**2)))
      
        if r< r_b6thick.value and r_b3thick.value>0 and r_b6thick.value>0:
            if r<1 and 1>r_b3thick.value:
                int_lolim_b3 = np.sqrt(1-r**2)
            else:
                int_lolim_b3 = np.sqrt(r_b3thick.value**2-r**2)
            if r<1 and r>r_b6thick.value:
                int_lolim_b6 = np.sqrt(1-r**2)
            else:
                int_lolim_b6 = np.sqrt(r_b6thick.value**2-r**2)
                
            if any((int_lolim_b3 > int_uplim,int_lolim_b6 > int_uplim)): #completely optically thick
                tau_b3 = 0
                tau_b6 = 0
                I_b3 = BBthick_b3 
                I_b6 = BBthick_b6 
            else:
                tau_b3 = (kappa_b3 * rho_0 * integrate.quad(lambda x: (np.sqrt(r**2+x**2)/r_0.value)**(-alpha), int_lolim_b3, int_uplim) *u.au).to(u.cm/u.cm)[0] 
                tau_b6 = (kappa_b6 * rho_0 * integrate.quad(lambda x: (np.sqrt(r**2+x**2)/r_0.value)**(-alpha), int_lolim_b6, int_uplim) *u.au).to(u.cm/u.cm)[0] 

                #tau_b3 = (kappa_b3 * rho_0 * r_0**alpha / (1-alpha) * ((int_uplim)**(1-alpha) - (int_lolim_b3)**(1-alpha))).to(u.cm/u.cm) 
                #tau_b6 = (kappa_b6 * rho_0 * r_0**alpha / (1-alpha) * ((int_uplim)**(1-alpha) - (int_lolim_b6)**(1-alpha))).to(u.cm/u.cm) 

                I_b3 = BBthick_b3 * np.exp(-tau_b3) + BBthin_b3 * (1-np.exp(-tau_b3))
                I_b6 = BBthick_b6 * np.exp(-tau_b6) + BBthin_b6 * (1-np.exp(-tau_b6))
            if verbose:
                print('case 1 ', tau_b3, tau_b6, I_b3, I_b6)
            
        elif r>=r_b6thick.value and r<r_b3thick.value and r_b3thick.value>0 and r_b6thick.value>0:
            if 1>r_b6thick.value:
                int_lolim_b6 = np.sqrt(1-r**2)
            else:
                int_lolim_b6 = 0
                
            if 1>r_b3thick.value:
                int_lolim_b3 = np.sqrt(1-r**2)
            else:
                int_lolim_b3 = np.sqrt(r_b3thick.value**2-r**2)
            
            
            
            if int_lolim_b3 > int_uplim:
                tau_b3=0
                I_b3=BBthick_b3
            else:
                tau_b3 = (kappa_b3 * rho_0 * integrate.quad(lambda x: (np.sqrt(r**2+x**2)/r_0.value)**(-alpha), int_lolim_b3, int_uplim) *u.au).to(u.cm/u.cm)[0]
                I_b3 = BBthick_b3 * np.exp(-tau_b3) + BBthin_b3 * (1-np.exp(-tau_b3))

            tau_b6 = (kappa_b6 * rho_0 * integrate.quad(lambda x: (np.sqrt(r**2+x**2)/r_0.value)**(-alpha), int_lolim_b6, int_uplim) *u.au).to(u.cm/u.cm)[0]
            I_b6 = BBthin_b6 * (1-np.exp(-tau_b6))
            if verbose:
                print('case 2 ', tau_b3, tau_b6, I_b3, I_b6)

     
        else:
            if 1>r_b3thick.value:
                int_lolim_b3 = np.sqrt(1-r**2)
            else:
                int_lolim_b3 = 0
            if 1>r_b6thick.value:
                int_lolim_b6 = np.sqrt(1-r**2)
            else:
                int_lolim_b6 = 0
                
            tau_b3 = (kappa_b3 * rho_0 * integrate.quad(lambda x: (np.sqrt(r**2+x**2)/r_0.value)**(-alpha), int_lolim_b3, int_uplim) *u.au).to(u.cm/u.cm)[0]
            tau_b6 = (kappa_b6 * rho_0 * integrate.quad(lambda x: (np.sqrt(r**2+x**2)/r_0.value)**(-alpha), int_lolim_b6, int_uplim) *u.au).to(u.cm/u.cm)[0]
            I_b3 = BBthin_b3 * (1-np.exp(-tau_b3))
            I_b6 = BBthin_b6 * (1-np.exp(-tau_b6))
            if verbose:
                print('case 3 ', tau_b3, tau_b6, I_b3, I_b6)

           
            
        flux_b3 = flux_b3 + 2 * np.pi * r*u.au * dr[i]*u.au * I_b3 /dist**2
        flux_b6 = flux_b6 + 2 * np.pi * r*u.au * dr[i]*u.au * I_b6 /dist**2
        if verbose:
            print('flux', flux_b3, flux_b6)
        if any((~np.isfinite(flux_b3),~np.isfinite(flux_b6))):
            print('r, r_b3thick, r_b6thick, flux_b3, flux_b6, I_b3, I_b6, tau_b3, tau_b6, rho_0, kappa_b3, kappa_b6, int_lolim_b3, int_lolim_b6, int_uplim')
            print(r, r_b3thick, r_b6thick, flux_b3, flux_b6, I_b3, I_b6, tau_b3, tau_b6, rho_0, kappa_b3, kappa_b6, int_lolim_b3, int_lolim_b6, int_uplim) 
            raise ValueError('nan in flux calc')
        #print('flux_b3',flux_b3)

        flux_b3_arr.append(flux_b3.to(u.Jy).value)
        flux_b6_arr.append(flux_b6.to(u.Jy).value)
        
        
    
    return flux_b3_arr, flux_b6_arr



def lnlike(theta, flux_b3, flux_b6, fluxerr_b3, fluxerr_b6, rarr,
           dist=5.41*u.kpc, freqb3=92982346121.91989*u.Hz, freqb6=226691598706.70853*u.Hz):
    model_b3, model_b6 = get_flux_aperture(theta, rarr, freqb3, freqb6, dist=dist)
    #print('hoho', flux_b3,model_b3)
    
    return -0.5* (np.sum((flux_b3-model_b3)**2/fluxerr_b3**2) + 
                  np.sum((flux_b6-model_b6)**2/fluxerr_b6**2))


def lnprior(theta):
    T1, T2, logrho_0, alpha = theta
    if 0 < T1 < 2000 and 0 < T2 < 2000 and -14 < logrho_0 < 5 and 1<alpha<10:
        return 0.0
    return -np.inf

def lnprob(theta, flux_b3, flux_b6, fluxerr_b3, fluxerr_b6, rarr):
    lp = lnprior(theta) 
    if not np.isfinite(lp):
        return - np.inf
    return lp + lnlike(theta, flux_b3, flux_b6, fluxerr_b3, fluxerr_b6, rarr)

def main(p0,nwalkers,niter,ndim,lnprob,data):
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)

    print("Running burn-in...")
    p0, _, _ = sampler.run_mcmc(p0, 100)
    sampler.reset()

    print("Running production...")
    pos, prob, state = sampler.run_mcmc(p0, niter)

    return sampler, pos, prob, state


def flux_gaussian_multiple_apertures(flux, major, minor, rarr): # major, minor in FWHM, au unit
    flux_rescaled = flux*major/minor # making a 2d symmetric Gaussian with FWHM = FWHM_major
    major_sigma = major/np.sqrt(8*np.log(2))
    peak_height = flux / 2 / np.pi/ major/ minor * 8 * np.log(2)
    flux_rarr =2*np.pi*major_sigma**2 * peak_height * (1 - np.exp(-rarr**2/2/major_sigma**2))
    return flux_rarr

def get_mass(rho_0, alpha):
    mass = rho_0 * u.g /u.cm**3 * 4* np.pi /(3-alpha) * (500**(3-alpha) - 1) * u.au **3
    return mass.to(u.Msun)






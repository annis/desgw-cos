import numpy as np
from scipy import interpolate
import time

#
# bayesan approach to measuring H_0 using LIGO distances
#

def posterior(h, q,  gw, gw_err, voxen) :
    from scipy.stats import norm
    from scipy.integrate import quad
    from scipy.integrate import romberg
    from scipy.integrate import simps

    gw_var = gw_err*gw_err

    # the limits to the integral are a simple model for
    # not used?
    dmin = 200.  ;# volume  and
    dmax = 600.  ;# sensitivity/ signal to noise falling

    #prior = norm.pdf(h, loc=70., scale=15.)
    prior = 1.0
    print "\t prior = ", prior,

    uid = voxen.voxelated_gals["pixels"]

    prob = 0
    for i in range(0,uid.size) :
        pix = uid[i]
        galaxies = voxen.voxelated_gals[pix]
        spatial = voxen.map_vals[pix]
        mean_v_z = voxen.voxelated_gals[pix,"mean_z"] ;# mean galaxy z inside voxel
        lum_dist = d_from_z (mean_v_z, h, q) 
            
        log_like = -((lum_dist-gw)**2)/(2*gw_var)
        like_vox = spatial * galaxies * np.e**(log_like)

        prob += like_vox.sum()
        
    post = prior * prob/np.sqrt(2*np.pi*gw_var)
    print " posterior = {:.6f}".format(post),

    return post

# luminosity distance in cosmographic terms, to second order in z
def d_from_z (z, h, q) :
    c = 3.0e5
    d = (1. + 0.5*(1.-q)*z) * z * c/h
    return d
# and its inverse, if we need it
def z_from_d (d, h, q) :
    c = 3.0e5
    z = (-1. + np.sqrt( 1. + 2.*(1.-q)*d*h/c))/(1.-q)
    return z
def testit( ) :
    z = 0.1; h=70.; q=-0.5
    d = d_from_z(z, h, q)
    z2 = z_from_d(d, h, q)
    print z, z2, z-z2
def dzdd (d, h, q) :
    c = 3.0e5
    pre = 1./(1.-q)
    chain1 = (1. + 2.*(1.-q)*d*h/c)**(-0.5)
    chain2 = 2.*(1.-q)*h/c
    deriv = pre * chain1*chain2
    return deriv


# =================================================
#
# now one wants to calculate the full posterior
#
# =================================================

def calc_posterior (
        q, gw, gwerr, 
        map_ra, map_dec, map_vals, 
        gal_ra, gal_dec, gal_zed, gal_zerr):
    import voxel
    timer= True

    linear_distance=7.36
    voxen = voxel.voxel(70., q, linear_distance)
    voxen.set_map(map_ra, map_dec, map_vals)
    if timer: 
        start = time.time()

    delH = 1; prior_Ho_low = 50; prior_Ho_high = 100
#    delH = 3
#    delH = 10;  prior_Ho_low = 65; prior_Ho_high = 85
    Ho, density =np.array([]),np.array([])
    for h in range(prior_Ho_low, prior_Ho_high+1, delH) :
        print h,
        voxen.reset(h, q, linear_distance)
        voxen.set_galcat(gal_ra, gal_dec, gal_zed, gal_zerr)
        #if (np.mod(h+1, 10) == 0) : print ""

        # evaluation of the probability of the data given h (eq 11)
        prob = posterior(h, q, gw, gwerr, voxen)

        # bookkeeping
        Ho = np.append(Ho, h)
        density = np.append(density, prob)
        if timer: print "\t {:.1f} sec".format(time.time()-start),
        print ""
    density = density/density.sum()
    if timer: print "total time: {:.1f} sec".format(time.time()-start)
    return Ho, density

def do_sdss(map_ra, map_dec, map_vals, gw=410., gwerr = 170., \
        target_ra=0, target_dec=0, inj_ra=0, inj_dec=0, q=-1) :
    import ho_prep_data

    ra, dec, zed, zerr = ho_prep_data.getRotatedSDSSCat( 
        target_ra=target_ra, target_dec=target_dec, 
        inj_ra=inj_ra, inj_dec=inj_dec)

    Ho, density = calc_posterior (
        q, gw, gwerr, map_ra, map_dec, map_vals, ra, dec, zed,zerr)
    return Ho, density



import numpy as np
import healpy as hp
from scipy import interpolate
#t_gal_zed, t_voxel_limits_z, t_unique_pix, t_pix = [],[],[],[]

# ====================================================================
#
# voxel technology
#
# ====================================================================

class voxel(object) :
    """
    Constuct the voxels.

    Algorithm:
        Given a linear size construct a volume from linear_distance**3.
        Assume a nsides=64 healpy map via assumeing fraction_sky=2e-5.
        Work out from z=0 to z=0.2 building constant volume voxels
        by using brentq to solve for the zero crossings of
            (d_vol_interp(z)-d_vol_last) - voxel_vol = 0
        where d_vol_last was the volume out to the end of the previous voxel
        and d_vol_interp is a interp1d interpolation object fit to
        the volume element (which assumes a spatially flat cosmology):
            frac_sky * (4pi/3) d_c**3
        and d_c is the line of sight comoving distance in a cosmographic expansion:
            dc = (1. - 0.5*(1.+q)*z) * z * c/h

        Prior to setting the galaxy catalog one must set the astropy map using set_map.
    
    Object methods:
        voxel.voxel(h,q,linear_distance, fraction_sky)
        voxel.build_voxels()
        voxel.set_map(ra, dec, vals)
        voxel.set_galcat(ra, dec, z, zerr) - this produces the voxelated map 
        voxel.reset(h,q,linear_distance)   - this resets the cosmology and rebuilds the voxels

    Usage:
        Once the voxels have been constructed by voxel.build_voxels(),
        the voxel edges are known in both redshift and comoving distance.
            voxel.voxel_limits_z
            voxel.voxel_limits_dc
        Once the galaxy map has been set, the dictionary "voxelated_gals" contains
        the voxels  as "[pix]" where pix is the pixel number from "[pixels]"-
        this is a vector of length voxel_limits_z.size and 
        whose outer edge is given by voxel_limits_z and voxel_limits_dc

    Notes:
        linear_distance = 6 for size of nsize=64 pixels
        linear_distance = 7.36 for pixels of 400 Mpc^3 and start z ~0.04 (0.0397)
    """
    def __init__(self, h, q, linear_distance=7.36, fraction_sky=2e-5, verbose=False):
        self.verbose      = verbose
        self.h            = h   ;# Ho for this voxelization
        self.q            = q   ;# qo for this voxelization
        self.linear_distance = linear_distance  ;# V=l^3, so l sets the voxel volume
        self.fraction_sky = fraction_sky

        self.voxel_z      = 0
        self.voxel_del_vol = 0   ;# these 3 are debugging vectors set by build_voxel
        self.voxel_del_d   = 0

        self.voxel_limits_z = 0 ;# these are build_voxel constructed voxel limits
        self.voxel_limits_dc= 0

        self.map_vals       = "" ;# ra,dec, vals of a astropy map
        self.map_ra         = ""
        self.map_dec        = ""

    def reset(self, h, q, linear_distance) :
        self.h            = h   ;# Ho for this voxelization
        self.q            = q   ;# qo for this voxelization
        self.linear_distance = linear_distance  ;# V=l^3, so l sets the voxel volume
        self.build_voxels()

    # one sets the map ra, dec
    def set_map (self, ra, dec, map) :
        self.map_vals = map
        self.map_ra = ra
        self.map_dec = dec

    # one sets the galaxy catalog
    def set_galcat (self, gal_ra, gal_dec, gal_zed, gal_zerr) :
        #global t_gal_zed
        #global t_voxel_limits_z
        #global t_unique_pix
        #global t_pix
        map_ra = self.map_ra
        map_dec = self.map_dec
        map_vals = self.map_vals
        voxel_limits_z = self.voxel_limits_z

        ix = gal_ra > 180
        gal_ra[ix] = gal_ra[ix]-360
        self.gal_ra  = gal_ra
        self.gal_dec  = gal_dec
        self.gal_zed = gal_zed
        self.gal_zerr = gal_zerr
        # will have to smooth these out, eventually

        nsides = hp.get_nside(map_ra)
    
        phi = gal_ra*2*np.pi/360.;
        theta = (90-gal_dec)*2*np.pi/360.
        pix = hp.ang2pix(nsides,theta,phi)

        ix = map_vals[pix] > 1e-9
        pix=pix[ix]
        gal_ra = gal_ra[ix]
        gal_dec = gal_dec[ix]
        gal_zed = gal_zed[ix]
        gal_zerr = gal_zerr[ix]
        unique_pix = np.unique(pix)

        voxelated_gals = dict()
        voxelated_gals["pixels"] = unique_pix

        print "{} gals over {} pixels each w/ {} voxels".format(
            gal_ra.size, unique_pix.size,voxel_limits_z.shape[0]),
        for i in range(0, unique_pix.size) :
            # find every galaxy in the given pixel
            upix = unique_pix[i]
            ix = upix == pix
            voxels = np.zeros(voxel_limits_z.shape[0])
            voxels_mean_z = np.zeros(voxel_limits_z.shape[0])
            for j in range(0,voxel_limits_z.shape[0]) :
                zed_1,zed_2 = voxel_limits_z[j]
                # fill out the voxels of this pixel
                ix2 = (gal_zed[ix]  >= zed_1) & (gal_zed[ix] < zed_2)
                voxels[j] = ( gal_zed[ix][ix2]).size
                if ( gal_zed[ix][ix2]).size > 0 :
                    voxels_mean_z[j] = np.median( gal_zed[ix][ix2] )
            # save the voxels labeled by pixel,z_bin
            voxelated_gals[upix] = voxels
            voxelated_gals[upix,"mean_z"] = voxels_mean_z

        self.voxelated_gals = voxelated_gals

        #print "I am here and starting mapping"
        #t_gal_zed, t_voxel_limits_z, t_unique_pix, t_pix = \
                #gal_zed, voxel_limits_z, unique_pix, pix
        #voxels = map(maploop, range(0,unique_pix.size))
        #print len(voxels),len(voxels[0])

        return 


    #
    # linear_distance = 6 for size of nsize=64 pixels
    # linear_distance = 7.36 for pixels of 400 Mpc^3 and start z ~0.04 (0.0397)
    #
    #   voxel_z carries the z of the outer edge of the shell
    #
    def build_voxels (self) :
        from scipy.optimize import brentq
        verbose = self.verbose
        h = self.h
        q = self.q
        linear_distance = self.linear_distance
        fraction_sky = self.fraction_sky

        voxel_vol = linear_distance**3
        del_vol_interp, first_z = self.volume_samples(h, q, fraction_sky, voxel_vol)
        print "\t ---voxel vol {:d} Mpc^3 => initial z={:.3f}  ".format( int(np.round(voxel_vol)), first_z),
    
        # testing outputs
        # e.g.,  z,v,d, vox_z,vox_dc =ho_measure.build_voxels()
        #        plt.plot(z,v), plt.scatter(z,d), plt.semilogy(z[:-1],z[1:]-z[:-1])
        voxel_z = np.array([0])
        voxel_del_vol=np.array([])
        voxel_del_d=np.array([])
    
        # voxel coordinates
        voxel_limits_z = np.array([]).reshape(0,2)
        voxel_limits_dc = np.array([]).reshape(0,2)
    
        d_vol_last = 0
        z_last = first_z
        zmax = 0.200
        zm_del = -0.01
        zp_del = 0.015
        while voxel_z[-1] < zmax :
            # print "a",self.delta_volume(z_last+zm_del, del_vol_interp, d_vol_last, voxel_vol)
            # print "b",self.delta_volume(z_last+zp_del, del_vol_interp, d_vol_last, voxel_vol)
            if verbose: print "="
            new_z = brentq(self.delta_volume, 
                z_last+zm_del, z_last+zp_del, args=(del_vol_interp, d_vol_last, voxel_vol),
                xtol=1e-12)
    
            voxel_limits_z = np.vstack([voxel_limits_z, [voxel_z[-1], new_z]])
            voxel_limits_dc = np.vstack([voxel_limits_dc,  
                [self.comoving_distance(z_last),self.comoving_distance(new_z) ]])
     
    
            voxel_z = np.append(voxel_z, new_z)
            voxel_del_vol = np.append(voxel_del_vol, 
                    del_vol_interp(new_z)-d_vol_last)
            voxel_del_d = np.append(voxel_del_d, 
                self.comoving_distance(new_z) - self.comoving_distance(z_last) )
    
            d_vol_last = del_vol_interp(new_z)
            z_last = new_z
    
        voxel_z = voxel_z[1:]
        voxel_limits_dc[0][0] = 0.0

        self.voxel_z = voxel_z
        self.voxel_del_vol = voxel_del_vol
        self.voxel_del_d = voxel_del_d

        self.voxel_limits_z = voxel_limits_z
        self.voxel_limits_dc = voxel_limits_dc
        return 
    
    # a function designed to compare the volume in a shell to a fiducial volume
    def delta_volume (self, z, d_vol_interp, d_vol_last, voxel_vol) :
        verbose = self.verbose
        if verbose: print "\t ",z,"vol=",d_vol_interp(z),
        value = (d_vol_interp(z)-d_vol_last) - voxel_vol
        if verbose: print "val=",value, "d=",self.comoving_distance(z)
        return value

    # build a volume interpreter and give its first valid z 
    def volume_samples (self, h, q, fraction_sky = 1.0, voxel_vol=0) :
        dz = 0.001
        zeds = np.arange(0.0001,0.251,dz)
        d_c = self.comoving_distance(zeds)
        d_vol = fraction_sky*(4*np.pi/3.)*(d_c**3 )
        d_vol_interp = interpolate.interp1d(zeds,d_vol)
    
        zeds = np.arange(0.0001,0.25,dz/10.)
        d_vol = d_vol_interp (zeds)
        first = np.argmax(d_vol>voxel_vol)
        first_z = zeds[first]
     
        return d_vol_interp, first_z

    # dc  distance_commoving
    def comoving_distance(self, z) :
        h= self.h
        q = self.q
        c = 3.0e5
        d = (1. - 0.5*(1.+q)*z) * z * c/h
        return d

# end voxel technology
# ====================================================================
    

#def maploop(i,unique_pix,pix,gal_zed,voxel_limits_z):
def maploop(i) :
    global t_gal_zed
    global t_voxel_limits_z
    global t_unique_pix
    global t_pix
    gal_zed, voxel_limits_z, unique_pix, pix = \
                t_gal_zed, t_voxel_limits_z, t_unique_pix, t_pix
    upix = unique_pix[i]
    ix = upix == pix
    voxels = np.zeros(voxel_limits_z.shape[0])
    for    j in range(0,voxel_limits_z.shape[0]):
        zed_1,zed_2 = voxel_limits_z[j]
        ix2 = (gal_zed[ix] >= zed_1) & (gal_zed[ix] < zed_2)
        voxels[j] = (gal_zed[ix][ix2]).size

    return voxels


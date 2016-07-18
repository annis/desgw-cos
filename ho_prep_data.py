import numpy as np
import fitsio

#
# bayesan approach to measuring H_0 using LIGO distances
#
# =====================================================================
#
# routines to work the first two year ligo maps
#
# =====================================================================
#
# These two routines are bypassed in ho_sim
#
# mra,mdec,map = ho_measure.getMap()
def getMap (map="LALInference_skymap.fits") :
    import hp2np
    dir = "/data/des30.a/data/annis/des-gw/lmc-event/512/"
    map = dir+map
    ra,dec,vals = hp2np.hp2np(map, degrade=False, fluxConservation=True)
    return ra,dec,vals

def getSimMap(sim=1087, bay=False, year=2016) :
    import hp2np
    if bay :
        if year == 2015:
            dir = "/data/des30.a/data/annis/des-gw/ligo/sims/"
        elif year == 2016:
            dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016/"
        file = "bayestar-{}.fits.gz".format(str(sim))
    else :
        dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-la/"
        file = "lalinference-{}.fits.gz".format(str(sim))
    map = dir+file
    ra,dec,vals = hp2np.hp2np(map, degrade=False, fluxConservation=True)
    return ra,dec,vals

#
# =====================================================================
#
#   routines to work the sdss catalog
#
# =====================================================================
#
def getSDSSCat (dir="./data/") :
    cat = "sdss.csv" 
    cat = dir+cat
    ra,dec,imag,gi,zed = np.genfromtxt(cat, unpack=True,comments="#", delimiter=",")
    ix = zed >= 0.02
    ra = ra[ix]; dec = dec[ix]
    imag = imag[ix]; gi = gi[ix]
    zed = zed[ix]
    dm = 5*np.log10(zed*3e5/70)+25
    abs = imag-dm
    ix = (abs < -21.8)
    ra = ra[ix]; dec = dec[ix]; imag=imag[ix]; zed=zed[ix];  gi=gi[ix]
    print "HACK ho_prep_data"
    ix=(ra>160)&(ra<180)&(dec>10)&(dec<30)  ;# constant surface density
    ix=(ra>100)&(ra<260)&(dec>-2)&(dec<70)  ;# the bulk of the north galactic cap
    ra = ra[ix]; dec = dec[ix]; imag=imag[ix]; zed=zed[ix];  gi=gi[ix]
    ix = ra > 180
    ra[ix] = ra[ix] -360.
    #return ra,dec,imag,gi,zed
    return ra,dec,zed

def getRotatedSDSSCat(target_ra, target_dec, inj_ra, inj_dec, 
        gal_cat_dir="./data/") :
    ra,dec,zed = getSDSSCat(dir=gal_cat_dir)
    zerr = 0.0033*np.ones(zed.size)
    ra, dec = rotate_galaxies( ra, dec, target_ra, target_dec, inj_ra, inj_dec) 
    return ra,dec,zed,zerr


#===========================
#
# rotate the galaxies onto a new place on the sky
#   say, at a LIGO sim map position
#
#===========================
#
# target_ra,dec is the chosen galaxy ra, dec
# inj_ra,dec is the sim map injection point
#
#   returns the new ra, dec of the galaxy catalog
#
def rotate_galaxies (vgal_ra, vgal_dec, target_ra, target_dec, inj_ra,  inj_dec) :
    import healpy as hp
    import rotate
    # rotate the galaxy coordinates to that place on the map
    # the way this works is that first the target galaxy is rotated onto the origin, bringing
    # its friends with it, then the origin is shifted onto the sim map position.
    ra, dec = rotate.rotateSkyToOrigin(vgal_ra, vgal_dec, target_ra, target_dec)
    ra, dec = rotate.rotateSkyAwayFromOrigin(ra, dec, inj_ra, inj_dec)
    return ra, dec

# ===========================================================
#
# obsolete
#


#===========================
#
# generic routines
#
#===========================

#
# target_ra,dec is the chosen galaxy ra, dec
#   we'll rotate the sim map to this coordinate-
#
#
def local_prob (vgal_ra, vgal_dec, map_ra, map_dec, map_val, target_ra, target_dec, inj_ra,  inj_dec) :
    import healpy as hp
    import rotate
    # rotate the galaxy coordinates to that place on the map
    # the way this works is that first the target galaxy is rotated onto the origin, bringing
    # its friends with it, then the origin is shifted onto the sim map position.
    ra, dec = rotate.rotateSkyToOrigin(vgal_ra, vgal_dec, target_ra, target_dec)
    ra, dec = rotate.rotateSkyAwayFromOrigin(ra, dec, inj_ra, inj_dec)
    # get the weights for that map. this works because the galaxy vector remains ordered.
    nsides = hp.get_nside(map_val)
    phi = ra*2*np.pi/360.;
    theta = (90-dec)*2*np.pi/360.
    pix = hp.ang2pix(nsides,theta,phi)
    local = map_val[pix]
    ix = np.invert(np.isnan(local)) & ( local > 3e-7)  # 99% inside GW map
    return local, ix

# test patten
def testit () :
    a=np.array([45.5,42,48]);b=np.array([45.5,48,42])
    xra,xdec=rotate.rotateSkyToOrigin(ra,dec,rot_ra[i],rot_dec[i]);
    xra,xdec=rotate.rotateSkyAwayFromOrigin(xra,xdec,a[0],b[0]);
    plt.clf();plt.hexbin(xra,xdec,vals);plt.scatter(a[0],b[0],c="w",s=70); plt.scatter(a[1:3],b[1:3],c="w")
    dist1=(xra-a[0])**2+(xdec-b[0])**2; dist2=(xra-a[1])**2+(xdec-b[1])**2; dist3=(xra-a[2])**2+(xdec-b[2])**2; 
    ix1=np.argmin(dist1);ix2=np.argmin(dist2);ix3=np.argmin(dist3);print vals[ix1],vals[ix2],vals[ix3]
    xra,xdec=rotate.rotateSkyToOrigin(a,b,a[0],b[0]); 
    xra,xdec=rotate.rotateSkyAwayFromOrigin(xra,xdec,rot_ra[i],rot_dec[i]); 
    plt.clf();plt.hexbin(ra,dec,vals);plt.scatter(xra,xdec,c="w")
    dist1=(ra-xra[0])**2+(dec-xdec[0])**2; dist2=(ra-xra[1])**2+(dec-xdec[1])**2; dist3=(ra-xra[2])**2+(dec-xdec[2])**2; 
    ix1=np.argmin(dist1);ix2=np.argmin(dist2);ix3=np.argmin(dist3);print vals[ix1],vals[ix2],vals[ix3]

#   if spatial_prob = true, use a random number to decide which
#   pixel in the spatial localization map the target_ra,dec corresponds to,
#   as opposed to using, say, the maximum probability pixel.
#   There is no "center", per se, for the spatial localization map.
# 
def local_prob_mod (map_ra, map_dec, map_val, target_ra, target_dec, spatial_prob=True) :
    # choose where the "center" of the sim map is
    spatial_prob = True
    sim_ra, sim_dec = whereInMap(map_ra, map_dec, map_val, spatial_prob=spatial_prob)
    # calculate angle to rotate galaxies to that sim map "center"
    rot_ra_angle = sim_ra-target_ra
    rot_dec_angle = sim_dec-target_dec
    if rot_ra_angle < 0 : rot_ra_angle += 360.
    if rot_ra_angle > 360 : rot_ra_angle -= 360.
    if rot_dec_angle < -90 : rot_dec_angle += 90
    if rot_dec_angle > 90 : rot_dec_angle -= 90
    # rotate tha galaxy coordinates to that center
    #       no -1 for alpha, beta (see rotate.py) as we're doing the inverse rotation, gals onto sim
    ra, dec = rotate.rotateSky(vgal_ra, vgal_dec, rot_ra_angle, rot_dec_angle)

def whereInMap(ra,dec,map,spatial_prob = True) :
    ix = np.argsort(map)[::-1] ;# descending order
    if not spatial_prob :
        out_ra = ra[ix][0]
        out_dec = dec[ix][0]
    else :
        which = np.random.uniform(0,1)
        wix = np.argmax(map[ix]/map.sum() >= which)
        out_ra = ra[ix][wix]
        out_dec = dec[ix][wix]
    return out_ra, out_dec




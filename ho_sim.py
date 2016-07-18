import numpy as np
import ho_prep_data
import ho_measure

#=================================================
#
# RUN the LIGO Simulations!!!
#
#=================================================

def veniVidi (
        output="maps-2016-la.txt", 
        pickleDir = "d_pdf/", 
        year=2016, 
        bay=False, 
        derr = 0.15, 
        targetHo=66.,
        targetQo=-1.,
        grid=False,
        data_catalog="sdss") :
    locationDict =  makeLocationDict (metaMapFile=output, pickleDir=pickleDir, 
        year=year, bay=bay, distance_error_scale=derr, targetHo=targetHo, targetQo=targetQo,
        grid=grid, data_catalog=data_catalog) 
    sims, mjds, distances, inj_ra, inj_dec =veni(locationDict)
    vidi(sims, mjds, distances, inj_ra, inj_dec, locationDict)
    return locationDict

def nowVici (
        output="maps-2016-la.txt", 
        pickleDir = "d_pdf/", 
        year=2016, 
        bay=False, 
        derr = 0.15, 
        targetHo=66.,
        targetQo=-1.,
        grid=False,
        data_catalog="sdss") :
    locationDict =  makeLocationDict (metaMapFile=output, pickleDir=pickleDir, 
        year=year, bay=bay, distance_error_scale=derr, targetHo=targetHo, targetQo=targetQo,
        grid=grid, data_catalog=data_catalog) 
    vici(locationDict)

def veniVidiVici (
        output="maps-2016-la.txt", 
        pickleDir = "d_pdf/", 
        year=2016, 
        bay=False, 
        derr = 0.15, 
        targetHo=66.,
        targetQo=-1.,
        grid=False,
        data_catalog="sdss") :
    import os

    veniVidi(output, pickleDir, year, bay, derr, targetHo, grid)
    locationDict =  makeLocationDict (metaMapFile=output, pickleDir=pickleDir, 
        year=year, bay=bay, distance_error_scale=derr, targetHo=targetHo, targetQo=targetQo,
        grid=grid, data_catalog=data_catalog) 
    vici(locationDict)

def makeLocationDict (metaMapFile="maps.txt", pickleDir="d_pdf/", 
        year=2016, bay=True, distance_error_scale=0.20, targetHo=66., targetQo=-1, grid=False,
        data_catalog="sdss") :
    import os
    pickleDir = os.path.join(pickleDir,"")

    bay15 = "sims/"
    bay16 = "sims-2016/"
    bay16la = "sims-2016-la/"
    masterDir  = "/data/des30.a/data/annis/des-gw/ligo/"
    if grid: 
        bay15 = "sims-2015/"
        masterDir  = "../"

    if year==2015 and not grid:
        data_dir = masterDir + bay15
        simsFile = masterDir + "2015_inj.txt"
    if year==2016 and bay and not grid:
        data_dir = masterDir + bay16
        simsFile = masterDir + "2016_inj.txt"
    if year==2016 and not bay and not grid:
        data_dir = masterDir + bay16la
        simsFile = masterDir + "2016_inj.txt"

    if year==2015 and grid:
        data_dir = masterDir + bay15
        simsFile = masterDir + bay15 + "2015_inj.txt"
    if year==2016 and bay and grid:
        data_dir = masterDir + bay16
        simsFile = masterDir + bay16 + "2016_inj.txt"
    if year==2016 and not bay and grid:
        data_dir = masterDir + bay16la
        simsFile = masterDir + bay16la + "2016_inj.txt"


    simFileFormat = "bayestar-{:d}.fits.gz"
    stripStringDir = data_dir + "bayestar-"
    if not bay :
        simFileFormat = "lalinference-{:d}.fits.gz"
        stripStringDir = data_dir + "lalinference-"

    data = dict()
    data["year"] = year
    data["simsFile"] = simsFile
    data["data_dir"] = data_dir
    data["simFileFormat"] = simFileFormat
    data["metaMapFile"] = metaMapFile
    data["stripStringDir"] = stripStringDir
    data["pickleDir"] = pickleDir
    data["distance_error_scale"] = distance_error_scale
    data["area_scale"] = 0.90
    data["targetHo"] = targetHo
    data["targetQo"] = targetQo
    data["catalog"] = data_catalog

    return data
#
# I need the list of simulations to do
#
def veni( locationDict ) :
    simsFile = locationDict["simsFile"]
    sims, mjds, distances, inj_ra, inj_dec = np.genfromtxt(simsFile, unpack=True, skip_header=40, usecols=(0,2,8,3,4))
    sims = sims.astype("int")
    return sims, mjds, distances, inj_ra, inj_dec
#
# I'll need to compute a gravitational wave event distance and error
# So, we rewrite a new master file that has the location of the sim file
#   and several metadata entries, like
#    file, distance, sigma, ra, dec, z, rot_ra, rot_dec
#
def vidi(sims, mjds, distances, inj_ras, inj_decs, locationDict, verbose=False) :
    import os.path

    pickleDir = locationDict["pickleDir"] 
    output = pickleDir + locationDict["metaMapFile"]

    data_dir = locationDict["data_dir"]
    area_scale = locationDict["area_scale"]
    distance_error_scale = locationDict["distance_error_scale"]
    file = locationDict["simFileFormat"]
    targetHo = locationDict["targetHo"] 
    targetQo = locationDict["targetQo"] 
    catalog = locationDict["catalog"] 

    if not os.path.exists(pickleDir): 
        print "=== making ", pickleDir
        os.makedirs(pickleDir)

    # what relation does catalog bear to getSDSSCat ? just a label?
    # yes. And how does sdss cat get to masterDir on the grid?
    # clearly getSDSSCat wants a real place to look for a catalog
    #   gal_ra,gal_dec,gal_zed = ho_prep_data.getSDSSCat("/home/s1/annis/daedalean/desgw-cos/data/")
    # the answer is that it is always in "data/"
    gal_ra,gal_dec,gal_zed = ho_prep_data.getSDSSCat()

    fd = open(output, "w")
    fd.write("# galaxy catalog= {}\n".format(catalog))
    fd.write("# file, distance, sigma, target_ra, target_dec, target_z, inj_ra, inj_dec \n")
    fd.close()
    counter = 0
    for sim, mjd, distance, inj_ra, inj_dec in zip(sims,mjds,distances,inj_ras,inj_decs) :

        simfile = file.format(sim)

        simfile = data_dir+simfile
        if not os.path.isfile(simfile) : 
            if verbose: print "skipping sim ","not existent  ================================"
            continue
        if os.path.getsize(simfile) == 0 : 
            if verbose: print "skipping sim ","empty ================================"
            continue
        print simfile, 
        print "\t\t ",

        ra, dec, zed = pickGalaxy(gal_ra, gal_dec, gal_zed)
        distance = ho_measure.d_from_z (zed, targetHo, targetQo) 
        distance_sigma = distance_error_scale*distance
        # reample with mean = dist and sigma = err
        distance2 = np.random.normal(distance, distance_sigma)
        if distance2 < 200 or distance2 > 600 :
            distance2 = np.random.normal(distance, distance_sigma)
        distance = distance2

        fd = open(output, "a")
        fd.write("{} {:.1f} {:.2f} {:.6f} {:.5f} {:.3f} {:.4f} {:.4f}\n".format(
            simfile, distance, distance_sigma, ra, dec, zed, inj_ra, inj_dec))
        fd.close()
        print "id, distance: {}    {:.1f} ".format( sim, distance)
        #print "skipping sim ", sim, distance, " ================================"
        counter +=1
    print "n= ",counter

#
# Now I need to measure H_o for each simmap
#
def vici (locationDict, run_qo=-1, start=0,end=0) :
    import hp2np
    import ho_measure
    import cPickle as pickle
    pickleDir = locationDict["pickleDir"] 
    pickle.dump(locationDict, open(pickleDir+"locationDict", "wb"))
    infile = pickleDir + locationDict["metaMapFile"]
    data_catalog = locationDict["catalog"]

    gw,gwerr = np.genfromtxt(infile, unpack=True, skip_header=1, usecols=(1,2))
    target_ra,target_dec = np.genfromtxt(infile, unpack=True, skip_header=1, usecols=(3,4))
    inj_ra,inj_dec = np.genfromtxt(infile, unpack=True, skip_header=1, usecols=(6,7))
    file = np.genfromtxt(infile, unpack=True, skip_header=1, dtype="str", usecols=0)
    if end == 0: end = gw.size
    for i in range(start, end) :
        print file[i], gw[i]
        repDir = locationDict["stripStringDir"] 
        simfile = file[i]
        sim = simfile.replace(repDir,"")
        sim = sim.replace(".fits.gz","")
        sim = int(sim)
        map_ra,map_dec,map_vals = hp2np.hp2np(simfile, degrade=64, fluxConservation=True)
        if data_catalog == "sdss" :
            h,d = ho_measure.do_sdss(map_ra, map_dec, map_vals, gw[i], gwerr[i], \
                target_ra=target_ra[i], target_dec=target_dec[i], \
                inj_ra=inj_ra[i], inj_dec=inj_dec[i], q=run_qo) 
        
        pickle.dump([h,d,sim,run_qo], open(pickleDir+str(i), "wb"))

        # now, for debugging purposes make it easy to rerun this without loading map
        if end-start == 1:
            print "ho_measure.do_sdss(map_ra,map_dec,map_vals,", gw[i], ",",gwerr[i], ",", target_ra[i],",", target_dec[i],",",  inj_ra[i],",",  inj_dec[i],",",  run_qo,")"
            return map_ra,map_dec,map_vals

def runOne(locationDict, i) :
    map_ra,map_dec,map_vals= vici(locationDict, start=i, end=i+1)
    return map_ra,map_dec,map_vals


# replace with marcelle code- it will likely take the galaxy catalog
# and return the index
def pickGalaxy (gra, gdec, gzed) :
    return 171.,20.,0.1

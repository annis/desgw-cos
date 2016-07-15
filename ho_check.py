import numpy as np
from scipy import interpolate
import ho_prep_data
import cPickle as pickle

# import os; import ho_check; import ho_measure; import ho_sim
# simNum=503473; cwd=os.getcwd();cwd2="/data/des30.a/data/annis/des-gw/cosmo/Runs"
# os.chdir(cwd); reload(ho_check); os.chdir(cwd2)
# locDic=ho_check.getLocationDict("d_pdf_2016la_sdss_66/")
# ho_check.findOne(locDic["pickleDir"], 50,simNum)
# ho_check.check(locDic,0,30,include=np.array([7,]))
# ho_check.check_one(locDic["pickleDir"],7)
# ho_check.rho_check(7,simNum, 65., locDic["pickleDir"])

def getLocationDict(dir) :
    locationDict = pickle.load(open(dir + "locationDict", "rb"))
    return locationDict
    #data["year"] = year
    #data["simsFile"] = simsFile
    #data["data_dir"] = data_dir
    #data["simFileFormat"] = simFileFormat
    #data["metaMapFile"] = metaMapFile
    #data["stripStringDir"] = stripStringDir
    #data["pickleDir"] = pickleDir
    #data["distance_error_scale"] = distance_error_scale
    #data["area_scale"] = 0.90
    #data["targetHo"] = targetHo

#   dist,sigma,ng,nw,area = ho_check.meta(locDict)
def meta(locDict) :
    file = locDict["pickleDir"] + locDict["metaMapFile"]
    dist,sigma,ng,nw,area  = np.genfromtxt(file, unpack=True, skiprows=1, usecols=(1,2,5,6,7))
    return dist,sigma,ng,nw,area

def checkPlot(dir, start, end, include="", exclude="") :
    import cPickle as pickle
    import matplotlib.pyplot as plt
    verbose = True
    locDict = pickle.load(open(dir+"locationDict","rb"))
    dir = locDict["pickleDir"]
    if exclude != "": skip = exclude.tolist()
    if include != "": get  = include.tolist()
    first = True
    h,sig,icount = np.array([]),np.array([]), np.array([])
    counter = 0
    for i in range(start,end+1) :
        if exclude != "": 
            if i in skip : continue
        if include != "": 
            if i not in get  : continue
        data=pickle.load(open("{}{}".format(dir, str(i)),"rb"))
        if first == True:
            pdf = np.ones(data[1].size)
            distance = data[0]
            first = False
        old_pdf = pdf
        pdf = pdf * data[1]
        pdfn = pdf/pdf.sum()
        pdf = pdf/pdf.sum()
        if np.any(np.isnan(pdf)) :
            pdf = old_pdf
            print i," has NANs, skipping"
            continue
        try :
            print i,
            amplitude,mean,sigma=fit_curve(distance,pdfn,0.16,70,50, verbose=verbose); 
        except :
            print "boom"
            pdf = old_pdf 
            continue
        h = np.append(h, mean)
        sig = np.append(sig, sigma)
        icount = np.append(icount, counter)
        counter += 1
    plt.clf()
    plt.plot(icount, h)
    plt.errorbar(icount,h,yerr=sig,fmt=None)
    xmin,xmax=icount.min(),icount.max()
    plt.ylim(50,100)
    plt.xlim(xmin,xmax)
    plt.plot([xmin,xmax],[h[-1]+2, h[-1]+2],c="0.33")
    plt.plot([xmin,xmax],[h[-1]-2, h[-1]-2],c="0.35")
    plt.title("gray lines at $\pm$ 2 km/s/Mpc.")
    xfractional = 100*np.abs(sig[-1]/h[-1])
    plt.annotate( '{:.0f} $\pm${:.1f} km/s/Mpc     a {:.0f}% measurement'.format( 
        np.abs(h[-1]),np.abs(sig[-1]), xfractional),
        xy=(0.25, 0.90), xycoords='axes fraction')
    return icount, h, sig


# n=46; ix = np.nonzero((dist <= 500) & (dist >= 200) & (i<=n))[0]; print ix.size; ho_measure.check(n+1,include=ix)
def check (locDict, start, end, include="", exclude="") :
    import matplotlib.pyplot as plt
    dir = locDict["pickleDir"]
    end = end+1
    if include == "" :
        if exclude == "":
            h,d=getall(dir, start, end); 
        else :
            h,d=getall(dir, start, end, exclude_ix=exclude); 
    elif include != "":
        h,d=getall(dir, start, end, include_ix=include); 
    plt.clf();plt.plot(h,d)
    try :
        a,b,c=fit_curve(h,d,0.16,70,50); 
        plt.plot(h,gaussian(h, a,b,c),c="green")
        plt.plot(h,gaussian(h, a,b,15.),c="orange")
    except :
        a,b,c = d.max(), 70., 30.
        plt.plot(h,gaussian(h, a,b,c),c="red")
        print "blammo"
    max=d.max(); ii=d==max; 
    print a,b,c,"  max=", h[ii][0]

# n=5; sim,max_h =ho_measure.check_one(n)
# ho_measure.rho_check(n,sim, max_h)
def check_one(dir, n) :
    import matplotlib.pyplot as plt
    print n, "{}{}".format(dir,n)
    h,d,sim = pickle.load(open("{}{}".format(dir,n),"rb"))
    if np.any(np.isnan(d)) : raise Exception ("d has NANs")
    max=d.max(); ii=d==max; max_h = h[ii][0]
    plt.clf();plt.plot(h,d)
    try :
        a,b,c=fit_curve(h,d,0.16,70,50); 
        plt.plot(h,gaussian(h, a,b,c),c="orange")
        best_fit_h=b
        print a,b,c,"  max=", max_h
    except:
        print "blew up,  max=", max_h
        best_fit_h=d.max()
    plt.title("fit in blue, gaussian with same parameters in orange")
    return sim, best_fit_h

def rho_check(n, simCheck, max_h, dir, bay=True, year=2016) :
    import matplotlib.pyplot as plt
    locDict = pickle.load(open(dir+"locationDict","rb"))
    pickleDir = locDict["pickleDir"] 
    mapFile = pickleDir + locDict["metaMapFile"]
    repDir = locDict["stripStringDir"] 
    dist,sigma,ng,nw,area = np.genfromtxt(mapFile,unpack=True,skiprows=1,usecols=(1,2,5,6,7))
    file = np.genfromtxt(mapFile,unpack=True,skiprows=1,usecols=(0), dtype="str")
    file = file[n].replace(repDir,"")
    file = file.replace(".fits.gz","")
    sim = int(file)
    #mra, mdec, map = ho_prep_data.getSimMap(sim);
    #return mra, mdec, map 
    #ra,dec,imag,z,zerr = ho_prep_data.get2mCat(dir="/home/s1/annis/daedalean/desgw-cos/data/")
    #ra,dec,imag,z,zerr = ho_prep_data.getSDSSCat(dir="/home/s1/annis/daedalean/desgw-cos/data/")
    #return ra, dec, z
    if sim != simCheck : raise Exception("how? {} {}".format(sim, simCheck))
    d,r = look_at_rho(sim, h=max_h, bay=bay, year=year)
    plt.clf()
    plt.plot(d,r); plt.title("distance = {} h={}".format(dist[n], max_h))
    ymin,ymax = plt.ylim(); 
    plt.plot([dist[n],dist[n]],[ymin,ymax], c= "r")
    print sim, dist[n], sigma[n]

def look_at_rho (sim, h=70, bay=True, year=2016) :
    import ho_measure
    mra, mdec, map = ho_prep_data.getSimMap(sim, bay=bay, year=year);
    zeds, pdf_c = ho_prep_data.get2mCatPdf( map, dir="/home/s1/annis/daedalean/desgw-cos/data/" ); 
    #zeds, pdf_c = ho_prep_data.getSDSSCatPdf( map, dir="/home/s1/annis/daedalean/desgw-cos/data/" ); 
    dndz, weight_dndz = ho_prep_data.getCatDndz( zeds ); 
    d=np.arange(150,600,1);
    r = ho_measure.rho(d,h, pdf_c, weight_dndz); 
    return d, r

# plt.clf();plt.scatter(h,d);plt.scatter(h1,d1,c="g"); 
# plt.plot(h1, gaussian(h1, 0.0155959433626, 72.1251199107, 30),c="orange"); 
#plt.xlabel("$H_0$");plt.ylabel("posterior probability density"); plt.title("2MASS XSC and GW150914");plt.text(45,0.005,"blue: probaility density with pr($H_0$)=1"); plt.text(45,0.004,"green: with pr($H_0$) of gaussian $H_0$=70, $\sigma$=30"); plt.text(45,0.003,"width of green = 28")


# h5,d5=allMaps.getall(77,exclude_ix=ix)
def getall(dir, start, end, exclude_ix="", include_ix="") :
    counter = 0; first = -1
    if exclude_ix != "" :
        #skip = np.nonzero(exclude_ix)[0].tolist()
        skip =exclude_ix.tolist()
    if include_ix != "" :
        #get = np.nonzero(include_ix)[0].tolist()
        get =include_ix.tolist()
    data = dict()
    for i in range(start,end) :
        if exclude_ix != "" :
            if i in skip : continue
        if include_ix != "" :
            if i not in get : continue
        print i, "{}{}".format(dir,i)
        data[i]=pickle.load(open("{}{}".format(dir, str(i)),"rb"))
        counter += 1
        if first == -1: first = i
    print "collected ",counter
    pdf=np.copy(data[first][1])
    distance = data[first][0]
    for i in range(start,end) :
        if exclude_ix != "" :
            if i in skip : continue
        if include_ix != "" :
            if i not in get : continue
        if i != first:
            pdf = pdf*data[i][1]
            pdf = pdf/pdf.max()
    return distance, pdf

def findOne (dir, max, simNum) :
    for i in range(0,max+1):
        d=pickle.load(open("{}{}".format(dir,i),"rb"))
        if d[2] == simNum: 
            print i

# ================================================================
# support routines
#

def gaussian(x, amp, mu, sig): 
    return amp*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
def fit_curve(h, d, guess_amp, guess_mean, guess_sigma, verbose=True) :
    from scipy.optimize import curve_fit
    amp = guess_amp
    mean = guess_mean
    sigma = guess_sigma
    bf, covar = curve_fit(gaussian, h, d, p0=[amp,mean,sigma])

    #interp_d = interpolate.interp1d(h,d,kind="cubic")
    #interp_h = np.arange(h.min(),h.max(),0.05)
    #bf, covar = curve_fit(gaussian, interp_h, interp_d(interp_h), p0=[amp,mean,sigma])
    #bf, covar = curve_fit(gaussian, h, d, p0=[amp,mean,sigma], bounds=[[0.,30.,.1],[1.,150.,100.,]])

    if verbose: print bf[0], bf[1], bf[2]
    return bf[0], bf[1], bf[2]

def sigmaClip (vector, sigma=3.0) :
    ix = sigmaclipIndex(vector,sigma)
    mean = vector[ix].mean()
    sigma = vector[ix].std()
    return mean, sigma
def sigmaclipIndex ( vector, sigma=3.0) :
    mean = vector.mean()
    std = vector.std()
    ix = np.nonzero( abs(vector-mean) < sigma*std)
    mean = vector[ix].mean()
    std = vector[ix].std()
    ix = np.nonzero( abs(vector-mean) < sigma*std)
    mean = vector[ix].mean()
    std = vector[ix].std()
    ix = np.nonzero( abs(vector-mean) < sigma*std)
    #print np.nonzero(abs(vector-mean) > sigma*std)
    return ix

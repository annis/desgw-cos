import numpy as np
# import os; import gridwrx; os.chdir("~/gridworks")
# os.chdir("/home/s1/annis/daedalean/desgw-cos/"); reload(gridwrx); os.chdir("/home/s1/annis/gridworks"); 
# gridwrx.prep(10)
# gridwrx.prep(11, "map-2015-74.txt", "d_pdf_2015_74/", targetHo=74.)

# i work in /home/s1/annis/gridworks -> /pnfs/des/scratch/gw/annis/

def prep(jobnum, output="maps-2015-66.txt", pickleDir = "d_pdf_2015_66/", year=2015, 
        bay=True, derr = 0.35, targetHo=66., grid=True, cat="2mpz")  :

    sim_tar_file = "sims-2015.tar"
    if year == 2016 and bay == True : sim_tar_file = "sims-2016.tar"
    if year == 2016 and bay == False : sim_tar_file = "sims-2016-la.tar"
    outdir = pickleDir.rstrip('/')

    write_shfile(jobnum=jobnum, sim_tar_file=sim_tar_file, outdir=outdir) 
    write_py1(jobnum=jobnum, output=output, pickleDir = pickleDir, year=year, 
        bay=bay, derr = derr, targetHo=targetHo, cat=cat, grid=grid)  
    write_py2(jobnum=jobnum, output=output, pickleDir = pickleDir, year=year, 
        bay=bay, derr = derr, targetHo=targetHo, cat=cat, grid=grid)  


def write_shfile(jobnum, sim_tar_file="sims-2015.tar",outdir="pdf_2015_66") :
    sh_file = "run-{}.sh".format(jobnum)
    py1_file = "py1-{}.py".format(jobnum)
    py2_file = "py2-{}.py".format(jobnum)
    print "writing ", sh_file
    fd = open(sh_file,"w")
    fd.write("source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob.sh\n")
    fd.write("ifdh cp -D /pnfs/des/scratch/gw/annis/staging_out/{} ./\n".format(sim_tar_file))
    fd.write("ifdh cp -D /pnfs/des/scratch/gw/annis/staging_out/desgw-cos.tar ./\n")
    fd.write("tar xvf {}\n".format(sim_tar_file))
    fd.write("tar xvf desgw-cos.tar\n")
    fd.write("cd desgw-cos\n")
    fd.write("setup Y2Nstack 1.0.6+14\n")
    fd.write("setup healpy\n")
    fd.write("ifdh cp -D /pnfs/des/scratch/gw/annis/{} ./\n".format(py1_file))
    fd.write("ifdh cp -D /pnfs/des/scratch/gw/annis/{} ./\n".format(py2_file))
    fd.write("python {} \n".format(py1_file))
    fd.write("python {} \n".format(py2_file))
    fd.write("\n")
    fd.write("tar -cvf {}.tar {}\n".format(outdir,outdir))
    fd.write("ifdh cp {}.tar /pnfs/des/scratch/gw/annis/staging_in/{}.tar\n".format(outdir,outdir))
    fd.close()

def write_py1(jobnum, output="maps-2015-66.txt", pickleDir = "d_pdf_2015_66/", year=2015, 
        bay=True, derr = 0.35, targetHo=66., cat="2mpz", grid=True)  :
    file = "py1-{}.py".format(jobnum)
    print "writing ", file
    fd = open(file,"w")
    fd.write("import ho_sim\n\n")
    fd.write("ho_sim.veniVidi (\n")
    fd.write("\t output=\"{}\", \n".format(output))
    fd.write("\t pickleDir = \"{}\", \n".format(pickleDir))
    fd.write("\t year={}, bay={}, derr = {:.2f}, targetHo={:.1f}, catalog={},  grid={})\n".format(
        year, bay,  derr, targetHo, cat, grid) )
    fd.close()
def write_py2(jobnum, output="maps-2015-66.txt", pickleDir = "d_pdf_2015_66/", year=2015, 
        bay=True, derr = 0.35, targetHo=66., cat="2mpz", grid=True)  :
    file = "py2-{}.py".format(jobnum)
    print "writing ", file
    fd = open(file,"w")
    fd.write("import ho_sim\n\n")
    fd.write("ho_sim.nowVici (\n")
    fd.write("\t output=\"{}\", \n".format(output))
    fd.write("\t pickleDir = \"{}\", \n".format(pickleDir))
    fd.write("\t year={}, bay={}, derr = {:.2f}, targetHo={:.1f},  catalog={}, grid={})\n".format(
        year, bay,  derr, targetHo, cat, grid) )
    fd.close()


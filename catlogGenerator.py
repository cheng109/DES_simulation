"""
  @package
  @file gridCatalog.py
  @brief generate stars on square grid

  @brief Created by:
  @author Jun Cheng (Purdue)

  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python gridCatalog.py

Notes: 1. You could set the 'magnitude', 'density'

"""
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import math
#from random import gauss
from astropy.io import fits
import subprocess
import conversion
import numpy as np
import commons
import measurepsf
import math
import random

CLASSIFICATION_THRESHOLD = 0.95
PIXEL_SCALE = 0.27
MAGNITUDE_THRESHOLD = 22.5

class SexSource():
    def __init__(self, chipID,  objDict):

        self.chipID = chipID
        self.xcenter = float(objDict["X_IMAGE"])
        self.ycenter = float(objDict["Y_IMAGE"])
        self.ra = float(objDict["X_WORLD"])
        self.dec = float(objDict["Y_WORLD"])
        self.xx = float(objDict["X2_IMAGE"])
        self.yy = float(objDict["Y2_IMAGE"])
        self.xy = float(objDict["XY_IMAGE"])
        self.e1 = (self.xx-self.yy)/(self.xx+self.yy)
        self.e2 = 2*self.xy/(self.xx+self.yy)
        self.angle = 0 #math.atan(self.e2/self.e1)
        self.e = np.sqrt(self.e1**2 + self.e2**2)
        self.new_ra = 0
        self.new_dec = 0

        self.id = objDict["NUMBER"]
        self.mag = float(objDict["MAG_AUTO"])
        self.gamma1 = 0
        self.gamma2= 0
        self.mu = 0
        self.delta_ra = 0
        self.delta_dec = 0
        self.z = 0
        self.type = "star"
        self.sed = "../sky/sed_flat.txt"

        self.dx = 0
        self.dy = 0
        self.matchMag = 0
        self.matchID = "none"
        self.matchFWHM = 0

        self.new_e1 = 0
        self.new_e2 = 0
        self.new_e = 0
        self.new_angle = 0
        self.rms = 0
        self.flux = 0
        self.converge = False
        self.alphax=0
        self.alphay=0

        self.fwhm = float(objDict["FWHM_IMAGE"])

        self.class_star = float(objDict["CLASS_STAR"])
        if  self.class_star < CLASSIFICATION_THRESHOLD or self.mag > MAGNITUDE_THRESHOLD :
            self.type = "gauss"
            self.gauss_sigma = self.fwhm*PIXEL_SCALE/2.3355   #convert degree to arcsec

    def writeOutCatalog(self, magCorrection):
        s = "object " \
            + str(self.id) + " "\
            + str(self.ra)  + " "\
            + str(self.dec) + " "\
            + str(float(self.mag)+magCorrection) +" "\
            + str(self.sed) + " "\
            + str(self.z)   + " "\
            + str(self.gamma1) + " "\
            + str(self.gamma2) + " "\
            + str(self.mu) + " "\
            + str(self.delta_ra) + " "\
            + str(self.delta_dec) + " "\
            + str(self.type) + " "
        if self.type == "star":
            s += "none none\n"
        # Assume Gaussion Galaxies:
        if self.type == "gauss":
            #s= ''
            s += str(self.gauss_sigma)+" none none\n"
        return s


    def updateEllipticity(self, data):
        #print self.xcenter, self.ycenter
        self.new_e1, self.new_e2, self.new_e, self.new_angle, self.flux, self.rms, self.converge, self.alphax, self.alphay = measurepsf.getEllipticity(data)
        #print self.new_angle



def readSexCatalog(chipID, sexCatalogFile):
    sexFile = open(sexCatalogFile, 'r')
    paramList = []
    objList = []
    objDict = {}
    for line in sexFile.readlines():
        temp = line.split();
        if line[0]=='#':
            paramList.append(temp[2])
        if line[0]!='#':
            for i in range(len(paramList)):
                objDict[paramList[i]]=temp[i]
            m = SexSource(chipID,objDict)
            #m.updateWCS(imageHeader)
            objList.append(m)
    sexFile.close()
    return objList


def writePhosimCatalog(phosimCatalogFile, objList, header, magCorrection):
    phosimFile = open(phosimCatalogFile, 'w')
    phosimFile.write(header)
    for source in objList:
        phosimFile.write(source.writeOutCatalog(magCorrection))
    phosimFile.close()

def writePhosimCommand(phosimCommandName, imageHeader, shift, tilt):
    coefArcsecToRad = 4.8481e-6
    file = open(phosimCommandName,'w')
    file.write("zenith_v 1000.0\n"
               + "raydensity 0.0\n"
               + "clearperturbations\n"
               #+ "well_depth 43796"
               + "overdepbias -40.0\n"
               + "siliconthickness 250\n"
               #+ "throughputfile 1\n"
               + "domeseeing 0.0\n"
               + "lascatprob 0.0\n"
               #+ "straylight 0.0\n"
               #+ "cleartracking\n"
               #+ "atmosphericdispersion 0\n"
               #+ "detectormode 0\n"
               #+ "telescopemode 0\n"
               + "body 0 0 " + str(tilt["phi"]  *coefArcsecToRad) + "\n"
               + "body 0 1 " + str(tilt["psi"]  *coefArcsecToRad) + "\n"
               + "body 0 2 " + str(tilt["theta"]*coefArcsecToRad) + "\n"
               + "body 0 3 " + str(shift["x"]) + "\n"
               + "body 0 4 " + str(shift["y"]) + "\n"
               + "body 0 5 " + str(shift["z"]) + "\n"
               )
    file.close()


def createSingleChipCatalog(chipID, dRa, dDec ,rotation, imageName, catalog , command, magCorrection, rawSeeing, shift, tilt, suffix):
    catalogHeader ="""Opsim_rottelpos 0
Opsim_moondec -90
Opsim_moonra 180
Opsim_expmjd 49552.3
Opsim_moonalt -90
Opsim_sunalt -90
Opsim_filter 1
Slalib_date 1994/7/19/0.298822999997
Opsim_dist2moon 180.0
Opsim_moonphase 10.0
Opsim_obshistid 99999999
SIM_SEED     1000
SIM_MINSOURCE 1
SIM_TELCONFIG 0
SIM_CAMCONFIG 1
SIM_NSNAP 1
"""
    Catalog = "dataCatalog"
    imageHeader = fits.getheader(imageName, 0)

    #rotskypos = 269.881 + 0.0108
    pointingDEC = conversion.convDMS(imageHeader["TELDEC"])

    #rawSeeing = str(imageHeader["DIMMSEE"])
    #rawSeeing = str(0.9)
    factor = np.cos(np.radians(pointingDEC))
    subprocess.call("sex -c data.sex " + imageName + " -CATALOG_NAME " + Catalog, shell=True)
    catalogHeader = "Unrefracted_RA_deg " +str(conversion.convHMS(imageHeader["TELRA"])+dRa*factor)+\
                    "\nUnrefracted_Dec_deg " + str(conversion.convDMS(imageHeader["TELDEC"])+dDec) +\
                    "\nUnrefracted_Azimuth " + str(imageHeader["AZ"]) +\
                    "\nUnrefracted_Altitude " + str(90-float(imageHeader["ZD"])) +\
                    "\nOpsim_rotskypos " + str(rotation) +\
                    "\nOpsim_rawseeing " + str(rawSeeing) +\
                    "\nSIM_VISTIME " + str(imageHeader["EXPTIME"]) +\
                    "\n" + catalogHeader
                     #"\nSlalib_date " + str(imageHeader["DATE-OBS"]) +\


    objList = readSexCatalog(chipID, Catalog)
    writePhosimCatalog(catalog, objList, catalogHeader, magCorrection)
    writePhosimCommand(command, imageHeader, shift, tilt)

def createRun(scriptName, chip, shift, tilt,  DestinationDir,  suffix):


    outputName = DestinationDir + "/" + "Images_" + chip \
                 + "_x_" + str(shift["x"])\
                 + "_y_" + str(shift["y"])\
                 + "_z_" + str(shift["z"])\
                 + "_phi_" + str(tilt["phi"])\
                 + "_psi_" + str(tilt["psi"])\
                 + "_theta_" + str(tilt["theta"])\
                 + ".fits"

    f=open(scriptName, 'w')
    runScript = """#!/bin/bash
pwd=$(pwd)
PHOSIM_PATH=/Users/cheng109/toberemoved/phosim/phosim_core/
cd $PHOSIM_PATH
"""

    workDirectory = chip + "_work_ID_" + suffix
    outputDirectory = chip + "_output_ID_" + suffix

    createDirectory = "mkdir "+ workDirectory +" " + outputDirectory
    deleteDirectory = "rm -rf $pwd/" +  workDirectory +"\n" \
                      + "rm -rf $pwd/"+ outputDirectory + "\n" \
                      + "rm -rf $pwd/" + chip + "_phosimCatalog_ID_"+suffix + "\n" \
                      + "rm -rf $pwd/" + chip + "_phosimCommand_ID_"+suffix

    subprocess.call(createDirectory, shell=True)
    runScript += "./phosim $pwd/" + chip +"_phosimCatalog_ID_" + suffix + " -c $pwd/" \
                     + chip +"_phosimCommand_ID_" + suffix+ " -e 0 -w $pwd/"+ workDirectory + " -o $pwd/" + outputDirectory+ " " \
                     +  "-i deCam " \
                     +  "-s " +chip \
                     +" && gunzip -f $pwd/" + outputDirectory + "/*"+chip +"*.gz \n"\
                     +"cp $pwd/" + outputDirectory+ "/*"+chip +"*.fits " + "$pwd/" + outputName + "\n" \
                    + deleteDirectory + "\n\n"





    f.write(runScript)
    f.close()

    subprocess.call("chmod +x " + scriptName, shell=True)

def catalogGenerator(chip,  coarseDec, coarseRa,coarseRotation, fineCorrect, magCorrection, rawSeeing, shift, tilt, suffix):

    Dec=(coarseDec+fineCorrect[chip]["dec"])*0.27/3600
    Ra=(coarseRa+fineCorrect[chip]["ra"])*0.27/3600
    rotation = coarseRotation + fineCorrect[chip]["rotation"]
    createSingleChipCatalog(chip, Ra, Dec, rotation, imageName = chip+ "_test_image.fits", catalog = chip+ "_phosimCatalog_ID_"+suffix, command= chip + "_phosimCommand_ID_" +suffix,
                            magCorrection=magCorrection, rawSeeing=rawSeeing, shift=shift, tilt=tilt, suffix = suffix)

    return 0


### For coarse correction reference:
### For one specific star with format (x,y) or (dec, ra)
# Data       =  (841,   1941),
# Sim        =  (845,   1937).
# Data - Sim =  (-4,    4   )
###     Then the coarse correction is:
#            =  (-4,    -8  )
####     Rotation slope :
####     Rotation slope correction:

def main():





    CHIPS = ["N1", "N4", "N7", "N22", "S22"]
    #CHIPS = ["S22"]
    shift = {"x": 0.0, "y": 0.0, "z": 0.0}
    tilt = {"phi": 0.0, "psi": 0.0, "theta": 0.0}    # degree
    #tilt = {"phi": 324000.0, "psi": 0.0, "theta": 0.0}


    fineCorrect={}
    rawSeeing = 0.9


    DestinationDir = "new_all_chips_results_tilt"

    combRUN = ""
    for chip in CHIPS:
        for angle in np.arange(-20,0, 200):
            suffix = str(random.randint(0, 999999999))
            tilt["theta"] = angle
            dDEC, dRA, dROTATION = commons.posCorrection(shift, tilt)


            confCoarseMap = commons.parseConfigure("conf.txt")
            DEC, RA, ROTATION = confCoarseMap[chip][0]+ dDEC, confCoarseMap[chip][1]+ dRA, 270.0

            fineCorrect[chip] = {"dec":confCoarseMap[chip][0], "ra":confCoarseMap[chip][1], "rotation":ROTATION}    #magCorrection = 2;
            catalogGenerator(chip,  coarseDec=dDEC, coarseRa=dRA, coarseRotation=0, fineCorrect=fineCorrect, magCorrection=2.0, rawSeeing=rawSeeing, shift=shift, tilt=tilt, suffix = suffix)

            # create bash script
            scriptName = "run_" + chip + "_phi_" + str(tilt["phi"])  + "_" + str(angle)
            createRun(scriptName, chip, shift,tilt,  DestinationDir, suffix=suffix)
            combRUN = combRUN + "./" + scriptName + " & \n"

    f = open("combRUN", 'w')
    f.write(combRUN)
    f.close()
    subprocess.call("chmod +x combRUN" , shell=True)

if __name__=='__main__':
    main()


# theta = -80 ;   dDEC =
# theta = -60 ;   dDEC =

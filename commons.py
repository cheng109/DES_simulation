__author__ = 'cheng109'
from astropy.io import fits

class ChipConfig:

    def __init__(self, confList):
        self.chip = confList[0]
        self.shiftX = float(confList[1])
        self.shiftY = float(confList[2])
        self.shiftZ = float(confList[3])
        self.fineDEC = float(confList[4])
        self.fineRA = float(confList[5])
        self.coarseDEC = 0
        self.coarseRA = 0
        self.DEC = 0
        self.RA = 0
    def updatePosition(self, confCourseMap):
        self.coarseDEC =float(confCourseMap[self.chip][0])
        self.coarseRA = float(confCourseMap[self.chip][1])
        self.DEC = self.coarseDEC + self.fineDEC
        self.RA  = self.coarseRA  + self.fineRA

def copyWCS(dataImageName,simulationImageName):
    dataHeader = fits.getheader(dataImageName, 0)
    data,  SimulationHeader= fits.getdata(simulationImageName, header=True)

    SimulationHeader["CRPIX1"] = dataHeader["CRPIX1"]
    SimulationHeader["CRPIX2"] = dataHeader["CRPIX2"]
    SimulationHeader["CRVAL1"] = dataHeader["CRVAL1"]
    SimulationHeader["CRVAL2"] = dataHeader["CRVAL2"]

    SimulationHeader["CD1_1"] = dataHeader["CD1_1"]
    SimulationHeader["CD1_2"] = dataHeader["CD1_2"]
    SimulationHeader["CD2_1"] = dataHeader["CD2_1"]
    SimulationHeader["CD2_2"] = dataHeader["CD2_2"]

    SimulationHeader["CTYPE1"] = dataHeader["CTYPE1"]
    SimulationHeader["CTYPE2"] = dataHeader["CTYPE2"]

    SimulationHeader["CUNIT1"] = dataHeader["CUNIT1"]
    SimulationHeader["CUNIT2"] = dataHeader["CUNIT2"]

    SimulationHeader["MJD-OBS"] = dataHeader["MJD-OBS"]

    SimulationHeader["PV1_1"] = dataHeader["PV1_1"]
    SimulationHeader["PV1_2"] = dataHeader["PV1_2"]
    SimulationHeader["PV1_3"] = dataHeader["PV1_3"]
    SimulationHeader["PV1_4"] = dataHeader["PV1_4"]
    SimulationHeader["PV1_5"] = dataHeader["PV1_5"]
    SimulationHeader["PV1_6"] = dataHeader["PV1_6"]
    SimulationHeader["PV1_7"] = dataHeader["PV1_7"]
    SimulationHeader["PV1_8"] = dataHeader["PV1_8"]
    SimulationHeader["PV1_9"] = dataHeader["PV1_9"]
    SimulationHeader["PV1_10"] = dataHeader["PV1_10"]

    SimulationHeader["PV1_0"] = dataHeader["PV1_0"]
    SimulationHeader["PV1_1"] = dataHeader["PV1_1"]
    SimulationHeader["PV1_2"] = dataHeader["PV1_2"]
    SimulationHeader["PV1_3"] = dataHeader["PV1_3"]
    SimulationHeader["PV1_4"] = dataHeader["PV1_4"]
    SimulationHeader["PV1_5"] = dataHeader["PV1_5"]
    SimulationHeader["PV1_6"] = dataHeader["PV1_6"]
    SimulationHeader["PV1_7"] = dataHeader["PV1_7"]
    SimulationHeader["PV1_8"] = dataHeader["PV1_8"]
    SimulationHeader["PV1_9"] = dataHeader["PV1_9"]
    SimulationHeader["PV1_10"] = dataHeader["PV1_10"]

    SimulationHeader["PV2_0"] = dataHeader["PV2_0"]
    SimulationHeader["PV2_1"] = dataHeader["PV2_1"]
    SimulationHeader["PV2_2"] = dataHeader["PV2_2"]
    SimulationHeader["PV2_3"] = dataHeader["PV2_3"]
    SimulationHeader["PV2_4"] = dataHeader["PV2_4"]
    SimulationHeader["PV2_5"] = dataHeader["PV2_5"]
    SimulationHeader["PV2_6"] = dataHeader["PV2_6"]
    SimulationHeader["PV2_7"] = dataHeader["PV2_7"]
    SimulationHeader["PV2_8"] = dataHeader["PV2_8"]
    SimulationHeader["PV2_9"] = dataHeader["PV2_9"]
    SimulationHeader["PV2_10"] = dataHeader["PV2_10"]

    fits.writeto(simulationImageName,data, header=SimulationHeader, clobber=True)


def rmBackground(icon, background):
    for i in range(icon.shape[0]):
        for j in range(icon.shape[1]):
            icon[i][j]= icon[i][j]-background
            if icon[i][j]<0:
                icon[i][j]= 0

    return icon

def plotClassStar(new_dataObjList):
    mag = []
    classStar = []
    for obj in new_dataObjList:
        mag.append(obj.mag)
        classStar.append(obj.class_star)
    plt.plot(mag, classStar, '.')
    plt.xlabel("Magnitude")
    plt.ylabel("Class_Star")
    plt.show()
    return 0


def writeStampe(dataDict):
    dirName = "stampes/"

    for key, data in dataDict.iteritems():
        hdu = fits.PrimaryHDU(data)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dirName+key+".fits", clobber=True)

def removeSubplotPanel(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9):

    ax1.axis('off')
    ax3.axis('off')
    ax7.axis('off')
    ax9.axis('off')
    ax2.axes.get_yaxis().set_visible(False)
    ax4.axes.get_yaxis().set_visible(False)
    ax5.axes.get_yaxis().set_visible(False)
    ax6.axes.get_yaxis().set_visible(False)
    ax8.axes.get_yaxis().set_visible(False)


def showLegend(axList) :

    for ax in axList:
        ax.legend()


def parseConfigure(confFileName):
    confCoarseMap = {}

    confFile = open(confFileName, 'r')
    for line in confFile.readlines():

        temp = line.split()
        length = len(temp)
        if line[0]!='#' and length >0 :
            if length==4:
                confCoarseMap[temp[0]] = (float(temp[2]), float(temp[3]))


    confFile.close()
    return confCoarseMap


def posCorrection(shift, tilt):
    #tilt = {"phi": 324000.0, "psi": 0.0, "theta": 0.0}
    dDEC2 = 0
    dRA2 = 0
    # Correction for SHIFT
    dDEC1 = -71.467*shift["x"]
    dRA1 = 190.91* shift["y"] + 3.889

    # Correction for TILT:
    if tilt["phi"]==90.0:
        dDEC2 =  -7.3744 *tilt["theta"] - 0.869
    else:
        dRA2 = -19.829*tilt["theta"] -12.786
    dROTATION = 0
    return dDEC1+dDEC2, dRA1+dRA2, dROTATION


def main():
    parseConfigure("conf.txt")




    return "Nothing"

if __name__=='__main__':

    main()


    #fineCorrect["N4"] = {"dec":-4.0, "ra":-8.0, "rotation":269.814}
    #catalogGenerator(chip, coarseRa=6.2, coarseDec=-0.74, coarseRotation=0.019, fineCorrect=fineCorrect, magCorrection=1.95)
    #fineCorrect["N7"] = {"dec":-43, "ra":-82, "rotation":270.0}   #magCorrection = 2;
    #fineCorrect[chip] = {"dec":-32, "ra":80, "rotation":270.0}    #magCorrection = 2;


# N22:  x =-0.5: (37, 0)             y=-0.5: (0, -91)
# N22:  x =-1.0: (71, 0)             y=-1.0: (0, -188)
# N22:  x =-1.5: (107,0)            y=-1.5: (0, -283)
# N22:  x =-2.0: (143, 0)            y=-2.0: (0, -376)
# N22:  x =0.5: (-37, 0)            y=0.5: (0, 100)
# N22:  x =1.0: (-71, 0)            y=1.0: (0, 196)
# N22:  x =1.5: (-107, 0)            y=1.5: (0,291 )
# N22:  x =2.0: (-143, 0)            y=2.0: (0,386 )


# S22:   y = -0.5:  (0,  -91)
# S22:   y = -1.0:  (0, -188)
# S22:   y = -1.5:  (0, -283)
# S22:   y = -2.0:  (0, -376)
# S22:   y = 0.0:   (0, 0)
# S22:   y = 0.5:  (0,  100)
# S22:   y = 1.0:  (0, 196)
# S22:   y = 1.5:  (0, 291)
# S22:   y = 2.0:  (0, 386)


# N1:   y = -0.5:  (0,  -91)
# N1:   y = -1.0:  (0, -188)
# N1:   y = -1.5:  (0, -283)
# N1:   y = -2.0:  (0, -376)

# N1:   y = 0.5:  (0,  100)
# N1:   y = 1.0:  (0, 196)
# N1:   y = 1.5:  (0, 291)
# N1:   y = 2.0:  (0, 386)



#globals()
# N4:  x =-0.5: (37, 0)             y=-0.5: (0, -91)
# N4:  x =-1.0: (71, 0)             y=-1.0: (0, -188)
# N4:  x =-1.5: (107,0)            y=-1.5: (0, -283)
# N4:  x =-2.0: (143, 0)            y=-2.0: (0, -376)
# N4:  x =0.5: (-37, 0)            y=0.5: (0, 100)
# N4:  x =1.0: (-71, 0)            y=1.0: (0, 196)
# N4:  x =1.5: (-107, 0)            y=1.5: (0,291 )
# N4:  x =2.0: (-143, 0)            y=2.0: (0,386 )


# N7:  x= -0.5: (
# N7:  x= -1.0: (71, 0)
# N7:  x=-1.5: ( 107, 0)           y = -2.0: (0, -283)
# N7:  x=-2.0: (143, 0)            y = -2.0: (0, -376)
# N7:  x=2.0: (-143, 0)            y = 2.0: (0, 386)


# S30: x=-2.0: (     , 0)           y =-2.0: (0, -371)
# S30: x=0.0:  (     , 0)           y = 0.0: (0, 0)
# S30: x=2.0:  (     , 0)           y = 2.0: (0, 386)
# S30: x=1.0:  (     , 0)           y = 1.0: (0, 196)
# S30: x=1.0:  (     , 0)           y = -1.0: (0, -188)



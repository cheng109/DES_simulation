
import numpy as np
import matplotlib.pyplot as plt
import catlogGenerator
import subprocess
import matchCatalog
import math

def getObjectList(CHIPS):
    writeStamps = False
    dataFileNameList = []
    simFileNameList = []

    for chip in CHIPS:
        dataFileNameList.append(chip + "_test_image.fits")
        simFileNameList.append("output/" +"deCam_e_99999999_f1_" + chip + "_E000.fits")
    dataObjDict = {}
    simulationObjDict = {}
    for i in range(len(CHIPS)):
        chip = CHIPS[i]
        dataSExName = "output/" + chip+"_data_sex_cat"
        simSExName = "output/" + chip+"_sim_sex_cat"
        subprocess.call("/opt/local/bin/sex -c data.sex " + dataFileNameList[i] + " -CATALOG_NAME " + dataSExName, shell=True)
        subprocess.call("/opt/local/bin/sex -c sim.sex " + simFileNameList[i] + " -CATALOG_NAME " + simSExName, shell=True)

        dataObj  = catlogGenerator.readSexCatalog(chip, dataSExName)
        dataObjDict[chip]  = matchCatalog.updateEllipticity(dataFileNameList[i],  dataObj, type='data', writeStamps=writeStamps)

        simObj = catlogGenerator.readSexCatalog(chip, simSExName)
        simulationObjDict[chip]= matchCatalog.updateEllipticity(simFileNameList[i], simObj, type='simulation', writeStamps=writeStamps)
    return dataObjDict, simulationObjDict

def combineChips(objDict):
    objList = []
    for key, value in objDict.iteritems():
        objList = objList + value
    return objList

def getCorrelation(objDict):
    ObjList = combineChips(objDict)
    l = len(ObjList)
    print "number of obj: ", l
    Distance = []
    Distance2 = []
    Exx = []
    counter = 0;
    for i in range(l):
        if ObjList[i].type =="star":
            starA = ObjList[i]

            for j in np.arange(i+1, l-i-1):
                if ObjList[j].type =="star":
                    counter = counter+1
                    starB = ObjList[j]
                    dRa = starA.ra - starB.ra
                    dDec = starA.dec - starB.dec
                    dx = starA.xcenter - starB.xcenter
                    dy = starA.ycenter - starB.ycenter
                    if dRa<1.0e-8:
                        continue
                    if dx<1.0e-6:
                        continue
                    angle = -2* math.atan(dDec/dRa)   # in radian
                    angle = -2* math.atan(dy/dx)   # in radian
                    XA_t=(starA.new_e1*np.cos(angle)-starA.new_e2*np.sin(angle))
                    XA_x=(starA.new_e1*np.sin(angle)+starA.new_e2*np.cos(angle))
                    XB_t=(starB.new_e1*np.cos(angle)-starB.new_e2*np.sin(angle))
                    XB_x=(starB.new_e1*np.sin(angle)+starB.new_e2*np.cos(angle))
                    #Distance.append( np.sqrt(dRa**2+dDec**2)*3600/0.27)
                    Distance2.append(np.sqrt(dx**2+dy**2))
                    print Distance2[-1]

                    Exx.append( XA_t*XB_t + XA_x*XB_x)
    print "number of stars conbinations: ", counter
    interval = 300
    Exxbins = [[] for i in range(int(10000/interval))]
    print len(Distance2)
    for i in range(len(Exx)):   # bins[0] means interval between [0, 50]

        Exxbins[int(Distance2[i]/interval)].append(Exx[i])
    newDistance = []
    newExx = []
    stdErr = []
    for i in range(len(Exxbins)):
        if len(Exxbins[i])>0:
            newDistance.append((i+0.5)*interval*0.27/60)
            newExx.append(np.mean(Exxbins[i]))
            stdErr.append(np.std(Exxbins[i])/np.sqrt(len(Exxbins[i])))
    return newDistance, newExx, stdErr


def plotCor(dataObjDict,simulationObjDict ):

    dataDistance, dataExx, dataStdErr = getCorrelation(dataObjDict)
    simDistance, simExx, simStdErr = getCorrelation(simulationObjDict)
    plt.errorbar(dataDistance, dataExx, dataStdErr, fmt='o', label="data")
    plt.errorbar(simDistance, simExx, simStdErr, fmt='o', label="simulation")
    plt.xlabel("Distance (arcmin)")
    plt.ylabel("Ellipticity correlation")
    plt.title("Elliptictiy correlation")
    plt.legend()

    plt.show()

def readinfo(fileName, CHIP):

    infoFile = open(fileName, 'r')
    centerX = []
    centerY = []
    e1 = []
    e2 = []
    e = []

    for line in infoFile.readlines():
        temp = line.split();
        centerX.append(float(temp[0]))
        if CHIP=="N7":
            centerY.append(float(temp[1])+12000)
        centerY.append(float(temp[1]))
        e1.append(float(temp[2]))
        e2.append(float(temp[3]))
        e.append(float(temp[4]))
        print temp

    infoFile.close()
    return centerX, centerY, e1, e2, e

def getCor(xcenter, ycenter, e1, e2, e):
    l = len(xcenter)
    Distance = []
    Exx = []
    for i in range(l):
        for j in np.arange(i+1, l-i-1):
            dx = xcenter[i] - xcenter[j]
            dy = ycenter[i] - ycenter[j]
            if dx<1.0e-6:
                continue
            angle = -2* math.atan(dy/dx)   # in radian
            XA_t=(e1[i]*np.cos(angle)-e2[i]*np.sin(angle))
            XA_x=(e1[i]*np.sin(angle)+e2[i]*np.cos(angle))
            XB_t=(e1[j]*np.cos(angle)-e2[j]*np.sin(angle))
            XB_x=(e1[j]*np.sin(angle)+e2[j]*np.cos(angle))
            Distance.append(np.sqrt(dx**2+dy**2))
            print Distance[-1]

            Exx.append( XA_t*XB_t + XA_x*XB_x)

    interval = 300
    Exxbins = [[] for i in range(int(18000/interval))]
    print len(Distance)
    for i in range(len(Exx)):   # bins[0] means interval between [0, 50]

        Exxbins[int(Distance[i]/interval)].append(Exx[i])
    newDistance = []
    newExx = []
    stdErr = []
    for i in range(len(Exxbins)):
        if len(Exxbins[i])>0:
            newDistance.append((i+0.5)*interval*0.27/60)
            newExx.append(np.mean(Exxbins[i]))
            stdErr.append(np.std(Exxbins[i])/np.sqrt(len(Exxbins[i])))



    return newDistance, newExx, stdErr



def main():

    CHIPS = ["N4","N7"]
    dataCenterX = []
    dataCenterY = []
    dataE1 = []
    dataE2 = []
    dataE = []

    simCenterX = []
    simCenterY = []
    simE1 = []
    simE2 = []
    simE = []
    for chip in CHIPS:
        dataFileName = "data_info_" + chip+".txt"
        simFileName = "sim_info_" + chip + ".txt"
        centerX, centerY, e1, e2, e = readinfo(dataFileName, chip)
        dataCenterX = dataCenterX + centerX
        dataCenterY = dataCenterY + centerY
        dataE1 = dataE1 + e1
        dataE2 = dataE2 + e2
        dataE  = dataE + e


    for chip in CHIPS:
        simFileName = "sim_info_" + chip + ".txt"
        centerX, centerY, e1, e2, e = readinfo(simFileName, chip)
        simCenterX = simCenterX + centerX
        simCenterY = simCenterY + centerY
        simE1 = simE1 + e1
        simE2 = simE2 + e2
        simE  = simE + e


    dataDistance, dataExx, dataStdErr = getCor(dataCenterX, dataCenterY, dataE1, dataE2, dataE)
    simDistance, simExx, simStdErr = getCor(simCenterX, simCenterY, simE1, simE2, simE)

    plt.errorbar(dataDistance, dataExx, dataStdErr, fmt='o', label="data")
    plt.errorbar(simDistance, simExx, simStdErr, fmt='o', label="simulation")
    plt.xlabel("Distance (arcmin)")
    plt.ylabel("Ellipticity correlation")
    plt.title("Elliptictiy correlation")
    plt.legend()

    plt.show()

    #dataObjDict, simulationObjDict = getObjectList(CHIPS)
    #plotCor(dataObjDict,simulationObjDict)


    return

if __name__=='__main__':
    main()


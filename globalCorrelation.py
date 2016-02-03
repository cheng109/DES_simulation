
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


def getCorrelation(e, e1, e2, RA, DEC):
    l = len(e)
    Distance = []
    Exx = []
    for i in range(l):
        for j in np.arange(i+1, l):
            dDEC = DEC[i] - DEC[j]
            dRA = RA[i] - RA[j]
            if dRA<1.0e-6:
                continue
            angle = -2* math.atan(dRA/dDEC)   # in radian
            XA_t=(e1[i]*np.cos(angle)-e2[i]*np.sin(angle))
            XA_x=(e1[i]*np.sin(angle)+e2[i]*np.cos(angle))
            XB_t=(e1[j]*np.cos(angle)-e2[j]*np.sin(angle))
            XB_x=(e1[j]*np.sin(angle)+e2[j]*np.cos(angle))
            Distance.append(np.sqrt(dDEC**2+dRA**2)*60)      # distance is in degree;


            Exx.append( XA_t*XB_t + XA_x*XB_x)

    print max(Distance)
    interval = 2.0   # in arcmin
    Exxbins = [[] for i in range(int(120/interval))]
    #print len(Distance)
    for i in range(len(Exx)):   # bins[0] means interval between [0, 50]

        Exxbins[int(Distance[i]/interval)].append(Exx[i])
    newDistance = []
    newExx = []
    stdErr = []
    for i in range(len(Exxbins)):
        if len(Exxbins[i])>0:
            newDistance.append((i+0.5)*interval)
            newExx.append(np.mean(Exxbins[i]))
            stdErr.append(np.std(Exxbins[i])/np.sqrt(len(Exxbins[i])))



    return newDistance, newExx, stdErr


def main1():

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


def combineFiles(CHIPS, shift, tilt,  type, DIR):
    # type = "_Data_"
    # type = "_Simu_"
    e = []
    e1 = []
    e2 = []
    RA = []
    DEC = []

    prefix = DIR + "output_" + "_x_" + str(shift["x"]) + "_y_" + str(shift["y"]) + "_z_" +str(shift["z"]) +  "_phi_" + str(tilt["phi"])+ "_psi_" + str(tilt["psi"]) + "_theta_" + str(tilt["theta"]) + "_"
    for chip in CHIPS:
        fileName = prefix + "ELL" + type + chip + ".txt"
        f = open(fileName, 'r')
        for line in f.readlines():
            temp = line.split()
            e.append(float(temp[0]))
            e1.append(float(temp[1]))
            e2.append(float(temp[2]))
            RA.append(float(temp[3]))
            DEC.append(float(temp[4]))
    return e, e1, e2, RA, DEC


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


def main():
    CHIPS = ["N1", "N4", "N7", "N22", "S22"]
    #CHIPS  = ["N1"]

    DIR = "new_all_chips_results_tilt/"

    shift = {"x": 0.0, "y": 0.0, "z": 0.0}
    tilt = {"phi": 0.0, "psi": 0.0, "theta": -20}

    data_e, data_e1, data_e2, data_RA, data_DEC = combineFiles(CHIPS, shift, tilt, type = "_Data_", DIR=DIR)
    simu_e, simu_e1, simu_e2, simu_RA, simu_DEC = combineFiles(CHIPS, shift, tilt, type = "_Simu_", DIR=DIR)

    data_Distance, data_Exx, data_StdErr = getCorrelation(data_e, data_e1, data_e2, data_RA, data_DEC)
    simu_Distance, simu_Exx, simu_StdErr = getCorrelation(simu_e, simu_e1, simu_e2, simu_RA, simu_DEC)

    plt.errorbar(data_Distance, data_Exx, data_StdErr, fmt='o', label="data")
    plt.errorbar(simu_Distance, simu_Exx, simu_StdErr, fmt='o', label="simulation")
    plt.xlabel("Distance (arcmin)")
    plt.ylabel("Ellipticity correlation")
    plt.title("Elliptictiy correlation")
    plt.legend()

    plt.show()



if __name__=='__main__':
    main()


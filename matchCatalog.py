import catlogGenerator
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from astropy.io import fits
import math
import imageMatch
import commons
import measurepsf
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from astropy.wcs import WCS

def matchCatalog(dataObjList, simulationObjList):
    distance_tolerance = 20
    new_dataObjList=[]
    new_simulationObjList = []
    for obj in dataObjList:
        corrected_X= float(obj.xcenter)
        corrected_Y = float(obj.ycenter)
        matchID = "none"
        minDistance = 100000
        minDx = 0
        minDy = 0
        for simulation in simulationObjList:
            dy = corrected_Y-float(simulation.ycenter)
            dx = corrected_X-float(simulation.xcenter)
            d = np.sqrt(dx**2 + dy**2)
            if d<minDistance:
                minDistance = d
                minDx = dx
                minDy = dy
                matchID = simulation.id
                matchMag =  simulation.mag
                matchFWHM = simulation.fwhm
        if minDistance < distance_tolerance:
            obj.matchID = matchID
            new_dataObjList.append(obj)
            obj.dx = minDx
            obj.dy = minDy
            obj.matchMag = matchMag
            obj.matchFWHM = matchFWHM

    for obj in new_dataObjList:
        for simObj in simulationObjList:
            if obj.matchID==simObj.id:
                simObj.type=obj.type
                new_simulationObjList.append(simObj)
    print len(new_dataObjList)
    print len(dataObjList)

    return new_dataObjList, new_simulationObjList


def getPeak(binArray, bins):
    start = 0
    maxDensity = 0
    maxBinVal = 0
    for i in range(len(binArray[0])):
        if binArray[0][i] > maxDensity:
            maxDensity = binArray[0][i]
            maxBinVal = i

    peak = bins[maxBinVal]

    if maxBinVal!=0 and maxBinVal!=(len(binArray[0])-1):
        peak = (bins[maxBinVal-1]*binArray[0][maxBinVal-1] + bins[maxBinVal+1]*binArray[0][maxBinVal] + bins[maxBinVal+1]*binArray[0][maxBinVal+1])/\
               (binArray[0][maxBinVal-1] + binArray[0][maxBinVal] + binArray[0][maxBinVal+1])

    return peak

def plotPositionDiff(chip, new_dataObjList, plot ):
    of=open("position.txt", 'w')
    count = 0
    if plot==True:
        f, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2) #, sharex='col', sharey='row')

    dxList = []
    dyList = []
    for obj in new_dataObjList:
        if obj.matchID!="none":
            dxList.append(obj.dx)
            dyList.append(obj.dy)
            if plot==True:
                ax1.plot(obj.xcenter, obj.dx, 'r.')
                ax2.plot(obj.xcenter, obj.dy, 'r.')
                ax3.plot(obj.ycenter, obj.dx, 'r.')
                ax4.plot(obj.ycenter, obj.dy, 'r.')
                of.write(str(obj.xcenter) + "\t" +str(obj.dx) + "\t" + str(obj.xcenter) + "\t" +str(obj.dy) +"\t"\
                        +str(obj.ycenter) + "\t" +str(obj.dx) + "\t" + str(obj.ycenter) + "\t" +str(obj.dy) + "\n")
            count += 1

    of.close()
    bins = [-20]
    dbins = 0.2
    for i in range(200):
        bins.append(bins[i]+dbins)
    retX =  np.histogram(dxList, bins=bins, density=True)
    retY =  np.histogram(dyList, bins=bins, density=True)

    dxList = []
    dyList = []
    xcenter = []
    ycenter = []
    for obj in new_dataObjList:
        dxList.append(obj.dx)
        dyList.append(obj.dy)
        xcenter.append(obj.xcenter)
        ycenter.append(obj.ycenter)
    slopeXX, interceptXX = np.polyfit(xcenter, dxList, 1)
    slopeXY, interceptXY = np.polyfit(ycenter, dxList, 1)
    slopeYX, interceptYX = np.polyfit(xcenter, dyList, 1)
    slopeYY, interceptYY = np.polyfit(ycenter, dyList, 1)

    print "Dx VS X function:   Dx = ", slopeXX, "*X + ", interceptXX
    print "Dx VS Y function:   Dx = ", slopeXY, "*Y + ", interceptXY
    print "Dy VS X function:   Dy = ", slopeYX, "*Y + ", interceptYX
    print "Dy VS Y function:   Dy = ", slopeYY, "*Y + ", interceptYY

    meanDx = np.mean(dxList)
    meanDy = np.mean(dyList)
    print "Mean(dx) = ", meanDx, " std(dx) = ", np.std(dxList)
    print "Mean(dy) = ", meanDy, " std(dy) = ", np.std(dyList)

    peakX = getPeak(retX, bins)
    peakY = getPeak(retY, bins)
    print "dxHist_Max = ", peakX
    print "dyHist_Max = ", peakY


    angle = math.degrees(math.atan(slopeXY))
    print "Rotation Angle: ", math.degrees(math.atan(slopeXY)), " degree"


    lim = 10
    ax1.set_ylim([-lim, lim])
    ax2.set_ylim([-lim, lim])
    ax3.set_ylim([-lim, lim])
    ax4.set_ylim([-lim, lim])
    if plot==True:
        ax1.set_xlabel("X (pixels)")
        ax1.set_ylabel("dX (pixels)")
        ax2.set_xlabel("X (pixels)")
        ax2.set_ylabel("dY (pixels)")
        ax3.set_xlabel("Y (pixels)")
        ax3.set_ylabel("dX (pixels)")
        ax4.set_xlabel("Y (pixels)")
        ax4.set_ylabel("dY (pixels)")

        plt.show()
        #f.savefig(chip + "position.pdf")

    return peakX, peakY, angle


def createMissMatchRegionFile(new_dataObjList, regionFileName):
    regionFile = open(regionFileName, 'w')
    regionFile.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nimage\n")
    for obj in new_dataObjList:
        if obj.matchID=="none":
            regionFile.write("circle("+str(obj.xcenter)+","+str(obj.ycenter)+",10)\n")

    regionFile.close()

def plotMagDiff(new_dataObjList, new_simulationObjList):
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    magDiff = []
    f2=open("mag.txt",'w')
    ymag = []
    xmag = []

    for obj in new_dataObjList:
        if obj.matchID!="none":
            ymag.append(obj.matchMag)
            xmag.append(obj.mag)
            ax1.plot(obj.mag, obj.matchMag, 'b.')

            ax2.plot(obj.mag, obj.mag-obj.matchMag, 'r.')

            ax3.plot(obj.mag, obj.class_star, 'b.')
            f2.write(str(obj.mag)+"\t" + str(obj.matchMag) + "\t" + str(obj.mag) + "\t" + str(obj.mag-obj.matchMag) + "\n")
            magDiff.append(obj.mag-obj.matchMag)
    slopeMag, interceptMag = np.polyfit(xmag, ymag,1)
    print "Mag_simulation VS Mag_data function:   Y = ", slopeMag, "*X + (", interceptMag, ")"
    #ax1.plot([14,14],[25,25],'k-')
    ax1.set_xlabel("Mag_data")
    ax1.set_ylabel("Mag_simuation")
    #ax2.plot([14,0],[25,0],'k-')
    ax2.set_xlabel("Mag_data")
    ax2.set_ylabel("Mag_data- Mag_simulation")
    ax3.set_xlabel("Magnitude")
    ax3.set_ylabel("Star_Class")
    f2.close()
    fwhmDiff = []

    data_fwhm = []
    simulation_fwhm = []
    f3 = open("PSF_size.txt", 'w')
    for obj in new_dataObjList:
        if obj.matchID!="none" and obj.type=="star": # and obj.fwhm<5 and obj.matchFWHM<5:
            #print obj.id, obj.fwhm, obj.matchID, obj.matchFWHM
            data_fwhm.append(obj.fwhm)
            simulation_fwhm.append(obj.matchFWHM)
            ax4.plot(obj.fwhm, obj.matchFWHM, 'b.')
            f3.write(str(obj.fwhm) + "\t" + str(obj.matchFWHM) + "\n" )
            fwhmDiff.append(obj.fwhm-obj.matchFWHM)
    print "PSF_data = ", np.mean(data_fwhm)
    print "PSF_simulation = ", np.mean(simulation_fwhm)
    f3.close()
    ax4.set_xlabel("STAR_fwhm_data")
    ax4.set_ylabel("STAR_fwhm_simuation")

    ax4.set_xlim([3.2,4.2])
    ax4.set_ylim([3.2,4.2])



    #plt.show()



def plotStarSizeDiff(new_dataObjList, new_simulationObjList):
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fwhmDiff = []

    data_fwhm = []
    simulation_fwhm = []
    for obj in new_dataObjList:
        if obj.matchID!="none" and obj.type=="star": # and obj.fwhm<5 and obj.matchFWHM<5:
            #print obj.id, obj.fwhm, obj.matchID, obj.matchFWHM
            data_fwhm.append(obj.fwhm)
            simulation_fwhm.append(obj.matchFWHM)
            ax1.plot(obj.fwhm, obj.matchFWHM, 'b.')
            ax2.plot(obj.fwhm, obj.fwhm-obj.matchFWHM, 'r.')
            #ax3.plot(obj.mag, obj.fwhm,'r.')
            #ax4.plot(obj.matchFWHM, obj.matchMag, 'r.')
            fwhmDiff.append(obj.fwhm-obj.matchFWHM)
    print "PSF_data = ", np.mean(data_fwhm)
    print "PSF_simulation = ", np.mean(simulation_fwhm)

    for obj in new_simulationObjList:
        if obj.type=="star":
            #if obj.fwhm>4 and obj.mag<15:
            #print obj.id, obj.fwhm, obj.mag
            ax4.plot(obj.mag, obj.fwhm, 'r.')

    for obj in new_dataObjList:
        if obj.type=="star":
            ax3.plot(obj.mag, obj.fwhm,'r.')

    ax1.set_xlabel("STAR_fwhm_data")
    ax1.set_ylabel("STAR_fwhm_simuation")
    #ax1.set_xlim([4,6])
    #ax1.set_ylim([3,4])
#

    #ax2.plot([14,0],[25,0],'k-')
    ax2.set_xlabel("STAR_fwhm_data")
    ax2.set_ylabel("STAR_fwhm_data- STAR_fwhm_simulation")
    #ax2.set_xlim([3.4, 4])
    #ax2.set_ylim([-0.6,0.2])

    ax3.set_ylabel("FWHM_data")
    ax3.set_xlabel("Mag_data")
    #ax3.set_xlim([14.5,18])
    #ax3.set_ylim([3.4, 4])

    ax4.set_ylabel("FWHM_simulation")
    ax4.set_xlabel("Mag_simulation")
    #ax4.set_xlim([14.5,18])
    #ax4.set_ylim([3.7, 4.3])
    #plt.show()

def getEllipticityDiff(new_dataObjList, new_simulationObjList, dataFileName, simFileName):
    fdataEll = open("dataEll.txt", 'w')
    fsimEll = open("simEll.txt", 'w')

    dataEllipse = []
    dataFWHM = []
    dataPA = []
    dataRA = []
    dataDEC = []

    simulationElipse = []
    simulationFWHM = []
    simulationPA = []
    simulationRA = []
    simulationDEC = []

    dataFlux = []
    simulationFlux = []

    eLimit = 0.35

    #wcsData = WCS(dataFileName)
    #wcsSim  = WCS(simFileName)
    for obj in new_dataObjList:
        if obj.type=="star" and obj.new_e < eLimit and obj.new_e > 0.000001:
            dataFlux.append(obj.flux)
            dataEllipse.append(obj.new_e)
            dataPA.append(obj.new_angle)
            dataFWHM.append(obj.fwhm)
            simulationFWHM.append(obj.matchFWHM)

            dataRa, dataDec = commons.pix2world(simFileName, obj.xcenter, obj.ycenter)

            dataRA.append(dataRa)
            dataDEC.append(dataDec)
            simulationRA.append(dataRa)
            simulationDEC.append(dataDec)


    for obj in new_simulationObjList:
        if obj.type=="star"and obj.new_e < eLimit and obj.new_e > 0.000001:
            simulationElipse.append(obj.new_e)
            simulationFlux.append(obj.flux)
            simulationPA.append(obj.new_angle)



    for e in dataEllipse:
        fdataEll.write(str(e) + "\n")
    for e in simulationElipse:
        fsimEll.write(str(e)+"\n")
    fdataEll.close()
    fsimEll.close()
    print "==================="
    print "meanDataE: ", np.mean(dataEllipse), np.std(dataEllipse)
    print "meanSimE: ", np.mean(simulationElipse), np.std(simulationElipse)
    print "==================="
    print "meanDataPA:", np.mean(dataPA), "\tSTD: ", np.std(dataPA)
    print "meanSimPA: ", np.mean(simulationPA), "\tSTD: ", np.std(simulationPA)
    print "==================="
    return dataEllipse, simulationElipse, dataFlux, simulationFlux, dataPA, simulationPA, dataFWHM, simulationFWHM, dataRA, dataDEC, simulationRA, simulationDEC



def getCorrelation(ObjList, outFileName, infoFileName):
    fCor = open(outFileName, 'w')
    multiFile = open(infoFileName, 'w')

    l = len(ObjList)
    Exx = []
    Distance = []

    #### write star ellipticity information to a file #####

    for i in range(l):
        if ObjList[i].type =="star":
            starA = ObjList[i]
            multiFile.write(str(starA.xcenter) + "\t" + str(starA.ycenter) + "\t" + str(starA.new_e1) + "\t" + str(starA.new_e2) + "\t" + str(starA.new_e) + "\t"+ "\n")


    multiFile.close()
    #########################################################
    for i in range(l):
        if ObjList[i].type =="star":
            starA = ObjList[i]
            for j in np.arange(i+1, l-i-1):
                if ObjList[j].type =="star":
                    starB = ObjList[j]
                    dx = starA.xcenter - starB.xcenter
                    dy = starA.ycenter - starB.ycenter
                    if dx<1.0e-6:
                        continue
                    angle = -2* math.atan(dy/dx)   # in radian

                    XA_t=(starA.new_e1*np.cos(angle)-starA.new_e2*np.sin(angle))
                    XA_x=(starA.new_e1*np.sin(angle)+starA.new_e2*np.cos(angle))
                    XB_t=(starB.new_e1*np.cos(angle)-starB.new_e2*np.sin(angle))
                    XB_x=(starB.new_e1*np.sin(angle)+starB.new_e2*np.cos(angle))
                    Distance.append( np.sqrt(dx**2+dy**2))
                    Exx.append( XA_t*XB_t + XA_x*XB_x)

                    #Exx.append( XA_t*XB_t + XA_x*XB_x)

    for i in range(len(Distance)) :
          fCor.write(str(Distance[i]) + "\t" + str(Exx[i]) + "\n")
    fCor.close()
    interval = 300
    Exxbins = [[] for i in range(int(10000/interval))]

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

def getWhisker(objList, R):
    start =[]
    end = []
    #length = []
    center = []
    counter = 0
    for obj in objList:
        if obj.type=="star" and obj.new_e<0.50 and obj.new_e>0.000001:
            #angle = math.atan(obj.e2/obj.e1)
            counter = counter +  1
            PA = obj.new_angle
            e = obj.new_e
            start.append((obj.xcenter- R*e*np.cos(PA),obj.ycenter-R*e* np.sin(PA)))
            end.append((obj.xcenter+ R*e*np.cos(PA),obj.ycenter+e* np.sin(PA)))
            center.append((obj.xcenter, obj.ycenter))
    print "counter" , counter
    return start, end, center

def plotEllipticityAndCorrelation(dataImageName, simulationImageName, new_dataObjList, new_simulationObjList, writeStamps):

    new_dataObjList         = updateEllipticity(dataImageName, new_dataObjList, type='data', writeStamps=writeStamps)
    new_simulationObjList   = updateEllipticity(simulationImageName, new_simulationObjList, type='simulation', writeStamps=writeStamps)

    dataDistance, dataExx, dataStdErr   = getCorrelation(new_dataObjList, "dataCor.txt", "data_info.txt")
    simDistance, simExx , simStdErr     = getCorrelation(new_simulationObjList,  "simCor.txt", "sim_info.txt")

    print simulationImageName
    #f = open(simulationImageName+".txt", 'w')
    #for i in range(len(simDistance)):
    #    f.write(str(simDistance[i]) +"\t" + str(simExx[i]) + "\t" + str(simStdErr[i]) + "\n")
    #f.close()

    dataEllipse, simulationElipse, dataFlux, simulationFlux, dataPA, simulationPA, dataFWHM, simulationFWHM,  dataRA, dataDEC, simulationRA, simulationDEC = getEllipticityDiff(new_dataObjList, new_simulationObjList, dataImageName, simulationImageName)
    #####  PLOT ELLIPTICITY VS  FLUX ################
    #plt.plot(dataFlux, dataEllipse, 'bo', label='Data' )
    #plt.plot(simulationFlux, simulationElipse, 'og', label='Simulation')
    #plt.legend()
    #plt.xlabel('Flux')
    #plt.ylabel('Ellipticity')
    #plt.xscale('log')
    #plt.show()
    #f, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2, 3)
    f2 = plt.figure()

    gs = gridspec.GridSpec(2, 3,
                       width_ratios=[1,1, 2],
                       height_ratios=[1,1, 1]
                       )

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = plt.subplot(gs[3])
    ax5 = plt.subplot(gs[4])
    ax6 = plt.subplot(gs[5])

    data_fwhm = []
    simulation_fwhm = []
    fwhmDiff = []
    for obj in new_dataObjList:
        if obj.matchID!="none" and obj.type=="star": # and obj.fwhm<5 and obj.matchFWHM<5:
            #print obj.id, obj.fwhm, obj.matchID, obj.matchFWHM
            data_fwhm.append(obj.fwhm)
            simulation_fwhm.append(obj.matchFWHM)
            ax5.plot(obj.fwhm, obj.matchFWHM, 'b.')

            fwhmDiff.append(obj.fwhm-obj.matchFWHM)
    print "PSF_data = ", np.mean(data_fwhm)
    print "PSF_simulation = ", np.mean(simulation_fwhm)

    ax5.set_xlabel("FWHM_data")
    ax5.set_ylabel("FWHM_sim")
    ax5.set_xlim([3.8,4.8])
    ax5.set_ylim([3.6,4.6])




    ax1.hist(dataEllipse, 30, normed=1, facecolor='blue', label='Data')
    ax1.hist(simulationElipse, 30, normed=1, facecolor='green', label='Simulation', alpha=0.5)
    ax1.set_xlabel("Ellipticity ")
    ax1.legend()

    ax4.hist(dataPA, 30, normed=1, facecolor='blue', label='Data')
    ax4.hist(simulationPA, 30, normed=1, facecolor='green', label='Simulation', alpha=0.5)
    ax4.set_xlabel("Position Angle (degree)")
    ax4.legend()



    #ax2.set_yscale('log')
    ax2.errorbar(dataDistance, dataExx, dataStdErr, fmt='o', label="data")
    ax2.errorbar(simDistance, simExx, simStdErr, fmt='o', label="simulation")
    ax2.set_xlabel("Distance (arcmin)")
    ax2.set_ylabel("Ellipticity correlation")
    ax2.set_title("Elliptictiy correlation")
    ax2.legend()

    #R = 1000
    R = 400
    dataStart, dataEnd, dataCenter = getWhisker(new_dataObjList, R)
    simStart, simEnd, simCenter    = getWhisker(new_simulationObjList, R)

    for i in range(len(dataStart)):
        ax6.plot( [dataStart[i][1], dataEnd[i][1]],[dataStart[i][0], dataEnd[i][0]], markersize=0, label="data", color='b', linewidth=2.0)
        #ax3.plot(dataCenter[i][1], dataCenter[i][0], '.', color="red")
    for i in range(len(simStart)):
        ax3.plot( [simStart[i][1], simEnd[i][1]],[simStart[i][0], simEnd[i][0]], markersize=0, label="simulation", color='g', linewidth=2.0)
        #ax4.plot(simCenter[i][1], simCenter[i][0], '.', color="green")
    ax6.set_xlabel("X (pixel)")
    ax6.set_ylabel("Y (pixel)")
    ax6.set_ylim([0,2100])
    ax6.set_xlim([0,4200])
    ax6.set_title("Data Ellipticity")
    #ax3.axes().set_aspect('equal')


    ax3.set_xlabel("X (pixel)")
    ax3.set_ylabel("Y (pixel)")
    ax3.set_ylim([0,2100])
    ax3.set_xlim([0,4200])
    ax3.set_title("Simulation Ellipticity")

    #plt.show()

    return new_dataObjList, new_simulationObjList



def updateEllipticity(imageName, ObjList, type, writeStamps):

    convergeCounter = 0
    data = measurepsf.readFitsImage(imageName)
    dirName = "stampes/"
    dataDict = {}
    winSizeX = 8
    winSizeY = 8
    new_ObjList = []

    alphaX = []
    alphaY = []

    print "Test_Test_Test"
    for obj in ObjList:
        #obj.type=="gauss"
        if obj.type=="star"  and  not (obj.xcenter <15 and obj.ycenter> 955 and obj.ycenter<980 ):#and obj.ycenter-winSizeY>0 and obj.ycenter+winSizeY < 4000 and obj.xcenter-winSizeX>0 and obj.xcenter+winSizeX < 2000:
            yrange_start = int(obj.ycenter-winSizeY)
            yrange_end = int(obj.ycenter+winSizeY)

            xrange_start = int(obj.xcenter-winSizeX)
            xrange_end = int(obj.xcenter+winSizeX)

            if obj.ycenter-winSizeY<=0:
                yrange_start = 0
            if obj.ycenter+winSizeY >= 4000:
                yrange_end = 4000
            if obj.xcenter-winSizeX<=0:
                xrange_start = 0
            if obj.xcenter+winSizeX >= 2000:
                xrange_end = 2000

            if yrange_end -  yrange_start >8 and xrange_end- xrange_start >8:

                icon = data[yrange_start: yrange_end,xrange_start: xrange_end]
                background = 0
                if type=="data":
                    background=90
                icon = commons.rmBackground(icon, background)
                obj.updateEllipticity(icon)
                # print obj.new_e, "\t", obj.alphax,"\t", obj.alphay,"\t", np.sqrt(obj.alphax**2+ obj.alphay**2)
                if(obj.converge == True):
                    convergeCounter +=1

                if type=='data':
                     filter = [] #[719,1094,212,204,555,599 ] #1123, 197, 1281, 364]
                if type=='simulation':
                     filter = [] #[1074 ]
                filter = [str(x) for x in filter]
                if obj.id not in filter:

                    # if(writeStamps==True):
                    #     hdu = fits.PrimaryHDU(icon)
                    #     hdulist = fits.HDUList([hdu])
                    #     hdulist.writeto(dirName+str(type + "_" +str(format(obj.new_e,'.3f')) +"_" +obj.id+".fits"), clobber=True)
                    new_ObjList.append(obj)


   # print "Ratio of convergence: ", convergeCounter/len(new_ObjList)
    return new_ObjList

def transformationMatrix(new_dataObjList):
    numStars = 20
    dataPoints = []
    simPoints = []
    counter = 0
    for obj in new_dataObjList:
        if obj.type=="star" and obj.mag<17 and counter < numStars:
            counter +=1
            dataPoints.append([obj.xcenter, obj.ycenter])
            simPoints.append([obj.xcenter-obj.dx, obj.ycenter-obj.dy])

    Rotation, Shift = imageMatch.rigid_transform_3D(np.matrix(dataPoints), np.matrix(simPoints))
    return Rotation, Shift


def twoPoints(A, B):
    [x0, y0] = A
    [x1, y1] = B
    # Assume Y = a*X + b
    a = (y1-y0)/(x1-x0)
    b = y1 - a*x1
    return -b/a




def singleMatch(chipID, dataImageName, simulationImageName, simulationCatalog, dataCatalog, plot, writeStamps):

    commons.copyWCS(dataImageName,simulationImageName)
    subprocess.call("/opt/local/bin/sex -c data.sex " + dataImageName + " -CATALOG_NAME " + dataCatalog, shell=True)
    subprocess.call("/opt/local/bin/sex -c sim.sex " + simulationImageName + " -CATALOG_NAME " + simulationCatalog, shell=True)



    dataObjList     = catlogGenerator.readSexCatalog(chipID, dataCatalog)
    simulationObjList   = catlogGenerator.readSexCatalog(chipID, simulationCatalog)

    print "Number of objects in data images: ",  len(dataObjList)
    print "Number of objects in simulation images: ",  len(simulationObjList)

    new_dataObjList, new_simulationObjList = matchCatalog(dataObjList, simulationObjList)
    #plotPositionDiff(chipID, new_dataObjList, plot=True )
    #plotMagDiff(new_dataObjList, new_simulationObjList)


    new_dataObjList, new_simulationObjList = plotEllipticityAndCorrelation(dataImageName, simulationImageName, new_dataObjList, new_simulationObjList, writeStamps)
    return new_dataObjList, new_simulationObjList





def matchCHIPS(CHIPS, simName, prefix):
    #prefix = "output_"
    for i in range(len(CHIPS)):

        data_PSF_file = open(prefix + "PSF_Data_" + CHIPS[i] + ".txt", 'w')
        data_ELL_file = open(prefix +"ELL_Data_"  + CHIPS[i] + ".txt", 'w')
        data_PA_file  = open(prefix +"PAA_Data_"   + CHIPS[i] + ".txt", 'w')

        sim_PSF_file  = open(prefix +"PSF_Simu_"   + CHIPS[i] + ".txt", 'w')
        sim_ELL_file  = open(prefix +"ELL_Simu_"   + CHIPS[i] + ".txt", 'w')
        sim_PA_file   = open(prefix +"PAA_Simu_"    + CHIPS[i] + ".txt", 'w')

        dataImageName = CHIPS[i] + "_test_image.fits"
        new_dataObjList, new_simulationObjList = singleMatch(CHIPS[i], dataImageName=dataImageName, simulationImageName=simName[i], simulationCatalog = CHIPS[i] +"_simCatalog", dataCatalog=CHIPS[i]+"_dataCatalog", plot=True, writeStamps=False)
        dataEllipse, simulationEllipse, dataFlux, simulationFlux, dataPA, simulationPA, dataFWHM, simulationFWHM,  dataRA, dataDEC, simulationRA, simulationDEC \
            = getEllipticityDiff(new_dataObjList, new_simulationObjList, dataImageName, simName[i])

        for i in range(len(simulationEllipse)):
            sim_ELL_file.write(str(simulationEllipse[i]) + "\t")
            sim_ELL_file.write(str(simulationRA[i]) + "\t")
            sim_ELL_file.write(str(simulationDEC[i]) + "\n")
        for i in range(len(simulationFWHM)):
            sim_PSF_file.write(str(simulationFWHM[i]) + "\n")
        for i in range(len(simulationPA)):
            sim_PA_file.write(str(simulationPA[i]) + "\n")


        for i in range(len(dataEllipse)):
            data_ELL_file.write(str(dataEllipse[i]) + "\t")
            data_ELL_file.write(str(dataRA[i]) + "\t")
            data_ELL_file.write(str(dataDEC[i]) + "\n")
        for i in range(len(dataFWHM)) :
            data_PSF_file.write(str(dataFWHM[i]) + "\n")
        for i in range(len(dataPA)):
            data_PA_file.write(str(dataPA[i]) + "\n")



        data_ELL_file.close()
        data_PSF_file.close()
        data_PA_file.close()

        sim_ELL_file.close()
        sim_PSF_file.close()
        sim_PA_file.close()


def txtToList(fileName):
    txtFile = open(fileName, 'r')
    list = []
    for line in txtFile.readlines():
        temp = line.split()
        list.append(float(temp[0]))

    return list

def plotMultCHIPS(CHIPS, prefix, PSF_plot , ELL_plot):
    #prefix  = "output_"
    numEllBins = 40
    numPSFBins = 40

    minEll = 0.0
    maxEll = 0.3
    minPSF = 3.0
    maxPSF = 6.0
    ellBins= []
    psfBins= []
    for i in range(numEllBins):
        ellBins.append(i*(maxEll-minEll)/numEllBins)
    for i in range(numPSFBins):
        psfBins.append(minPSF + i*(maxPSF-minPSF)/numPSFBins)

    ## ax1:   Ellipticity
    ## ax11:  PSF_FWHM

    if ELL_plot==True:
        f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)
        commons.removeSubplotPanel(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)

    if PSF_plot==True:
        f1, ((ax11, ax22, ax33), (ax44, ax55, ax66), (ax77, ax88, ax99)) = plt.subplots(3, 3)
        commons.removeSubplotPanel(ax11, ax22, ax33, ax44, ax55, ax66, ax77, ax88, ax99)

    for chip in CHIPS :

        data_PSF_file = prefix + "PSF_Data_" + chip + ".txt"
        data_ELL_file = prefix +"ELL_Data_"  + chip + ".txt"
        sim_PSF_file  = prefix +"PSF_Simu_"   + chip + ".txt"
        sim_ELL_file  = prefix +"ELL_Simu_"   + chip + ".txt"

        psfData = txtToList(data_PSF_file)
        ellData = txtToList(data_ELL_file)
        psfSim  = txtToList(sim_PSF_file)
        ellSim  = txtToList(sim_ELL_file)

        if chip=="N22":
            if ELL_plot==True:
                ax2.hist(ellData, bins = ellBins,  normed=1, color='blue',  label = "Ell_Data",  alpha=1.0)
                ax2.hist(ellSim,  bins = ellBins,  normed=1, color='green', label = "Ell_Sim",  alpha=0.5)
                ax2.set_title(chip)
            if PSF_plot==True:
                ax22.hist(psfData, bins = psfBins,  normed=1, color='blue',  label = "PSF_Data",  alpha=1.0)
                ax22.hist(psfSim,  bins = psfBins,  normed=1, color='green', label = "PSF_Sim",  alpha=0.5)
                ax22.set_title(chip)

        if chip=="N1":
            if ELL_plot==True:
                ax4.hist(ellData, bins = ellBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax4.hist(ellSim,  bins = ellBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax4.set_title(chip)
            if PSF_plot==True:
                ax44.hist(psfData, bins = psfBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax44.hist(psfSim,  bins = psfBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax44.set_title(chip)

        if chip=="N4":
            if ELL_plot==True:
                ax5.hist(ellData, bins = ellBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax5.hist(ellSim,  bins = ellBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax5.set_title(chip)

            if PSF_plot==True:
                ax55.hist(psfData, bins = psfBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax55.hist(psfSim,  bins = psfBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax55.set_title(chip)

        if chip=="N7":
            if ELL_plot==True:
                ax6.hist(ellData, bins = ellBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax6.hist(ellSim,  bins = ellBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax6.set_title(chip)

            if PSF_plot==True:
                ax66.hist(psfData, bins = psfBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax66.hist(psfSim,  bins = psfBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax66.set_title(chip)

        if chip=="S22":
            if ELL_plot==True:
                ax8.hist(ellData, bins = ellBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax8.hist(ellSim,  bins = ellBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax8.set_title(chip)

            if PSF_plot==True:
                ax88.hist(psfData, bins = psfBins,  normed=1, color='blue',  label = chip,  alpha=1.0)
                ax88.hist(psfSim,  bins = psfBins,  normed=1, color='green', label = chip,  alpha=0.5)
                ax88.set_title(chip)


        if ELL_plot==True:
            commons.showLegend([ax2])
        if PSF_plot==True:
            commons.showLegend([ax22])

        ax1.set_xlim([minPSF, maxPSF])
        ax3.set_xlim([minPSF, maxPSF])
        ax2.set_xlim([minEll,maxEll])
        ax4.set_xlim([minEll,maxEll])


    plt.show()

def main():
    ### preface #####
    DIR  = "new_all_chips_results_tilt/"
    CHIPS = []
    simNameList = []
    shift = {"x": 0.0, "y": 0.0, "z": 0.0}
    tilt = {"phi": 0.0, "psi": 0.0, "theta": -20}
    #tilt = {"phi": 324000.0, "psi": 0.0, "theta": 40}
    CHIPS.append("N1")
    #CHIPS.append("N4")
    #CHIPS.append("N7")
    #CHIPS.append("N22")
    #CHIPS.append("S22")

    CHIPS = ["N1", "N4", "N7", "N22", "S22"]

    for chip in CHIPS :
        simName = DIR + "Images_" + chip + "_x_" + str(shift["x"]) + "_y_" + str(shift["y"]) + "_z_" + str(shift["z"]) + "_phi_" + str(tilt["phi"])+ "_psi_" + str(tilt["psi"]) + "_theta_" + str(tilt["theta"])+ ".fits"
        simNameList.append(simName)


    ### Operation ######
    prefix = DIR + "output_" + "_x_" + str(shift["x"]) + "_y_" + str(shift["y"]) + "_z_" +str(shift["z"]) +  "_phi_" + str(tilt["phi"])+ "_psi_" + str(tilt["psi"]) + "_theta_" + str(tilt["theta"]) + "_"
    matchCHIPS(CHIPS, simNameList, prefix = prefix)
    plotMultCHIPS(CHIPS, prefix = prefix, PSF_plot=False, ELL_plot=True)


if __name__=='__main__':
    main()

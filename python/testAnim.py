#
# testAnim.py
#

# 2018-09-09 WIC - get familiar with matplotlib's animation
# interface. It may be easier than my old practice of dumping every
# frame to an image and using other software to stitch together the
# movie...

# On my laptop, ffmpeg creates blocky output but is stil uesful for
# checking the output. For production, Quicktime 7's open image
# sequence functionality should work well.
#
#

import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
from matplotlib import colors

# for importing data
from astropy.table import Table

# for trend interpolation
from scipy.interpolate import interp1d

# for coordinates
from astropy.coordinates import SkyCoord
from astropy import units as u

# for importing the trend
import cPickle as pickle

# for drawing the convex hull around the objects
import polyUtils as pu

import os, sys
import subprocess
import glob

# 2018-09-09 - the backend save doesn't seem to work. Recurrent error
# "I/O on closed file". So we try a different approach...

# 2018-09-15 - updated to use real data

plt.ion()

class Points(object):

    """Points object for animation"""

    def __init__(self, nSteps=50, tStep=0.1, nObj=25, figNum=1, \
                     cScatt='b', dirFrames='./TMPframes', \
                     frameByFrame=True, \
                     frameStem='frame', \
                     frameTyp='png', \
                     filMov='TEST.mp4', \
                     framerate=25, \
                     cleanFrames=False, \
                     labelX='X', \
                     labelY='Y', \
                     xLims = np.array([]), \
                     yLims = np.array([]), \
                     distEffects=True, \
                     figSz=[6,6], \
                     doGrid=True, \
                     labelTimes=True, \
                     fszLabel=20, \
                     pathData='', \
                     pathTrend='', \
                     showHull=False):

        # timestep
        self.tStep = np.copy(tStep)
        self.nSteps = np.copy(nSteps)
        
        # the figure object
        self.fig = None
        self.ax = None
        self.figNum = figNum
        self.figSz=np.copy(figSz)

        # label information
        self.labelX = labelX[:]
        self.labelY = labelY[:]

        # show convex hull?
        self.showHull = showHull

        # do grid?
        self.doGrid = doGrid

        # path for input data
        self.pathData = pathData[:]
        self.tData = Table()

        # path for trends
        self.pathTrend = pathTrend[:]
        self.tTrend = Table()

        # interpolation object for trend
        self.fTrend = None

        # label the times?
        self.labelTimes = labelTimes
        
        # fontsize for labels
        self.fszLabel = fszLabel

        # scatter object to update
        self.scatt = None

        # plot color for scatter
        self.cScatt = cScatt[:]

        # switch - are we doing this frame-by-frame?
        self.frameByFrame = frameByFrame
        self.frameName = 'TEST_FRAME.jpg'
        self.frameStem = frameStem[:]
        self.frameTyp = frameTyp[:]

        # list of frames generated
        self.lFrames = []

        # clean the intermediate frames?
        self.cleanFrames = cleanFrames

        # movie file name
        self.filMov=filMov[:]

        # for data plotting
        self.xLims = np.copy(xLims)
        self.yLims = np.copy(yLims)

        # add distance effects to the scatterplot?
        self.distEffects = distEffects
        self.szMin = 2.
        self.szMax = 49.
        self.alphaMin = 0.2
        self.alphaMax = 0.8

        # for data generation
        self.nObj = np.copy(nObj)
        self.xMax = 4096.
        self.xMin = 0.
        self.yMax = 4096.
        self.yMin = 0.

        self.vXmin = -10.
        self.vXmax = 10.
        self.vYmin = -2.
        self.vYmax = +2.
        self.vXstd = 2.0
        self.vYstd = 0.5

        # distances (to generate fake rotation curve)
        self.dists = np.array([])
        self.dMin = 7.
        self.dMax = 11.

        # starting points and velocities
        self.xStart = np.array([])
        self.yStart = np.array([])
        self.vX = np.array([])
        self.vY = np.array([])

        # output for frames
        self.dirFrames = dirFrames[:]

        # command for the SYSTEM ffmpeg 
        self.FFMPEG = 'ffmpeg'
        self.fmtFFMPG = '%s_%03d.%s' 
        self.framerate = framerate
        self.nZeros = -1

    def genStartingPoints(self):

        """Generates starting points"""
        
        self.xStart = np.random.uniform(self.xMin, self.xMax, self.nObj)
        self.yStart = np.random.uniform(self.yMin, self.yMax, self.nObj)
        
        # for generation of velocities
        self.dists = np.random.uniform(self.dMin, self.dMax, self.nObj)

    def genVelocities(self):

        """Generates velocities"""

        # random component
        randX = np.random.normal(0., self.vXstd, self.nObj)
        randY = np.random.normal(0., self.vYstd, self.nObj)
        
        # linear component
        gradientX = (self.vXmax - self.vXmin) / (self.dMax - self.dMin)
        gradientY = (self.vYmax - self.vYmin) / (self.dMax - self.dMin)

        dMid = 0.5*(self.dMin + self.dMax)
        linrX = gradientX * (self.dists - dMid)
        linrY = gradientY * (self.dists - dMid)

        # add them together
        self.vX = linrX + randX
        self.vY = linrY + randY 

    def loadData(self):

        """Loads the data from file"""

        if len(self.pathData) < 3:
            return

        if not os.access(self.pathData, os.R_OK):
            print("Points.dataFrompath WARN - cannot read path %s" \
                % (self.pathData))
            return

        self.tData = Table.read(self.pathData)
        
    def dataFromPath(self):

        """Loads data and slots into the expected variables"""

        self.loadData()
        if len(self.tData) < 5:
            return

        # look for "good" objects
        bGood = np.abs(self.tData['mul_SWEEPS'] ) < 50.

        self.tData = self.tData[bGood]

        ra = self.tData['RA']
        de = self.tData['DEC']

        # convert to galactics
        c = SkyCoord(ra = np.asarray(ra)*u.deg, \
                         dec=np.asarray(de)*u.deg, frame='fk5')
        l = c.galactic.l.degree
        b = c.galactic.b.degree

        #self.xStart = self.tData['RA']
        #self.yStart = self.tData['DEC']
        self.xStart = np.copy(l)
        self.yStart = np.copy(b)
        self.vX = self.tData['mul_SWEEPS']
        self.vY = self.tData['mub_SWEEPS']

        self.labelX=r"Galactic longitude $l$, degrees"
        self.labelY=r"Galactic latitude $b$, degrees"

        # if we're using mas/yr, then we have a scale factor to apply
        self.vX /= 3.6e7 # to convert from mas/yr to degrees per year
        self.vY /= 3.6e7
        
        distCol = 'd_FromLoI'
        if self.pathData.find('ich') > -1:
            distCol = 'd_FromHiI'

        self.dists = self.tData[distCol]

        # limits
        self.xMin = np.min(self.xStart)
        self.xMax = np.max(self.xStart)
        self.yMin = np.min(self.yStart)
        self.yMax = np.max(self.yStart)

        # some hard limits for zoomed plot
        self.xMin = 1.24
        self.xMax = 1.27
        self.yMin = -2.67
        self.yMax = -2.64

        # some size limits for the plot
        self.szMin=2
        self.szMax=36
        
    def loadTrend(self):

        """Load the trend data from disk"""

        if len(self.pathTrend) < 3:
            return

        if not os.access(self.pathTrend, os.R_OK):
            print("Points.trendFromData WARN - cannot read path %s" \
                      % (self.pathTrend))
            return

        DDum = pickle.load(open(self.pathTrend, 'r'))
        if self.pathData.find('ich') > -1:
            self.tTrend = DDum['Metal-rich']
        else:
            self.tTrend = DDum['Metal-poor']
                 
        # convert the trends into degrees per year
        self.tTrend['muL'] /= 3.6e7
        self.tTrend['muB'] /= 3.6e7
             
        # correct the midpoint offsets
        if self.pathData.find('oor') > -1:
            self.tTrend['muL'] += 0.5e-8
        else:
            self.tTrend['muL'] += 1.0e-8

    def replaceVelWithTrend(self):

        """Finds the trend and replaces velocities with the trend values"""

        lTrend = np.argsort(self.tTrend['dMod'])

        if self.pathData.find('oor') > -1:
            bVis = self.tTrend['dMod'][lTrend] < 1.0
            self.tTrend['muL'] += 1.0e-8
        else:
            bVis = np.isfinite(self.tTrend['muL'])

        self.fTrendL = interp1d(self.tTrend['dMod'][lTrend][bVis], \
                                    self.tTrend['muL'][lTrend][bVis], \
                                    fill_value='extrapolate', \
                                    kind='slinear')

        self.fTrendB = interp1d(self.tTrend['dMod'][lTrend], \
                                    self.tTrend['muB'][lTrend], \
                                    fill_value='extrapolate', \
                                    kind='slinear')

        self.vX = self.fTrendL(self.dists)
        self.vY = self.fTrendB(self.dists)

    def setHull(self):

        """Sets the convex hull around the dataset"""

        self.xHull, self.yHull = pu.GetHull(self.xStart, self.yStart)

    def makeTimes(self):

        """Generates the times"""
        
        tMax = self.tStep * self.nSteps
        self.times = np.arange(0., tMax, self.tStep)
        
    def makeFigure(self):
        
        """Makes the figure"""

        # NOTE - when working from data, will need to ensure
        # self.xMin, etc. are set from the data.

        # what limits to use?
        if np.size(self.xLims) < 2:
            xLimits = (self.xMin, self.xMax)
        else:
            xLimits = np.copy(self.xLims)

        if np.size(self.yLims) < 2:
            yLimits = (self.yMin, self.yMax)
        else:
            yLimits = np.copy(self.yLims)

        
        self.fig = plt.figure(self.figNum)
        self.fig.set_size_inches(self.figSz, forward=True)
        self.fig.clf()
        self.ax = self.fig.add_subplot(111, \
                                           xlim=xLimits,\
                                           ylim=yLimits)

    def plotFirstScatt(self):

        """Plots the first scatter frame"""

        # do distance effects?
        sz = np.repeat(4., np.size(self.xStart))
        co = self.cScatt[:]
        if self.distEffects:
            distRange = np.max(self.dists) - np.min(self.dists)
            distMin = np.min(self.dists)
            distScal = (self.dists - distMin) / distRange
            szRange = self.szMax - self.szMin
            sz = self.szMin + distScal * szRange

            # try something similar with the colors
            rgba = colors.to_rgba(self.cScatt[:])
            co = np.zeros((np.size(self.xStart), 4))
            for iCol in range(3):
                co[:,iCol] = rgba[iCol]

            # now assign the alphas
            alphas = self.alphaMax - \
                (1.0 - distScal) * (self.alphaMax - self.alphaMin)
            co[:,3] = alphas


            print np.min(sz), np.max(sz)

        self.scatt = self.ax.scatter(self.xStart, self.yStart, \
                                         c=co, \
                                         s=sz, \
                                         edgecolor='None')

        if self.doGrid:
            self.ax.grid(which='both', alpha=0.6)

    def labelAxes(self):

        """Labels the axes"""
        
        self.ax.set_xlabel(self.labelX, size=self.fszLabel)
        self.ax.set_ylabel(self.labelY, size=self.fszLabel)

    def initScatterOffsets(self):

        """Initialize the scatter offsets"""

        # This returns an object because that's how the animation
        # object appears to work.

        scattInit = np.transpose(np.vstack((self.xStart, self.yStart)))
        
        self.scatt.set_offsets(scattInit)

        return [self.scatt]

    def getFmtForFFMPEG(self):

        """Produces file format string for ffmpeg to search"""

        if self.nZeros < 0:
            self.getNzeros()
        self.fmtFFMPG = '%s/%s_%%0%id.%s' \
            % (self.dirFrames, self.frameStem, self.nZeros, self.frameTyp)

    def movieFromFrames(self, doRun=True):

        """Makes mp4 from frames"""

#        sCommand = '%s -framerate %i -q:v 2 -i %s %s' % \
#            (self.FFMPEG, self.framerate, self.fmtFFMPG, \
#                 self.filMov)

        sCommand = '%s -framerate %i -q:v 2 -i %s %s' % \
            (self.FFMPEG, self.framerate, self.fmtFFMPG, \
                 self.filMov)

        
        if not doRun:
            return

        if os.access(self.filMov, os.R_OK):
            os.remove(self.filMov)
        sOut = subprocess.check_output(sCommand, shell=True)

    def getNzeros(self):

        """Sets the number of zeros for the frame filenames"""

        self.nZeros = int(np.log10(self.nSteps)+1)

    def makeIthFilename(self, i=0):

        """Generates the i'th frame filename"""

        if self.nZeros < 0:
            self.getNzeros()

        sIndex = str(i).zfill(self.nZeros)
        
        self.frameName = '%s_%s.%s' % (self.frameStem, sIndex,self.frameTyp)
        self.framePath = '%s/%s' % (self.dirFrames, self.frameName)
        self.lFrames.append(self.framePath)

    def updateScatterPos(self, i=0, fig=None, scatt=None):

        """Updates the scatter plot with offsets"""

        newX = np.float(i)*self.tStep*self.vX + self.xStart
        newY = np.float(i)*self.tStep*self.vY + self.yStart

        # we want an Nx2 array...
        newPosns = np.transpose(np.vstack(( newX, newY )) )

        scatt.set_offsets(newPosns)

        # updates the title
        if self.labelTimes:
            self.ax.set_title('%i yr' % (self.times[i]), loc='right', \
                                  color=self.cScatt, \
                                  fontsize=self.fszLabel)


        return [scatt]

    def writeFrame(self):

        """Writes the current scatterplot to disk"""

        if not os.access(self.dirFrames, os.R_OK):
            os.makedirs(self.dirFrames)

        self.fig.savefig(self.framePath, quality=100)

    def wipeFrames(self):

        """Cleans the intermediate frames to save disk space"""

        if not self.cleanFrames:
            return

        if len(self.lFrames) < 1:
            return

        for sFrame in self.lFrames:
            if os.access(sFrame, os.W_OK):
                os.remove(sFrame)

class SliderPlot(object):

    """Plots timeseries with slider for visualization"""

    def __init__(self, nGen=50, figNum=1, \
                     labelX='MJD', \
                     labelY='Scale factor', \
                     figsz=(8,3), \
                     nFrames = 100, \
                     tMin=-1., tMax=-1., \
                     tBuf = 50., \
                     xMin = 1e6, \
                     xMax = -1e6, \
                     sliderLW = 4, \
                     sliderColor='k',\
                     sliderAlpha=0.25, \
                     sliderZorder=10, \
                     dirFrames='./tmpFramesSlider', \
                     frameTyp='jpg', \
                     hlSym='o', \
                     hlSz=9, \
                     hlColo='g', \
                     scattColo='b', \
                     scattSz=4, \
                     scattAlpha = 0.6, \
                     clumpInterval = 30., \
                     clumpAlpha=0.25, \
                     clumpColor='darkslateblue', \
                     clumpLW=1):

        # time, vertical data
        self.t = np.array([])
        self.x = np.array([])

        # slider positions
        self.nFrames = nFrames
        self.tSlider = np.array([])
        self.iData = -1

        # slider plot parameters
        self.sliderAlpha = sliderAlpha
        self.sliderColor = sliderColor
        self.sliderLW = sliderLW
        self.sliderZorder = sliderZorder

        # scatter plot parameters
        self.scattColo = scattColo
        self.scattSz = scattSz
        self.scattAlpha = scattAlpha

        # for plot limits
        self.tMin = np.copy(tMin)
        self.tMax = np.copy(tMax)
        self.tBuf = np.copy(tBuf)

        self.xMin = np.copy(xMin)
        self.xMax = np.copy(xMax)

        # for joining the dots between clumps
        self.clumpInterval = clumpInterval
        self.clumpAlpha = clumpAlpha
        self.clumpColor = clumpColor
        self.clumpLW = clumpLW

        # figure
        self.fig = None
        self.figNum = figNum
        self.ax = None
        self.figSz = figsz

        # some more frame parameters
        self.dirFrames = dirFrames[:]
        self.frameStem = 'slider'
        self.frameTyp = frameTyp[:]
        self.nZeros = -1
        self.lFrames = []

        # slider vertical values
        self.ySlider = np.array([])

        # highlight plot information
        self.hlSym = hlSym
        self.hlSz = hlSz
        self.hlColo = hlColo

        # plot objects to update
        self.scatt = None
        self.slider = None
        self.hilite = None

        self.labelX = labelX[:]
        self.labelY = labelY[:]

        # when generating data
        self.nGen = nGen

    def makePoints(self):

        """For testing: makes fake data"""

        tMin = 50000.
        tMax = 51000.

        self.t = np.random.uniform(size=self.nGen) * \
            (tMax - tMin) + tMin
        self.t = self.t[np.argsort(self.t)]
        
        self.x = np.random.normal(size=self.nGen)

    def makeFigure(self):

        """Makes the figure"""

        self.setPlotLims()
        xLimits = np.array([self.tMin, self.tMax])

        # (copy from the above object)
        yLimits = np.array([np.min(self.x)*0.95, np.max(self.x)*1.05])

        self.fig = plt.figure(self.figNum)
        self.fig.set_size_inches(self.figSz, forward=True)
        self.fig.clf()
        self.ax = self.fig.add_subplot(111, \
                                           xlim=xLimits,\
                                           ylim=yLimits)
        
        # ensure sufficient room at the bottom to fit the label
        self.fig.subplots_adjust(bottom=0.2, left=0.2)

        # label the axes right away
        self.ax.set_xlabel(self.labelX)
        self.ax.set_ylabel(self.labelY)

    def setPlotLims(self):

        """Sets the plot limits"""

        if self.tMin < 0:
            self.tMin = np.min(self.t) - self.tBuf
        if self.tMax < 0:
            self.tMax = np.max(self.t) + self.tBuf

    def getPlotLims(self):

        """Gets the (vertical) plot limits having plotted once"""

        yLims = np.copy(self.ax.get_ylim())
        self.xMin = np.copy(yLims[0])
        self.xMax = np.copy(yLims[1])

    def adjustVerticalLimits(self, factor=0.2):

        """Can adjust vertical plot limits"""

        yMed =  0.5*(self.xMax + self.xMin)
        yDist = 0.5*(self.xMax - self.xMin)
        
        self.xMax = yMed + yDist * (1.0 + factor)
        self.xMin = yMed - yDist * (1.0 + factor)
        


    def makeSliderTimes(self):

        """Generates fine-grid of slider times"""

        self.tSlider = np.linspace(self.tMin, self.tMax, \
                                       self.nFrames, endpoint=True)

    def scatterData(self):

        """Adds the scatterplot for data to the frame"""

        self.scatt = self.ax.scatter(self.t, self.x, zorder=5, \
                                         c=self.scattColo, \
                                         s=self.scattSz, \
                                         alpha=self.scattAlpha)

    def drawClumpLines(self):

        """Draws lines connecting clumps"""

        # partition the times into clumps
        dt = self.t - np.roll(self.t, 1)
        gFirstClump = np.where(dt > self.clumpInterval)[0]
        gLastClump = np.roll(gFirstClump, 1)
        gLastClump[0] = 0 

        for iClump in range(np.size(gLastClump)):
            iLo = gLastClump[iClump]
            iHi = gFirstClump[iClump]

            dum = self.ax.plot(self.t[iLo:iHi], \
                                   self.x[iLo:iHi], \
                                   color=self.clumpColor, \
                                   alpha=self.clumpAlpha, \
                                   lw=self.clumpLW)

    def addSliderToPlot(self, iSlider=0):

        """Add the slider to the plot"""

        if np.size(self.ySlider) <> 2:
            self.getPlotLims()
            self.adjustVerticalLimits()
            self.ax.set_ylim(self.xMin, self.xMax)
            self.ySlider = np.array([self.xMin, self.xMax])

        # times for the slider (as a vector so we can plot)
        tSlider = self.tSlider[iSlider]*np.array([1., 1.])

        # select the datapoint below the slider
        bPos = self.t <= tSlider[0]
        if np.sum(bPos) < 1:
            self.iData = -1
        else:
            gPos = np.where(bPos)[0]
            iMostRecent = np.argmin(tSlider[0] - self.t[gPos])
            self.iData = gPos[iMostRecent]
        
        # now add the slider
        if iSlider < 1:
            self.slider,  = self.ax.plot(tSlider, self.ySlider, \
                                             alpha=self.sliderAlpha, \
                                             color=self.sliderColor, \
                                             lw=self.sliderLW, \
                                             zorder = self.sliderZorder)
        else:
            self.slider.set_data(tSlider, self.ySlider)

    def hilightCurrentScatter(self):

        """Highlights the current scatter object"""

        # don't plot anything if we haven't got to a datapoint yet
        if self.iData < 0:
            return

        tThis = self.t[self.iData]
        xThis = self.x[self.iData]
        if not self.hilite:
            self.hilite,   = self.ax.plot(tThis, xThis, \
                                              c=self.hlColo, \
                                              ms=self.hlSz, \
                                              marker=self.hlSym, \
                                              ls='None', \
                                              zorder=15)            
            return

        else:
            self.hilite.set_data(tThis, xThis)

    # some familiar stuff for filenames
    def getNzeros(self):

        """Sets the number of zeros for the frame filenames"""

        self.nZeros = int(np.log10(self.nFrames)+1)

    def makeIthFilename(self, i=0):

        """Generates the i'th frame filename"""

        if self.nZeros < 0:
            self.getNzeros()

        sIndex = str(i).zfill(self.nZeros)
        
        self.frameName = '%s_%s.%s' % (self.frameStem, sIndex,self.frameTyp)
        self.framePath = '%s/%s' % (self.dirFrames, self.frameName)
        self.lFrames.append(self.framePath)

    def precleanFrames(self):

        """Pre-cleans the frames"""

        if not os.access(self.dirFrames, os.R_OK):
            return

        # don't remove from the current directory
        if len(self.dirFrames) < 4:
            return

        lPre = glob.glob("%s/%s_*.%s" % \
                             (self.dirFrames, self.frameStem, self.frameTyp))

        if len(lPre) < 1:
            return

        for sPre in lPre:
            if os.access(sPre, os.W_OK):
                os.remove(sPre)

    def writePlot(self):

        """Writes the current plot to disk"""

        if not os.access(self.dirFrames, os.R_OK):
            os.makedirs(self.dirFrames)
            
        self.fig.savefig(self.framePath, quality=100)


def TestAnim(nPts=1000, nFrames=200, tStep=0.5):

    PP = Points(nPts, tStep=tStep)
    PP.genStartingPoints()
    PP.genVelocities()
    PP.makeTimes()

    # now try making and animating the figure
    PP.makeFigure()
    PP.plotFirstScatt()
    dum = PP.initScatterOffsets()
    
    # try a scatter that stays
    blah = PP.ax.scatter(PP.xStart, PP.yStart, s=3, c='g')

    anim = animation.FuncAnimation(\
        PP.fig, PP.updateScatterPos, \
            init_func=PP.initScatterOffsets,\
            fargs=(PP.fig, PP.scatt), \
            interval=1000, frames=nFrames, blit=True, repeat=True)

    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, bitrate=1800)

    PP.fig.show()

    anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    #anim.save('test.mp4', writer=writer)


def TestFrameByFrame(nPts=10000, tStep=0.5, nFrames=200, framerate=50, \
                         subwindow=True, plotColor='b'):

    """Test frame by frame"""

    PP = Points(nSteps=nFrames, nObj=nPts, tStep=tStep, framerate=framerate, \
                    cScatt=plotColor)

    if subwindow:
        PP.xLims = [1000., 2000.]
        PP.yLims = [1000., 2000.]

    PP.genStartingPoints()
    PP.genVelocities()
    PP.makeTimes()

    # now try making and animating the figure
    PP.makeFigure()
    PP.plotFirstScatt()
    PP.labelAxes()

    for iFrm in range(PP.nSteps):
        PP.makeIthFilename(iFrm)
        sys.stdout.write("\r %s" % (PP.framePath))
        sys.stdout.flush()

        PP.updateScatterPos(iFrm, PP.fig, PP.scatt)
        PP.writeFrame()
        

    # now make the movie
    PP.getFmtForFFMPEG()
    PP.movieFromFrames()

    PP.wipeFrames()

def TestWithData(tStep=200.0, nFrames=200, framerate=50, \
                     subwindow=True, plotColor='r', \
                     pathData='TEST_bothDists_metalRich.fits', \
                     dirFrames='TMP_MR_frames', distCol='dist_FromLoI', \
                     useTrend=False, \
                     showHull=True, MP=False, \
                     dbgPlot=False):


    """Tests building an animation from actual datapoints"""

    # example:
    #
    # testAnim.TestWithData(nFrames=501, useTrend=True, tStep=1000., showHull=True, MP=False)

    if MP:
        pathData='TEST_bothDists_metalPoor.fits'
        dirFrames='TMP_MP_frames'
        plotColor='b'

    PD = Points(nFrames, tStep=tStep, pathData = pathData, \
                    dirFrames=dirFrames[:], \
                    cScatt=plotColor, \
                    pathTrend='rotnCurves_250pBin_1000trials.pickle')

    PD.fszLabel = 16

    PD.dataFromPath()


    if useTrend:
        PD.loadTrend()
        #if not dbgPlot:
        #    
        PD.replaceVelWithTrend()

    # let's try interpolation
    if dbgPlot:
        lTrend = np.argsort(PD.tTrend['dMod'])
        interp = interp1d(PD.tTrend['dMod'][lTrend], \
                              PD.tTrend['muB'][lTrend], \
                              fill_value='extrapolate',\
                              kind='slinear')

        fig2 = plt.figure(2)
        fig2.clf()
        ax2 = fig2.add_subplot(111)

        dumScatt = ax2.scatter(PD.dists, PD.vY)
        dumTrend = ax2.plot(PD.tTrend['dMod'], PD.tTrend['muB'], 'k-')
    
        dumEval = interp(PD.dists)
        dum3 = ax2.plot(PD.dists, dumEval, 'k.', zorder=25)

        print np.min(dumEval), np.max(dumEval)
    
        ax2.set_ylim(-3e-7, 3e-7)

        return

    PD.makeTimes()

    # now try making and animating the figure
    PD.makeFigure()
    PD.plotFirstScatt()
    PD.labelAxes()

    PD.fig.subplots_adjust(bottom=0.18, left=0.15)

    if showHull:
        PD.setHull()
        dumHull = PD.ax.plot(PD.xHull, PD.yHull, 'k-', lw=3, alpha=0.25, \
                                 zorder=30)

    for iFrm in range(PD.nSteps):
        PD.makeIthFilename(iFrm)
        sys.stdout.write("\r %s" % (PD.framePath))
        sys.stdout.flush()

        PD.updateScatterPos(iFrm, PD.fig, PD.scatt)
        PD.writeFrame()
        

    PD.wipeFrames()


def TestSlider(nFrames=10):

    """Test our slider method"""

    SP = SliderPlot(nFrames=nFrames)
    SP.makePoints()
    
    SP.makeFigure()
    SP.scatterData()
    SP.drawClumpLines()

    SP.makeSliderTimes()
    
    # pre-clean the frames
    SP.precleanFrames()

    for iFrame in range(SP.nFrames):
        SP.addSliderToPlot(iFrame)
        SP.hilightCurrentScatter()
        SP.makeIthFilename(iFrame)

        sys.stdout.write("\r %s" % (SP.framePath))
        sys.stdout.flush()

        SP.writePlot()

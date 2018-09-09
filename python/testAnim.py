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

import os, sys
import subprocess

# 2018-09-09 - the backend save doesn't seem to work. Recurrent error
# "I/O on closed file". So we try a different approach...

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
                     fszLabel=20):

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

        # do grid?
        self.doGrid = doGrid

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
            self.ax.set_title('%.2f yr' % (self.times[i]), loc='right', \
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



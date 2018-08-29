#
# accumPosns.py
#

# 2018-08-27 WIC - having matched positions into a big FITS table,
# accumulate the average positions to produce a reference list for
# matching objects

from astropy.table import Table, Column
from astropy import units as u
import os
import numpy as np

from scipy.optimize import leastsq

import matplotlib.pylab as plt
plt.ion()

class PosnSet(object):

    def __init__(self, pathAccum='TEST_matchedMulti_51-150.fits'):

        self.pathAccum = pathAccum[:]
        self.tBig = Table()

        self.pathTransf = 'TEST_transf.fits'
        self.pathMed = 'TEST_medians.fits'
        self.setTransFile()
        self.setMedList()
        
        # arrays for X, Y, M, Q
        self.aX = np.array([])
        self.aY = np.array([])
        self.aM = np.array([])
        self.aQ = np.array([])

        # table of transformation parameters
        self.tPars = Table()

        # housekeeping: array sizes
        self.nSets = 0
        self.nRows = 0

        # apparent magnitude selection (not that it matters much)
        self.faintM = -11.
        
        # array-column pairs
        self.DCols = {'aX':'X', 'aY':'Y', 'aM':'M', 'aQ':'Q'}

        # choice of reference array (for shifting lists before fitting
        # cross-transformations)
        self.iRef = 0

        # a few auxiliary arrays we'll need later
        self.dXmed = np.array([])
        self.dYmed = np.array([])

        # co-ords shifted onto reference frame
        self.aXshifted = np.array([])
        self.aYshifted = np.array([])

        # separate array for the positions re-centered to boresight
        self.aXrecen = np.array([])
        self.aYrecen = np.array([])

        # objects mapped onto the common frame
        self.aXmapped = np.array([])
        self.aYmapped = np.array([])
        
        # boresight for the corrected-frame
        self.cenX = 2096.
        self.cenY = 2096.

        # Median transformed quantities
        self.tMed = Table()
        
    def setTransFile(self):

        """Sets the filename for the transformation parameters"""

        stem = os.path.split(self.pathAccum)[-1]
        self.pathTransf = 'transf_%s' % (stem)

    def setMedList(self):

        """Sets the filename for the median list of positions"""

        stem = os.path.split(self.pathAccum)[-1]
        self.pathMed = 'mednPosn_%s' % (stem)

    def clearBigTable(self):

        """Re-initializes the big table"""

        self.tBig = Table()

    def initArrays(self):

        """Initializes X, Y, M, Q arrays"""

        for sAttr in ['aX', 'aY', 'aM', 'aQ']:
            setattr(self, sAttr, np.zeros((self.nSets, self.nRows)) )
        
    def loadAccum(self):

        """Loads the accumulated positions"""

        self.tBig = Table.read(self.pathAccum, 'fits')


    def extractArrays(self, freemem=True):

        """Extracts separate arrays for X, Y, M, Q"""

        if len(self.tBig) < 1:
            return
        
        # how many datasets went into this?
        self.nSets = int(self.tBig.colnames[-1].split('_')[-1])
        self.nRows = len(self.tBig)

        # initialize the arrays
        self.initArrays()

        # now go through the sets and slot in the values
        for iSet in range(self.nSets):
            for sAttr in self.DCols.keys():
                sCol = '%s_%i' % (self.DCols[sAttr], iSet+1)

                thisArr = getattr(self, sAttr)
                thisArr[iSet] = np.asarray(self.tBig[sCol])

        if freemem:
            self.clearBigTable()

    def estOffsets(self):
        
        """Estimates differences from the median in each case"""

        self.medX = np.median(self.aX, axis=0)
        self.medY = np.median(self.aY, axis=0)
        self.medM = np.median(self.aM, axis=0)
        self.medQ = np.median(self.aQ, axis=0)

        # replicate these all to 2D arrays
        aMedX = self.aX * 0. + self.medX
        aMedY = self.aY * 0. + self.medY
        aMedM = self.aM * 0. + self.medM
        aMedQ = self.aQ * 0. + self.medQ

        self.dX = self.aX - aMedX
        self.dY = self.aY - aMedY
        self.dM = self.aM - aMedM
        self.dQ = self.aQ - aMedQ

        self.dXmed = np.median(self.dX, axis=1)
        self.dYmed = np.median(self.dY, axis=1)
        self.dMmed = np.median(self.dM, axis=1)
        self.dQmed = np.median(self.dQ, axis=1)

    def pickRefFrame(self):

        """Pick the reference frame"""

        if np.size(self.dXmed) < 1 or np.size(self.dYmed) < 1:
            return
        
        drMed = self.dXmed**2 + self.dYmed**2
        self.iRef = np.argmin(drMed)

    def shiftOntoMedian(self):

        """Shifts each list onto the median position"""

        shiftX = self.aX * 0.
        shiftY = np.copy(shiftX)
        for iSet in range(np.size(self.dXmed)):
            shiftX[iSet, :] = self.dXmed[iSet]
            shiftY[iSet, :] = self.dYmed[iSet]

        self.aXshifted = self.aX - shiftX
        self.aYshifted = self.aY - shiftY
        
        # do the boresight now
        self.shiftToBoresight()
        
    def shiftToBoresight(self):

        """Shifts the positions onto the boresight"""

        self.aXrecen = self.aXshifted - self.cenX
        self.aYrecen = self.aYshifted - self.cenY

    def medianTransformed(self):

        """Finds the median of the transformed quantities. Does not do any
selection by measurement."""

        # initialise the table
        self.tMed = Table()

        self.tMed['X'] = np.median(self.aXmapped, axis=0)
        self.tMed['Y'] = np.median(self.aYmapped, axis=0)
        self.tMed['eX'] = np.std(self.aXmapped, axis=0)
        self.tMed['eY'] = np.std(self.aYmapped, axis=0)
        
        # now for the other quantities:
        for sCol in ['M','Q']:
            thisAttr = 'a%s' % (sCol)
            aThis = getattr(self, thisAttr)
            eCol = 'e%s' % (sCol)
            self.tMed[sCol] = np.median(aThis, axis=0)
            self.tMed[eCol] = np.std(aThis, axis=0)
            
        
    def initTransformTable(self):

        """Initialises the table of transformation parameters"""

        zInt = np.zeros(self.nSets, 'int')
        zFlt = np.zeros(self.nSets, 'float')

        # Column names. Note that the floating-point array names
        # should match the attributes in the TransfPair object.
        lInts = ['set']
        lFlts = ['dX', 'dY', 'sX', 'sY', 'rot', 'skew', 'a', 'b', 'c', 'd', 'e', 'f', 'stdX', 'stdY']

        # let's add units too.
        DUnits = {'dX':u.pix, 'dY':u.pix, \
                  'sX':u.Unit(' '), 'sY':u.Unit(' '), \
                  'rot':u.radian, 'skew':u.radian, \
                  'a':u.pix, 'd':u.pix, \
                  'b':u.Unit(' '), \
                  'c':u.Unit(' '), \
                  'e':u.Unit(' '), \
                  'f':u.Unit(' '), \
                  'stdX':u.pix, 'stdY':u.pix}
        
        self.tPars = Table()
        for iInt in range(len(lInts)):
            self.tPars[lInts[iInt]] = np.copy(zInt)

            self.tPars[lInts[iInt]].unit=u.Unit(' ')
            
        for iFlt in range(len(lFlts)):
            sCol = lFlts[iFlt]
            self.tPars[sCol] = np.copy(zFlt)

            if sCol in DUnits.keys():
                self.tPars[sCol].unit = DUnits[sCol]
            
    def fitAllTransforms(self):

        """Wrapper - goes through all the translated sets and fits the
transformations."""

        if np.size(self.aXmapped) < 2:
            self.aXmapped = self.aX * 0.
            self.aYmapped = self.aY * 0.
        
        for iSet in range(np.shape(self.aXrecen)[0]):
            self.fitTransform(iSet)
        
    def fitTransform(self, iSet=0):

        """Fits the transformation for the iSet'th dataset onto the median
reference frame. Also maps onto the common frame for convenient
accumulation of median positions.

        """

        # let's try selecting objects
        bUse = self.aM[iSet] < self.faintM

        xRef = np.median(self.aXrecen, axis=0)
        yRef = np.median(self.aYrecen, axis=0)

        xInp = self.aXrecen[iSet]
        yInp = self.aYrecen[iSet]
        
        TP = TransfPair(xRef[bUse], yRef[bUse], \
                        xInp[bUse], yInp[bUse], True)

        # Initialise the params table if it's not already set
        if len(self.tPars.colnames) < 1:
            self.initTransformTable()

        if not 'faintM' in self.tPars.meta.keys():
            self.tPars.meta['faintM'] = np.float(self.faintM)
            
        if 'set' in self.tPars.colnames:
            self.tPars['set'][iSet] = np.copy(iSet)
            
        # now drop the attributes into the table.
        for sCol in self.tPars.colnames:
            if hasattr(TP, sCol):
                self.tPars[sCol][iSet] = np.copy(getattr(TP, sCol))

        # now apply the transformation to the mapped array
        parsX = np.array([TP.a, TP.b, TP.c])
        parsY = np.array([TP.d, TP.e, TP.f])
        
        self.aXmapped[iSet] = fLinear(parsX, xInp, yInp)
        self.aYmapped[iSet] = fLinear(parsY, xInp, yInp)

        # just for interest, show the histogram of residuals
        if iSet < 2:
            fig2 = plt.figure(2)
            fig2.clf()
            ax2 = fig2.add_subplot(111)
            diffsX = err_leastsq(parsX, xInp, yInp, xRef)
            diffsY = err_leastsq(parsY, xInp, yInp, yRef)
            dum = ax2.scatter(diffsX[bUse][TP.bUse], diffsY[bUse][TP.bUse], \
                              alpha=0.15, \
                              s=4)

            ax2.set_xlabel(r'$\delta X$ (pix)')
            ax2.set_ylabel(r'$\delta Y$ (pix)')

            #print("fitTransform INFO: %.2e, %.2e" % \
            #      (np.std(diffsX[bUse][TP.bUse]), \
            #       np.std(diffsY[bUse][TP.bUse]) ) )

    def writeTransfTable(self):

        """Writes the transformation table to disk"""

        if len(self.tPars) < 1:
            return
        
        self.tPars.write(self.pathTransf, overwrite=True)

    def writeMedList(self):

        """Writes the median position list to disk"""

        if len(self.tMed) < 1:
            return
        
        self.tMed.write(self.pathMed, overwrite=True)
        
    def showPosns(self):

        """Debug routine - shows the positions"""

        #lPosn = np.arange(np.size(self.dXmed))

        # let's pick random integers for our set
        lPosn = np.random.random_integers(0, self.nRows-1, 500)
        lTimes = np.arange(np.size(self.dXmed))
        
        fig1 = plt.figure(1)
        fig1.clf()

        ax1 = fig1.add_subplot(111)
        dumAll = ax1.scatter(self.aXrecen[0], self.aYrecen[0], \
                             color='0.7', s=2)
        for iShow in range(500):
            xThis = self.aXrecen[:,lPosn[iShow]]
            yThis = self.aYrecen[:,lPosn[iShow]]
            dum = ax1.scatter(xThis, yThis, c=lTimes)

        ax1.set_xlabel(r'$\Delta X$ (pixels)')
        ax1.set_ylabel(r'$\Delta Y$ (pixels)')

class TransfPair(object):

    def __init__(self, xRef=np.array([]), yRef = np.array([]), \
                 xIn = np.array([]), yIn=np.array([]), \
                 nSigm=4., nClip=3, \
                 runOnInit=True):

        """Methods to find the transformation between two position-lists via
direct application of the method of least squares

        """
        
        self.xRef = np.copy(xRef)
        self.yRef = np.copy(yRef)
        self.xIn = np.copy(xIn)
        self.yIn = np.copy(yIn)

        # objects to use
        self.bUse = np.isfinite(self.xIn)

        # clipping factors
        self.nSigm = np.float(nSigm)
        self.nClip = np.int(nClip)
        
        
        # parameters in two forms
        self.a = 0.
        self.b = 1.
        self.c = 0.
        self.d = 0.
        self.e = 0.
        self.f = 1.

        self.dX = 0.
        self.dY = 0.
        self.sX = 1.
        self.sY = 1.
        self.rot = 0.
        self.skew = 0.

        # some statistics of interest
        self.stdX = 0.
        self.stdY = 0.
        
        if runOnInit:
            self.solveLeastsq()

    def solveLeastsq(self):

        """Solves the parameters by calling an optimizer"""

        guessX = np.array([0., 1., 0.])
        guessY = np.array([0., 0., 1.])

        for iFit in range(self.nClip):
            parsX, statusX = leastsq(err_leastsq, guessX, \
                                     args=(self.xIn[self.bUse], \
                                           self.yIn[self.bUse], \
                                           self.xRef[self.bUse]))
            parsY, statusY  = leastsq(err_leastsq, guessY, \
                                      args=(self.xIn[self.bUse], \
                                            self.yIn[self.bUse], \
                                            self.yRef[self.bUse]))

            deltsX = err_leastsq(parsX, self.xIn, \
                                 self.yIn, \
                                 self.xRef)
            deltsY = err_leastsq(parsY, self.xIn, \
                                 self.yIn, \
                                 self.yRef)

            deltsR = np.sqrt((deltsX - np.median(deltsX))**2 \
                        + (deltsY - np.median(deltsY))**2)
            dStd = np.std(deltsR)
            bGood = deltsR / dStd < np.float(self.nSigm)

            self.bUse = (self.bUse) & (bGood)

            guessX = np.copy(parsX)
            guessY = np.copy(parsY)

            self.stdX = np.std(deltsX[self.bUse])
            self.stdY = np.std(deltsY[self.bUse])
        
        self.a = parsX[0]
        self.b = parsX[1]
        self.c = parsX[2]
        self.d = parsY[0]
        self.e = parsY[1]
        self.f = parsY[2]

        self.translateParams()
        
    def translateParams(self):

        """Given the four parameters of the linear transformation, recast them
as human-readable parameters"""

        self.sX = np.sqrt(self.b**2 + self.e**2)
        self.sY = np.sqrt(self.c**2 + self.f**2)

        xRot = np.arctan2(self.e, self.b)
        yRot = np.arctan2(self.c, self.f)

        self.rot = 0.5 * (xRot + yRot)
        self.skew = yRot - xRot

        self.dX = self.a
        self.dY = self.d

def fLinear(pars, x=np.array([]), y=np.array([]) ):

    """Returns transformed positions"""

    return x*pars[1] + y*pars[2] + pars[0]
    
def err_leastsq(pars, x=np.array([]), y=np.array([]), z=np.array([]) ):

    """Returns sqrt(chisq) for linear model"""

    return z - fLinear(pars, x, y)
    
def TestLoad(pathAccum='TEST_matchedMulti_51-150.fits', iShow = 10000):

    PS = PosnSet(pathAccum)
    PS.loadAccum()
    PS.extractArrays()
    PS.estOffsets()
    PS.pickRefFrame()
    PS.shiftOntoMedian()
    PS.initTransformTable()
    
    PS.fitAllTransforms()
    PS.writeTransfTable()

    PS.medianTransformed()
    PS.writeMedList()

    # Some syntax to plot the results. Could move into a method.
    fig3 = plt.figure(3)
    fig3.clf()
    ax3 = fig3.add_subplot(211)
    ax4 = fig3.add_subplot(212, sharex=ax3, sharey=ax3)
    
    dum = ax3.scatter(PS.tMed['M'], PS.tMed['eX'], alpha=0.25, s=4)
    ax3.set_xlabel('Instrumental mag')
    ax3.set_ylabel(r'$\sigma(X)$, pix')
    ax3.grid(which='both', visible=True, zorder=1)

    dum4 = ax4.scatter(PS.tMed['M'], PS.tMed['eY'], alpha=0.25, s=4)
    ax4.set_xlabel('Instrumental mag')
    ax4.set_ylabel(r'$\sigma(Y)$, pix')
    ax4.grid(which='both', visible=True, zorder=1)

    for ax in [ax3, ax4]:
        ax.set_yscale('log')
        ax.set_ylim(1e-3,1e-1)
        ax.set_xlim(-13.9, -11.9) 

#
# gatherTranstables.py
#

# 2018-09-04 WIC - gather pointing information from various places

from astropy.table import Table, vstack, Column, join
from astropy.io import fits
import os
import glob
import numpy as np

class Pointings(object):

    def __init__(self, \
                     pathVF='vafactors.txt', \
                     srchTrans='./transf_j*fits', \
                     dirTrans='./transf', \
                     dirHeaders='./', \
                     tailHdr='_flc_xymc_noAst_trim.fits', \
                     pathFused='./pointings.fits'):

        self.pathVF = pathVF[:]
        self.dirTrans = dirTrans[:]
        self.srchTrans = srchTrans[:]
        self.dirHeaders = dirHeaders
        self.pathFused = pathFused[:]

        # tail for header
        self.tailHdr = tailHdr[:]

        self.rowRead = 1

        # list of transformations
        self.LTrans = []

        # vafactors
        self.tVF = Table()
        self.DVF = {}
        
        # fused table
        self.tFused = Table()

    def findTransfs(self):

        """Searches for transformation files"""
        
        self.LTrans = glob.glob('%s/%s' % (self.dirTrans, self.srchTrans))

    def loadVAFACTOR(self):

        """loads the file-vafactor set"""

        if not os.access(self.pathVF, os.R_OK):
            print("loadVAFACTOR WARN - cannot read VAFACTORS: %s" \
                % (self.pathVF))
            return

        self.tVF = Table.read(self.pathVF, format='ascii')
        self.tVF['col1'].name='flt'
        self.tVF['col2'].name='VAFACTOR'
        self.translateVAFACTOR()

    def translateVAFACTOR(self):

        """Translates VAFACTOR to dictionary for convenience"""

        self.DVF = {}
        for iRow in range(len(self.tVF)):
            thisRow = self.tVF[iRow]
            thisStem = thisRow['flt'].split('_flt')[0]
            thisVF = thisRow['VAFACTOR']

            self.DVF[thisStem] = thisVF

    def accumTransfs(self):

        """Accumulates all the transformation tables"""

        for iTransf in range(len(self.LTrans)):
            self.accumTransf(iTransf)

    def accumTransf(self, iEntry=0):

        """Loads the iEntry'th transformation file"""

        sFil = self.LTrans[iEntry]
        sStem = os.path.split(sFil)[-1].split('transf_')[-1].split('_')[0]
        if not os.access(sFil, os.R_OK):
            return

        # header file
        sHeader = '%s%s' % (sStem, self.tailHdr)
        thisMJD = 0.
        thisV3 = 0.
        try:
            hdr = fits.getheader(sHeader, 1)
            thisMJD = hdr['MJD'] 
            thisV3  = hdr['PA_V3']
        except:
            thisMJD = 0.

        tThis = Table.read(sFil)
        # existing rows
        LNames = tThis.colnames

        # vafactor
        thisVF = 0.
        if sStem in self.DVF.keys():
            thisVF = self.DVF[sStem]
            
        tThis.add_column(Column([thisVF, thisVF], name='VAFACTOR'))

        # column
        tThis.add_column(Column([sStem, sStem], name='stem'))

        tThis.add_column(Column([thisMJD, thisMJD], name='MJD'))
        tThis.add_column(Column([thisV3, thisV3], name='PA_V3'))

        # reorder the table rows
        LReorder = ['MJD', 'stem', 'VAFACTOR', 'PA_V3'] + LNames
        tThis = tThis[LReorder]

        # we only want one row. Add a column which gives the filename
        # stem.
        tThis = tThis[self.rowRead]
        
        # ensure the master table exists with the right columns
        if len(self.tFused) < 1:
            self.tFused = tThis
        else:
            self.tFused = vstack([self.tFused, tThis])
        
    def writeFused(self):

        """Writes fused table to disk"""

        if len(self.tFused) < 1:
            return

        if len(self.pathFused) < 3:
            return

        self.tFused.write(self.pathFused, overwrite=True)

class Transf(object):

    def __init__(self, pathTransf=''):
        
        self.pathTransf = pathTransf[:]
        self.tTransf = Table()

    def loadThisTransf(self):

        """Loads the transformation"""

        if len(self.pathTransf) < 3:
            return

        if not os.access(self.pathTransf, os.R_OK):
            return

        self.tTransf = Table.read(self.pathTransf)

def TestGather():

    # example use:
    
    # cd  /media/datadrive0/Data/HST/9750/PROC/f606w/gathered/xymc

    # !cp -p /home/wiclarks/Projects/properVisualization/tmpIntermed/vafactors.txt .

    # gatherTranstables.TestGather()

    # takes about 5-10s for 250 images
    
    POINT = Pointings()
    POINT.findTransfs()
    POINT.loadVAFACTOR()
    POINT.accumTransfs()
    POINT.writeFused()

def joinPointings(set1='pointings_f814w.fits', \
                  set2='pointings_f606w.fits', \
                  setOut='pointingsAll.fits'):
    
    """Utility - joins the two pointings"""

    # cd /media/datadrive0/Data/HST/9750/PROC/collected
    
    for tThis in [set1, set2]:
        if not os.access(tThis, os.R_OK):
            print("joinPointings WARN - cannot read input table: %s" \
                  % (tThis))
            return

    tOne = Table.read(set1)
    tTwo = Table.read(set2)

    # add a column for filter
    tOne['FILTER'] = tOne['set'] * 0 + 1
    tTwo['FILTER'] = tTwo['set'] * 0 + 2
    
    tBoth = join(tOne, tTwo, join_type='outer')

    tBoth.meta['file1'] = set1[:]
    tBoth.meta['file2'] = set2[:]
    
    # it looks like the default behavior is to sort on the first
    # column. Nice!
    
    tBoth.write(setOut, overwrite=True)
    
def standardisePointings(pathIn='pointingsAll.fits', \
                         pathOut='pointingsStan.fits'):

    """Standardizes the pointings to the median vafactor"""

    if not os.access(pathIn, os.R_OK):
        print("standardisePointings WARN - cannot read path %s" \
              % (pathIn))
        return

    tRaw = Table.read(pathIn)

    medVF = np.median(tRaw['VAFACTOR'])

    # create new columns
    tRaw['ssRaw'] = 0.5*(tRaw['sX']+tRaw['sY'])
    tRaw['ssScaled'] = tRaw['ssRaw']*1.0
    
    lFilters = np.unique(np.sort(tRaw['FILTER']))
    for iFilt in lFilters:
        bThis = tRaw['FILTER'] == iFilt

        scaleThis = np.median(tRaw['VAFACTOR'][bThis] / \
                              tRaw['ssRaw'][bThis])

        tRaw.meta['mult%i' % (iFilt)] = np.float(scaleThis)

        tRaw['ssScaled'][bThis] *= scaleThis
        

    tRaw.write(pathOut, overwrite=True)

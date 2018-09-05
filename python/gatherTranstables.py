#
# gatherTranstables.py
#

# 2018-09-04 WIC - gather pointing information from various places

from astropy.table import Table, vstack, Column
from astropy.io import fits
import os
import glob
import numpy

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

    POINT = Pointings()
    POINT.findTransfs()
    POINT.loadVAFACTOR()
    POINT.accumTransfs()
    POINT.writeFused()


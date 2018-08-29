#
# getPosns.py
#

# Python wrapper for Jay Anderson's img2xym fortran routine. Also does
# some filename carpentry.

import shutil
import os, sys, time
import subprocess
from astropy.io import fits
import glob

import numpy as np

from astropy.table import Table
from astropy.time import Time

class PositionSet(object):

    def __init__(self, pathImg='', doTest=True, Verbose=True):

        self.Verbose = Verbose # control variable
        
        # initialize the paths
        self.initPaths()
        self.pathImg=pathImg[:]

        # command and arguments for image finding
        self.i2xCOMMAND = '/home/wiclarks/progs/AKWFC/img2xym_WFC_09x10'
        self.i2xHMIN = 5
        self.i2xFMIN = 50
        self.i2xPMAX = 999999
        self.filEPSF = 'PSFEFF.F750L.fits' # deliberately wrong to check
        self.i2xPERT = 'PERT'
        self.i2xSUBT = 'SUBT'
        self.i2xQSEL = ''

        self.magMin = -13.7
        self.magMax =  -12.0
        self.qMin = 0.
        self.qMax = 0.5

        # output directory
        self.dirOut = 'OUTPUT'

        # position-file (useful when merging output with headers)
        self.filPos = 'NONE'
        self.filTmp = 'tmp.xym'

        # for fusing positions and header
        self.tPosn = Table()
        self.filFused = 'NONE'
        
        # file to hold screen output from the routine
        self.filScreen='stdout.txt'
        
        # filter
        self.sFiltr='CLEAR'
        self.sExp='0s'
        
        # directory holding PSFEFF files
        self.dirEPSF = '/home/wiclarks/progs/AKWFC/PSFLIB/'
        
        # image header
        self.hdr = []
        
        if doTest:
            self.setupTest()

        self.parseImgPath()
        
    def initPaths(self):

        """Initializes the paths"""

        self.dirImg='NONE'
        self.filImg='NONE'
        self.filUse='NONE'
        self.compressTail=''
        self.cmdUnpack = 'NONE'
        self.cmdPack = 'NONE'
        self.filEPSF = 'PSFEFF.F850L.fits'
        self.dirOut = 'OUTPUT'
        
    def setupTest(self):

        """Sets some variables for testing"""

        if self.Verbose:
            print("PositionSet.setupTest INFO: setting test conditions")
        
        #self.pathImg='j8q623yvq_flc.fits.fz'

        self.pathImg='/media/datadrive0/Data/HST/9750/flc/f814w/j8q616krq_flc.fits.fz'
        
    def parseImgPath(self):

        """Parses the image path"""

        dirImg, filImg = os.path.split(self.pathImg)

        if len(dirImg) < 1:
            dirImg = os.getcwd()

        # whatever we do, we want to know the original pieces of the
        # path. So, pass these up to the instance.
        self.dirImg = dirImg[:]
        self.filImg = filImg[:]
        self.filUse = filImg[:]
            
        # is compressed?
        self.compresstail = ''
        if filImg[-1].find('z') > -1:
            #self.filUse, self.compressTail = os.path.splitext(filImg)[-1]
            lSplit = os.path.splitext(filImg)
            self.filUse = lSplit[0]
            self.compressTail = lSplit[-1]

        self.chooseUnpacker()
            
    def chooseUnpacker(self):

        """Chooses the appropriate routine to unpack the image file"""

        DUnpack = {'gz':'gunzip', \
                   'bz':'bunzip2', \
                   'fz':'funpack'}
        DPack = {'gz':'gzip', \
                 'bz':'bzip2', \
                 'fz':'fpack'}

        if len(self.compressTail) < 1:
            return
        
        sKey = self.compressTail[1::]
        self.cmdUnpack = DUnpack[sKey]
        self.cmdPack = DPack[sKey]

    def unpackImg(self):

        """Unpacks the image file"""

        if self.filUse.find(self.filImg) > -1:
            if self.Verbose:
                print("PositionSet.unpackImg WARN - zipped and unzipped files are identical: %s and %s" % (self.filUse, self.filImg))
            return
        
        if self.cmdUnpack.find('NONE') > -1:
            if self.Verbose:
                print("PositionSet.unpackImg WARN - unpack command %s. Not running." % (self.cmdUnpack) )
            return

        # The source file needs to be in the current directory
        # (possibly copied over from another location)
        if not os.access(self.filImg, os.R_OK):
            if self.Verbose:
                print("PositionSet.upackImg WARN - file %s is not in the current working directory)" % (self.filImg))
            return
        
        
        # the unpacking should move harmlessly on if the output file
        # already exists. However, the output string is annoying. So
        # we include the conditional on existence here.
        if not os.access(self.filUse, os.R_OK):
            if self.Verbose:
                print("PositionSet.unpackImg INFO - unpacking compressed image %s" % (self.filImg))
            subprocess.call([self.cmdUnpack, self.filImg])


    def removeUnpacked(self):

        """Removes the unpacked image, leaving the packed image"""

        # do nothing if the packed and unpacked filenames are identical:
        if self.filUse.find(self.filImg) > -1:
            if self.Verbose:
                print "PositionSet.removeUnpacked WARN - packed and unpacked files are identical. Not removing."
            return

        # don't do anything if the unpack command is not set
        if self.cmdPack.find('NONE') > -1:
            if self.Verbose:
                print "PositionSet.removeUnpacked WARN - no repack command set."
            return

        # now if our unpack command removed the original packed file,
        # we repack it here.
        if not os.access(self.filImg, os.R_OK):
            subprocess.call(self.cmdPack, self.filUse)

        # if we get here, and only if BOTH the unpacked and repacked
        # file exist, we remove the unpacked file.
        if not os.access(self.filImg, os.R_OK):
            if self.Verbose:
                print "PositionSet.removeUnpacked WARN - packed file not readable. Not removing unpacked file."

            return

        if os.access(self.filUse, os.R_OK):
            os.remove(self.filUse)

    def removeLocalImg(self, override=False):

        """Removes the local copy of the image"""

        # don't remove the only copy
        if self.dirImgIsWD():
            if not override:
                print("PositionSet.removeLocalImg WARN - source directory is local directory and override=False. Not removing file %s" % (self.filImg))
            return
                
        if os.access(self.filImg, os.R_OK):
            os.remove(self.filImg)
        
    def copyImgHere(self):

        """Copies the image from storage to working directory"""

        dirSrc = self.dirImg
        dirDes = os.getcwd()

        if self.dirImgIsWD():        
            if self.Verbose:
                print("PositionSet.copyImgHere WARN - source and destination directories identical. Not copying.")
            return

        if self.Verbose:
            print("PositionSet.copyImgHere INFO - copying in file %s" % (self.pathImg))
        shutil.copy(self.pathImg, dirDes)

    def getHeader(self):

        """Load the FITS header"""

        if not os.access(self.filUse, os.R_OK):
            if self.Verbose:
                print("PositionSet.getHeader WARN - image filr not readable: %s" % (self.filUse))
            return

        self.hdr=fits.getheader(self.filUse, verify='silentfix')
        
    def parseFilter(self):

        """Use the FITS header to find the filter"""

        if not os.access(self.filUse, os.R_OK):
            if self.Verbose:
                print("PositionSet.parseFilter WARN - image filr not readable: %s" % (self.filUse))
            return

        self.getHeader()

        # currently we only care about f814w or f606w

        self.sFiltr='CLEAR'
        for sKey in ['FILTER1', 'FILTER2']:
            try:
                sThis = self.hdr[sKey]
            except:
                sThis = ''

            if sThis.find('CLEAR') > -1:
                continue

            self.sFiltr = sThis[:]

        # get the exposure time too (might need to rename this method)
        expTime = self.hdr['EXPTIME']
        self.sExp='%.1f' % (expTime)
        self.sExp = self.sExp.replace('.','p')
            

    def chooseEPSF(self):

        """Finds EPSF file given filter"""

        self.pathEPSF = 'NONE'
        srchEPSF = '%s/*%s.fits' % (self.dirEPSF, self.sFiltr)

        LEPSF = glob.glob(srchEPSF)

        if len(LEPSF) > 0:
            self.pathEPSF = LEPSF[0]
            self.filEPSF = os.path.split(self.pathEPSF)[-1]
            
    def bringEPSF(self):

        """Copy EPSF file to local disk"""

        self.filEPSF = os.path.split(self.pathEPSF)[-1]

        if not os.access(self.filEPSF, os.R_OK):
            shutil.copy(self.pathEPSF, os.getcwd())
        
        
    def runFind(self):

        """Runs the finder routine"""

        if not os.access(self.filUse, os.R_OK):
            if self.Verbose:
                print("PositionSet.runFind WARN - image file not readable: %s" % (self.filUse))
            return

        if not os.access(self.filEPSF, os.R_OK):
            if self.Verbose:
                print("PositionSet.runFind WARN - EPSF file not found: %s" \
                      % (self.filEPSF))
        
        lCommand = [self.i2xCOMMAND, str(self.i2xHMIN), \
                    str(self.i2xFMIN), \
                    str(self.i2xPMAX), str(self.filEPSF), \
                    str(self.i2xPERT), \
                    str(self.i2xSUBT), \
                    str(self.i2xQSEL), \
                    self.filUse]

        sCommand='%s %i %i %i %s %s %s %s %s ' \
            % (self.i2xCOMMAND, self.i2xHMIN, \
               self.i2xFMIN, self.i2xPMAX, \
               self.filEPSF, \
               self.i2xPERT, \
               self.i2xSUBT, \
               self.i2xQSEL, \
               self.filUse)

        # subprocess.call(sCommand, shell=True)

        # UPDATE - try capturing the output in a string...
        tZero = time.time()
        sys.stdout.write("Command: %s ... " % (sCommand))
        sys.stdout.flush()
        sOut = subprocess.check_output(sCommand, shell=True)

        sys.stdout.write("Done in %.2f min\n" % ((time.time()-tZero)/60.))
        
        # ... and writing to file
        fOut = open(self.filScreen, 'w')
        fOut.write(sOut)
        fOut.close()
        
    def buildOutputDir(self):

        """Sets up output directory"""

        self.dirOut='./%s/%s' % (self.sFiltr, self.sExp)

        if not os.access(self.dirOut, os.R_OK):
            os.makedirs(self.dirOut)
        
    def organizeOutput(self):

        """Moves the various outputs into results directories"""

        self.buildOutputDir()
        
        # identify the filename stem
        fStem = os.path.splitext(self.filUse)[0]

        # build a list of files to move
        lMove = []
        
        # compress the subtracted fits file
        fSub = '%s_SUB.fits' % (fStem)
        if os.access(fSub, os.R_OK):
            fSubFz='%s.fz' % (fSub)
            if not os.access(fSubFz, os.R_OK):
                subprocess.call('fpack %s' % (fSub), shell=True)

            # remove original if the compressed output was produced.
            if os.access(fSubFz, os.R_OK):
                lMove.append(fSubFz)
                os.remove(fSub)

        # now the position files...
        for sExt in ['.xym', '.xymu', '.xymc']:
            fThis = '%s%s' % (fStem, sExt)
            if os.access(fThis, os.R_OK):
                lMove.append(fThis)

        # now the logfiles
        for sLog in ['LOG.pert.fits', \
                     'LOG.perts.fits', \
                     'LOG.peaksat_CEN.reg', \
                     'LOG.PSFPERT.OUT', \
                     'fort.98', 'fort.92', \
                     self.filScreen]:

            sMoved = '%s_%s' % (fStem, sLog)

            if os.access(sLog, os.R_OK):
                os.rename(sLog, sMoved)                        
            lMove.append(sMoved)

        # now move to destination directory
        for sMove in lMove:
            if not os.access(sMove, os.R_OK):
                continue

            os.rename(sMove, '%s/%s' % (self.dirOut, sMove))
                
    def makeTmpNoAst(self, noComment=True):

        """Copies the position file into a local temporary file without
        asterisks

        """

        if not os.access(self.filPos, os.R_OK):
            if self.Verbose:
                print("getPosn.makeTmpNoAst WARN - position file not readable: %s" % (self.filPos))
            return

        # Loops through the input, moving it to tmp
        with open(self.filPos, 'r') as rObj:
            with open(self.filTmp, 'w') as wObj:
                for line in rObj:
                    if line.find('**') > -1:
                        continue

                    if noComment:
                        if line.find('#') > -1:
                            continue
                    
                    wObj.write(line)

    def loadTmpfil(self):

        """Loads the temporary position file without asterisks"""

        if not os.access(self.filTmp, os.R_OK):
            if self.Verbose:
                print("PosnSet.loadTmpFil WARN - tmpfile not readable: %s" \
                      % (self.filTmp))

            return

        self.tPosn = Table.read(self.filTmp, format='ascii')

    def renameTableColnames(self):

        """Sets the column names for the table. Hardcoded to AKWFC."""

        colNames = self.tPosn.colnames[:]

        if len(colNames) < 4:
            return

        # the initial three are the same no matter what
        namesThree = ['X', 'Y', 'M']
        for iFirst in range(len(namesThree)):
            self.tPosn[colNames[iFirst]].name = namesThree[iFirst][:]

        # the final column will also be 'Q'
        sLast = colNames[-1]
        self.tPosn[sLast].name = 'Q'

        # Those are the columns we need... if there are any more, we
        # could modify this later.

    def trimTable(self):

        """Trims the photometry table by apparent magnitude and q"""

        if len(self.tPosn) < 1:
            return

        bGood = np.repeat(True, len(self.tPosn))
        bMag = np.copy(bGood)
        bQ = np.copy(bGood)

        if 'M' in self.tPosn.colnames:
            mag = self.tPosn['M']
            bMag = (mag >= self.magMin) & (mag < self.magMax)
            self.tPosn.meta['magMin'] = np.float(self.magMin)
            self.tPosn.meta['magMax'] = np.float(self.magMax)

        if 'Q' in self.tPosn.colnames:
            qFac = self.tPosn['Q']
            bQ = (qFac >= self.qMin) & (qFac < self.qMax)
            self.tPosn.meta['qMin'] = np.float(self.qMin)
            self.tPosn.meta['qMax'] = np.float(self.qMax)

        # combine the conditions...
        bGood = (bGood) & (bMag) & (bQ)

        # ... and trim the table
        self.tPosn = self.tPosn[bGood]

    def appendHeaderToTableMeta(self):

        """Appends selected FITS metadata to the table metadata"""

        # I used the following command to make this easier:
        #
        # dfits j8q642dqq_flc_SUB.fits | awk '{print $1}' | sed "s/=/,/g" | sed "s/,//g" | awk '{print "\""$1"\","}'

        lKeys = ["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CTYPE1", \
                 "CTYPE2", "CD1_1", "CD1_2", "CD2_1", "CD2_2", \
                 "ORIENTAT", "PA_APER", "PA_V3", "DATE-OBS", \
                 "TIME-OBS", "EXPTIME", \
                 "RA_TARG", "DEC_TARG", "PROPOSID", "FILTER1", \
                 "FILTER2", "CCDGAIN"]

        # (TARGNAME and ROOTNAME appear to be unparseable, verify.fix
        # does not seem to fix this...)
        
        for sKey in lKeys:
            self.tPosn.meta[sKey] = self.hdr[sKey]

        # let's also calculate the JD if we have both DATE-OBS and
        # TIME-OBS
        tStr = '%sT%s' \
               % (self.hdr['DATE-OBS'], self.hdr['TIME-OBS'])
        tt = Time(tStr)
        self.tPosn.meta['MJD'] = tt.mjd
        self.tPosn.meta['JD'] = tt.jd
        
    def dirImgIsWD(self):

        """Utility - Returns True if the image directory is the working
directory, False otherwise.

        """

        dirSrc = self.dirImg[:]
        dirDes = os.getcwd()

        dirsSame = False
        if dirSrc.find(dirDes) > -1:
            if dirDes.find(dirSrc) > -1:
                dirsSame = True

        return dirsSame

    def wrapSetupForRun(self):

        """Wrapper that gets set up for the run"""

        self.copyImgHere()
        self.unpackImg()
        self.parseFilter()
        self.chooseEPSF()
        self.bringEPSF()

    def wrapOrganizeOutput(self):

        """Wrapper to organize the output"""

        self.organizeOutput()
        self.removeUnpacked()
        self.removeLocalImg()

    def wrapFusePosnHeader(self, doTrim=True):

        """Wrapper: loads a position file and fuses it with the header from
the subtracted fits file. Then writes out to fits position file for
rapid input.

        """

        # ENSURE filenames are correctly set
        filHead = os.path.splitext(self.filPos)[0]
        filStem = os.path.split(filHead)[-1]

        # tail of original file for output filename
        filExten = self.filPos.split('.')[-1]
        
        # filename flag for trimmed
        sTrim = ''
        if doTrim:
            sTrim = '_trim'

        self.filImg = '%s_SUB.fits.fz' % (filStem)
        self.filUse = self.filImg.split('.fz')[0]
        self.filFused = '%s_%s_noAst%s.fits' % (filStem, filExten, sTrim)

        # now remove the asterisks from the image file, putting them
        # into the tmp file. To be sure there is no confusion, remove
        # the tmp file first.
        if os.access(self.filTmp, os.R_OK):
            os.remove(self.filTmp)

        self.makeTmpNoAst()
        self.loadTmpfil()

        # do the renaming and trimming
        self.renameTableColnames()
        if doTrim:
            self.trimTable()

        # Uncompress the SUB.fits, load the header, and remove the
        # uncompressed version IF the compressed version is still on
        # disk
        self.chooseUnpacker()
        self.unpackImg()
        self.getHeader()
        self.removeUnpacked()

        # add the image header to the table. Because we want to steer
        # clear of format problems, we select the keywords we want.
        self.appendHeaderToTableMeta()

        # now write the resulting table to fits. Include the output
        # directory in the path so it all goes to one place
        if not os.access(self.dirOut, os.R_OK):
            os.makedirs(self.dirOut)

        pathFused = '%s/%s' % (self.dirOut, self.filFused)
        
        self.tPosn.write(pathFused, overwrite=True, format='fits')
        
        #print self.tPosn.meta
        
        #print("PositionSet.wrapFusePosnHeader INFO: %s, %s, %s" \
        #      % (self.filPos, self.filUse, self.filFused))
        
def TestOne():

    PS = PositionSet(doTest=True)

    PS.copyImgHere()
    PS.unpackImg()

    PS.parseFilter()
    PS.chooseEPSF()
    PS.bringEPSF()

    PS.runFind()

    PS.organizeOutput()
    
    PS.removeUnpacked()
    PS.removeLocalImg()
    

def TestMany(dirSrc='/media/datadrive0/Data/HST/9750/flc/f814w/', \
             iMin=0, iMax=2):

    """Run two find instances"""

    lImgs = glob.glob('%s/*.fits.fz' % (dirSrc))

    if iMin < 0:
        iStart = 0
    else:
        iStart = iMin
    
    if iMax < 0:
        nRows = len(lImgs)
    else:
        nRows = iMax
        
    tStart = time.time()

    print("TestMany INFO: starting set of (%i-%i) = %i images at %s ..." % \
          (iStart, nRows, nRows - iStart, time.strftime('%X %x %Z')))
    
    for iImg in range(iStart, nRows):
        PS = PositionSet(lImgs[iImg], doTest=False, Verbose=False)
        PS.wrapSetupForRun()
        PS.runFind()
        PS.wrapOrganizeOutput()

    dt = time.time()-tStart
    print("TestMany INFO: did %i images in %.2f minutes" \
          % (nRows, dt/60.))


def TestFuseOne(filPos='j8q642dqq_flc.xymu'):

    """Test fusing a single file with its header"""

    PF = PositionSet(Verbose=True)
    PF.filPos = filPos[:]
    PF.dirOut = '/media/datadrive0/Data/HST/9750/PROC/f814w/gathered'

    # now let's see where we are
    PF.wrapFusePosnHeader()

def FuseMany(srchPos='xymu', filtr='f814w'):

    """Fuses all the position-files with their appropriate headers,
removing bad-readings (***) along the way.

    """

    lPosns = glob.glob('./*.%s' % (srchPos))

    if len(lPosns) < 1:
        return

    tZer = time.time()
    print("FuseMany INFO - about to fuse %i .%s files:" \
          % (len(lPosns), srchPos))

    for iFil in range(len(lPosns)):
        PF = PositionSet(Verbose=False)
        PF.filPos = lPosns[iFil]
        PF.dirOut = '/media/datadrive0/Data/HST/9750/PROC/%s/gathered/%s' \
            % (filtr, srchPos)
        print("FuseMany INFO: %i of %i: %s" \
              % (iFil, len(lPosns), lPosns[iFil]))
        PF.wrapFusePosnHeader()

    tDone = time.time() - tZer
    print("FuseMany INFO: did %i files in %.2f minutes" % \
          (iFil+1, tDone/60.))
        

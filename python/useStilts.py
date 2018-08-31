#
# useStarlink.py
#

#
# 2018-08-27 WIC - calls various Starlink routines. On my system, the
# "starlink" command must be called BEFORE launching IPython, to make
# the various packages available to these methods.
#

import os
import glob
import subprocess
import time

class MatchSet(object):

    def __init__(self, dirTables='', outFil = 'TEST_matched.fits', \
                 srchString = '*_flc_xymc_noAst_trim.fits', \
                 iStart=0, iEnd=-1, \
                 Verbose=True):

        self.Verbose = Verbose
        
        self.dirTables = ''
        self.defaultDir()

        # output directory
        self.dirOut='matches'
        
        # match criteria
        self.ifmt1 = "fits"
        self.ifmt2 = "fits"
        self.find = "best"
        self.join = "1and2"
        self.matcher = "2d"
        self.values1 = "X Y"
        self.values2 = "X Y"
        self.maxDist = 2.
        self.omode = 'out'
        self.outFil = 'TEST_matchedMulti.fits'

        # similar arguments for N-table matching
        self.ifmtN = "fits"
        self.valuesN = "X Y"
        self.joinN = 'match'
        
        # for the moment, we'll use a pair of objects for our
        # matching.
        self.inList1 = 'j8q627j2q_flc_xymc_noAst_trim.fits'
        self.inList2 = 'j8q627mbq_flc_xymc_noAst_trim.fits'

        # search string when tracking down lists of positions
        self.srchString = srchString[:]
        self.lFiles = []
        
        # start and end positions in the list-matching
        self.iStart = iStart
        self.iEnd = iEnd
        
        # match command
        self.cmdMatch = 'NONE'

        # file to capture screen output
        self.filStdout = 'matchSTDOUT.txt'
        
    def defaultDir(self, Clobber=False):

        """Sets default directory"""

        if len(self.dirTables) > 0 and not Clobber:
            return

        self.dirTables='/media/datadrive0/Data/HST/9750/PROC/f814w/testMatch'

    def clobberOutputFile(self):

        """Removes the output file if it exists already"""

        if len(self.outFil) < 3:
            if self.Verbose:
                print("MatchSet.clobberOutput INFO - refusing to remove short-pathname %s" % (self.outFil))
            return

        # nothing to remove if the output file doesn't actually exist
        if not os.access(self.outFil, os.R_OK):
            return

        if os.access(self.outFil, os.W_OK):
            os.remove(self.outFil)

    def findLists(self):

        """Finds lists of positions"""

        self.lFiles = glob.glob(self.srchString)

    def buildOutNamesTwo(self):

        """Constructs the output filenames when being used in a pair"""

        stemInp = os.path.splitext(self.inList2)[0].split('_xy')[0]
        stemRef = 'refpos'
        stemOut = '%s_to_%s' % (stemInp, stemRef)

        self.filStdout = '%s_stdout.txt' % (stemOut)
        self.outFil = '%s.fits' % (stemOut)

        if len(self.dirOut) < 2:
            return

        if not os.access(self.dirOut, os.R_OK):
            os.makedirs(self.dirOut)

        self.filStdout = '%s/%s' % (self.dirOut, self.filStdout)
        self.outFil = '%s/%s' % (self.dirOut, self.outFil)
        
    def buildCommandTwo(self):

        """Builds the match command for a pair of files"""

        # we could do this by looking up named attributes. However,
        # stilts appears to want some of these delineated as strings
        # on the cmdline, so we'll hardcode them term-by-term for the
        # moment.
        
        self.cmdMatch = "stilts tmatch2"
        self.cmdMatch = "%s ifmt1='%s' ifmt2='%s'" % \
                        (self.cmdMatch, self.ifmt1, self.ifmt2)        

        self.cmdMatch = "%s in1='%s' in2='%s'" % \
                        (self.cmdMatch, self.inList1, self.inList2)
        
        self.cmdMatch = "%s find=%s join=%s" % \
                        (self.cmdMatch, self.find, self.join)

        self.cmdMatch = "%s matcher=%s Values1='%s' values2='%s'" % \
                        (self.cmdMatch, self.matcher, self.values1, \
                         self.values2)

        self.cmdMatch = "%s params=%.2f omode=%s out='%s'" % \
                        (self.cmdMatch, self.maxDist, \
                         self.omode, self.outFil)
        
    def buildCommandMulti(self):

        """Builds command for generic crossmatching"""

        if len(self.lFiles) < 2:
            if self.Verbose:
                print("MatchSet.buildCommandMulti WARN - file list has <2 elements.")
            return

        # How many files?
        if self.iEnd < 0 or self.iEnd >= len(self.lFiles):
            iLast = len(self.lFiles)
        else:
            iLast = self.iEnd*1
        
        
        # update the output file with indices
        sOfil = self.outFil.split('.fits')[0]
        self.outFil = '%s_%i-%i.fits' % (sOfil, self.iStart, iLast)

        sStdout = self.filStdout.split('.txt')[0]
        self.filStdout = '%s_%i-%i.txt' % (sStdout, self.iStart, iLast)
        
        self.cmdMatch = "stilts tmatchn"
        #self.cmdMatch = "%s find=%s join=%s" % \
        #                (self.cmdMatch, self.find, self.join)

        self.cmdMatch = "%s matcher=%s" % \
                        (self.cmdMatch, self.matcher)
        
        self.cmdMatch = "%s params=%.2f omode=%s out='%s'" % \
                        (self.cmdMatch, self.maxDist, \
                         self.omode, self.outFil)

        self.cmdMatch = "%s multimode=group" % (self.cmdMatch)
        
        # OK now we have to loop through our list of files, adding one
        # set of format and filename for each. This might take a
        # while...
        inCount = 0
        for iList in range(self.iStart, iLast):
            sList = self.lFiles[iList]
            if not os.access(self.lFiles[iList], os.R_OK):
                continue

            inCount = inCount + 1
            
            sThis = "ifmt%i='%s' values%i='%s' in%i='%s' join%i=%s" % \
                    (inCount, self.ifmtN, \
                     inCount, self.valuesN, \
                     inCount, sList, \
                     inCount, self.joinN)

            self.cmdMatch = "%s %s" % (self.cmdMatch, sThis)

        # ensure the command knows how many counts
        self.cmdMatch = "%s nin=%i" % (self.cmdMatch, inCount)
            
    def runMatcher(self):

        """Performs the matching"""

        if len(self.cmdMatch) < 1 or self.cmdMatch.find('NONE') > -1:
            if self.Verbose:
                print("MatchSet.runMatcher WARN - command not yet set")
            return

        tStarted = time.time()
        sOut = subprocess.check_output(self.cmdMatch, shell=True, \
                                       stderr=subprocess.STDOUT)
        tEnded = time.time()
        
        # write output to disk
        if len(self.filStdout) > 3:
            tElapsed = (tEnded - tStarted)/1.0
            with open(self.filStdout, 'w') as wObj:
                wObj.write('### %s\n' % (time.strftime('%X %x %Z')))
                wObj.write('### Command was: %s\n' % (self.cmdMatch))
                wObj.write('### Time elapsed: %.2e seconds\n' % (tElapsed))

                
                
                if os.access(self.outFil, os.R_OK):
                    wObj.write("### Output file %s\n" \
                      % (self.outFil) )

                # now dump screen output
                wObj.write(sOut)

def TestMatchTwo(reflist='', maxDist=2.):

    """Tests matching two star-lists"""
    
    MS = MatchSet()

    # if a reference list was supplied, use that instead
    if len(reflist) > 1:
        if reflist.find('fits') > -1:
            MS.inList1=reflist[:]
            MS.buildOutNamesTwo()
    MS.maxDist=maxDist
    
    MS.buildCommandTwo()
    MS.clobberOutputFile()
    MS.runMatcher()

    
def TestMatchMulti(iStart=0, iEnd=10):

    MS = MatchSet(iStart=iStart, iEnd=iEnd)    
    MS.findLists()
    MS.buildCommandMulti()
    MS.clobberOutputFile()
    print MS.cmdMatch

    MS.runMatcher()

def TestMatchManyToRef(refList='./try2/mednPosn_TEST_matchedMulti_0-125.fits', \
                       srchString='*_flc_xymc_noAst_trim.fits', \
                       iLo=0, iHi=-1, \
                       maxDist=2., Verbose=True):

    """Matches many files to reference list"""

    # This might be slightly easier using the MatchSet's own glob
    # functionality.  But this way of doing things keeps things
    # separate.

    if not os.access(refList, os.R_OK):
        print("TestMatchManyToRef WARN - cannot read refList: %s" % (refList))
        return

    
    # list of frames to match against the reference
    lFrames = glob.glob(srchString)

    if len(lFrames) < 1:
        print("TestMatchManyToRef WARN - no frames match %s" % (srchString))
        return

    if iHi < 0 or iHi > len(lFrames):
        iHi = len(lFrames)

    for iFrm in range(iLo, iHi):
        inpList = lFrames[iFrm]
        MS = MatchSet()
        MS.Verbose = Verbose
        MS.maxDist = maxDist
        
        MS.inList1 = refList[:]
        MS.inList2 = inpList[:]
        MS.buildOutNamesTwo()

        if Verbose:
            print("TestMatchManyToRef INFO: %4i of %4i: %s" % \
                  (iFrm, iHi-iLo, MS.outFil))
        
        MS.buildCommandTwo()
        MS.clobberOutputFile()
        MS.runMatcher()

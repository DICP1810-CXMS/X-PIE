
#
# read structure files, fit to first one and write out a .edm file
#

from xplor import command
import re
import sys
import os
import protocol
import cdsList
import cdsVector
import vec3
from xplor import select

from xplorPot import XplorPot
from rdcPotTools import *
from pdbTool import *
from atomAction import *
from selectTools import *




sim = xplor.simulation
protocol.initParams('./parallhdg_new.pro')
protocol.initStruct(['./complex-3conf.psf'],erase=False)
protocol.initCoords(['./Random4/Calc_0.pdb'])





#
# names of the structure files are in structures.list
#

structures=open("Random/list").readlines()
    
#
# first structure in list is reference structure
#
#print "reference structure: %s\n" % structures[0]

command("coor disp=comp @%s" % structures[0])

atomPosList = []

for structure in structures:
    sys.stdout.write("reading %s..." % structure)
    command("""
    coor init end
    coor @%s
    """ %structure)

    # fit each structure to reference structure
    command("coor select ((segid ALT0 or segid BLT0) and name ca) fit end")

    atomPosList.append( sim.atomPosArr() )
    print
    pass

from atomProb import AtomProb
from atomSel import AtomSel

#
# make an AtomProb object with structures in atomPosList. The prob. map
# is made up of those atoms listed in the selection string
#
map = AtomProb( AtomSel(""" segid DLT1
"""), atomPosList)


#
# to set grid values by hand
#
#map.gridVals.xmin=-55.63
#map.gridVals.xmax=41.56
#map.gridVals.ymin=-33.61
#map.gridVals.ymax=65.40
#map.gridVals.zmin=-52.44
#map.gridVals.zmax=44.58
#map.gridVals.xdelta=2.0
#map.gridVals.ydelta=2.0
#map.gridVals.zdelta=2.0
#map.gridVals.cushion=100
#map.setAtomRadius(1.0)

#
# other tweaking:
map.setScaleType("flat")
map.setDistType("gaussian")
map.setAtomRadius(1.2)

#print "verbose: " , map.verbose()

map.setVerbose(1)

print "verbose: " , map.verbose()
#
# do the work and write out an edm file
#
map.calc()
map.writeEDM("complex.xplor")
#
#

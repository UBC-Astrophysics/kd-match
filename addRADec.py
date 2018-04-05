from astropy.wcs import WCS
import sys
import numpy as np


w = WCS(sys.argv[1])
data = np.loadtxt(sys.argv[2],unpack=True,skiprows=3)
PXX, PXY = w.all_pix2world(data[1],data[2],0)
data = np.vstack((PXX,PXY,data))

# is there a magnitude column to process
if len(sys.argv)>5:
    magcol=int(sys.argv[3])
    magoffset=float(sys.argv[4])

    
np.savetxt(sys.argv[3], np.transpose(data),fmt=' %14.8f')

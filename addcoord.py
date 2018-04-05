from astropy.wcs import WCS
import sys
import numpy as np

if len(sys.argv)<3:
    print("""Format:

python addcoord.py _FITS_file_ _kd_file_

addcoord.py will add the RA and Dec to the list of stars in the
kd-file using the information in the FITS header.
""")
else:
    w = WCS(sys.argv[1])
    data = np.loadtxt(sys.argv[2],unpack=True)
    RA, Dec = w.all_pix2world(data[1],data[2],0)
    data = np.vstack((data, RA, Dec))
    np.savetxt(sys.argv[2],np.transpose(data),fmt='%7d %8.3f %8.3f %9.3f %12.3f %6.2f %5.2f %8.4f %13.9f %13.9f')



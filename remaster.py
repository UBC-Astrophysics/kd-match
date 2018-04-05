import numpy as np
import sys
from scipy.optimize import leastsq
from scipy.spatial import KDTree

def f(xdata,*params):
    xdata[0]*param[0]+xdata[1]*param[1]+xdata[2]*param[2]+xdata[3]*param[3]+param[4]

def calcxyz(ra,dec):
    ra=np.radians(ra)
    dec=np.radians(dec)
    cosra=np.cos(ra)
    sinra=np.sin(ra)
    sindec=np.sin(dec)
    cosdec=np.sqrt(1-sindec*sindec)
    return sinra*cosdec, cosra*cosdec,sindec

if (len(sys.argv)<5):
    print("""
Format:

python remaster.py coordinate_file master_coordinate_file delta_output_file new_coordinate_file

""")
    exit(-1)
print('hello 1')
data1 = np.loadtxt(sys.argv[1],unpack=True)
ra1=data1[0]
dec1=data1[1]

print('hello 2')
data2 = np.loadtxt(sys.argv[2],unpack=True)
ra2=data2[0]
dec2=data2[1]
x2,y2,z2 = calcxyz(ra2,dec2)
print('hello 3')
tree=KDTree(zip(x2,y2,z2))
print('hello 4')

for dlimit in [8e-7,4e-7,2e-7,1e-7,0.5e-7]:
    x1,y1,z1 = calcxyz(ra1,dec1)
    d,pos=tree.query(zip(x1,y1,z1))
    print('hello 5 %g' % dlimit)

    ra2_match=ra2[pos]
    dec2_match=dec2[pos]
    keep=d<dlimit

    x1f=ra1[keep]
    y1f=dec1[keep]
    x2f=ra2_match[keep]
    y2f=dec2_match[keep]

    print('Number of Matches %d' % len(x1f))
    xmean=np.mean(x1f)
    ymean=np.mean(y1f)
    xstd=np.std(x1f)
    ystd=np.std(y1f)

    xnorm=(x1f-xmean)/xstd
    ynorm=(y1f-ymean)/ystd

    funcQuad=lambda tpl,xnorm,ynorm : tpl[0] + tpl[1]*xnorm + tpl[2]*ynorm + tpl[3]*xnorm*xnorm+tpl[4]*xnorm*ynorm+tpl[5]*ynorm*ynorm
    ErrorFunc=lambda tpl,xnorm,ynorm,Delta: funcQuad(tpl,xnorm,ynorm)-Delta

    tplInitial1=(1e-5,1e-5,1e-5,1e-5,1e-5,1e-5)
    tplx,success =  leastsq(ErrorFunc,tplInitial1[:],args=(xnorm,ynorm,(x2f-x1f)/xstd))
    tply,success =  leastsq(ErrorFunc,tplInitial1[:],args=(xnorm,ynorm,(y2f-y1f)/ystd))
    print("quadratic fit of x1 and x2" ,tplx,success)
    print("quadratic fit of y1 and y2" ,tply,success)

    xnorm=(ra1-xmean)/xstd
    ynorm=(dec1-ymean)/ystd
    ra1=(ra1+funcQuad(tplx,xnorm,ynorm)*xstd)
    dec1=(dec1+funcQuad(tply,xnorm,ynorm)*ystd)
    
deltax=ra2_match-ra1
deltay=dec2_match-dec1
np.savetxt(sys.argv[3],np.transpose([deltax,deltay,d,ra2_match-ra1,dec2_match-dec1,ra1,dec1]))
data = np.vstack((ra1,dec1,data1[2:]))
np.savetxt(sys.argv[4],np.transpose(data))

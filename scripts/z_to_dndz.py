import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt

work_dir = '.'

cat_name = work_dir + '/lensfit_goldshape_2022v1.fits'
nz_hdu1 = work_dir + '/blind_nz_cfis_lensfit_goldshape_2022v1.fits'

hdu = fits.open(cat_name)
hdu.info()
data = hdu1[1].data

blinds = ['A','B','C']


for blind in blinds:
    hdu = fits.open(nz_hdu1)
    z1 = hdu[1].data['Z_%s' %blind]

    (n,bins,_)= plt.hist(z1, bins=200, range=(0,5.0), density=True, histtype='step', weights=None,label='LensFit blind_%s' %blind)

    plt.xlabel('Redshifts')
    plt.ylabel('n(z)')
    print("zmin = ",min(z1))
    print("zmax = ",max(z1))
    plt.legend(fontsize=20)

    plt.show()

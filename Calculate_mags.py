#Last edit 16/02/2022 MF


import numpy as np
import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as mp
import random, os
from scipy.stats import norm
from scipy.interpolate import interp1d
try:
    import pyregion as pyreg
except:
    print('Cannot import pyregion')
    exit()

class PyMask(object):

    """ Class to handle masks in fits files """

    def __init__(self, x, y, regfile, header=None):

        """
        Parse the ds9 region file here. Also save the number 
	of regions and the x y dimensions of the data
        """
        
        if(header):
            self.reg = pyreg.open(regfile).as_imagecoord(header=header)
        else:
            self.reg = pyreg.open(regfile)


        self.filter = self.reg.get_filter()
        self.x = x
        self.y = y
        self.nreg = len(self.reg)
    
    def fillmask(self, ind):
        
        """ Given a region object and an index fill the mask"""

        if (self.reg[ind]).coord_format == 'image':
	
            self.mask = self.filter[ind].mask((self.y, self.x))
            self.maskname = (self.reg[ind]).name
            self.maskcomm = (self.reg[ind]).comment

        else:
            raise TypeError('The region coordinates are not in image units')
	
    def write(self,outname):
        """Write the mask output"""

        from astropy.io import fits
        
        #PrimaryHDU object to encapsulate the data:
	#Mask is a boolean --> convert on the fly
        hdu = fits.PrimaryHDU(self.mask*1)
        hdu.header['REGTYPE'] = self.maskname
        hdu.header['REGNAME'] = self.maskcomm
	
	#write
        hdu.writeto(outname,clobber=True)

def MW_extinction(lam):
  
  anchor_points = np.array((1111,1176,1250,1316,1393,1490,1600,1701,1799,1901,2000,2101,2188,2299,2398,2500,2740,3436,4000,4405,5495,6993,9009,12500,16390,22200))
  anchor_values = np.array((11.53,10.53,9.63,8.99,8.44,8.07,7.81,7.68,7.73,7.98,8.60,9.31,9.65,8.85,7.98,7.27,6.18,4.88,4.38,4.08,3.08,2.30,1.48,0.85,0.50,0.32))
  interpolator = interp1d(anchor_points, anchor_values, kind='cubic')
  
  return interpolator(lam)
  
print('---------------------------------')
print('  Aperture photometry code       ')
print('---------------------------------')

p = argparse.ArgumentParser()
p.add_argument("IMAGE")
p.add_argument("Regfile")
p.add_argument("--StarsReg", default=None)
p.add_argument("--AvoidExtReg", default=None)
#p.add_argument("--ebv", default = 0.04639)
p.add_argument("--ebv", default = 0.0092)
p.add_argument("--outdir",   default = 'pyphotom/')
p.add_argument("--majaxoff", default = 0, type=float)
p.add_argument("--maxbkgoff", default=5, help="Max offset in multiples of the region size")


args = p.parse_args()

hdu1 = fits.open(args.IMAGE)

data = hdu1[0].data


OUTMAG = 0
MONTEERR = 0

#Prepare special variables
if "UVIT_FUV" in args.IMAGE:
  Wave_pivot = 1534 #uvit website
  AG = args.ebv*MW_extinction(Wave_pivot)  #8.24
  band = 'UVIT FUV'
  Exptime = hdu1[0].header['RDCDTIME']
  ZP = 17.771
  ZP_err = 0.01
  OUTMAG = 1
  data /= Exptime
if "GALEX_FUV" in args.IMAGE:
  Wave_pivot = 1524 #galex website
  AG = args.ebv*MW_extinction(Wave_pivot)  #8.24
  band = 'Galex FUV'
  Exptime = hdu1[0].header['EXPTIME']
  ZP = 18.82
  ZP_err = 0.05
  OUTMAG = 1
if "GALEX_NUV" in os.path.basename(args.IMAGE):
  
  Wave_pivot = 2297 #galex website
  AG = args.ebv*MW_extinction(Wave_pivot) #8.2 
  band = 'Galex NUV'
  Exptime = hdu1[0].header['EXPTIME']
  ZP = 20.08
  ZP_err = 0.05
  OUTMAG = 1
# if "NGVS_u" in args.IMAGE:
#   Wave_pivot = 3799
#   AG = args.ebv*MW_extinction(Wave_pivot) #4.682
#   band = "NGVS_u'"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = 30
#   ZP_err = 0.07
#   OUTMAG = 1
# if "NGVS_g" in args.IMAGE:
#   Wave_pivot = 4846.4
#   AG = args.ebv*3.636
#   band = "NGVS_g'"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = 30
#   ZP_err = 0.035
#   OUTMAG = 1
# if "NGVS_i" in args.IMAGE:
#   Wave_pivot = 7543
#   AG = args.ebv*MW_extinction(Wave_pivot) #1.864
#   band = "NGVS_i'"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = 30
#   ZP_err = 0.035
#   OUTMAG = 1
# if "NGVS_z" in args.IMAGE:
#   Wave_pivot = 8823
#   AG = args.ebv*MW_extinction(Wave_pivot) #1.409
#   band = "NGVS_z'"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = 30
#   ZP_err = 0.035
#   OUTMAG = 1
# if "CFHT_r" in args.IMAGE:
#   Wave_pivot = 6369
#   AG = args.ebv*MW_extinction(Wave_pivot) #*2.545
#   band = "CFHT_r'"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = hdu1[0].header['ZP_STACK']+0.04547
#   ZP_err = 0.035
#   OUTMAG = 1
# if "VESTIGE_Ha" in args.IMAGE:
#   Wave_pivot = 6590
#   AG = args.ebv*MW_extinction(Wave_pivot) #*2.545
#   band = "CFHT_Ha (NB)'"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = 30
#   ZP_err = 0.035
#   OUTMAG = 1
# if "NET" in args.IMAGE :
#   Wave_pivot = 6590
#   AG = 0
#   band = "CFHT_Ha NET"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = 1E-18
#   ZP_err = 0.0
#   OUTMAG = -1 #Results in exponential notation
# if "VESTIGE_r" in args.IMAGE:
#   Wave_pivot = 6369
#   AG = args.ebv*MW_extinction(Wave_pivot) #*2.545
#   band = "CFHT_r'"
#   Exptime = hdu1[0].header['EXPTIME']
#   ZP = 30
#   ZP_err = 0.035
#   OUTMAG = 1
if "SDSS_g" in args.IMAGE:
  Wave_pivot = 4702 #From ezgal
  AG = args.ebv*MW_extinction(Wave_pivot) #3.636
  band = "SDSS_g"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['MAGZP']
  ZP_err = 0.00
  OUTMAG = 1
if "SDSS_r" in args.IMAGE:
  Wave_pivot = 6175
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "SDSS_r"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['MAGZP']
  ZP_err = 0.00
  OUTMAG = 1
if "VST_u" in args.IMAGE:
  Wave_pivot = 3500
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "VST_u"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['PHOTZP']
  ZP_err = 0.05
  OUTMAG = 1  
if "VST_g" in args.IMAGE:
  Wave_pivot = 4800
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "VST_g"
  Exptime = 100#Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['PHOTZP']
  ZP_err = 0.05
  OUTMAG = 1  
if "VST_r" in args.IMAGE:
  Wave_pivot = 6250
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "VST_r"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['PHOTZP']
  ZP_err = 0.05
  OUTMAG = 1  
if "VST_i" in args.IMAGE:
  Wave_pivot = 7700
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "VST_i"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['PHOTZP']
  ZP_err = 0.05
  OUTMAG = 1  
if "VST_z" in args.IMAGE:
  Wave_pivot = 9100
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "VST_z"
  Exptime = 90#Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['PHOTZP']
  ZP_err = 0.05
  OUTMAG = 1  
if "SDSS_i" in args.IMAGE:
  Wave_pivot = 7489
  AG = args.ebv*MW_extinction(Wave_pivot) #1.864
  band = "SDSS_i"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = hdu1[0].header['MAGZP']
  ZP_err = 0.00
  OUTMAG = 1

if "CAlto_H" in args.IMAGE:
  Wave_pivot = 16448
  AG = args.ebv*MW_extinction(Wave_pivot) #1.864
  band = "CAlto H (2MASS calib)"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = 24.60 +1.365 #The second term is to go from 2mass vega to AB mags
  ZP_err = 0.05
  OUTMAG = 1

if "VISTA_K" in args.IMAGE:
  Wave_pivot = 21569
  AG = args.ebv*MW_extinction(Wave_pivot) #1.864
  band = "VISTA Ks (2MASS calib)"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = 25.1379 +1.839 #The second term is to go from 2mass vega to AB mags
  ZP_err = 0.05
  OUTMAG = 1

if "IRAC_1" in args.IMAGE:
  AG = args.ebv*0.00
  band = "IRAC chan 1"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = ((0.75/3600)**2)*3.046118*(10**-4)*10**9 #Result in mJy
  ZP_err = 0.00
if "IRAC_2" in args.IMAGE:
  AG = args.ebv*0.00
  band = "IRAC chan 2"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = ((0.75/3600)**2)*3.046118*(10**-4)*10**9 #Result in mJy
  ZP_err = 0.00
if "IRAC_3" in args.IMAGE:
  AG = args.ebv*0.00
  band = "IRAC chan 3"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = ((0.75/3600)**2)*3.046118*(10**-4)*10**9 #Result in mJy
  ZP_err = 0.00
if "IRAC_4" in args.IMAGE:
  AG = args.ebv*0.00
  band = "IRAC chan 4"
  Exptime = hdu1[0].header['EXPTIME']
  ZP = ((0.75/3600)**2)*3.046118*(10**-4)*10**9 #Result in mJy
  ZP_err = 0.00
if "JWST_f090w" in args.IMAGE:
  Wave_pivot = 0.901*1e4
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "JWST_f090w"
  data = hdu1[1].data
  Exptime = hdu1[0].header['TEXPTIME']
  ZP = 28
  ZP_err = 0.05
  OUTMAG = 1
if "JWST_f150w" in args.IMAGE:
  Wave_pivot = 1.501*1e4
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "JWST_f150w"
  Exptime = hdu1[1].header['TEXPTIME']
  ZP = 28
  ZP_err = 0.05
  OUTMAG = 1
if "JWST_f200w" in args.IMAGE:
  Wave_pivot = 1.990*1e4
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "JWST_f200w"
  Exptime = hdu1[1].header['TEXPTIME']
  ZP = 28
  ZP_err = 0.05
  OUTMAG = 1
if "JWST_f277w" in args.IMAGE:
  Wave_pivot = 2.786*1e4
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "JWST_f277w"
  Exptime = hdu1[1].header['TEXPTIME']
  ZP = 26.47
  ZP_err = 0.05
  OUTMAG = 1
if "JWST_f356w" in args.IMAGE:
  Wave_pivot = 3.563*1e4
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "JWST_f356w"
  Exptime = hdu1[1].header['TEXPTIME']
  ZP = 26.47
  ZP_err = 0.05
  OUTMAG = 1
if "JWST_f444w" in args.IMAGE:
  Wave_pivot = 4.421*1e4
  AG = args.ebv*MW_extinction(Wave_pivot) #2.545
  band = "JWST_f444w"
  Exptime = hdu1[1].header['TEXPTIME']
  ZP = 26.47
  ZP_err = 0.05
  OUTMAG = 1
# if "WISE_1" in args.IMAGE:
#   AG = args.ebv*0.00
#   band = "WISE chan 1"
#   Exptime = 1
#   ZP = (1.935E-3) #Result in mJy (DN to mJy conversion)
#   ZP_err = 0.00
# if "WISE_2" in args.IMAGE:
#   AG = args.ebv*0.00
#   band = "WISE chan 2"
#   Exptime = 1
#   ZP = (2.7048E-3) #Result in mJy (DN to mJy conversion)
#   ZP_err = 0.00
# if "WISE_3" in args.IMAGE:
#   AG = args.ebv*0.00
#   band = "WISE chan 3"
#   Exptime = 1
#   ZP = (1.8326E-3) #Result in mJy (DN to mJy conversion)
#   ZP_err = 0.00
#   Gain = 1
# if "WISE_4" in args.IMAGE:
#   AG = args.ebv*0.00
#   band = "WISE chan 4"
#   Exptime = 1
#   ZP = (5.2269E-2) #Result in mJy (DN to mJy conversion)
#   ZP_err = 0.00
# if "Spitzer_24" or "MIPS_24" in args.IMAGE:
#   AG = args.ebv*0.00
#   band = "SPITZER MIPS 24um"
#   Exptime = 1
#   ZP = ((1.5/3600)**2)*3.046118*(10**-4)*10**9 #Result in mJy 
#   ZP_err = 0.00
# if "Spitzer_70" in args.IMAGE:
#   AG = args.ebv*0.00
#   band = "SPITZER MIPS 70um"
#   Exptime = 1
#   ZP = ((4.5/3600)**2)*3.046118*(10**-4)*10**9 #Result in mJy 
#   ZP_err = 0.00
if "PACS_100" in args.IMAGE:
  AG = args.ebv*0.00
  band = "Herschel PACS 100um"
  Exptime = 1
  ZP = 1000 #Result in mJy from Jy
  ZP_err = 0.00
  data= (hdu1[0].data)[0,:,:]
  error = (hdu1[0].data)[1,:,:]
  MONTEERR = 1
if "PACS_160" in args.IMAGE:
  AG = args.ebv*0.00
  band = "Herschel PACS 160um"
  Exptime = 1
  ZP = 1000 #Result in mJy from Jy
  ZP_err = 0.00
  data = (hdu1[0].data)[0,:,:]
  error = (hdu1[0].data)[1,:,:]
  MONTEERR = 1


size = data.shape   

print('Image loaded: {0}'.format(args.IMAGE))
print('Detected band: {0}'.format(band))
print('  Image is {0} x {1} pixels'.format(size[1], size[0]))

#Prepare quantities for aperture photometry

print('  Exposure time: {0}'.format(Exptime))
print('  Photo Zero point: {0}'.format(ZP))
print('  Galactic extinction: {0}'.format(AG))

#Prepare masks
if args.Regfile[-4:] == '.reg':

  Mask_str = PyMask(size[1], size[0], args.Regfile)
  Nregions = Mask_str.nreg+1

  print('Science Region file loaded: {0}'.format(args.Regfile))
  print('  {0} regions found'.format(Nregions))

  #Change position of regions if offset along the major axis is required
  if args.majaxoff !=0:
     centralreg = int(Nregions/2.)
     DeltaX = np.abs(Mask_str.reg[centralreg].coord_list[0] - Mask_str.reg[centralreg-1].coord_list[0])*args.majaxoff
     DeltaY = -1.*np.abs(Mask_str.reg[centralreg].coord_list[1] - Mask_str.reg[centralreg-1].coord_list[1])*args.majaxoff
     #Redefine shape region
     for ii in range(Nregions-1):
       coords = Mask_str.reg[ii].coord_list
       regionstr = 'image;box({},{},{},{},{})'.format(coords[0]+DeltaX,coords[1]+DeltaY,coords[2], coords[3], coords[4])
       Mask_str.reg[ii] = pyreg.pyreg.parse(regionstr)[0]
     Mask_str.filter = Mask_str.reg.get_filter()

  Masks = np.zeros((Nregions, size[0], size[1]))

  for ii in range(int(Nregions)-1):
    Mask_str.fillmask(ii)
    Masks[ii,:,:] = Mask_str.mask*1

  Mask_total = np.sum(Masks, axis=0)
  overlap = (Mask_total>1)
  Mask_total[overlap]=1

  Masks[-1,:,:] = Mask_total


elif args.Regfile[-5:] == '.fits':
  reghdu = fits.open(args.Regfile)
  immask = reghdu[0].data
  
  valbin = np.unique(immask[immask>0])
  Nregions = len(valbin)+1
  
  Masks = np.zeros((Nregions, size[0], size[1]))
  
  for bb in range(Nregions-1):
      tmpmask = np.zeros((size[0], size[1]))
      tmpmask[immask==valbin[bb]] = 1
      
      Masks[bb,...] = tmpmask

  Mask_total = np.zeros((size[0], size[1]))
  Mask_total[immask>0] = 1
  
  Masks[-1,:,:] = Mask_total


else:
  print("Regfile type not understood aborting")
  exit()

#Prepare extra masks if they are required
Mask_bad = np.zeros((size[0], size[1]))
Mask_ext = np.zeros((size[0], size[1]))

if args.StarsReg:
  Mask_str_bad = PyMask(size[1], size[0], args.StarsReg)
  
  print('Star Region file loaded: {0}'.format(args.StarsReg))
  print('  {0} regions found'.format(Mask_str_bad.nreg))
  
  for ii in range(Mask_str_bad.nreg):
     Mask_str_bad.fillmask(ii)
     Mask_bad +=  Mask_str_bad.mask*1

overlap = (Mask_bad>1)
Mask_bad[overlap]=1

if args.AvoidExtReg:
  
  Mask_extra = PyMask(size[1], size[0], args.AvoidExtReg)

  print('Extra Region file loaded: {0}'.format(args.AvoidExtReg))
  print('  {0} regions found'.format(Mask_extra.nreg))
  
  for ii in range(Mask_extra.nreg):
     Mask_extra.fillmask(ii)
     Mask_ext +=  Mask_extra.mask*1

overlap = (Mask_ext>1)
Mask_ext[overlap]=1

TotBadMask = Mask_ext+Mask_bad

hduout = fits.PrimaryHDU(TotBadMask)
hduout.writeto('temp.fits', overwrite=True)

print('---------------------------------')
print('Estimating background in regions')

Background     = np.zeros((Nregions))
Background_sig = np.zeros((Nregions))
Success        = np.ones((Nregions))

#Now for each mask compute background levels and sigma by moving the mask around
for ii in range(int(Nregions)):
   
   #Remove stars from apertures if any
   before = Masks[ii,:,:].sum()
   
   tmp = Masks[ii,:,:] - Mask_bad
   negative = (tmp <0)
   tmp[negative] = 0
   Masks[ii,:,:] = tmp
   
   after = Masks[ii,:,:].sum()
   if after<before:
     print('Removed {0} pixels from region {1}'.format(before-after,ii+1))
   
   #Find indices
   indices_ok = np.where(Masks[ii,:,:] == 1) 
   edges_ok = [np.nanmin(indices_ok[0])-1, np.nanmax(indices_ok[0])+1, np.nanmin(indices_ok[1])-1, np.nanmax(indices_ok[1])+1]
   tmpmask = Masks[ii,edges_ok[0]:edges_ok[1],edges_ok[2]:edges_ok[3]]
   
   if tmpmask.sum() != after:
      print('WARNING Region {}: mask appears to have been trimmed'.format(ii))
       
   cnt = 0
   
   roll_size = int(0.4*np.max(size))
   roll_axis0 = np.random.uniform(-1*roll_size,roll_size,500000).astype(int)
   roll_axis1 = np.random.uniform(-1*roll_size,roll_size,500000).astype(int)
   
   #Apply cut above minimum shift
   #minroll = 200
   #roll_axis0 = roll_axis0[np.abs(roll_axis0)>minroll]
   #roll_axis1 = roll_axis1[np.abs(roll_axis1)>minroll]
   
   #if len(roll_axis0) < 100000 or len(roll_axis1) < 100000:
   #  print 'WARNING WARNING min roll size'
   
   tmp_back = np.zeros(5000)
   tmpbkgmask = np.copy(tmpmask)
   Nbkgpix_orig  = np.sum(tmpbkgmask*1.)
   Nbkgpix_final = np.sum(tmpbkgmask*1.)
   
   FAIL=0
   
   for jj in range(500000):
          
     if cnt>=5000:
       break
     
     if jj>5000 and cnt < 5 and FAIL==0:
        
        print('FAIL')
        FAIL = 1
        tmp_back = np.zeros(5000)
        cnt=0
        
        tmpbkgmask *= 0
        Nbkgsize = int((np.sqrt(Nbkgpix_orig))/2)
        
        tmpbkgmask[0:Nbkgsize,0:Nbkgsize] = 1

        Nbkgpix_final = np.sum(tmpbkgmask*1.)
        
     if (edges_ok[0]+roll_axis0[jj]<5) or (edges_ok[1]+roll_axis0[jj]>=size[0]-5) or (edges_ok[2]+roll_axis1[jj]<5) or (edges_ok[3]+roll_axis1[jj]>=size[1]-5):
       continue
	     
     tmpdata = data[edges_ok[0]+roll_axis0[jj]:edges_ok[1]+roll_axis0[jj],edges_ok[2]+roll_axis1[jj]:edges_ok[3]+roll_axis1[jj]]
     tmptotmask = Mask_total[edges_ok[0]+roll_axis0[jj]:edges_ok[1]+roll_axis0[jj],edges_ok[2]+roll_axis1[jj]:edges_ok[3]+roll_axis1[jj]]
     tmpbadmask = TotBadMask[edges_ok[0]+roll_axis0[jj]:edges_ok[1]+roll_axis0[jj],edges_ok[2]+roll_axis1[jj]:edges_ok[3]+roll_axis1[jj]]
     
     #Verify this mask does not overlap with total mask or with bad stars or that there are nans
     if ((tmptotmask+tmpbadmask)*tmpbkgmask).sum() >0 or np.logical_not(np.isfinite(tmpdata[tmpbkgmask==1])).sum() > 0:
       continue
     
     tmp_back[cnt] = np.nansum(tmpbkgmask*tmpdata)
     cnt +=1
     print(jj, cnt)
   
   if cnt < 5000:
     Success[ii] = 0
     print("WARNING {0} {1}".format(ii+1, cnt))
   
   #Now apply scaling if necessary
   tmp_back = tmp_back/Nbkgpix_final*Nbkgpix_orig
   
   #Now clean the signal
   
   for rep in range(3):
     sigma = np.nanstd(tmp_back)
     mean = np.nanmean(tmp_back)
     
     ok = (tmp_back>mean-5.*sigma) & (tmp_back<mean+5.*sigma) & (tmp_back!=0)
   
     tmp_back = tmp_back[ok]
     
   Background[ii] = mean
   Background_sig[ii] = sigma
   
   if ii == Nregions-1:
     n, bins, patches = mp.hist(tmp_back, 30, density=True)
     import matplotlib.mlab as mlab
     #mp.plot(bins, mlab.normpdf(bins, mean, sigma), 'r-', linewidth=2)
     #mp.show()

Flux = np.zeros(Nregions)
Flux_err = np.zeros(Nregions)

#Now do proper aperture photometry
if MONTEERR:
  print('Estimating Montecarlo errors on the source flux')

print('---------------------------------')

print('Source SourceERR  Background BackgroundERR')
for ii in range(Nregions):
  
  #Total flux in aperture
  Flux[ii] = np.nansum(Masks[ii,:,:]*data)-Background[ii]

  
  if MONTEERR:
    #Montecarlo errors are requested data and error must exists
    Monteflux = np.zeros(100)
    for monte in range(100):
      Montedata = data+np.random.normal(size=(size[0],size[1]))*error
      Monteflux[monte] = np.nansum(Masks[ii,:,:]*Montedata)
    Err_source = np.std(Monteflux)
      
  else:  
    if Flux[ii]+Background[ii]>0:
      Err_source = np.sqrt((Flux[ii])*Exptime)/(Exptime)
    else:
      Err_source = 0
    if Exptime==1:
      Err_source = 0
    
  Err_back = Background_sig[ii]
  
  print('{0:6.5f} {1:6.5f} {2:6.5f} {3:6.5f}'.format(Flux[ii], Err_source, Background[ii], Err_back))
  
  Flux_err[ii] = np.sqrt(Err_source**2+Err_back**2)

  if ii == Nregions-1 and Success[ii] ==0:
     Flux[ii] = np.sum(Flux[:ii])

#Are we writing the output in mJy or in magnitudes?  
if OUTMAG==1:
  magAB = (-2.5*np.log10(Flux)+ZP)-AG
  magAB_err = 1.086*Flux_err/Flux
  magAB_err = np.sqrt(magAB_err**2+(ZP_err**2))
  
  fl_mjy = 3.631E6*10**(-0.4*magAB)
  fl_mjy_err = fl_mjy*0.4*np.log(10)*magAB_err
  
else:
  fl_mjy = Flux*ZP
  fl_mjy_err = Flux_err*ZP

#Output
imfile = os.path.basename(args.IMAGE)
scr=open(args.outdir+'/'+imfile.replace('.fits','.mag'),"w")

print('-----------------------------')

if OUTMAG==1:
  print('Region   FlmJy sigFlmJy  SN  MagAB  sigMag ')
  scr.write("#Region  FlmJy sigFlmJy  SN  MagAB  sigMag \n") 
else:
  print('Region  FlmJy  sigFlmJy SN ')
  scr.write("#Region FlmJy  sigFlmJy SN \n") 

for ii in range(Nregions):	       
  if ii == Nregions-1:		       
    reg = 'TOT'
  else:
    reg = '{:02d}'.format(ii+1)        
      				       
  if OUTMAG==1:
    #Output in mJy and Magnitudes
    print(    '{0}     {1:5.4f}  {2:5.4f}  {3:5.4f}  {4:5.4f}  {5:5.4f}     '.format(reg, fl_mjy[ii], fl_mjy_err[ii], fl_mjy[ii]/fl_mjy_err[ii], magAB[ii], magAB_err[ii])) 
    scr.write('{0}     {1:5.4f}  {2:5.4f}  {3:5.4f}  {4:5.4f}  {5:5.4f}   \n'.format(reg, fl_mjy[ii], fl_mjy_err[ii], fl_mjy[ii]/fl_mjy_err[ii], magAB[ii], magAB_err[ii]))
  elif OUTMAG==0:
    #Output in mJy
    print(    '{0}     {1:5.4f}  {2:5.4f}  {3:5.4f}   '.format(reg, fl_mjy[ii], fl_mjy_err[ii], fl_mjy[ii]/fl_mjy_err[ii]))
    scr.write('{0}     {1:5.4f}  {2:5.4f}  {3:5.4f} \n'.format(reg, fl_mjy[ii], fl_mjy_err[ii], fl_mjy[ii]/fl_mjy_err[ii]))
  elif OUTMAG==-1:
    #Output in erg/cm2/s exponential notation
    print(    '{0}     {1:.3e}  {2:.3e}  {3:5.4f}   '.format(reg, fl_mjy[ii], fl_mjy_err[ii], fl_mjy[ii]/fl_mjy_err[ii]))
    scr.write('{0}     {1:.3e}  {2:.3e}  {3:5.4f} \n'.format(reg, fl_mjy[ii], fl_mjy_err[ii], fl_mjy[ii]/fl_mjy_err[ii]))

   
scr.close()





# script for imaging the calibrated data
# run on CASA 5.6.1

# Obs Date: 01-Sep-2019 22:14:28.8~00:53:45.9

# J1625-2527: phase calibration
# J1751+0939: parallactic angle calibration
# J1517-2422: flux/bandpass calibration

# see other file (tclean) for details about obs

contvis='calibrated_final_irs63.ms'
contimagename='irs63_c561'

# From the other clean we found that some antenna in the middle are:
# DA49, DA60, DA42 and DV19
# We tested DA49 and it wasn't great so we will use DA42

spwmap = [0,0,0,0]
cell = '0.018arcsec'
imsize = [1024,1024]
refant = 'DA42'
### May want to give multiple reference antenna

# save initial flags in case you don't like the final
# self-calibration. The task applycal will flag data that doesn't have
# solutions.
flagmanager(vis=contvis,mode='save',versionname='before_selfcal',merge='replace')

# Get rid of any models that might be hanging around in the image header
clearcal(contvis)
delmod(vis=contvis,field='',otf=True)



# flag data that have looked problematic and see if that improves the clean
# these flags were determined with prior runs of selfcal, 
# where some spws/antenna had scatter in phase
# in some cases, it may have affected the amplitudes

# If you are re-doing your self-cal, uncomment the next line to reset
# your corrected data column back to its original state.
#clearcal(vis=contvis)

# Strategy for self-cal
# 1) Do a shallow clean with a mask only toward the brightest real emission
# 2) Check S/N ratio and rms, if SNR > 20 and rms is higher than expected, then
#    selfcal
# 3) do gaincal/apply cal steps and then re-clean going deeper, but still be
#    conservative with the masking - it is better to miss emission than to 
#    include noise
# 4) re-image and check S/N ratio and rms, if both have improved, 
#    go back to step 3
# 5) re-image and check S/N ratio and rms, if improvement is minor, then try
#    amplitude calibration
# 6) do a final clean, going very deep and including all the emission you think
#    could be real

# Self calibration, round 1
# using tclean this time

tclean(vis = contvis,
         imagename = contimagename+'_p0',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mtmfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.00mJy',
         interactive = True,
         savemodel='modelcolumn',
         usepointing=False,
         pbcor=False,
         pblimit = 0.2,
         mask = ''
         )
# Number of iterations: 3X50 (focus on the brigh emission at the phase center)
# max residual 
# rms , peak ~ , SNR ~ 
#### NOTE: with nterms = 2 extra files are made.  The image file is the one ending in *.tt0



rmtables('pcal1_c561_irs63*')
gaincal(vis=contvis,
          caltable='pcal1_561_irs63',
          gaintype='T', #G does correction for each pol but can get spurious changes with elevation
          refant=refant,
          calmode='p',
          combine='', #changed from spw
          solint='inf',
          minsnr=3.0,
          minblperant=4)  
# some flagged data
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 
# spw 1: 
# spw 2: 
# spw 3: 

# Check the solution
plotcal(caltable='pcal1_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-90,90])
# should be flat

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
           field='',
           spwmap='',
           gaintable=['pcal1_561_irs63'],
           gainfield='',
           flagbackup=False,
           interp='linear',
           applymode='calonly') #if there are a lot of flags, calonly could be an issue

flagmanager(vis=contvis,mode='save',versionname='after_pcal1')

# split out the data
split(vis=contvis,outputvis='calibrated_final_irs63_c561_p1.ms',datacolumn='corrected')

contvis='calibrated_final_irs63_c561_p1.ms'

flagmanager(vis=contvis,mode='save',versionname='before_pcal2')

# Second round of selfcal

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p1'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p1',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mtmfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         mask=contimagename+'_p0.mask',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.0mJy',
         interactive = True,
         pbcor = False,
	 pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False
         )
# Number of iterations 1X50, 4X100
# rms ~ , peak ~ , SNR ~ 


# Have 2 sec intervals

rmtables('pcal2_561_irs63*')
gaincal(vis=contvis,
          caltable='pcal2_561_irs63',
          gaintype='T', 
          refant=refant,
          calmode='p',
          combine='', #changed from spw
          solint='30.3s', #30s gets you 15 integrations
          minsnr=3.0,
          minblperant=4)  #changed from 6
# some flagged data
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 
# spw 1: 
# spw 2:
# spw 3:


# Check the solution
plotcal(caltable='pcal2_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-90,90])
# still a fair bit of scatter

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
           field='',
           spwmap='',
           gaintable=['pcal2_561_irs63'],
           gainfield='',
           flagbackup=False,
           interp='linear',
           applymode='calonly')

flagmanager(vis=contvis,mode='save',versionname='after_pcal2')

# split out the data
os.system('rm -rf calibrated_final_irs63_c561_p2.ms*')
split(vis=contvis,outputvis='calibrated_final_irs63_c561_p2.ms',datacolumn='corrected')

contvis='calibrated_final_irs63_c561_p2.ms'

flagmanager(vis=contvis,mode='save',versionname='before_pcal3')


# Third round of selfcal

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p2'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p2',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mtmfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.0mJy',
         interactive = True,
         mask=contimagename+'_p0.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False
         )
# Number of iterations 8X100
# rms ~ , peak ~ , SNR ~ 

# try shorter solutions

rmtables('pcal3_561_irs63*')
gaincal(vis=contvis,
          caltable='pcal3_561_irs63',
          gaintype='T', #from T so correction for each pol
          refant=refant,
          calmode='p',
          combine='', 
          solint='20.2s', #10s gets you 5 integrations
          minsnr=3.0,
          minblperant=4)  #changed from 6
# lots more flagged data
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 
# spw 1: 
# spw 2:
# spw 3:
# 

# Check the solution
plotcal(caltable='pcal3_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-30,30])
# fairly flat (mostly within 10 deg)

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
           field='',
           spwmap='',
           gaintable=['pcal3_561_irs63'],
           gainfield='',
           flagbackup=False,
           interp='linear',
           applymode='calonly')

flagmanager(vis=contvis,mode='save',versionname='after_pcal3')

# split out the data
os.system('rm -rf calibrated_final_irs63_c561_p3.ms*')
split(vis=contvis,outputvis='calibrated_final_irs63_c561_p3.ms',datacolumn='corrected')

contvis='calibrated_final_irs63_c561_p3.ms'

flagmanager(vis=contvis,mode='save',versionname='before_pacal1')









# Fourth round of selfcal

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p3'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p3',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mtmfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.00mJy',
         interactive = True,
         mask=contimagename+'_p0.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False,
         )
# Number of iterations 1X100, 4X200
# rms ~ , peak ~ , SNR ~ 
# still has ripples

# no change, try amp selfcal

# try amplitude selfcal
rmtables('pacal1_561_irs63*')
gaincal(vis=contvis,
          caltable='pacal1_561_irs63',
          gaintype='T', 
          refant=refant,
          calmode='ap',
          combine='', #changed from spw
          solint='inf', 
          minsnr=3.0,
          minblperant=6)
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 
# spw 1: 
# spw 2: 
# spw 3: 
# WARNING: amplitude selfcal could change results (and you don't want it to)
# could mitigate this by setting solnorm=True in gaincal

# Check the solution
plotcal(caltable='pacal1_561_irs63',
          xaxis='time',
          yaxis='amp',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,0.5,1.5])
# solutions look fairly flat

plotcal(caltable='pacal1_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-10,10])
# some phase ripples of about 5 deg

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
           field='',
           spwmap='',
           gaintable=['pacal1_561_irs63'],
           gainfield='',
           flagbackup=False,
           applymode='calonly',
           interp='linear')

flagmanager(vis=contvis,mode='save',versionname='after_ampcal1')

# split out the data
os.system('rm -rf calibrated_final_irs63_c561_p4.ms*')
split(vis=contvis,outputvis='calibrated_final_irs63_c561_p4.ms',datacolumn='corrected')

contvis='calibrated_final_irs63_c561_p4.ms'

flagmanager(vis=contvis,mode='save',versionname='before_pacal2')


# now clean deeply

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p4'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p4',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mtmfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.00mJy',
         interactive = True,
         mask=contimagename+'_p0.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False,
         )
# Number of iterations 1X100, 7X200
# rms , peak , SNR 
# should be much better

# if ripples persist, could do a second round of selfcal with 60.6 sec intervals













# This is where i am going to mock up how to add another clean
# --------------------------------------------------------------------
# --------------------------------------------------------------------
# try amplitude selfcal again
rmtables('pacal2_561_irs63*')
gaincal(vis=contvis,
          caltable='pacal2_561_irs63',
          gaintype='T', 
          refant=refant,
          calmode='ap',
          combine='', #changed from spw
          solint='100s', # previously it said infinity and now we will try 100s
          minsnr=3.0,
          minblperant=6,
          solnorm=True)  # Normalize amplitudes to avoid artificial flux changes
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 
# spw 1: 
# spw 2: 
# spw 3: 

# Check the solution
plotcal(caltable='pacal2_561_irs63',
          xaxis='time',
          yaxis='amp',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,0.5,1.5])
# solutions look 
plotcal(caltable='pacal2_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-10,10])
# solutions look

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
           field='',
           spwmap='',
           gaintable=['pacal2_561_irs63'],
           gainfield='',
           flagbackup=False,
           applymode='calonly',
           interp='linear')

flagmanager(vis=contvis,mode='save',versionname='after_ampcal2')

# split out the data
os.system('rm -rf calibrated_final_irs63_c561_p5.ms*')
split(vis=contvis,outputvis='calibrated_final_irs63_c561_p5.ms',datacolumn='corrected')

contvis='calibrated_final_irs63_c561_p5.ms'

flagmanager(vis=contvis,mode='save',versionname='before_pacal2')


# now clean deeply

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p5'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p5',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mtmfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.00mJy',
         interactive = True,
         mask=contimagename+'_p0.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False,
         )
# Number of iterations 
# rms , peak , SNR 
# --------------------------------------------------------------------
# --------------------------------------------------------------------





# --------------------------------------------------------------------
tclean(vis = contvis,
         imagename = contimagename+'_p4_stokesQ',
         stokes = 'Q',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 1,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.02mJy',
         interactive = False,
         mask=contimagename+'_p4.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False
         )
# max residual is 
# rms , min  and max , SNR 

tclean(vis = contvis,
         imagename = contimagename+'_p4_stokesU',
         stokes = 'U',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 1,
         imsize = imsize,
         cell = cell,
         deconvolver = 'hogbom',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.02mJy',
         interactive = False,
         mask=contimagename+'_p4.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False
         )
# max residual is
# rms , min  and max , SNR 

exportfits(imagename='irs63_c561_p4.image',fitsimage='IRS63_StokesI_clean_nopbcorr.fits')
exportfits(imagename='irs63_c561_p4_stokesQ.image',fitsimage='IRS63_StokesQ_clean_nopbcorr.fits')
exportfits(imagename='irs63_c561_p4_stokesU.image',fitsimage='IRS63_StokesU_clean_nopbcorr.fits')
# --------------------------------------------------------------------


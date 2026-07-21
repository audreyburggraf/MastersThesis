# Things got a bit messed up with the time between cleans and shutting stuff 
# down. I was not that far along so I will start fresh

# Start date: Monday January 5th 2026 @ 12:23 PM

# script for imaging the calibrated data
# run on CASA 5.6.1


# Obs Date: 01-Sep-2019 22:14:28.8~00:53:45.9

# J1625-2527: phase calibration
# J1751+0939: parallactic angle calibration
# J1517-2422: flux/bandpass calibration


# Here are the things that I need to load in
#-----------------------------------------------------------------
contvis='calibrated_final_irs63_c561.ms'
contimagename='irs63_c561'
#-----------------------------------------------------------------
plotants(vis = contvis)
# Some in the middle are:
# DV17, DV02, DV20
# From the other clean we found that some antenna in the middle are:

spwmap = [0,0,0,0]
cell = '0.018arcsec'
imsize = [1024,1024]
refant = 'DV17'
#------------------------------------------------------------------



# Save initial flags
#-----------------------------------------------------------------
# in case you don't like the final
# self-calibration. The task applycal will flag data that doesn't have
# solutions.
flagmanager(vis=contvis,mode='save',versionname='before_selfcal',merge='replace')
#----------------------------------------------------------------------------


# If redoing self-cal:
#----------------------------------------------------------------------
# Get rid of any models that might be hanging around in the image header
clearcal(contvis)
delmod(vis=contvis,field='',otf=True)
#----------------------------------------------------------------------



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




# ROUND 1 OF SELF CALIBRATION
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p0'+ ext)

# Self calibration, round 1
#------------------------------------
# using tclean this time
tclean(vis = contvis,
         imagename = contimagename+'_p0',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mfs', 
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.00mJy',
         interactive = True,
         savemodel='modelcolumn',
         usepointing=False,
         pbcor=False,
         pblimit = -0.2,
         mask = ''
         )
#-----------------------------------------------------------
# Number of iterations: 3/3X50

# Max Residual: approx 2 mJy
# rms ~0.035 - 0.046 mJy, peak ~26  mJy, SNR ~565-750 mJy
#### NOTE: with nterms = 2 extra files are made.  The image file is the one ending in *.tt0
#-------------------------------------------------------------


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
#-----------------------------------------------------------
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 17/17/17
# spw 1: 17/17/17
# spw 2: 17/17/17
# spw 3: 17/17/17
#-----------------------------------------------------------

# Check the solution
plotcal(caltable='pcal1_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-90,90])
# should be flat
# results: not all plots were flat

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
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------





# ROUND 2 OF SELF CALIBRATION
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p1'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p1',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
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
#-----------------------------------------------------------
# Number of iterations: 1/1X50, 4/4X100
# Max Residual: 0.44 mJy
# rms ~0.014-0.019 mJy, peak ~26 mJy, SNR ~1370-1860 mJy
#-------------------------------------------------------------

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
#-----------------------------------------------------------
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 76/76/76
# spw 1: 76/76/76
# spw 2: 76/76/76
# spw 3: 76/76/76
#-----------------------------------------------------------


# Check the solution
plotcal(caltable='pcal2_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-90,90])
# still a fair bit of scatter (this text was already here but holds true)

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
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------




# ROUND 3 OF SELF CALIBRATION
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
	rmtables(contimagename + '_p2'+ ext)

tclean(vis = contvis,
	imagename = contimagename+'_p2',
	stokes = 'I',
	outframe = 'LSRK',
	specmode = 'mfs', 
	nterms = 2,
	imsize = imsize,
	cell = cell,
	deconvolver = 'mtmfs',
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
#-----------------------------------------------------------
# Number of iterations: 8/8X100

# left off here on monday january 5th at 6:27 pm 


# start back up tuesday january 6th at 3:39 pm

# Max Residual: 0.145 mJy
# rms ~0.011-0.012 mJy, peak ~26 mJy, SNR ~2170-2370 mJy
#-------------------------------------------------------------

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

#-----------------------------------------------------------
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 116/116/116
# spw 1: 116/116/116
# spw 2: 116/116/116
# spw 3: 116/116/116
#-----------------------------------------------------------

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

flagmanager(vis=contvis,mode='save',versionname='before_pacal4')
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------





# ROUND 4 OF SELF CALIBRATION
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p3'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p3',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mfs', 
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
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
#-----------------------------------------------------------
# Number of iterations: 1/1X100, 4/4X200
# Max Residual: 0.0113
# rms ~0.0099 - 0.0116 mJy, peak ~26 mJy, SNR ~2240 - 2600 mJy
# Last clean S/N was approx. 2170-2370

# if no change, try amp selfcal
#-------------------------------------------------------------

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
#-----------------------------------------------------------
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 17/17/17 
# spw 1: 17/17/17
# spw 2: 17/17/17
# spw 3: 17/17/17
#-----------------------------------------------------------

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
# solutions look fairly flat (this text was here before but holds)

plotcal(caltable='pacal1_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-10,10])
# some phase ripples of about 5 deg (holds)

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
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------





# ROUND 5 OF SELF CALIBRATION
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p4'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p4',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
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
#-----------------------------------------------------------
# Number of iterations: 1/1X100, 7/7X200

# leaving off here tuesday january 6th at 6:30 pm 
# starting back up on wednesday january 7th at 12:20 pm

# Max Residual: 0.038 mJy
# rms ~0.0095-0.0101 mJy, peak ~26 mJy, SNR ~2575-2735 mJy
# reminder last clean was S/N of 2240 - 2600 


# if ripples persist, could do a second round of selfcal with 60.6 sec intervals
#-------------------------------------------------------------

rmtables('pacal2_561_irs63*')
gaincal(vis=contvis,
          caltable='pacal2_561_irs63',
          gaintype='T', 
          refant=refant,
          calmode='ap',
          combine='', #changed from spw
          solint='100s', # changed to 100s 
          minsnr=3.0,
          minblperant=6,
          solnorm=True)  # solnorm=True keeps flux scale stable
#-----------------------------------------------------------
# calibration solve statistics per spw: (exp/attempt/suc):
# spw 0: 29/29/29
# spw 1: 29/29/29
# spw 2: 29/29/29
# spw 3: 29/29/29
#-----------------------------------------------------------

# Check the solution
plotcal(caltable='pacal2_561_irs63',
          xaxis='time',
          yaxis='amp',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,0.8,1.2])
# most solutions look very flat
# some have some scatter but within 0.1

plotcal(caltable='pacal2_561_irs63',
          xaxis='time',
          yaxis='phase',
          timerange='',
          iteration='antenna',
          subplot=331,
          plotrange=[0,0,-5,5])
# solutions look: flat with a little spread within 4 ish degrees

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

flagmanager(vis=contvis,mode='save',versionname='before_pacal3')
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------





# ROUND 6 OF SELF CALIBRATION
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
      rmtables(contimagename + '_p5'+ ext)

tclean(vis = contvis,
         imagename = contimagename+'_p5',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mfs', 
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
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
#-----------------------------------------------------------
# I am not sure if i need to clean this much. I will do it and clean Q and U wrt p5 and not p5
# I can always come back and re-clean Q and U wrt p4 
# Update: xterm crashed during the clean anyways

# Number of iterations: 1/1X100, 0/7X200
# Max Residual:
# rms ~ mJy, peak ~ mJy, SNR ~ mJy
#-------------------------------------------------------------
#-------------------------------------------------------------



# CLEAN STOKES Q
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
tclean(vis = contvis,
         imagename = contimagename+'_p4_stokesQ',
         stokes = 'Q',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 1,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
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
#-----------------------------------------------------------
# Max Residual:
# rms ~ mJy, peak ~ mJy, SNR ~ mJy
#-------------------------------------------------------------
#-------------------------------------------------------------




# CLEAN STOKES U
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
tclean(vis = contvis,
         imagename = contimagename+'_p4_stokesU',
         stokes = 'U',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 1,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
         niter = 10000,
         weighting = 'briggs',
         robust = 0.5,
         gridder = 'standard',
         threshold = '0.02mJy',
         interactive = False,
         mask=contimagename+'_4.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False
         )
#-----------------------------------------------------------
# Max Residual:
# rms ~ mJy, peak ~ mJy, SNR ~ mJy
#-------------------------------------------------------------
#-------------------------------------------------------------




# Save data
# ----------------------------------------------------------------------------------------------------
exportfits(imagename='irs63_c561_p4.image.tt0',fitsimage='IRS63_BAND4_nterms2_StokesI_clean_nopbcorr.fits')
exportfits(imagename='irs63_c561_p4_stokesQ.image.tt0',fitsimage='IRS63_BAND4_nterms2_StokesQ_clean_nopbcorr.fits')
exportfits(imagename='irs63_c561_p4_stokesU.image.tt0',fitsimage='IRS63_BAND4_nterms2_StokesU_clean_nopbcorr.fits')
# -----------------------------------------------------------------------------------------------------


# Monday March 9th 
# We want to smooth the data w.r.t Band 5 beam size 

# Smooth Stokes I
imsmooth(
imagename="irs63_c561_p4.image.tt0",
outfile="irs63_c561_p4_smooth_wrt_b5.image.tt0",
major="0.329219arcsec",
minor="0.268545arcsec",
pa="-88.931deg",
targetres=True
)

# Smooth Stokes Q
imsmooth(
imagename="irs63_c561_p4_stokesQ.image.tt0",
outfile="irs63_c561_p4_StokesQ_smooth_wrt_b5.image.tt0",
major="0.329219arcsec",
minor="0.268545arcsec",
pa="-88.931deg",
targetres=True
)


# Smooth Stokes U
imsmooth(
imagename="irs63_c561_p4_stokesU.image.tt0",
outfile="irs63_c561_p4_StokesU_smooth_wrt_b5.image.tt0",
major="0.329219arcsec",
minor="0.268545arcsec",
pa="-88.931deg",
targetres=True
)




# -----------------------------------------------------------------------------------------------------
# Save the smoothed data
# -----------------------------------------------------------------------------------------------------
exportfits(imagename='irs63_c561_p4_smooth_wrt_b5.image.tt0',fitsimage='B4_SI_smooth_nt2_wrt_b5.fits')
exportfits(imagename='irs63_c561_p4_StokesQ_smooth_wrt_b5.image.tt0',fitsimage='B4_SQ_smooth_nt2_wrt_b5.fits')
exportfits(imagename='irs63_c561_p4_StokesU_smooth_wrt_b5.image.tt0',fitsimage='B4_SU_smooth_nt2_wrt_b5.fits')
# -----------------------------------------------------------------------------------------------------


# Wednesday March 11th
# Second smooth
# Smooth Band 4 data to match the Band 6 beam (targetres = true)


# Smooth Stokes I
imsmooth(
imagename="irs63_c561_p4.image.tt0",
outfile="BAND4_StokesI_smooth_wrt_b6.image",
major="0.265288arcsec",
minor="0.203718arcsec",
pa="-68.5533deg",
targetres=True
)


# Smooth Stokes Q
imsmooth(
imagename="irs63_c561_p4_stokesQ.image.tt0",
outfile="BAND4_StokesQ_smooth_wrt_b6.image",
major="0.265288arcsec",
minor="0.203718arcsec",
pa="-68.5533deg",
targetres=True
)


# Smooth Stokes U
imsmooth(
imagename="irs63_c561_p4_stokesU.image.tt0",
outfile="BAND4_StokesU_smooth_wrt_b6.image",
major="0.265288arcsec",
minor="0.203718arcsec",
pa="-68.5533deg",
targetres=True
)




# -----------------------------------------------------------------------------------------------------
# Save the smoothed data
# -----------------------------------------------------------------------------------------------------
exportfits(imagename='BAND4_StokesI_smooth_wrt_b6.image',fitsimage='B4_SI_smooth_wrt_B6.fits')
exportfits(imagename='BAND4_StokesQ_smooth_wrt_b6.image',fitsimage='B4_SQ_smooth_wrt_B6.fits')
exportfits(imagename='BAND4_StokesU_smooth_wrt_b6.image',fitsimage='B4_SU_smooth_wrt_B6.fits')
# -----------------------------------------------------------------------------------------------------



# Wednesday March 11th
# Third smooth
# Now, smooth the (Band 4 smoothed by Band 6 beam) with the Band 7 beam



# Smooth Stokes I
imsmooth(
imagename="BAND4_StokesI_smooth_wrt_b6.image",
outfile="BAND4_StokesI_smooth_wrt_B6_B7.image",
major="0.280187arcsec",
minor="0.172577arcsec",
pa="-81.229965deg",
targetres=False
)


# Smooth Stokes Q
imsmooth(
imagename="BAND4_StokesQ_smooth_wrt_b6.image",
outfile="BAND4_StokesQ_smooth_wrt_B6_B7.image",
major="0.280187arcsec",
minor="0.172577arcsec",
pa="-81.229965deg",
targetres=False
)


# Smooth Stokes U
imsmooth(
imagename="BAND4_stokesU_smooth_wrt_b6.image",
outfile="BAND4_StokesU_smooth_wrt_B6_B7.image",
major="0.280187arcsec",
minor="0.172577arcsec",
pa="-81.229965deg",
targetres=False
)




# -----------------------------------------------------------------------------------------------------
# Save the smoothed data
# -----------------------------------------------------------------------------------------------------
exportfits(imagename='BAND4_StokesI_smooth_wrt_B6_B7.image',fitsimage='B4_SI_smooth_wrt_B6_B7.fits')
exportfits(imagename='BAND4_StokesQ_smooth_wrt_B6_B7.image',fitsimage='B4_SQ_smooth_wrt_B6_B7.fits')
exportfits(imagename='BAND4_StokesU_smooth_wrt_B6_B7.image',fitsimage='B4_SU_smooth_wrt_B6_B7.fits')
# -----------------------------------------------------------------------------------------------------


# 
# Thursday March 19th at 11:38 AM
# We are going to test robust -1 on the Band 4 data 


# reputting this down here to make it easy 
contvis = 'calibrated_final_irs63_c561_p4.ms'

# REDEO ROUND 5 OF SELF CALIBRATION with robust = -1
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------



tclean(vis = contvis,
         imagename = contimagename+'_p4_robust_minus1',
         stokes = 'I',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
         niter = 10000,
         weighting = 'briggs',
         robust = -1,
         gridder = 'standard',
         threshold = '0.00mJy',
         interactive = True,
         mask=contimagename+'_p0.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False,
         )

#-----------------------------------------------------------
# Number of iterations: 1/1X100, 7/7X200
# Max Residual:
# rms ~ mJy, peak ~ mJy, SNR ~ mJy

#-------------------------------------------------------------



# CLEAN STOKES Q with robust = -1
# ----------------------------------------------------------->
# ----------------------------------------------------------->


tclean(vis = contvis,
         imagename = contimagename+'_p4_robust_minus1_stokesQ',
         stokes = 'Q',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
         niter = 10000,
         weighting = 'briggs',
         robust = -1,
         gridder = 'standard',
         threshold = '0.02mJy',
         interactive = False,
         mask=contimagename+'_p4.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False
         )



#-----------------------------------------------------------
# Max Residual:
# rms ~0.0283-0.0291 mJy, peak ~18.2 mJy, SNR ~ mJy
#-------------------------------------------------------------
#-------------------------------------------------------------




# CLEAN STOKES U with robust = -1
# ---------------------------------------------------------------------
# --------------------------------------------------------------------


tclean(vis = contvis,
         imagename = contimagename+'_p4_robust_minus1_stokesU',
         stokes = 'U',
         outframe = 'LSRK',
         specmode = 'mfs',
         nterms = 2,
         imsize = imsize,
         cell = cell,
         deconvolver = 'mtmfs',
         niter = 10000,
         weighting = 'briggs',
         robust = -1,
         gridder = 'standard',
         threshold = '0.02mJy',
         interactive = False,
         mask=contimagename+'_p4.mask',
         pbcor = False,
         pblimit = -0.2,
         savemodel='modelcolumn',
         usepointing=False
         )




#-----------------------------------------------------------
# Max Residual:
# rms ~ mJy, peak ~ mJy, SNR ~ mJy
#-------------------------------------------------------------
#-------------------------------------------------------------









# --------------------------------------------------------
# Save the robust = -1 data 
# --------------------------------------


# -----------------------------------------------------------------------------------------------------
# Save the robust = -1 data 
# -----------------------------------------------------------------------------------------------------
exportfits(imagename='irs63_c561_p4_robust_minus1.image.tt0',fitsimage='B4_SI_nter2_rob_min1.fits')
exportfits(imagename='irs63_c561_p4_robust_minus1_stokesQ.image.tt0',fitsimage='B4_SQ_nter2_rob_min1.fits')
exportfits(imagename='irs63_c561_p4_robust_minus1_stokesU.image.tt0',fitsimage='B4_SU_nter2_rob_min1.fits')
# -----------------------------------------------------------------------------------------------------




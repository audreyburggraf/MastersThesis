import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm

# File path to your FITS image
fits_image_filename = '/Users/pnozari/Library/Mobile Documents/com~apple~CloudDocs/Documents/MSc Thesis/noema/ANALYSIS/NW167/NW167-LSB-I-sc.fits'

# Open the FITS file and load the primary HDU (Header Data Unit)
hdu = fits.open(fits_image_filename)[0]

# Define a function to flatten the HDU by removing unwanted axes
def flatten_alma_hdu(hduX):
    # Remove 3rd and 4th axes (e.g., Stokes and Frequency/Velocity axes) from the header
    del hduX.header['NAXIS3']
    del hduX.header['CTYPE3']
    del hduX.header['CDELT3']
    del hduX.header['CRPIX3']
    del hduX.header['CROTA3']
    del hduX.header['CRVAL3']
    del hduX.header['CUNIT3']
    
    # Update the header to indicate a 2D image
    hduX.header['NAXIS'] = 2
    
    # Squeeze the data to remove extra dimensions
    hduX.data = np.squeeze(hduX.data)
    return hduX

# Flatten the HDU to work with a 2D image
hdu = flatten_alma_hdu(hdu)

# Create a WCS (World Coordinate System) object from the 2D header
hdr = WCS(hdu.header)

# Convert the data units from Jy/beam to mJy/beam for easier visualization
hdu.data = hdu.data * 1000

# Prepare the figure for plotting with a 6x6 inch size and black edge color
fig = plt.figure(figsize=(6, 6), edgecolor='black')

# Set up the plot with WCS coordinates and adjust spacing within the figure
ax = fig.add_axes([0.17, 0.17, 0.8, 0.8], projection=hdr)

# Remove grid lines for a cleaner look
ax.grid(alpha=0.0)

# Customize tick parameters for width, length, and direction; set font size for tick labels
ax.tick_params(width=3, length=8, direction='in', labelsize=15)

# Display the image data, setting the color map to 'magma' and adjusting intensity scaling
plt.imshow(hdu.data, origin='lower', cmap='magma_r', vmin=hdu.data.min(), vmax=hdu.data.max())

# Set axis labels with font size and padding
ax.set_ylabel("Declination", fontsize=15, labelpad=0.0)
ax.set_xlabel("Right Ascension", fontsize=15, labelpad=1)

# Add a colorbar to show the intensity scale, with label size and orientation adjustments
cbar = plt.colorbar(shrink=0.8)
cbar.ax.tick_params(labelsize=15)
cbar.ax.set_ylabel('Intensity (mJy/beam)', fontsize=15, rotation=270, labelpad=30)

# Add a title to the plot, indicating the source name, with specified position and font size
ax.text(0.05, 0.95, 'NW167', transform=ax.transAxes, fontsize=20, color='black', va='top', ha='left')

# Display the plot
plt.show()

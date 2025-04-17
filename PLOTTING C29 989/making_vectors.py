import run_all_intro_stuff



max_length_pix = 400  # Maximum length of the vector in pixels for 100% polarization
reference_fraction = 0.03
nx, ny = StokesI_data_2d_mJy.shape
step = 7 



# List to store vector data
vector_data_actual = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10): 
            
            # Get the polarization angle at this pixel
            angle_rad = calculated_polarization_angle_rad[y, x] + np.pi / 2
            
            # Get the polarization fraction at this pixel
            polarization_fraction = calc_polarized_frac[y, x]  # Add your data source
            
            # Scale vector length by polarization fraction
            vector_length_pix = max_length_pix * polarization_fraction
            
            # Compute the vector components
            dx = vector_length_pix * np.cos(angle_rad)
            dy = vector_length_pix * np.sin(angle_rad)
            
            # Store vector data in a row
            vector_data_actual.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
            
            
            
# List to store vector data
vector_data_actual_same_length = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10): 
            
            # Get the polarization angle at this pixel
            angle_rad = calculated_polarization_angle_rad[y, x] + np.pi / 2
            
            # Compute the vector components
            dx = vector_length_pix_const * np.cos(angle_rad)
            dy = vector_length_pix_const * np.sin(angle_rad)
            
            # Store vector data in a row
            vector_data_actual_same_length.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
            
            
            
            
# Set the angle 
uniform_angle_rad = minor_angle_rad + np.pi/2
#-np.pi/4

# List to store vector data
vector_data_uniform = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10): 
            
            # Compute the vector components
            dx = vector_length_pix_const * np.cos(uniform_angle_rad)
            dy = vector_length_pix_const * np.sin(uniform_angle_rad)
            
            # Store vector data in a row
            vector_data_uniform.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
            
            
# List to store vector data for azimuthal polarization
vector_data_azimuthal = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10): 
            
            
            # Calculate the azimuthal angle (in radians) based on the position
            dx_center = x - RA_centre_pix  # Center of the disk (or image)
            dy_center = y - Dec_centre_pix  # Center of the disk (or image)
            
            # I got chatgpt to help me with this angle but i added the np.pi/2
            azimuthal_angle_rad = np.arctan2(dy_center, dx_center) + np.pi/2 # Calculate the angle

            # Compute the vector components using the azimuthal angle
            dx = vector_length_pix_const * np.cos(azimuthal_angle_rad)
            dy = vector_length_pix_const * np.sin(azimuthal_angle_rad)
            
            # Store vector data in a row
            vector_data_azimuthal.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
            
            
            
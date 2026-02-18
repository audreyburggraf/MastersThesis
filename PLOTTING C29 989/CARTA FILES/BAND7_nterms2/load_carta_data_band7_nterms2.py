import numpy as np
import pandas as pd 

loc = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7_nterms2/"




# ------------------------------------------------------------------------------------------------------
carta_major_df_StokesI = pd.read_csv(loc + "IRS63_StokesI_CARTA_MAJOR_AXIS_BAND7_nterms2.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_major_offset_StokesI = np.array(carta_major_df_StokesI["Offset (arcsec)"])
carta_major_data_mJy_StokesI   = np.array(carta_major_df_StokesI["Value (Jy/beam)"]) * 1000
# ------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------
carta_minor_df_StokesI = pd.read_csv(loc + "IRS63_StokesI_CARTA_MINOR_AXIS_BAND7_nterms2.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_minor_offset_StokesI = np.array(carta_minor_df_StokesI["Offset (arcsec)"])
carta_minor_data_mJy_StokesI   = np.array(carta_minor_df_StokesI["Value (Jy/beam)"]) * 1000
# ------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------
carta_major_df_POLI = pd.read_csv(loc + "IRS63_POLI_CARTA_MAJOR_AXIS_BAND7_nterms2_debiased.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_major_offset_POLI = np.array(carta_major_df_POLI["Offset (arcsec)"])
carta_major_data_mJy_POLI   = np.array(carta_major_df_POLI["Value (Jy/beam)"]) 
# ------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------
carta_minor_df_POLI = pd.read_csv(loc + "IRS63_POLI_CARTA_MINOR_AXIS_BAND7_nterms2_debiased.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_minor_offset_POLI = np.array(carta_minor_df_POLI["Offset (arcsec)"])
carta_minor_data_mJy_POLI = np.array(carta_minor_df_POLI["Value (Jy/beam)"]) 
# ------------------------------------------------------------------------------------------------------







# # ------------------------------------------------------------------------------------------------------
# carta_major_df_POLF = pd.read_csv(carta_image_path_band5 + "c2d_989_POLF_233GHz_CARTA_MAJOR_AXIS_BAND5.csv", 
#                        delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
#                        header=None, 
#                        skiprows = 6,
#                        names=["Offset (arcsec)", "Value"])  

# carta_major_offset_POLF = np.array(carta_major_df_POLF["Offset (arcsec)"])
# carta_major_data_POLF   = np.array(carta_major_df_POLF["Value"]) 
# # ------------------------------------------------------------------------------------------------------




# # ------------------------------------------------------------------------------------------------------
# carta_minor_df_POLF = pd.read_csv(carta_image_path_band5 + "c2d_989_POLF_233GHz_CARTA_MINOR_AXIS_BAND5.csv", 
#                        delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
#                        header=None, 
#                        skiprows = 6,
#                        names=["Offset (arcsec)", "Value"])  

# carta_minor_offset_POLF = np.array(carta_minor_df_POLF["Offset (arcsec)"])
# carta_minor_data_POLF = np.array(carta_minor_df_POLF["Value"])
# # ------------------------------------------------------------------------------------------------------


# # Write down some information to check if my slice and the CARTA slice are different 
# CARTA_band7_nterms2_info = 
# {
# centre_ra_pix: 515.045
# centre_dec_pix: 504.278

# }

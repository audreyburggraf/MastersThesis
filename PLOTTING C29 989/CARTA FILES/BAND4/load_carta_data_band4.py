import numpy as np
import pandas as pd 

carta_image_path_band4 = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND4/"




# ------------------------------------------------------------------------------------------------------
carta_major_df_StokesI = pd.read_csv(carta_image_path_band4 + "c2d_989_StokesI_233GHz_CARTA_MAJOR_AXIS_BAND4.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_major_offset_StokesI = np.array(carta_major_df_StokesI["Offset (arcsec)"])
carta_major_data_mJy_StokesI   = np.array(carta_major_df_StokesI["Value (Jy/beam)"]) * 1000
# ------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------
carta_minor_df_StokesI = pd.read_csv(carta_image_path_band4 + "c2d_989_StokesI_233GHz_CARTA_MINOR_AXIS_BAND4.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_minor_offset_StokesI = np.array(carta_minor_df_StokesI["Offset (arcsec)"])
carta_minor_data_mJy_StokesI   = np.array(carta_minor_df_StokesI["Value (Jy/beam)"]) * 1000
# ------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------
carta_major_df_POLI = pd.read_csv(carta_image_path_band4 + "c2d_989_POLI_233GHz_CARTA_MAJOR_AXIS_BAND4.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_major_offset_POLI = np.array(carta_major_df_POLI["Offset (arcsec)"])
carta_major_data_mJy_POLI   = np.array(carta_major_df_POLI["Value (Jy/beam)"]) * 1000
# ------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------
carta_minor_df_POLI = pd.read_csv(carta_image_path_band4 + "c2d_989_POLI_233GHz_CARTA_MINOR_AXIS_BAND4.csv", 
                       delim_whitespace=True,  # Auto-detects spaces/tabs as delimiter
                       header=None, 
                       skiprows = 6,
                       names=["Offset (arcsec)", "Value (Jy/beam)"])  

carta_minor_offset_POLI = np.array(carta_minor_df_POLI["Offset (arcsec)"])
carta_minor_data_mJy_POLI = np.array(carta_minor_df_POLI["Value (Jy/beam)"]) * 1000
# ------------------------------------------------------------------------------------------------------

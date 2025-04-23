import sys

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")

# Now you can import constants.py
import constants

# Use the variable from constants.py
functions_folder_path = constants.functions_folder_path


from FITS_Image_Functions import *
from PlottingWithFunction import * 
from custom_colormap import *
from PolarizationFunctions import *
from DataAnalysisFunctions import *
from GaussianFunctions import *
from RatioModelFunctions import *
from MakingGridFunctions import * 
from SlicesFunctions import * 
from IntroductionFunctions import *
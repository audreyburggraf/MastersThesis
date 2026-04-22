# ------------------------------------------------------------
# Import nexessary things
# ------------------------------------------------------------
import numpy as np
# ------------------------------------------------------------
#
#
#
## ------------------------------------------------------------
# Function to convert mm to cm
# ------------------------------------------------------------
def micron_to_cm(value_micron):
    value_micron = np.array(value_micron)  # convert input to numpy array
    return value_micron * 1e-4  # 1 micron = 1e-4 cm
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# Function to convert micron to cm
# ------------------------------------------------------------
def mm_to_cm(value_mm):
    value_mm = np.array(value_mm)  # convert input to numpy array
    return value_mm * 1e-1  # 1 mm = 0.1 cm = 1e-1 cm
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# Function to convert cm to micron
# ------------------------------------------------------------
def cm_to_micron(value_cm):
    value_cm = np.array(value_cm)  # convert input to numpy array
    return value_cm * 10000
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# Function to convert mm to micron
# ------------------------------------------------------------
def mm_to_micron(value_mm):
    value_mm = np.array(value_mm)  # convert input to numpy array
    return value_mm * 1000
# ------------------------------------------------------------
#
#
#
#
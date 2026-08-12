import sys
import matplotlib.pyplot as plt

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")

# Now you can import constants.py
import constants

# Use the variable from constants.py
testing_ratios = constants.testing_ratios


# Import Functions
from FITS_Image_Functions import *
from PolarizationFunctions import *
from DataAnalysisFunctions import *



def generate_ratio_key(ratio1, ratio2):
    return f"{int(ratio1 * 100)}{int(ratio2 * 100)}"


def unpack_results(results, ratios):
    unpacked_data = {}
    for r1, r2 in ratios:
        key = generate_ratio_key(r1, r2)
        unpacked_data[key] = {
            "PA_grid": results[f"PA_grid_{key}"],
            "StokesQ_grid": results[f"StokesQ_grid_{key}"],
            "StokesU_grid": results[f"StokesU_grid_{key}"],
            "vectors_data": results[f"vectors_data_{key}"],
            "vectors_angle": results[f"vectors_angle_{key}"]
        }
    return unpacked_data


# def ratio_model_band6(StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
#                       StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
#                       ny, nx,
#                       vector_length_pix_const,
#                       StokesI_mJy, StokesI_err_mJy, 
#                       POLI_mJy, POLI_err_mJy,
#                       PA_err_deg,
#                       step = None):
    
#     if step is None:
#         step = constants.step_band6

#     results = {}

#     for ratio1, ratio2 in testing_ratios:
#         PA_grid, StokesQ_grid, StokesU_grid, vectors_data, vectors_angle = mix_StokesQU_and_generate_vectors_band6(
#             ratio1, ratio2,
#             StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
#             StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
#             ny, nx,
#             StokesI_mJy, StokesI_err_mJy, 
#             POLI_mJy, POLI_err_mJy,
#             PA_err_deg,
#             step)

#         key = generate_ratio_key(ratio1, ratio2)

#         results[f"PA_grid_{key}"]       = PA_grid
#         results[f"StokesQ_grid_{key}"]  = StokesQ_grid
#         results[f"StokesU_grid_{key}"]  = StokesU_grid
#         results[f"vectors_data_{key}"]  = vectors_data
#         results[f"vectors_angle_{key}"] = vectors_angle

#     unpacked = unpack_results(results, testing_ratios)

#     # Create clearly named vector data variables
#     vectors_data_100U_0A_cartesian  = unpacked["1000"]["vectors_data"]
#     vectors_data_90U_10A_cartesian  = unpacked["9010"]["vectors_data"]
#     vectors_data_80U_20A_cartesian  = unpacked["8020"]["vectors_data"]
#     vectors_data_70U_30A_cartesian  = unpacked["7030"]["vectors_data"]
#     vectors_data_60U_40A_cartesian  = unpacked["6040"]["vectors_data"]
#     vectors_data_50U_50A_cartesian  = unpacked["5050"]["vectors_data"]
#     vectors_data_40U_60A_cartesian  = unpacked["4060"]["vectors_data"]
#     vectors_data_30U_70A_cartesian  = unpacked["3070"]["vectors_data"]
#     vectors_data_20U_80A_cartesian  = unpacked["2080"]["vectors_data"]
#     vectors_data_10U_90A_cartesian  = unpacked["1090"]["vectors_data"]
#     vectors_data_0U_100A_cartesian  = unpacked["0100"]["vectors_data"]

#     # Save in list for return
#     vector_data_plotting_grid = [
#         vectors_data_100U_0A_cartesian, vectors_data_0U_100A_cartesian,
#         vectors_data_90U_10A_cartesian, vectors_data_10U_90A_cartesian,
#         vectors_data_80U_20A_cartesian, vectors_data_20U_80A_cartesian,
#         vectors_data_70U_30A_cartesian, vectors_data_30U_70A_cartesian,
#         vectors_data_60U_40A_cartesian, vectors_data_40U_60A_cartesian,
#         vectors_data_50U_50A_cartesian, vectors_data_50U_50A_cartesian  # repeated for consistency
#     ]

#     vector_data_list_100U_to_100A = [
#         vectors_data_100U_0A_cartesian, 
#         vectors_data_90U_10A_cartesian, 
#         vectors_data_80U_20A_cartesian, 
#         vectors_data_70U_30A_cartesian, 
#         vectors_data_60U_40A_cartesian, 
#         vectors_data_50U_50A_cartesian, 
#         vectors_data_40U_60A_cartesian,
#         vectors_data_30U_70A_cartesian,
#         vectors_data_20U_80A_cartesian,
#         vectors_data_10U_90A_cartesian,
#         vectors_data_0U_100A_cartesian,
#     ]
    
    
#     # Create clearly named vector data variables
#     vectors_angle_100U_0A_rad_astronomy  = unpacked["1000"]["vectors_angle"]
#     vectors_angle_90U_10A_rad_astronomy  = unpacked["9010"]["vectors_angle"]
#     vectors_angle_80U_20A_rad_astronomy  = unpacked["8020"]["vectors_angle"]
#     vectors_angle_70U_30A_rad_astronomy  = unpacked["7030"]["vectors_angle"]
#     vectors_angle_60U_40A_rad_astronomy  = unpacked["6040"]["vectors_angle"]
#     vectors_angle_50U_50A_rad_astronomy  = unpacked["5050"]["vectors_angle"]
#     vectors_angle_40U_60A_rad_astronomy  = unpacked["4060"]["vectors_angle"]
#     vectors_angle_30U_70A_rad_astronomy  = unpacked["3070"]["vectors_angle"]
#     vectors_angle_20U_80A_rad_astronomy  = unpacked["2080"]["vectors_angle"]
#     vectors_angle_10U_90A_rad_astronomy  = unpacked["1090"]["vectors_angle"]
#     vectors_angle_0U_100A_rad_astronomy  = unpacked["0100"]["vectors_angle"]

#     testing_vector_angles = [
#     np.array(vectors_angle_100U_0A_rad_astronomy),  
#     np.array(vectors_angle_90U_10A_rad_astronomy),
#     np.array(vectors_angle_80U_20A_rad_astronomy), 
#     np.array(vectors_angle_70U_30A_rad_astronomy), 
#     np.array(vectors_angle_60U_40A_rad_astronomy), 
#     np.array(vectors_angle_50U_50A_rad_astronomy), 
#     np.array(vectors_angle_40U_60A_rad_astronomy), 
#     np.array(vectors_angle_30U_70A_rad_astronomy), 
#     np.array(vectors_angle_20U_80A_rad_astronomy),
#     np.array(vectors_angle_10U_90A_rad_astronomy),
#     np.array(vectors_angle_0U_100A_rad_astronomy)
#     ]

#     return vector_data_plotting_grid, vector_data_list_100U_to_100A, testing_vector_angles





# def ratio_model_band47(StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
#                        StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
#                        ny, nx,
#                        vector_length_pix_const, 
#                        POLI_mJy, POLI_err_mJy,
#                        PA_err_deg,
#                        band, 
#                        step = None):
#     if step is None:
#         if band == 'band 4':
#             step = constants.step_band4
#         elif band == 'band 5':
#              step = constants.step_band5
#         elif band == 'band 7':
#              step = constants.step_band7
#         else:
#             raise ValueError("Unsupported band. Only Band 4, Band 5 and Band 7 are currently supported.")
        

#     results = {}

#     for ratio1, ratio2 in testing_ratios:
#         PA_grid, StokesQ_grid, StokesU_grid, vectors_data, vectors_angle = mix_StokesQU_and_generate_vectors_band47(
#             ratio1, ratio2,
#             StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
#             StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
#             ny, nx,
#             POLI_mJy, POLI_err_mJy,
#             PA_err_deg,
#             band, step)

#         key = generate_ratio_key(ratio1, ratio2)

#         results[f"PA_grid_{key}"]       = PA_grid
#         results[f"StokesQ_grid_{key}"]  = StokesQ_grid
#         results[f"StokesU_grid_{key}"]  = StokesU_grid
#         results[f"vectors_data_{key}"]  = vectors_data
#         results[f"vectors_angle_{key}"] = vectors_angle

#     unpacked = unpack_results(results, testing_ratios)
    
#        # Create clearly named vector data variables
#     PA_grid_100U_0A_rad_sky  = unpacked["1000"]["PA_grid"]
#     PA_grid_90U_10A_rad_sky  = unpacked["9010"]["PA_grid"]
#     PA_grid_80U_20A_rad_sky  = unpacked["8020"]["PA_grid"]
#     PA_grid_70U_30A_rad_sky  = unpacked["7030"]["PA_grid"]
#     PA_grid_60U_40A_rad_sky  = unpacked["6040"]["PA_grid"]
#     PA_grid_50U_50A_rad_sky  = unpacked["5050"]["PA_grid"]
#     PA_grid_40U_60A_rad_sky  = unpacked["4060"]["PA_grid"]
#     PA_grid_30U_70A_rad_sky  = unpacked["3070"]["PA_grid"]
#     PA_grid_20U_80A_rad_sky  = unpacked["2080"]["PA_grid"]
#     PA_grid_10U_90A_rad_sky  = unpacked["1090"]["PA_grid"]
#     PA_grid_0U_100A_rad_sky  = unpacked["0100"]["PA_grid"]

#     PA_mixed_grid_ratio_rad_sky = [
#         PA_grid_100U_0A_rad_sky, 
#         PA_grid_90U_10A_rad_sky, 
#         PA_grid_80U_20A_rad_sky, 
#         PA_grid_70U_30A_rad_sky, 
#         PA_grid_60U_40A_rad_sky, 
#         PA_grid_50U_50A_rad_sky, 
#         PA_grid_40U_60A_rad_sky,
#         PA_grid_30U_70A_rad_sky,
#         PA_grid_20U_80A_rad_sky,
#         PA_grid_10U_90A_rad_sky,
#         PA_grid_0U_100A_rad_sky,
#     ]
    
    
    

#     # Create clearly named vector data variables
#     vectors_data_100U_0A_cartesian  = unpacked["1000"]["vectors_data"]
#     vectors_data_90U_10A_cartesian  = unpacked["9010"]["vectors_data"]
#     vectors_data_80U_20A_cartesian  = unpacked["8020"]["vectors_data"]
#     vectors_data_70U_30A_cartesian  = unpacked["7030"]["vectors_data"]
#     vectors_data_60U_40A_cartesian  = unpacked["6040"]["vectors_data"]
#     vectors_data_50U_50A_cartesian  = unpacked["5050"]["vectors_data"]
#     vectors_data_40U_60A_cartesian  = unpacked["4060"]["vectors_data"]
#     vectors_data_30U_70A_cartesian  = unpacked["3070"]["vectors_data"]
#     vectors_data_20U_80A_cartesian  = unpacked["2080"]["vectors_data"]
#     vectors_data_10U_90A_cartesian  = unpacked["1090"]["vectors_data"]
#     vectors_data_0U_100A_cartesian  = unpacked["0100"]["vectors_data"]

#     # Save in list for return
#     vector_data_plotting_grid = [
#         vectors_data_100U_0A_cartesian, vectors_data_0U_100A_cartesian,
#         vectors_data_90U_10A_cartesian, vectors_data_10U_90A_cartesian,
#         vectors_data_80U_20A_cartesian, vectors_data_20U_80A_cartesian,
#         vectors_data_70U_30A_cartesian, vectors_data_30U_70A_cartesian,
#         vectors_data_60U_40A_cartesian, vectors_data_40U_60A_cartesian,
#         vectors_data_50U_50A_cartesian, vectors_data_50U_50A_cartesian  # repeated for consistency
#     ]

#     vector_data_list_100U_to_100A = [
#         vectors_data_100U_0A_cartesian, 
#         vectors_data_90U_10A_cartesian, 
#         vectors_data_80U_20A_cartesian, 
#         vectors_data_70U_30A_cartesian, 
#         vectors_data_60U_40A_cartesian, 
#         vectors_data_50U_50A_cartesian, 
#         vectors_data_40U_60A_cartesian,
#         vectors_data_30U_70A_cartesian,
#         vectors_data_20U_80A_cartesian,
#         vectors_data_10U_90A_cartesian,
#         vectors_data_0U_100A_cartesian,
#     ]
    
    
#     # Create clearly named vector data variables
#     vectors_angle_100U_0A_rad_astronomy  = unpacked["1000"]["vectors_angle"]
#     vectors_angle_90U_10A_rad_astronomy  = unpacked["9010"]["vectors_angle"]
#     vectors_angle_80U_20A_rad_astronomy  = unpacked["8020"]["vectors_angle"]
#     vectors_angle_70U_30A_rad_astronomy  = unpacked["7030"]["vectors_angle"]
#     vectors_angle_60U_40A_rad_astronomy  = unpacked["6040"]["vectors_angle"]
#     vectors_angle_50U_50A_rad_astronomy  = unpacked["5050"]["vectors_angle"]
#     vectors_angle_40U_60A_rad_astronomy  = unpacked["4060"]["vectors_angle"]
#     vectors_angle_30U_70A_rad_astronomy  = unpacked["3070"]["vectors_angle"]
#     vectors_angle_20U_80A_rad_astronomy  = unpacked["2080"]["vectors_angle"]
#     vectors_angle_10U_90A_rad_astronomy  = unpacked["1090"]["vectors_angle"]
#     vectors_angle_0U_100A_rad_astronomy  = unpacked["0100"]["vectors_angle"]

#     testing_vector_angles = [
#     np.array(vectors_angle_100U_0A_rad_astronomy),  
#     np.array(vectors_angle_90U_10A_rad_astronomy),
#     np.array(vectors_angle_80U_20A_rad_astronomy), 
#     np.array(vectors_angle_70U_30A_rad_astronomy), 
#     np.array(vectors_angle_60U_40A_rad_astronomy), 
#     np.array(vectors_angle_50U_50A_rad_astronomy), 
#     np.array(vectors_angle_40U_60A_rad_astronomy), 
#     np.array(vectors_angle_30U_70A_rad_astronomy), 
#     np.array(vectors_angle_20U_80A_rad_astronomy),
#     np.array(vectors_angle_10U_90A_rad_astronomy),
#     np.array(vectors_angle_0U_100A_rad_astronomy)
#     ]

#     return PA_mixed_grid_ratio_rad_sky, vector_data_plotting_grid, vector_data_list_100U_to_100A, testing_vector_angles
    

    
    

    
def ratio_model(StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
                StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
                ny, nx,
                vector_length_pix_const, 
                StokesI_mJy, StokesI_err_mJy, 
                POLI_mJy, POLI_err_mJy,
                PA_err_deg,
                band, 
                step = None):

    step_defaults = {
        'Band 4': constants.step_band4,
        'Band 4 nterms2': constants.step_band4_nterms2,
        'Band 4 nterms2 robust -1': constants.step_band4_nterms2_robust_minus1,
        'Band 4 nterms2 smooth': constants.step_band4_nterms2_smooth,
        'Band 4 nterms2 smooth B6': constants.step_band4_nterms2_smooth_B6,
        'Band 4 nterms2 smooth B6 B7': constants.step_band4_nterms2_smooth_B6_B7,

        'Band 5': constants.step_band5,
        'Band 5 v0': constants.step_band5_v0,
        'Band 5 robust -2': constants.step_band5_robust_minus2,
        'Band 5 robust -1': constants.step_band5_robust_minus1,

        'Band 6': constants.step_band6,
        'Band 6 smooth': constants.step_band6_smooth,
        'Band 6 smooth B7': constants.step_band6_smooth_B7,

        'Band 7 nterms2': constants.step_band7_nterms2,
        'Band 7 nterms2 smooth': constants.step_band7_nterms2_smooth,
        'Band 7 nterms2 smooth B6': constants.step_band7_nterms2_smooth_B6,
    }

    if step is None:
        if band not in step_defaults:
            raise ValueError(
                "Unsupported band. ")

        step = step_defaults[band]

        if band.startswith('Band 5'):
            print("All Band 5 using same step")
        

    results = {}

    for ratio1, ratio2 in testing_ratios:
        PA_grid, StokesQ_grid, StokesU_grid, vectors_data, vectors_angle, _ = mix_StokesQU_and_generate_vectors(
            ratio1, ratio2,
            StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
            StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
            ny, nx,
            StokesI_mJy, StokesI_err_mJy, 
            POLI_mJy, POLI_err_mJy,
            PA_err_deg,
            band, step)

        key = generate_ratio_key(ratio1, ratio2)

        results[f"PA_grid_{key}"]       = PA_grid
        results[f"StokesQ_grid_{key}"]  = StokesQ_grid
        results[f"StokesU_grid_{key}"]  = StokesU_grid
        results[f"vectors_data_{key}"]  = vectors_data
        results[f"vectors_angle_{key}"] = vectors_angle

    unpacked = unpack_results(results, testing_ratios)
    
       # Create clearly named vector data variables
    PA_grid_100U_0A_rad_sky  = unpacked["1000"]["PA_grid"]
    PA_grid_90U_10A_rad_sky  = unpacked["9010"]["PA_grid"]
    PA_grid_80U_20A_rad_sky  = unpacked["8020"]["PA_grid"]
    PA_grid_70U_30A_rad_sky  = unpacked["7030"]["PA_grid"]
    PA_grid_60U_40A_rad_sky  = unpacked["6040"]["PA_grid"]
    PA_grid_50U_50A_rad_sky  = unpacked["5050"]["PA_grid"]
    PA_grid_40U_60A_rad_sky  = unpacked["4060"]["PA_grid"]
    PA_grid_30U_70A_rad_sky  = unpacked["3070"]["PA_grid"]
    PA_grid_20U_80A_rad_sky  = unpacked["2080"]["PA_grid"]
    PA_grid_10U_90A_rad_sky  = unpacked["1090"]["PA_grid"]
    PA_grid_0U_100A_rad_sky  = unpacked["0100"]["PA_grid"]

    PA_mixed_grid_ratio_rad_sky = [
        PA_grid_100U_0A_rad_sky, 
        PA_grid_90U_10A_rad_sky, 
        PA_grid_80U_20A_rad_sky, 
        PA_grid_70U_30A_rad_sky, 
        PA_grid_60U_40A_rad_sky, 
        PA_grid_50U_50A_rad_sky, 
        PA_grid_40U_60A_rad_sky,
        PA_grid_30U_70A_rad_sky,
        PA_grid_20U_80A_rad_sky,
        PA_grid_10U_90A_rad_sky,
        PA_grid_0U_100A_rad_sky,
    ]
    
    
    

    # Create clearly named vector data variables
    vectors_data_100U_0A_cartesian  = unpacked["1000"]["vectors_data"]
    vectors_data_90U_10A_cartesian  = unpacked["9010"]["vectors_data"]
    vectors_data_80U_20A_cartesian  = unpacked["8020"]["vectors_data"]
    vectors_data_70U_30A_cartesian  = unpacked["7030"]["vectors_data"]
    vectors_data_60U_40A_cartesian  = unpacked["6040"]["vectors_data"]
    vectors_data_50U_50A_cartesian  = unpacked["5050"]["vectors_data"]
    vectors_data_40U_60A_cartesian  = unpacked["4060"]["vectors_data"]
    vectors_data_30U_70A_cartesian  = unpacked["3070"]["vectors_data"]
    vectors_data_20U_80A_cartesian  = unpacked["2080"]["vectors_data"]
    vectors_data_10U_90A_cartesian  = unpacked["1090"]["vectors_data"]
    vectors_data_0U_100A_cartesian  = unpacked["0100"]["vectors_data"]

    # Save in list for return
    vector_data_plotting_grid = [
        vectors_data_100U_0A_cartesian, vectors_data_0U_100A_cartesian,
        vectors_data_90U_10A_cartesian, vectors_data_10U_90A_cartesian,
        vectors_data_80U_20A_cartesian, vectors_data_20U_80A_cartesian,
        vectors_data_70U_30A_cartesian, vectors_data_30U_70A_cartesian,
        vectors_data_60U_40A_cartesian, vectors_data_40U_60A_cartesian,
        vectors_data_50U_50A_cartesian, vectors_data_50U_50A_cartesian  # repeated for consistency
    ]

    vector_data_list_100U_to_100A = [
        vectors_data_100U_0A_cartesian, 
        vectors_data_90U_10A_cartesian, 
        vectors_data_80U_20A_cartesian, 
        vectors_data_70U_30A_cartesian, 
        vectors_data_60U_40A_cartesian, 
        vectors_data_50U_50A_cartesian, 
        vectors_data_40U_60A_cartesian,
        vectors_data_30U_70A_cartesian,
        vectors_data_20U_80A_cartesian,
        vectors_data_10U_90A_cartesian,
        vectors_data_0U_100A_cartesian,
    ]
    
    
    # Create clearly named vector data variables
    vectors_angle_100U_0A_rad_astronomy  = unpacked["1000"]["vectors_angle"]
    vectors_angle_90U_10A_rad_astronomy  = unpacked["9010"]["vectors_angle"]
    vectors_angle_80U_20A_rad_astronomy  = unpacked["8020"]["vectors_angle"]
    vectors_angle_70U_30A_rad_astronomy  = unpacked["7030"]["vectors_angle"]
    vectors_angle_60U_40A_rad_astronomy  = unpacked["6040"]["vectors_angle"]
    vectors_angle_50U_50A_rad_astronomy  = unpacked["5050"]["vectors_angle"]
    vectors_angle_40U_60A_rad_astronomy  = unpacked["4060"]["vectors_angle"]
    vectors_angle_30U_70A_rad_astronomy  = unpacked["3070"]["vectors_angle"]
    vectors_angle_20U_80A_rad_astronomy  = unpacked["2080"]["vectors_angle"]
    vectors_angle_10U_90A_rad_astronomy  = unpacked["1090"]["vectors_angle"]
    vectors_angle_0U_100A_rad_astronomy  = unpacked["0100"]["vectors_angle"]

    testing_vector_angles = [
    np.array(vectors_angle_100U_0A_rad_astronomy),  
    np.array(vectors_angle_90U_10A_rad_astronomy),
    np.array(vectors_angle_80U_20A_rad_astronomy), 
    np.array(vectors_angle_70U_30A_rad_astronomy), 
    np.array(vectors_angle_60U_40A_rad_astronomy), 
    np.array(vectors_angle_50U_50A_rad_astronomy), 
    np.array(vectors_angle_40U_60A_rad_astronomy), 
    np.array(vectors_angle_30U_70A_rad_astronomy), 
    np.array(vectors_angle_20U_80A_rad_astronomy),
    np.array(vectors_angle_10U_90A_rad_astronomy),
    np.array(vectors_angle_0U_100A_rad_astronomy)
    ]

    return PA_mixed_grid_ratio_rad_sky, vector_data_plotting_grid, vector_data_list_100U_to_100A, testing_vector_angles
    

    
    
    


def find_best_fit_ratio_model(expected_angles, observed_angle_list, vector_PA_errors, print_results=True):
    """
    Identifies the best-fitting ratio model for a set of observed vector angles
    by minimizing the Chi-squared statistic compared to expected angles.

    This function is specific to ratio-based models (e.g., U/A mixtures like 100U/0A to 0U/100A).

    Parameters:
    - expected_angles (array-like): The expected vector angles (e.g., from a pure or reference model).
    - observed_angle_list (list of array-like): List of observed angle arrays from ratio-based models.
    - print_results (bool): If True, prints the chi-squared values for all ratio models and the best fit.

    Returns:
    - int: Index of the ratio model with the lowest Chi-squared value.
    """

    chi_squared_values = []

    for observed_angles in observed_angle_list:
        chi_squared = calculate_chi_squared_reduced(observed_angles, expected_angles, vector_PA_errors, 1)
#         chi_squared = calculate_chi_squared_v2(observed_angles, expected_angles)
        chi_squared_values.append(chi_squared)
       
    print('the length of chi_squared_values is :', len(chi_squared_values))
    min_index = chi_squared_values.index(min(chi_squared_values))

    if print_results:
        labels = [
            "100 U 0 A", "90 U 10 A", "80 U 20 A", "70 U 30 A", "60 U 40 A",
            "50 U 50 A", "40 U 60 A", "30 U 70 A", "20 U 80 A", "10 U 90 A", "0 U 100 A"
        ]

        print("Chi-squared values for ratio models:")
        for label, chi in zip(labels, chi_squared_values):
            print(f'  {label}: χ² = {chi * 10:.2f}')

        print(f'\nBest-fit ratio model: {labels[min_index]} (χ² = {chi_squared_values[min_index] * 10:.2f})')

    return min_index







def PlotRatioExample(band, StokesI_wcs, vector_data_plotting_grid,
                     xmin, xmax, ymin, ymax,
                    fig_x_size = 15, 
                    fig_y_size = 6.5,
                    xy_axes = True):
    
#     POLI_mJy, 
#                     StokesI_wcs, soft_colormap_v2, 
#                     
# reference_length_pix, reference_length_AU,
#                     BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
#                     max_length_pix, reference_fraction,
#                     vector_data_actual_cartesian, vector_data_plotting_grid, 
#                     fs_scale = 0.5):
    
    # Define the figure and subplots
    fig, axes = plt.subplots(1, 3, figsize=(fig_x_size, fig_y_size), 
                             #constrained_layout=True, 
                             subplot_kw={'projection': StokesI_wcs},
#                              gridspec_kw={'wspace': -1},
                             sharey = True)


    axes[0].text(0.05, 0.8, f"{constants.lambda_mm[band]} mm", transform=axes[0].transAxes, fontsize=constants.text_fs)
    
    plotting_indicies = [0, 10, 1]
    N = [100, 50, 0]
    titles = ['100Uni + 0 Azi', '50Uni + 50Azi', '0Uni + 100Azi']
    # Loop through subplots
    for i, ax in enumerate(axes.flat):


    

        # Add vector plots
        idx = plotting_indicies[i]
        #ax.set_title(ratio_grid_plot_titles[idx], fontsize=20)
        
        #ax.text(0.05, 0.90, f"$N$ = {N[i]}", transform=ax.transAxes, fontsize=constants.text_fs)
        ax.text(0.05, 0.90, f"{titles[i]}", transform=ax.transAxes, fontsize=constants.text_fs)
        
        
        for row in vector_data_plotting_grid[idx]:
            ax.plot([row[0], row[1]], [row[2], row[3]], color='black', lw = 3, label = 'Model')


        # ---------------------------------------------------------------------
        # Set y-axis labels and ticks 
        # ---------------------------------------------------------------------
        if xy_axes:
            if i == 0:
                ax.set_ylabel('Declination (ICRS)', fontsize=constants.axis_label_fs)
                ax.tick_params(axis="y", which="both", left=True, labelleft=True)
            else:
                ax.tick_params(axis="y", which="both", left=True, labelleft=False)
            # ---------------------------------------------------------------------
            # Set x-axis labels and ticks
            # ---------------------------------------------------------------------
            ax.set_xlabel('Right Ascension (ICRS)', fontsize=constants.axis_label_fs)
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
            # ---------------------------------------------------------------------
        else:
            ax.tick_params(axis="y", which="both", left=True, labelleft=False)
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            ax.set_xlabel(' ', fontsize = 0)
            


        # Set x and y limits
        # ---------------------------------------------------------------------
        diff = constants.plot_zoom_pixels - 10
        ax.set_xlim(xmin + diff, xmax - diff)
        ax.set_ylim(ymin + diff, ymax - diff)
        print(rf'({xmin}, {xmax}), ({ymin}, {ymax})')
        # ---------------------------------------------------------------------
        
        
        ax.minorticks_on()
        ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=constants.axis_num_fs - 2)
        ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=constants.axis_num_fs - 2)

        # Minor axis ticks
        # --------------------------------------------------------
#         ra = ax.coords['ra']

#         # --- Minor ticks ---
#         ra.display_minor_ticks(True)
#         ra.set_ticks_position(('b', 't'))
#         ra.tick_params(which='minor', length=4)
        # --------------------------------------------------------
    
    
    return fig, axes
# ----------------------------------------------------------------------------------------------------------------------
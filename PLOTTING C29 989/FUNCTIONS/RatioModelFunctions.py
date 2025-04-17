import sys

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


def ratio_model_band6(StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
                      StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
                      ny, nx,
                      step, vector_length_pix_const, 
                      StokesI_data_2d_mJy, StokesIerr_data_2d_mJy,
                      calculated_polarized_intensity, PolarizedIntensity_err_data_2d_mJy,
                      PolarizationAngle_err_data_2d_deg):

    results = {}

    # Loop through the ratios
    for ratio1, ratio2 in testing_ratios:

        # Run the function with the current ratio
        PA_grid, StokesQ_grid, StokesU_grid, vectors_data, vectors_angle = mix_StokesQU_and_generate_vectors(
            ratio1, ratio2, 
            StokesQ_grid_100Uniform, 
            StokesU_grid_100Uniform,
            StokesQ_grid_100Azimuthal, 
            StokesU_grid_100Azimuthal,
            ny, nx,
            step, vector_length_pix_const, 
            StokesI_data_2d_mJy, StokesIerr_data_2d_mJy,
            calculated_polarized_intensity, PolarizedIntensity_err_data_2d_mJy,
            PolarizationAngle_err_data_2d_deg)

        # Create a key for the dictionary based on the ratio
        ratio_str = f"{int(ratio1 * 100)}{int(ratio2 * 100)}"  # e.g., "5050" for (0.5, 0.5)

        # Save the results in the dictionary
        results[f"PA_grid_{ratio_str}"] = PA_grid
        results[f"StokesQ_grid_{ratio_str}"] = StokesQ_grid
        results[f"StokesU_grid_{ratio_str}"] = StokesU_grid
        results[f"vectors_data_{ratio_str}"] = vectors_data
        results[f"vectors_angle_{ratio_str}"] = vectors_angle
    
    
    
        PA_grid_100U_0A, StokesQ_grid_100U_0A, StokesU_grid_100U_0A, vectors_data_100U_0A_cartesian, vectors_angle_100U_0A_rad_astronomy = (results["PA_grid_1000"], results["StokesQ_grid_1000"], results["StokesU_grid_1000"], results["vectors_data_1000"], results["vectors_angle_1000"])
    PA_grid_90U_10A, StokesQ_grid_90U_10A, StokesU_grid_90U_10A, vectors_data_90U_10A_cartesian, vectors_angle_90U_10A_rad_astronomy = (results["PA_grid_9010"], results["StokesQ_grid_9010"], results["StokesU_grid_9010"], results["vectors_data_9010"], results["vectors_angle_9010"])
    PA_grid_80U_20A, StokesQ_grid_80U_20A, StokesU_grid_80U_20A, vectors_data_80U_20A_cartesian, vectors_angle_80U_20A_rad_astronomy = (results["PA_grid_8020"], results["StokesQ_grid_8020"], results["StokesU_grid_8020"], results["vectors_data_8020"], results["vectors_angle_8020"])
    PA_grid_70U_30A, StokesQ_grid_70U_30A, StokesU_grid_703U_0A, vectors_data_70U_30A_cartesian, vectors_angle_70U_30A_rad_astronomy = (results["PA_grid_7030"], results["StokesQ_grid_7030"], results["StokesU_grid_7030"], results["vectors_data_7030"], results["vectors_angle_7030"])
    PA_grid_60U_40A, StokesQ_grid_60U_40A, StokesU_grid_60U_40A, vectors_data_60U_40A_cartesian, vectors_angle_60U_40A_rad_astronomy = (results["PA_grid_6040"], results["StokesQ_grid_6040"], results["StokesU_grid_6040"], results["vectors_data_6040"], results["vectors_angle_6040"])
    PA_grid_50U_50A, StokesQ_grid_50U_50A, StokesU_grid_50U_50A, vectors_data_50U_50A_cartesian, vectors_angle_50U_50A_rad_astronomy = (results["PA_grid_5050"], results["StokesQ_grid_5050"], results["StokesU_grid_5050"], results["vectors_data_5050"], results["vectors_angle_5050"])
    PA_grid_40U_60A, StokesQ_grid_40U_60A, StokesU_grid_40U_60A, vectors_data_40U_60A_cartesian, vectors_angle_40U_60A_rad_astronomy = (results["PA_grid_4060"], results["StokesQ_grid_4060"], results["StokesU_grid_4060"], results["vectors_data_4060"], results["vectors_angle_4060"])
    PA_grid_30U_70A, StokesQ_grid_30U_70A, StokesU_grid_30U_70A, vectors_data_30U_70A_cartesian, vectors_angle_30U_70A_rad_astronomy = (results["PA_grid_3070"], results["StokesQ_grid_3070"], results["StokesU_grid_3070"], results["vectors_data_3070"], results["vectors_angle_3070"])
    PA_grid_20U_80A, StokesQ_grid_20U_80A, StokesU_grid_20U_80A, vectors_data_20U_80A_cartesian, vectors_angle_20U_80A_rad_astronomy = (results["PA_grid_2080"], results["StokesQ_grid_2080"], results["StokesU_grid_2080"], results["vectors_data_2080"], results["vectors_angle_2080"])
    PA_grid_10U_90A, StokesQ_grid_10U_90A, StokesU_grid_10U_90A, vectors_data_10U_90A_cartesian, vectors_angle_10U_90A_rad_astronomy = (results["PA_grid_1090"], results["StokesQ_grid_1090"], results["StokesU_grid_1090"], results["vectors_data_1090"], results["vectors_angle_1090"])
    PA_grid_0U_100A, StokesQ_grid_0U_100A, StokesU_grid_0U_100A, vectors_data_0U_100A_cartesian, vectors_angle_0U_100A_rad_astronomy = (results["PA_grid_0100"], results["StokesQ_grid_0100"], results["StokesU_grid_0100"], results["vectors_data_0100"], results["vectors_angle_0100"])




    vector_data_list = [
    vectors_data_100U_0A_cartesian, vectors_data_0U_100A_cartesian,
    vectors_data_90U_10A_cartesian, vectors_data_10U_90A_cartesian,
    vectors_data_80U_20A_cartesian, vectors_data_20U_80A_cartesian,
    vectors_data_70U_30A_cartesian, vectors_data_30U_70A_cartesian,
    vectors_data_60U_40A_cartesian, vectors_data_40U_60A_cartesian,
    vectors_data_50U_50A_cartesian, vectors_data_50U_50A_cartesian
]

    
    

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

    return vector_data_list, testing_vector_angles





def ratio_model_band4(StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
                      StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
                      ny, nx,
                      step, vector_length_pix_const, 
                      POLI_mJy, POLI_err_mJy,
                      PA_err_deg):

    results = {}

    for ratio1, ratio2 in testing_ratios:
        PA_grid, StokesQ_grid, StokesU_grid, vectors_data, vectors_angle = mix_StokesQU_and_generate_vectors_band4(
            ratio1, ratio2,
            StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
            StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
            ny, nx,
            POLI_mJy, POLI_err_mJy,
            PA_err_deg)

        key = generate_ratio_key(ratio1, ratio2)

        results[f"PA_grid_{key}"]       = PA_grid
        results[f"StokesQ_grid_{key}"]  = StokesQ_grid
        results[f"StokesU_grid_{key}"]  = StokesU_grid
        results[f"vectors_data_{key}"]  = vectors_data
        results[f"vectors_angle_{key}"] = vectors_angle

    unpacked = unpack_results(results, testing_ratios)

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

    return vector_data_plotting_grid, vector_data_list_100U_to_100A, testing_vector_angles
    

    
    



def find_best_fit_ratio_model(expected_angles, observed_angle_list, print_results=True):
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
        chi_squared = calculate_chi_squared_v2(observed_angles, expected_angles)
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
            print(f'  {label}: χ² = {chi:.3f}')

        print(f'\nBest-fit ratio model: {labels[min_index]} (χ² = {chi_squared_values[min_index]:.3f})')

    return min_index
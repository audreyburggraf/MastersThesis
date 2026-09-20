from matplotlib.colors import LinearSegmentedColormap

# Hex codes for the colors
hex_colors = [
    "#fefefe",  # White
    "#fff4f8",  # Very Light Pink
    "#fde5f3",  # Light Pink
    "#ffd8fb",  # Light Lavender
    "#f8ccff",  # Light Purple
    "#e3bfff",  # Light Lilac
    "#c5b3fd",  # Soft Purple
    "#a8aaff",  # Soft Blue
    "#9cbaff",  # Light Blue
    "#8ddafe",  # Sky Blue
    "#80f7ff",  # Pastel Cyan
    "#74fedd",  # Light Turquoise
    "#66ffaa",  # Lime Green
    "#59fe73",  # Light Green
    "#6aff4e",  # Bright Green
    "#8efe41",  # Lime Yellow
    "#ccff32",  # Bright Lime Yellow
    "#fef324",  # Soft Yellow
    "#ffaf1b",  # Soft Orange
    "#ff6a0f",  # Bright Orange
    "#ff1900",  # Red
]

# Convert hex codes to RGB values (normalized to 0-1)
rgb_colors = [
    tuple(int(color[i:i+2], 16) / 255 for i in (1, 3, 5)) for color in hex_colors
]

# Create the colormap using the pastel colors
soft_colormap_v2 = LinearSegmentedColormap.from_list("SoftColorMapV2", rgb_colors, N=1000)






# No Red
# ------------------------------------------------------------------------------------------------
# Hex codes for the colors
hex_colors_no_red = [
    "#fefefe",  # White
    "#fff4f8",  # Very Light Pink
    "#fde5f3",  # Light Pink
    "#ffd8fb",  # Light Lavender
    "#f8ccff",  # Light Purple
    "#e3bfff",  # Light Lilac
    "#c5b3fd",  # Soft Purple
    "#a8aaff",  # Soft Blue
    "#9cbaff",  # Light Blue
    "#8ddafe",  # Sky Blue
    "#80f7ff",  # Pastel Cyan
    "#74fedd",  # Light Turquoise
    "#66ffaa",  # Lime Green
    "#59fe73",  # Light Green
    "#6aff4e",  # Bright Green
    "#8efe41",  # Lime Yellow
    "#ccff32",  # Bright Lime Yellow
    "#fef324",  # Soft Yellow
    "#ffaf1b",  # Soft Orange
    "#ff6a0f",  # Bright Orange
#     "#ff1900",  # Red
]

# Convert hex codes to RGB values (normalized to 0-1)
rgb_colors_no_red = [
    tuple(int(color[i:i+2], 16) / 255 for i in (1, 3, 5)) for color in hex_colors_no_red
]

# Create the colormap using the pastel colors
soft_colormap_no_red = LinearSegmentedColormap.from_list("SoftColorMapNoRed", rgb_colors_no_red, N=1000)


# ------------------------------------------------------------------------------------------------


# Hex codes for the colors (first 10 from the list)
hex_colors_10 = [
    "#fefefe",  # White
    # "#fff4f8",  # Very Light Pink
    "#fde5f3",  # Light Pink
    # "#ffd8fb",  # Light Lavender
    "#f8ccff",  # Light Purple
    # "#e3bfff",  # Light Lilac
    "#c5b3fd",  # Soft Purple
    # "#a8aaff",  # Soft Blue
    "#9cbaff",  # Light Blue
    # "#8ddafe",  # Sky Blue
    "#80f7ff",  # Pastel Cyan
    # "#74fedd",  # Light Turquoise
    "#66ffaa",  # Lime Green
    # "#59fe73",  # Light Green
    "#6aff4e",  # Bright Green
    # "#8efe41",  # Lime Yellow
    "#ccff32",  # Bright Lime Yellow
    # "#fef324",  # Soft Yellow
    "#ffaf1b",  # Soft Orange
    # "#ff6a0f",  # Bright Orange
    "#ff1900",  # Red
]

# Convert hex codes to RGB values (normalized to 0-1)
rgb_colors_10 = [
    tuple(int(color[i:i+2], 16) / 255 for i in (1, 3, 5)) for color in hex_colors_10
]

# Create the colormap using the 10 selected colors
soft_colormap_10 = LinearSegmentedColormap.from_list("SoftColorMap10", rgb_colors_10, N=10)

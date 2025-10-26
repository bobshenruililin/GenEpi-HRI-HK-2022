# --- Input Data ---

if (!exists("unique_hospitals")) {
  stop("unique_hospitals not found. Please load the data first.")
}

if (!exists("unique_hospitals_official")) {
  stop("unique_hospitals_official not found. Please load the data first.")
}

if (!exists("unique_hospital_wards")) {
  stop("unique_hospital_wards not found. Please load the data first.")
}

if (!exists("unique_hospital_wards_official")) {
  stop("unique_hospital_wards_official not found. Please load the data first.")
}

# --- Color Assignment Logic ---

# Define the color for non-official / excluded wards
excluded_color <- "#808080" # A neutral grey

official_palette <- c(
  "#003366", # Very Dark Blue (Navy)
  "#006400", # Dark Green
  "#8B0000", # Dark Red
  "#4B0082", # Indigo (Dark Purple/Blue)
  "#8B4513", # Saddle Brown
  "#008080", # Teal
  "#B8860B", # Dark Goldenrod (Mustard/Dark Yellow)
  "#556B2F", # Dark Olive Green
  "#4682B4", # Steel Blue (Muted Medium Blue)
  "#800080", # Purple
  "#2E8B57", # Sea Green (Darker Bluish-Green)
  "#A0522D", # Sienna (Reddish Brown)
  "#CC5500", # Burnt Orange
  "#483D8B", # Dark Slate Blue (More Purple than Navy)
  "#008B8B", # Dark Cyan (Greener Teal)
  "#DC143C", # Crimson (Strong Red)
  "#2F4F4F", # Dark Slate Gray
  "#800000"  # Maroon
)
n_palette <- length(official_palette)
stopifnot(n_palette >= length(unique_hospital_wards_official))

# Initialize the named vector to store colors for each ward
ward_colors <- character(length(unique_hospital_wards))
names(ward_colors) <- unique_hospital_wards

current_index <- c()

# Assign colors
for (ward in sort(unique_hospital_wards)) {
  # Check if the ward is in the official list
  is_official <- ward %in% unique_hospital_wards_official

  if (is_official) {
    # Get the current color index for this hospital, initialize if first time
    if (is.null(current_index)) {
      current_index <- 1
    }

    # Assign color from the palette, cycling if necessary
    color_to_assign <- official_palette[((current_index - 1) %% n_palette) + 1]
    ward_colors[ward] <- color_to_assign

    # Increment the index for the next ward in this hospital
    current_index <- current_index + 1

  } else {
    # Ward is not in the official list, assign the excluded color
    ward_colors[ward] <- excluded_color
  }
}

# --- Output ---

# Display the generated colors
print(ward_colors)

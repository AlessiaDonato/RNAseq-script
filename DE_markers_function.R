library(Seurat)

ht <- readRDS("HT_cDC2_0207.Rds")

# Set the active identity class for the Seurat object to "Cluster_def"
# This determines how cells are grouped for downstream analysis
Idents(ht) <- "Cluster_def"
# Display the number of cells in each cluster based on the active identity class
table(Idents(ht))

# Identify differentially expressed genes using the MAST (Model-based Analysis 
# of Single-cell Transcriptomics) test

# Subset the Seurat object to include only cells with tumor and healthy condition
ht_T <- subset(ht, subset = cond == "T")
ht_H <- subset(ht, subset = cond == "H")

# Alternative way to split Seurat Object into one with lists
# ht_split <- SplitObject(ht, split.by = "cond")

# 1. Between cell populations of same tissue in tumor condition (looped) #######

# Define the nested folder structure
base_folder <- "DE_markers_csv"
subfolder_1 <- file.path(base_folder, "tissue_populations_csv")
subfolder_2 <- file.path(subfolder_1, "tumor_tissue_pop_csv")

# Create the nested folders if they don't already exist
dir.create(subfolder_2, recursive = TRUE, showWarnings = FALSE)

# Define the comparisons for intra-tumor cells
intra_tumor_comparisons <- list(
  c("DC2_Ifitm1", "DC3_Ifitm1")
)

# Define the comparisons for circulating cells
circulating_comparisons <- list(
  c("preDC2", "proDC3"),
  c("preDC2", "DC3 circ"),
  c("proDC3", "DC3 circ")
)

# Define the comparisons for lung cells
lung_comparisons <- list(
  c("DC2B_Ltb", "DC2B_Epcam"),
  c("DC2B_Ltb", "DC3 tiss"),
  c("DC2B_Ltb", "DC3 DCSIGN"),
  c("DC2B_Epcam", "DC3 tiss"),
  c("DC2B_Epcam", "DC3 DCSIGN"),
  c("DC3 tiss", "DC3 DCSIGN")
)

# Define a function to run FindMarkers and save the results
run_find_markers_and_save <- function(comparisons, subfolder) {
  final_folder <- file.path(subfolder_2, subfolder)
  dir.create(final_folder, recursive = TRUE, showWarnings = FALSE)
  
  for (comparison in comparisons) {
    ident1 <- comparison[1]
    ident2 <- comparison[2]
    marker_name <- paste0("T_", gsub(" ", "", ident1), "_vs_T_", gsub(" ", "", ident2), ".markers")
    
    # Perform the differential expression analysis using FindMarkers
    markers <- FindMarkers(ht_T,
                           ident.1 = ident1,
                           ident.2 = ident2,
                           test.use = "MAST")
    
    # Assign the result to a variable with the constructed name
    assign(marker_name, markers, envir = .GlobalEnv)
    
    # Construct the file name by appending '.csv' to the marker name and prepending the folder path
    file_name <- file.path(final_folder, paste0(marker_name, ".csv"))
    
    # Write the data frame to a CSV file using write.csv2
    write.csv2(markers, file = file_name, row.names = TRUE)
  }
}

# Run the function for each set of comparisons
run_find_markers_and_save(intra_tumor_comparisons, "intra_tumor_cells_csv")
run_find_markers_and_save(circulating_comparisons, "circulating_cells_csv")
run_find_markers_and_save(lung_comparisons, "lung_cells_csv")

# 1. Between cell populations of same tissue in healthy condition (looped) #####

# Define the nested folder structure
base_folder <- "DE_markers_csv"
subfolder_1 <- file.path(base_folder, "tissue_populations_csv")
subfolder_2 <- file.path(subfolder_1, "healthy_tissue_pop_csv")

# Create the nested folders if they don't already exist
dir.create(subfolder_2, recursive = TRUE, showWarnings = FALSE)

# Define the comparisons for circulating cells
circulating_comparisons <- list(
  c("preDC2", "proDC3"),
  c("preDC2", "DC3 circ"),
  c("proDC3", "DC3 circ")
)

# Define the comparisons for lung cells
lung_comparisons <- list(
  c("DC2B_Ltb", "DC2B_Epcam"),
  c("DC2B_Ltb", "DC3 tiss"),
  c("DC2B_Ltb", "DC3 DCSIGN"),
  c("DC2B_Epcam", "DC3 tiss"),
  c("DC2B_Epcam", "DC3 DCSIGN"),
  c("DC3 tiss", "DC3 DCSIGN")
)

# Define a function to run FindMarkers and save the results
run_find_markers_and_save <- function(comparisons, subfolder) {
  final_folder <- file.path(subfolder_2, subfolder)
  dir.create(final_folder, recursive = TRUE, showWarnings = FALSE)
  
  for (comparison in comparisons) {
    ident1 <- comparison[1]
    ident2 <- comparison[2]
    marker_name <- paste0("H_", gsub(" ", "", ident1), "_vs_H_", gsub(" ", "", ident2), ".markers")
    
    # Perform the differential expression analysis using FindMarkers
    markers <- FindMarkers(ht_H,
                           ident.1 = ident1,
                           ident.2 = ident2,
                           test.use = "MAST")
    
    # Assign the result to a variable with the constructed name
    assign(marker_name, markers, envir = .GlobalEnv)
    
    # Construct the file name by appending '.csv' to the marker name and prepending the folder path
    file_name <- file.path(final_folder, paste0(marker_name, ".csv"))
    
    # Write the data frame to a CSV file using write.csv2
    write.csv2(markers, file = file_name, row.names = TRUE)
  }
}

# Run the function for each set of comparisons
run_find_markers_and_save(circulating_comparisons, "circulating_cells_csv")
run_find_markers_and_save(lung_comparisons, "lung_cells_csv")

# 2. Between differentiated cells in tumor and healthy condition (looped) ######

# Define the nested folder structure
base_folder <- "DE_markers_csv"
subfolder_1 <- file.path(base_folder, "differentiated_cells_csv")
subfolder_healthy <- file.path(subfolder_1, "healthy_cells_csv")
subfolder_tumor <- file.path(subfolder_1, "tumor_cells_csv")

# Create the nested folders if they don't already exist
dir.create(subfolder_healthy, recursive = TRUE, showWarnings = FALSE)
dir.create(subfolder_tumor, recursive = TRUE, showWarnings = FALSE)

# Define the comparisons for healthy cells
healthy_comparisons <- list(
  c("DC3 circ", "DC3 tiss"),
  c("DC3 circ", "DC3 DCSIGN"),
  c("preDC2", "DC2A"),
  c("preDC2", "DC2B_Epcam"),
  c("preDC2", "DC2B_Ltb"),
  c("DC2B_Epcam", "DC_CD207"),
  c("DC2B_Epcam", "DC2B_Gm42418")
)

# Define the comparisons for tumor cells
tumor_comparisons <- list(
  c("DC3 circ", "DC3 tiss"),
  c("DC3 circ", "DC3 DCSIGN"),
  c("DC3 tiss", "DC3_Ifitm1"),
  c("DC3 DCSIGN", "DC3_Ifitm1"),
  c("preDC2", "DC2A"),
  c("preDC2", "DC2B_Epcam"),
  c("preDC2", "DC2B_Ltb"),
  c("DC2B_Ltb", "DC2_Ifitm1"),
  c("DC2B_Epcam", "DC_CD207"),
  c("DC2B_Epcam", "DC2B_Gm42418"),
  c("DC2B_Epcam", "DC2_Ifitm1")
)

# Define a function to run FindMarkers and save the results
run_find_markers_and_save <- function(comparisons, ht_object, subfolder) {
  for (comparison in comparisons) {
    ident1 <- comparison[1]
    ident2 <- comparison[2]
    marker_name <- paste0(gsub(" ", "", ident1), "_vs_", gsub(" ", "", ident2), ".markers")
    
    # Perform the differential expression analysis using FindMarkers
    markers <- FindMarkers(ht_object,
                           ident.1 = ident1,
                           ident.2 = ident2,
                           test.use = "MAST")
    
    # Assign the result to a variable with the constructed name
    assign(marker_name, markers, envir = .GlobalEnv)
    
    # Construct the file name by appending '.csv' to the marker name and prepending the folder path
    file_name <- file.path(subfolder, paste0(marker_name, ".csv"))
    
    # Write the data frame to a CSV file using write.csv2
    write.csv2(markers, file = file_name, row.names = TRUE)
  }
}

# Run the function for each set of comparisons
run_find_markers_and_save(healthy_comparisons, ht_H, subfolder_healthy)
run_find_markers_and_save(tumor_comparisons, ht_T, subfolder_tumor)

# 3. Between cell types in healthy and tumor condition (looped) ########

# Define the nested folder structure
base_folder <- "DE_markers_csv"
final_folder <- file.path(base_folder, "HvsT_csv")

# Create the nested folders if they don't already exist
dir.create(final_folder, recursive = TRUE, showWarnings = FALSE)

# Define the subset.ident values to iterate over
subset_idents <- c("preDC2", "proDC3", "DC3 circ", "DC2A", "DC3 DCSIGN", "DC3 tiss", 
                   "DC2B_Ltb", "DC2B_Epcam", "DC_CD207", "DC2B_Gm42418")

# Define a function to construct marker names
construct_marker_name <- function(subset_ident) {
  paste0("H_", gsub(" ", "", subset_ident), "_vs_T_", gsub(" ", "", subset_ident), ".markers")
}

# Loop through each subset.ident value and perform FindMarkers
for (subset_ident in subset_idents) {
  # Construct the marker name dynamically
  marker_name <- construct_marker_name(subset_ident)
  
  # Perform the differential expression analysis using FindMarkers
  markers <- FindMarkers(ht,
                         ident.1 = "H",
                         group.by = 'cond',
                         subset.ident = subset_ident,
                         test.use = "MAST")
  
  # Assign the result to a variable with the constructed name
  assign(marker_name, markers)
  
  # Construct the file name by appending '.csv' to the marker name and prepending the folder path
  file_name <- file.path(final_folder, paste0(marker_name, ".csv"))
  
  # Write the data frame to a CSV file using write.csv2
  write.csv2(markers, file = file_name, row.names = TRUE)
}




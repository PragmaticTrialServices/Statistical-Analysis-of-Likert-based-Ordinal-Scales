load_libraries <- function() {
  packages <- c("rje", "rms", "MASS", "winprob", "dplyr", "lmtest", "sandwich", "VGAM", "beepr")
  
  # Loop through the list of required packages
  for (pkg in packages) {
    # Install the package if it is not available
    if (!requireNamespace(pkg, quietly = TRUE)) {
      suppressMessages(install.packages(pkg, dependencies = TRUE))
    }
    # Load the package without showing startup messages
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  message("Libraries successfully loaded")
}

# Run the library-loading function
load_libraries()


# Open a connection to the requirements.txt file in write mode
file_conn <- file("requirements.txt", "w")

# Write the R version to the file
writeLines(paste("R version:", getRversion()), file_conn)

# Retrieve information about all loaded packages
loaded_packages <- sessionInfo()$otherPkgs

# Iterate over each loaded package and write its name and version to the file
for (pkg in loaded_packages) {
  writeLines(paste(pkg$Package, pkg$Version, sep="=="), file_conn)
}

# Close the file connection
close(file_conn)

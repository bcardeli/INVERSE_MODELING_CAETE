###                      INVERSE MODELING ALGORITHM                          ###
# ---------------------------------------------------------------------------- #
#      A new modeling approach to evaluate the effects of climate change       # 
#           on plant functional diversity in the Amazon rainforest             #
# ---------------------------------------------------------------------------- #

# This algorithm automates and integrates plant traits optimization (SLA, WD, G1) to 
# identify communities that sustain ecosystem processes and services  
# under climate change scenarios similar to stable climates.
# MODEL INTEGRATED: CAETÊ-DVM (Rius et al., 2023 - allom2light [branch] Cardeli et al., 2021)

#ATTENTION: Modify the paths to the folders according to where you are running: (1) ON SERVER or (2) MY MACHINE

#AUTHOR: Bárbara R. Cardeli (barbara.r.cardeli (at) gmail.com)

# Let's go!!!!! ;-)

# Load required libraries
library(jsonlite) # To create important config files 
library(rBayesianOptimization) # To Bayesian Optimization
library(ncdf4) # To read NetCDF files
library(fs) # For file and directory management
library(beepr) # For fun sounds :)

####### FUNCTIONS #######

# Constant for invalid score
INVALID_SCORE_PENALTY <- 1e6

# Function to generate values within the traits range [OPTIMIZATION]
generate_values <- function(npls, min_val, max_val) {
  if (min_val >= max_val) {
    warning("Minimum value is greater than or equal to the maximum value. Returning NA.")
    return(rep(NA, npls))
  }
  return(runif(npls, min_val, max_val))
}

# Define parameter bounds
bounds <- list(
  dwood_min = c(0.4, 0.5),
  dwood_max = c(0.9, 1.0),
  sla_min = c(0.006, 0.009),
  sla_max = c(0.040, 0.050),
  g1_min = c(0.1, 0.5),
  g1_max = c(19.0, 20.0)
)

objective_function_BO <- function(dwood_min, dwood_max, sla_min, sla_max, g1_min, g1_max) {
  
  # Generate values within the bounds
  dwood <- generate_values(npls, dwood_min, dwood_max)
  sla_var <- generate_values(npls, sla_min, sla_max)
  g1 <- generate_values(npls, g1_min, g1_max)
  
  # Check for NA values in any of the trait vectors
  if (any(is.na(c(dwood, sla_var, g1)))) {
    message("Invalid solution detected (NA values). Penalizing.")
    return(list(Score = INVALID_SCORE_PENALTY))  # Penalize invalid solutions
  }
  
  # Calculate the objective value
  obj_value <- sum(dwood) + sum(sla_var) + sum(g1)
  
  # Return as a list with the expected structure
  return(list(Score = as.numeric(obj_value)))
}

# Function to calculate the average of a variable in a NetCDF file [VERIFICATION]
calculate_mean <- function(nc_file, var_name) {
  nc <- nc_open(nc_file)
  var_data <- ncvar_get(nc, var_name)
  nc_close(nc)
  mean_value <- mean(var_data, na.rm = TRUE)
  return(mean_value)
}

# Function to calculate the sum and average of multiple variables in NetCDF files [VERIFICATION]
calculate_sum_and_mean <- function(directory, files, var_names) {
  var_list <- vector("list", length(var_names))
  
  for (i in seq_along(files)) {
    nc <- nc_open(file.path(directory, files[i]))
    var_list[[i]] <- ncvar_get(nc, var_names[i])
    nc_close(nc)
  }
  
  var_sum <- Reduce(`+`, var_list)
  var_mean <- mean(var_sum, na.rm = TRUE)
  
  return(var_mean)
}

# Function for performing calculations in an experiment (base or test) [VERIFICATION]
calculate_experiment <- function(directory, et_file, gpp_file, sum_files, var_names) {
  et_mean_value <- calculate_mean(file.path(directory, et_file), "et")
  gpp_mean_value <- calculate_mean(file.path(directory, gpp_file), "gpp")
  sum_mean_value <- calculate_sum_and_mean(directory, sum_files, var_names)
  
  return(list(et_mean_value = et_mean_value, gpp_mean_value = gpp_mean_value, sum_mean_value = sum_mean_value))
}

# Function to delete the output file in case of “NO MATCH”
delete_folder <- function(folder_path) {
  if (dir.exists(folder_path)) {
    unlink(folder_path, recursive = TRUE, force = TRUE)
    cat("Folder and all contents deleted successfully:", folder_path, "\n")
  } else {
    cat("Directory does not exist:", folder_path, "\n")
  }
}

####### PRE-PROCESSING #######

# Define base experiment results [VERIFICATION]
# (1) ON SERVER:
base_directory <- "/dmz/home/bcardeli/CAETE_INV_MODEL/RUN_BASE/nclim_base/nc_outputs"

base_et_file <- "evapm_20150101-20161231.nc4"
base_gpp_file <- "photo_20150101-20161231.nc4"
base_sum_files <- c("cleaf_20150101-20161231.nc4", "cawood_20150101-20161231.nc4", "cfroot_20150101-20161231.nc4")
base_var_names <- c("cleaf", "cawood", "cfroot")

# Calculate base experiment averages and sums
base_results <- calculate_experiment(base_directory, base_et_file, base_gpp_file, base_sum_files, base_var_names)

# Read the PLS number from the datafile
if (file.exists("npls.txt")) {
  npls_content <- readLines("npls.txt")
  if (length(npls_content) > 0) {
    npls <- as.integer(npls_content)
  } else {
    stop("The file npls.txt is empty.")
  }
} else {
  stop("Not possible to obtain PLS number.")
}

# Define a counter to track the number of successful matches
match_counter <- 0

# Function to copy the content of a folder
copy_directory <- function(from, to) {
  if (file.exists(from)) {
    dir.create(to, recursive = TRUE)
    file.copy(list.files(from, full.names = TRUE), to, recursive = TRUE)
  } else {
    stop(paste("Source directory does not exist:", from))
  }
}

iterations <- 200 # Loop number

# OPTIMIZATION LOOP
for (i in 1:iterations) {
  
  cat("Starting optimization...\n")
  
  # Bayesian Optimization
  opt_result_BO <- BayesianOptimization(
    FUN = objective_function_BO,
    bounds = bounds, #Parameters (traits) boundaries.
    init_points = 12, #number of random samples initially evaluated before starting the optimization.
    n_iter = 50, #Refers to the number of optimization iterations to balance exploration and refinement.
    acq = "ei", #Expected Improvement (EI): Searches for regions that have a high expectation of improving the solution.
    verbose = TRUE
  )
  
  # Check if the optimization result has the expected structure
  print(str(opt_result_BO))
  
  # Extract the best parameters
  best_params <- opt_result_BO$Best_Par
  
  # Ensure best_params is not NULL
  if (!is.null(best_params)) {
    cat("Optimized parameters:\n")
    print(best_params)
    
    # Saving optimized parameters to a JSON file
    params_to_save <- list(
      dwood_min = as.numeric(best_params["dwood_min"]),
      dwood_max = as.numeric(best_params["dwood_max"]),
      sla_min = as.numeric(best_params["sla_min"]),
      sla_max = as.numeric(best_params["sla_max"]),
      g1_min = as.numeric(best_params["g1_min"]),
      g1_max = as.numeric(best_params["g1_max"])
    )
    
    # Save parameters to JSON file
    
    # (1) ON SERVER:
    write_json(params_to_save, file.path("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/src/params.json"))
    
    if (file.size("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/src/params.json") == 0) {
      stop("params.json is empty. Cannot run the model.")
    }
  } else {
    stop("No optimized parameters found.")
  }
  
  # Call Python script to run the CAETÊ-DVM model
  
  cat("Calling the model...\n")
  system(paste("python3 /dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/src/model_driver.py"), wait = TRUE)
  cat("Model call completed.\n")
  
  # Read folder name from "run_name.txt" file
  run_name <-"/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/src/run_name.txt"
  result_folder <- readLines(run_name)
  
  # Debug: Show folder name read from file
  cat("Folder name, read from run_name.txt:", result_folder, "\n")
  
  ####### PROCESSING & VERIFICATION #######
  
  # Set the experiment directory
  test_directory <- file.path("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/outputs", result_folder, "nc_outputs")
  
  # Debug: Show experiment directory path
  cat("Experiment directory path:", test_directory, "\n")
  
  test_et_file <- "evapm_20150101-20161231.nc4"
  test_gpp_file <- "photo_20150101-20161231.nc4"
  test_sum_files <- c("cleaf_20150101-20161231.nc4", "cawood_20150101-20161231.nc4", "cfroot_20150101-20161231.nc4")
  test_var_names <- c("cleaf", "cawood", "cfroot")
  
  # Calculate experiment averages and sums
  test_results <- calculate_experiment(test_directory, test_et_file, test_gpp_file, test_sum_files, test_var_names)
  
  # Debug: Show the results to comparison
  cat("Comparing results:\n")
  cat("Base:\n")
  cat("Evap - mean:", base_results$et_mean_value, "\n")
  cat("GPP - mean:", base_results$gpp_mean_value, "\n")
  cat("Biomass - mean:", base_results$sum_mean_value, "\n")
  
  cat("Experiment:\n")
  cat("Evap - mean:", test_results$et_mean_value, "\n")
  cat("GPP - mean:", test_results$gpp_mean_value, "\n")
  cat("Biomass - mean:", test_results$sum_mean_value, "\n")
  
  # Comparing results
  if (test_results$et_mean_value >= base_results$et_mean_value &&
      test_results$gpp_mean_value >= base_results$gpp_mean_value &&
      test_results$sum_mean_value >= base_results$sum_mean_value) {
    match_counter <- match_counter + 1
    print("COMMUNITY FOUND - MATCH")
    beepr::beep(3)
    
    # Save optimized parameters to the new result folder
    
    # Output file path for optimization results
    output_file <- file.path("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/outputs", result_folder, "resultados_otimizacao.txt")
    
    # Range optimized
    traits_optim <- sprintf("Iteração %d: dwood_min = %.2f, dwood_max = %.2f, sla_min = %.4f, sla_max = %.4f, g1_min = %.2f, g1_max = %.2f\n",
                            i, best_params["dwood_min"], best_params["dwood_max"], best_params["sla_min"], best_params["sla_max"], best_params["g1_min"], best_params["g1_max"])
    
    # Save the optimized range
    write(traits_optim, file = output_file, append = TRUE)
    
    # Define a new folder for each successful match
    new_result_folder <- paste0("MATCH_", match_counter)
    
    # (1) ON SERVER:
    new_result_folder_path <- file.path("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/outputs", new_result_folder)
    
    # Copy the entire contents of the original output folder to the new folder
    
    #(1) ON SERVER:
    original_output_folder <- file.path("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/outputs", result_folder)
    
    copy_directory(original_output_folder, new_result_folder_path)
    
    # Check that the copy was successful (check that the “nc_outputs” folder exists)
    nc_outputs_path <- file.path(new_result_folder_path, "nc_outputs")
    if (dir.exists(nc_outputs_path)) {
      cat("Data copied successfully to", new_result_folder_path, "\n")
      
      # Delete original folder after successful copy
      delete_folder(original_output_folder)
      
      cat("Original output folder deleted:", original_output_folder, "\n")
    } else {
      stop("Error: 'nc_outputs' folder not found. The original folder will not be deleted.")
    }
    
    # Continue searching for more matches
    if (match_counter >= 30) {
      cat("30 COMMUNITIES FOUND. TERMINATING...\n")
      break #our work is done! ;-)
    }
    
  } else {
    print("NO MATCH!")
    print("Deleting the contents of the outputs folder and start it again...")
    
    # Deleting "NO MATCH!" outputs
    
    #(1) ON SERVER:
    output_base_path <- "/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/outputs"
    
    out_folder <- result_folder
    folder_path <- file.path(output_base_path, result_folder)
    print(folder_path)
    
    # Bye-bye...
    delete_folder(folder_path)
  }
  
  cat(sprintf("Iteração %d: dwood_min = %.2f, dwood_max = %.2f, sla_min = %.4f, sla_max = %.4f, g1_min = %.2f, g1_max = %.2f\n",
              i, best_params["dwood_min"], best_params["dwood_max"], best_params["sla_min"], best_params["sla_max"], best_params["g1_min"], best_params["g1_max"]))
  
  # The sound of a warning that we must begin again.
  system("echo -e '\a'")
}
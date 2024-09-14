# Load required libraries
library(jsonlite) #To create important config files 
library(ncdf4) #To read NetCDF files
library(fs) # For file and directory manage
library(beepr) #To do it funny - Some sounds :)

####### FUNCTIONS #######

# Function to generate values within the traits range [OPTMIZATION]
generate_values <- function(npls, min_val, max_val) {
  if (min_val >= max_val) {
    return(rep(NA, npls))
  }
  return(runif(npls, min_val, max_val))
}

# Objective function based on the parameters 'dwood', 'sla' and 'g1', 
# generating values for each of them within the minimum and maximum limits provided [OPTMIZATION]
objective_function <- function(params, npls) {
  dwood_min <- params[1]
  dwood_max <- params[2]
  sla_min <- params[3]
  sla_max <- params[4]
  g1_min <- params[5]
  g1_max <- params[6]
  
  dwood <- generate_values(npls, dwood_min, dwood_max)
  sla_var <- generate_values(npls, sla_min, sla_max)
  g1 <- generate_values(npls, g1_min, g1_max)
  
  if (any(is.na(c(dwood, sla_var, g1)))) {
    return(1e6)
  }
  
  return(sum(dwood) + sum(sla_var) + sum(g1))
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
# special case to biomass (sum of cleaf, cwoood and cfroot)
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

# In case of “NO MATCH” (see below) — Function to delete the output file
delete_folder <- function(folder_path) {
  if (dir.exists(folder_path)) {  # Check if the directory exists
    unlink(folder_path, recursive = TRUE, force = TRUE)
    cat("Pasta e todo o conteúdo deletados com sucesso:", folder_path, "\n")
  } else {
    cat("O diretório não existe:", folder_path, "\n")
  }
}

####### PRE-PROCESSING #######

# Define base experiment results [VERIFICATION]
base_directory <- "//home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/outputs/RUN_BASE/nc_outputs"
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

# Initial parameters dos traits (WD, SLA and G1, respectively - see pls_gen.py)
initial_params <- c(0.5, 0.9, 0.009, 0.040, 0.1, 19.0)
iterations <- 10 #Loop number

# OPTIMIZATION LOOP
for (i in 1:iterations) {
  
  cat("Inicializand optmizaion...\n")
  
  # Call Python script to run the CAETÊ model
  system(paste("python3 /home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/src/model_driver.py"), wait = TRUE)
  
  # Optimization to find the best distribution parameters
  opt_result <- optim(par = initial_params, 
                      fn = objective_function, 
                      npls = npls,
                      method = "L-BFGS-B",
                      lower = c(0.5, 0.5, 0.009, 0.009, 0.1, 0.1), 
                      upper = c(0.9, 0.9, 0.040, 0.040, 19.0, 19.0))
  
  # Checking if the optimization was successful
  if (opt_result$convergence != 0) {
    cat("Optimization didn't converge in iteration", i, "\n")
    next
  }
  
  # Saving optimized parameters to a JSON file
  params_to_save <- list(
    dwood_min = opt_result$par[1],
    dwood_max = opt_result$par[2],
    sla_min = opt_result$par[3],
    sla_max = opt_result$par[4],
    g1_min = opt_result$par[5],
    g1_max = opt_result$par[6]
  )
  
  # Define the current iteration's results folder name
  run_name <- "/home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/src/run_name.txt"
  
  # Read folder name from "run_name.txt" file
  result_folder <- readLines(run_name)
    
  # Debug: Show folder name read from file
  cat("Folder name, read from run_name.txt:", result_folder, "\n")
    
  # Debug: View the automatically generated folder name
  cat("Folder name generated for iteration:", result_folder, "\n")
  
  # Debug: View the full path to the results folder
  cat("Full path to the results folder:", file.path("/home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/outputs", result_folder), "\n")
  
  write_json(params_to_save, file.path("/home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/outputs", result_folder, "params.json"))
  
  ####### PROCESSING & VERIFICATION #######
  
  # Set the experiment directory (output file)
  test_directory <- file.path("/home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/outputs", result_folder, "nc_outputs")
  
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
  
  cat("Experiment (inv_model:\n")
  cat("Evap - mean:", test_results$et_mean_value, "\n")
  cat("GPP - mean:", test_results$gpp_mean_value, "\n")
  cat("Biomass - mean:", test_results$sum_mean_value, "\n")
  
  # Comparing results
  if (test_results$et_mean_value >= base_results$et_mean_value &&
      test_results$gpp_mean_value >= base_results$gpp_mean_value &&
      test_results$sum_mean_value >= base_results$sum_mean_value) {
    print("COMMUNITY FINDED - MATCH")
    beepr::beep(3) #tan-tan-tan-tan beat the drums!!!!!
    
    # Output file patch
    output_file <- file.path("/home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/outputs", result_folder, "resultados_otimizacao.txt")
    
    # Range optimazed
    traits_optim <- sprintf("Iteração %d: dwood_min = %.2f, dwood_max = %.2f, sla_min = %.4f, sla_max = %.4f, g1_min = %.2f, g1_max = %.2f\n",
                        i, opt_result$par[1], opt_result$par[2], opt_result$par[3], opt_result$par[4], opt_result$par[5], opt_result$par[6])
    
    # Save the optimazed range
    write(traits_optim, file = output_file, append = TRUE)
    
    break  #our work is done! ;-)
    
  } else {
    print("NO MATCH!")
    print("Deleting the contents of the outputs folder and start it again...")
    
    # Deleting "NO MATCH!" outputs
    output_base_path <- "/home/barbara/Documentos/CAETE-DVM_Branch/CAETE-DVM/outputs"
    out_folder <- result_folder
    folder_path <- file.path(output_base_path, result_folder)
    print(folder_path)
    
    # Bye-bye...
    delete_folder(folder_path)
  }
  
  cat(sprintf("Iteração %d: dwood_min = %.2f, dwood_max = %.2f, sla_min = %.4f, sla_max = %.4f, g1_min = %.2f, g1_max = %.2f\n",
              i, opt_result$par[1], opt_result$par[2], opt_result$par[3], opt_result$par[4], opt_result$par[5], opt_result$par[6]))
  
  # The sound of a warning that we must begin again.
  system("echo -e '\a'")
}

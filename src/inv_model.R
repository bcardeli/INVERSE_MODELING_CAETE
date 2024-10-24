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
  dwood_min <- min(params[1], params[2])
  dwood_max <- max(params[1], params[2])
  sla_min <- min(params[3], params[4])
  sla_max <- max(params[3], params[4])
  g1_min <- min(params[5], params[6])
  g1_max <- max(params[5], params[6])
  
  # Penalizar relações inconsistentes
  if (dwood_min >= dwood_max || sla_min >= sla_max || g1_min >= g1_max) {
    return(1e6)  # Penalidade alta para evitar estas configurações
  }
  
  dwood <- generate_values(npls, dwood_min, dwood_max)
  sla_var <- generate_values(npls, sla_min, sla_max)
  g1 <- generate_values(npls, g1_min, g1_max)
  
  if (any(is.na(c(dwood, sla_var, g1)))) {
    return(1e6)
  }
  
  # Adding a randomness factor to force more radical changes
  random_penalty <- runif(1, -2, 2) # range for higher variability
  
  return(-1 * (sum(dwood) + sum(sla_var) + sum(g1)) + random_penalty) #maximization (*-1)
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
    cat("Folder and all content successfully deleted:", folder_path, "\n")
  } else {
    cat("The directory doesn't exist:", folder_path, "\n")
  }
}

####### PRE-PROCESSING #######

# Define base experiment results [VERIFICATION]
# If you are running this algorithm for the first time, you should ask for this data. 

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

# Initial parameters dos traits (WD, SLA and G1, respectively - see pls_gen.py)
initial_params <- c(0.5, 0.9, 0.009, 0.040, 0.1, 19.0)
iterations <- 500 #Loop number

# Definir limites inferiores e superiores para cada parâmetro
lower_bounds <- c(0.5, 0.5,    # dwood_min, dwood_max
                  0.009, 0.009, # sla_min, sla_max
                  0.1, 0.1)     # g1_min, g1_max

upper_bounds <- c(0.9, 0.9,    # dwood_min, dwood_max
                  0.040, 0.040,  # sla_min, sla_max
                  19.0, 19.0)  # g1_min, g1_max


# Define a counter to track the number of successful matches
match_counter <- 0

# Função para copiar o conteúdo de uma pasta
copy_directory <- function(from, to) {
  if (file.exists(from)) {
    dir.create(to, recursive = TRUE)
    file.copy(list.files(from, full.names = TRUE), to, recursive = TRUE)
  } else {
    stop(paste("Source directory does not exist:", from))
  }
}

# OPTIMIZATION LOOP
for (i in 1:iterations) {
  
  cat("Starting optmization...\n")
  
  # Definir ponto de partida aleatório dentro dos limites para cada parâmetro em cada iteração
  initial_params <- runif(length(lower_bounds), lower_bounds, upper_bounds)
  
  # Optimization to find the best distribution parameters
  opt_result <- optim(par = initial_params, 
                      fn = objective_function, 
                      npls = npls,
                      method = "L-BFGS-B",
                      lower = lower_bounds, 
                      upper = upper_bounds)
  
  # Checking if the optimization was successful
  if (opt_result$convergence != 0) {
    cat("Optimization didn't converge in iteration", i, "\n")
    next
  }
  
  # Ajustar os resultados para garantir consistência antes de salvar
  dwood_min <- min(opt_result$par[1], opt_result$par[2])
  dwood_max <- max(opt_result$par[1], opt_result$par[2])
  sla_min <- min(opt_result$par[3], opt_result$par[4])
  sla_max <- max(opt_result$par[3], opt_result$par[4])
  g1_min <- min(opt_result$par[5], opt_result$par[6])
  g1_max <- max(opt_result$par[5], opt_result$par[6])
  
  # Salvar os parâmetros consistentes
  params_to_save <- list(
    dwood_min = dwood_min,
    dwood_max = dwood_max,
    sla_min = sla_min,
    sla_max = sla_max,
    g1_min = g1_min,
    g1_max = g1_max
  )
  
  # (1) ON SERVER:
  write_json(params_to_save, file.path("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/src/params.json"))
  
  # Call Python script to run the CAETÊ-DVM model
  
  # (1) ON SERVER:
  system(paste("python3 /dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/src/model_driver.py"), wait = TRUE)
  
  # Define the current iteration's results folder name
  
  # (1) ON SERVER:
  run_name <-"/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/src/run_name.txt"
  
  # Read folder name from "run_name.txt" file
  result_folder <- readLines(run_name)
    
  # Debug: Show folder name read from file
  cat("Folder name, read from run_name.txt:", result_folder, "\n")
    
  # Debug: View the automatically generated folder name
  cat("Folder name generated for iteration:", result_folder, "\n")
  
  # Debug: View the full path to the results folder
  
  ####### PROCESSING & VERIFICATION #######
  
  # Set the experiment directory (output file)
  
  # (1) ON SERVER:
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
  
  cat("Experiment (inv_model:\n")
  cat("Evap - mean:", test_results$et_mean_value, "\n")
  cat("GPP - mean:", test_results$gpp_mean_value, "\n")
  cat("Biomass - mean:", test_results$sum_mean_value, "\n")
  
  # Comparing results
  if (test_results$et_mean_value >= base_results$et_mean_value &&
      test_results$gpp_mean_value >= base_results$gpp_mean_value &&
      test_results$sum_mean_value >= base_results$sum_mean_value) {
    match_counter <- match_counter + 1  # Increment the match counter
    print("COMMUNITY FINDED - MATCH")
    beepr::beep(3) #tan-tan-tan-tan beat the drums!!!!!
    
    # Save optimized parameters to the new result folder
    
    # Output file path for optimization results
    output_file <- file.path("/dmz/home/bcardeli/CAETE_INV_MODEL/INVERSE_MODELING_CAETE/outputs", result_folder, "resultados_otimizacao.txt")
    
    # Range optimized
    traits_optim <- sprintf("Iteração %d: dwood_min = %.2f, dwood_max = %.2f, sla_min = %.4f, sla_max = %.4f, g1_min = %.2f, g1_max = %.2f\n",
                            i, opt_result$par[1], opt_result$par[2], opt_result$par[3], opt_result$par[4], opt_result$par[5], opt_result$par[6])
    
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
    
    # Verificar se a cópia foi bem-sucedida (verificar se a pasta "nc_outputs" existe)
    nc_outputs_path <- file.path(new_result_folder_path, "nc_outputs")
    if (dir.exists(nc_outputs_path)) {
      cat("Data copied successfully to", new_result_folder_path, "\n")
      
      # Excluir a pasta original após cópia bem-sucedida
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
              i, opt_result$par[1], opt_result$par[2], opt_result$par[3], opt_result$par[4], opt_result$par[5], opt_result$par[6]))
  
  # The sound of a warning that we must begin again.
  system("echo -e '\a'")
}

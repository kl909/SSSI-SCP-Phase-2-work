library(googledrive)
# authenticase (only need once per session)
drive_auth()
library(here)

# ------------------------------
# find the shared google folder
my_folder <- drive_find(
  pattern = "CritchlowEtAlDistributions10k",
  type = "folder",
  shared_drive = NULL # looks in shared with me area
)

# store folder ID
folder_id <- my_folder$id[1]

# list all files inside
all_files <- drive_ls(as_id(folder_id))

# files to exclude
exclude_list <-c("UK examples_july2016.zip", "uncertainty_maps.zip")

# filter out the excluded files
files_to_process <- all_files[!(all_files$name %in% exclude_list), ]

# empty list for file names
all_species_files <- list()

for (i in 1:nrow(files_to_process)) {
  current_zip_name <- files_to_process$name[i]
  cat("Downloading:", current_zip_name, "...\n")
  
  temp_zip <- tempfile(fileext = ".zip")
  drive_download(as_id(files_to_process$id[i]), path = temp_zip, overwrite = TRUE)
  
  # Just get the RAW names and store them
  all_species_files[[i]] <- data.frame(
    zip_source = current_zip_name,
    raw_name = unzip(temp_zip, list = TRUE)$Name
  )
  
  unlink(temp_zip)
}

# Combine everything into one big Master Table
master_species_list <- do.call(rbind, all_species_files)

# View the result
View(master_species_list)

# clean names
working_names <- master_species_list$raw_name
# 1. strip the folder name
clean_step <- basename(working_names)
# 2. remove everything after the species name
clean_step <- sub("_.*$", "", clean_step)
# 3. remove file extentions
clean_step <- sub("\\.[a-zA-Z0-9]+$", "", clean_step)

# put into dataframe
species_df <- data.frame(species = clean_step, stringsAsFactors = FALSE)

# remove "GTiffs", "Thumbs", and any empty strings
species_df <- species_df[!(species_df$species %in% c("GTiffs", "Thumbs", "")), , drop = FALSE]

# and remove rows where file extention is left over from species
species_df <- species_df[!grepl("db|txt|ini", species_df$species, ignore.case = TRUE), , drop = FALSE]

# get only unique species names
final_species_df <- data.frame(species = unique(species_df$species))

View(final_species_df)

write.csv(final_species_df, "critchlow_species.csv", row.names = FALSE)

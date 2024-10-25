library(data.table)
library(stringr)

data <- fread("C:/Users/mrooney/Desktop/WFO_Backbone/classification.csv", sep = "\t", encoding = "UTF-8")

str(data)

library(data.table)
library(stringr)

# Clean scientificName
data[, Taxon := str_trim(scientificName)]
data[, Taxon := str_replace_all(Taxon, '^"|"$', '')]

# Count of taxa by rank
rank_counts <- data[, .N, by = taxonRank]
setorder(rank_counts, -N)
print(rank_counts)

# Count of taxa by family
family_counts <- data[, .N, by = family]
setorder(family_counts, -N)
print(head(family_counts, 10))  # Top 10 families

# Count of genera
genus_count <- data[, uniqueN(genus)]
print(paste("Number of unique genera:", genus_count))

# Taxa with uncertain status
uncertain_taxa <- data[taxonomicStatus == "Unchecked", .N]
print(paste("Number of taxa with uncertain status:", uncertain_taxa))

# Most recent modifications
data[, modified := as.Date(modified)]
recent_mods <- data[order(-modified)][1:10, .(Taxon, modified)]
print(recent_mods)

# Count of missing values per column
missing_counts <- data[, lapply(.SD, function(x) sum(is.na(x) | x == "")), .SDcols = names(data)]
print(missing_counts)

# Subset of accepted names
accepted_names <- data[taxonomicStatus == "Accepted"]

# Summary of taxa by major group
major_group_summary <- data[, .(count = .N, 
                                genera = uniqueN(genus), 
                                families = uniqueN(family)), 
                            by = majorGroup]
print(major_group_summary)


# Export summaries
fwrite(rank_counts, "rank_counts.csv")
fwrite(family_counts, "family_counts.csv")
fwrite(major_group_summary, "major_group_summary.csv")

library(stringi)

# Identify rows with non-ASCII characters
non_ascii_rows <- which(stri_enc_isascii(data$Taxon) == FALSE)

# Print some examples of non-ASCII taxon names
print(data$Taxon[non_ascii_rows[1:10]])

library(data.table)
library(stringi)

# Function to replace hybrid symbols
replace_hybrid_symbol <- function(x) {
  # Replace common hybrid symbol variations with "×"
  x <- gsub("\\s[xX]\\s", " × ", x)  # Replace " x " or " X " with " × "
  x <- gsub("^[xX]\\s", "× ", x)     # Replace "x " or "X " at the start with "× "
  x <- gsub("\\s[xX]$", " ×", x)     # Replace " x" or " X" at the end with " ×"
  
  # Replace Unicode hybrid symbol (U+00D7) if present
  x <- gsub("\u00D7", "×", x)
  
  return(x)
}

# Apply the function to the Taxon column
data[, Taxon := replace_hybrid_symbol(Taxon)]

# Check for changes
hybrid_rows <- data[grep("×", Taxon)]
print(head(hybrid_rows$Taxon, 10))

library(data.table)
library(stringi)

clean_taxon <- function(x) {
  # First, preserve hybrid symbols by replacing them temporarily
  x <- gsub("×", "HYBRIDPLACEHOLDER", x)
  
  # Transliterate non-ASCII characters, preserving common diacritical marks
  x <- iconv(x, to = "ASCII//TRANSLIT", sub = "")
  
  # Keep only allowed characters: letters, numbers, spaces, hyphens, periods, parentheses, and specific diacritical marks
  x <- gsub("[^a-zA-Z0-9 .æœëöäüïÿ()'-]", "", x)
  
  # Ensure "subsp." and similar abbreviations keep their periods and following punctuation
  x <- gsub("subsp\\. ?", "subsp. ", x)
  x <- gsub("var\\. ?", "var. ", x)
  x <- gsub("f\\. ?", "f. ", x)
  x <- gsub("nothovar\\. ?", "nothovar. ", x)
  x <- gsub("nothosubsp\\. ?", "nothosubsp. ", x)
  
  # Ensure there's only one space between words and no leading/trailing spaces
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  
  # Now, replace the placeholder back with the hybrid symbol
  x <- gsub("HYBRIDPLACEHOLDER", "×", x)
  
  # Standardize spacing around hybrid symbol
  x <- gsub("\\s*×\\s*", " × ", x)
  x <- gsub("^× ", "×", x)  # Remove space after × if it's at the start of the string
  
  return(x)
}

# Apply the cleaning function to the Taxon column
data[, clean_Taxon := clean_taxon(Taxon)]

# Check the results
print(data[1:10, .(Taxon, clean_Taxon)])

# Count how many taxon names were modified
modified_count <- sum(data$Taxon != data$clean_Taxon)
print(paste("Number of modified taxon names:", modified_count))

# Show some examples of modifications
modified_examples <- data[Taxon != clean_Taxon, .(Taxon, clean_Taxon)][1:10]
print(modified_examples)

# Check for hybrid names
hybrid_examples <- data[grep("×", Taxon), .(Taxon, clean_Taxon)][1:20]
print("Sample of hybrid names:")
print(hybrid_examples)

# Check for any changes in subspecies, varieties, or forms
subsp_changes <- data[grep("subsp\\.|var\\.|f\\.|nothovar\\.|nothosubsp\\.", Taxon), .(Taxon, clean_Taxon)]
print("Sample of changes in names with subspecies, varieties, or forms:")
print(head(subsp_changes[Taxon != clean_Taxon], 10))

#Replace Taxon column with the cleaned column
data[, Taxon := clean_Taxon]
data[, clean_Taxon := NULL]  # Remove the temporary column

#Export clean csv
write.csv(data, "C:/Users/mrooney/Desktop/WFO_Backbone/cleaned_data.csv", row.names = FALSE)

library(data.table)

# Load the original data
tblTaxon_DurationSummary <- fread("C:/Users/mrooney/Desktop/WFO_Backbone/tblTaxon_DurationSummary.csv")
tblTaxon_HabitSummary <- fread("C:/Users/mrooney/Desktop/WFO_Backbone/tblTaxon_HabitSummary.csv")

# Inspect original Taxon column
print(unique(tblTaxon_DurationSummary$Taxon))
print(unique(tblTaxon_HabitSummary$Taxon))

# Load necessary library
library(data.table)

install.packages("readr")
install.packages("stringi")

library(stringi)
library(readr)
# Detect encoding
file_path <- "C:/Users/mrooney/Desktop/WFO_Backbone/tblTaxon_HabitSummary.csv"
encoding_info <- stri_enc_detect(file_path)
print(encoding_info)

# Specify the detected encoding (ISO-8859-1)
tblTaxon_HabitSummary <- read_csv(file_path, locale = locale(encoding = "ISO-8859-1"))

# Inspect unique values
print("Unique Values:")
print(unique(tblTaxon_HabitSummary$Taxon)[1:5])  # Print only first 5 unique values

# Load necessary libraries
library(readr)
library(data.table)

# Define file paths
habit_file_path <- "C:/Users/mrooney/Desktop/WFO_Backbone/tblTaxon_HabitSummary.csv"
duration_file_path <- "C:/Users/mrooney/Desktop/WFO_Backbone/tblTaxon_DurationSummary.csv"

# Load the datasets with specified encoding
tblTaxon_HabitSummary <- read_csv(habit_file_path, locale = locale(encoding = "ISO-8859-1"))
tblTaxon_DurationSummary <- read_csv(duration_file_path, locale = locale(encoding = "ISO-8859-1"))

# Define a cleaning function
clean_taxon <- function(x) {
  # Standardize spacing around hybrid symbol
  x <- gsub("\\s*×\\s*", " × ", x)  # Ensure consistent spacing
  
  # Keep valid characters: letters, numbers, spaces, hyphens, periods, parentheses, and specific diacritical marks
  x <- gsub("[^a-zA-Z0-9 .×æœëöäüïÿ()'-]", "", x)
  
  # Trim whitespace
  x <- trimws(x)
  
  return(x)
}

# Apply cleaning function to Taxon column in Habit Summary
tblTaxon_HabitSummary$Taxon <- clean_taxon(tblTaxon_HabitSummary$Taxon)

# Apply cleaning function to Taxon column in Duration Summary
tblTaxon_DurationSummary$Taxon <- clean_taxon(tblTaxon_DurationSummary$Taxon)

# Remove rows with NA in Taxon column for Habit Summary
tblTaxon_HabitSummary <- tblTaxon_HabitSummary[!is.na(tblTaxon_HabitSummary$Taxon), ]

# Remove rows with NA in Taxon column for Duration Summary
tblTaxon_DurationSummary <- tblTaxon_DurationSummary[!is.na(tblTaxon_DurationSummary$Taxon), ]

# Inspect cleaned values from both datasets
print("Cleaned Unique Values from Habit Summary:")
print(unique(tblTaxon_HabitSummary$Taxon)[1:5])

print("Cleaned Unique Values from Duration Summary:")
print(unique(tblTaxon_DurationSummary$Taxon)[1:5])

# Export cleaned data to CSV files
write_csv(tblTaxon_HabitSummary, "C:/Users/mrooney/Desktop/WFO_Backbone/cleaned_tblTaxon_HabitSummary.csv")
write_csv(tblTaxon_DurationSummary, "C:/Users/mrooney/Desktop/WFO_Backbone/cleaned_tblTaxon_DurationSummary.csv")

# Load necessary libraries
library(dplyr)
library(readr)

# Assuming classification data is already loaded as 'data'

# Load Habit and Duration datasets if not already done
# tblTaxon_HabitSummary <- read_csv("C:/Users/mrooney/Desktop/WFO_Backbone/cleaned_tblTaxon_HabitSummary.csv")
# tblTaxon_DurationSummary <- read_csv("C:/Users/mrooney/Desktop/WFO_Backbone/cleaned_tblTaxon_DurationSummary.csv")

# Use distinct() to remove duplicates based on Taxon
tblTaxon_HabitSummary_distinct <- tblTaxon_HabitSummary %>% distinct(Taxon, Habit)
tblTaxon_DurationSummary_distinct <- tblTaxon_DurationSummary %>% distinct(Taxon, Duration)

# Merge Habit data with multiple handling
classification_with_habit <- data %>%
  left_join(tblTaxon_HabitSummary_distinct, by = "Taxon", relationship = "many-to-many")

# Merge Duration data with multiple handling
classification_with_habit_duration <- classification_with_habit %>%
  left_join(tblTaxon_DurationSummary_distinct, by = "Taxon", relationship = "many-to-many")

# Inspect the resulting dataset
print("Resulting Dataset with Habit and Duration:")
head(classification_with_habit_duration)

# Optionally, save the resulting dataset to a CSV file
write_csv(classification_with_habit_duration, "C:/Users/mrooney/Desktop/WFO_Backbone/classification_with_habit_duration.csv")

# Load necessary library
library(dplyr)

# Assuming 'data' is your classification dataset
# Load Habit and Duration datasets if not already done
# tblTaxon_HabitSummary <- read_csv("C:/Users/mrooney/Desktop/WFO_Backbone/cleaned_tblTaxon_HabitSummary.csv")
# tblTaxon_DurationSummary <- read_csv("C:/Users/mrooney/Desktop/WFO_Backbone/cleaned_tblTaxon_DurationSummary.csv")

# Total number of unique taxa in the classification dataset
total_taxa <- n_distinct(data$Taxon)

# Count taxa with habit assigned
taxa_with_habit <- data %>%
  filter(Taxon %in% tblTaxon_HabitSummary$Taxon) %>%
  summarise(Count = n_distinct(Taxon)) %>%
  pull(Count)

# Count taxa with duration assigned
taxa_with_duration <- data %>%
  filter(Taxon %in% tblTaxon_DurationSummary$Taxon) %>%
  summarise(Count = n_distinct(Taxon)) %>%
  pull(Count)

# Count taxa with both assigned
taxa_with_both <- data %>%
  filter(Taxon %in% tblTaxon_HabitSummary$Taxon & Taxon %in% tblTaxon_DurationSummary$Taxon) %>%
  summarise(Count = n_distinct(Taxon)) %>%
  pull(Count)

# Count taxa with only habit assigned
taxa_only_habit <- data %>%
  filter(Taxon %in% tblTaxon_HabitSummary$Taxon & !(Taxon %in% tblTaxon_DurationSummary$Taxon)) %>%
  summarise(Count = n_distinct(Taxon)) %>%
  pull(Count)

# Count taxa with only duration assigned
taxa_only_duration <- data %>%
  filter(Taxon %in% tblTaxon_DurationSummary$Taxon & !(Taxon %in% tblTaxon_HabitSummary$Taxon)) %>%
  summarise(Count = n_distinct(Taxon)) %>%
  pull(Count)

# Count taxa with neither assigned
taxa_neither <- total_taxa - (taxa_with_habit + taxa_with_duration - taxa_with_both)

# Print summary results
summary_results <- tibble(
  TotalTaxa = total_taxa,
  TaxaWithHabit = taxa_with_habit,
  TaxaWithDuration = taxa_with_duration,
  TaxaWithBoth = taxa_with_both,
  TaxaOnlyHabit = taxa_only_habit,
  TaxaOnlyDuration = taxa_only_duration,
  TaxaNeither = taxa_neither
)

print("Summary of Taxa Assignments:")
print(summary_results)

# Print unique values for Habit and Duration
unique_habits <- unique(tblTaxon_HabitSummary$Habit)
unique_durations <- unique(tblTaxon_DurationSummary$Duration)

cat("\nUnique Values for Habit:\n")
print(unique_habits)

cat("\nUnique Values for Duration:\n")
print(unique_durations)
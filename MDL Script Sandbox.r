#Functions
#mdl                Calculates MDL for single xanalyte
#mdl_general        Calculates MDL for
#mdlbfunction       Calculates the MDLb using a single rule
#calc_mdlb          Calculates the MDLb applying rules 1, 2A, 2B, and 3 from EPA document
#clean_names_base   ???

#%>% means "and then"
#commas are used for steps that go together

#----------------------------------------------
#Setup environment

setwd('/Users/patrickheraghty/Downloads/2025 MDL Sandbox/Example for Collin') #Set working directory for remainder of code

#Define variables
#read.csv imports data from our working directory to R
#Once variables are defined, variable names are used in place of sheet names

LOQ <- read.csv("LOQs 2023-2024.csv")                         #Define LOQ values
Spikeold <- read.csv("TO-15 MDLs Data Init.csv")              #Define initial spike values
Spikecurrent <- read.csv("TO-15 MDLs Data Ongoing.csv")       #Define ongoing spike values
mdlb <- read.csv("TO-15 MDLb Data.csv")                       #Define blank values

#Call in packages to R
library(EnvStats)     #Official environmental stats package
library(readxl)       #Excel reader package
library(dplyr)        #Data manipulation package (part of tidyverse)
library(janitor)      #Cleaning variable naming conventions

#----------------------------------------------
# >>>CHANGE<<< Global rounding controls (display/export only; EPA math stays full precision)
ROUND_DIGITS <- 4                     # set precision to match verified values
ROUND_MODE   <- "half_up"             # "half_up" (Excel), "half_even" (bankers), or "half_odd"

round_half_up_scalar <- function(x, digits) {
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- floor(z + 0.5)
  posneg * z / 10^digits
}
round_even_odd <- function(x, digits = ROUND_DIGITS, mode = c("half_up","half_even","half_odd")) {
  mode <- match.arg(mode)
  if (mode == "half_up")   return(round_half_up_scalar(x, digits))
  if (mode == "half_even") return(round(x, digits))  # base R does ties-to-even
  s <- 10^digits   # ties-to-odd
  z <- x * s
  flo <- floor(z); frac <- z - flo; up <- flo + 1
  tol <- .Machine$double.eps^0.5
  is_tie <- abs(frac - 0.5) <= tol
  out <- ifelse(frac < 0.5, flo, ifelse(frac > 0.5, up, NA_real_))
  pick_odd <- ifelse((flo %% 2) != 0, flo, up)
  out[is_tie] <- pick_odd[is_tie]
  out / s
}
round_df_mode <- function(df, digits = ROUND_DIGITS, mode = ROUND_MODE) {
  dplyr::mutate(df, dplyr::across(where(is.numeric), ~ round_even_odd(., digits = digits, mode = mode)))
}

#----------------------------------------------
#Define MDLs function for Initial Spikes (complete dataset)

#mdl_general is a one imput function
#"%>%" translates to "and then"

#Define mdl_general function
mdl_general <- function(values)
{
  values <- as.numeric(values)       #Make sure all inputs are treated as numbers
  values <- values[!is.na(values)]   #Remove missing NA values, if needed
  n <- length(values)                #Number of spikes
  if (n < 7) return(NA)              #Ensures min 7 spikes are present
  s <- sd(values)                    #Standard deviation
  t_val <- qt(0.99, df = n - 1)      #Student T value (one-tailed 0.99)
  return(t_val * s)                  #Returns function result
}

#Values is set in the future of the code

#----------------------------------------------
#Call MDLs function for Initial Spikes (complete dataset)

Spikeold <- Spikeold %>%                        #Failsafe
  rowwise() %>%                                 #Tells r to process one row at a time (for parsing through xanayltes)
  mutate(                                       #Mutates data set. Adds extra column in our future dataset.
    MDL = mdl_general(c_across("X.12":"X.19"))  #Apply our Mdl_General function to rows 12-19
  ) %>%
  ungroup()                                     #Ungroup undoes the "rowwise()" modification we no longer need.

# View(Spikeold)        # >>>CHANGE<<< commented out to show only final results view

#----------------------------------------------
#Calculating MDL for LOQ data

LOQ_with_MDL <- LOQ %>%                        #Look at LOQ data, and then
  group_by(xlodsetid, xanalyte, xunits) %>%    #Group columns xlodsetid, xanaylte, and xunits
  mutate(
    n_reps = n(),                              #Count how many replicate rows per group
    sd_val = sd(xqcspmes, na.rm = TRUE),       #Calculate SD of measured values
    MDL = mdl_general(xqcspmes)                # >>>CHANGE<<< ensure we use this defined function
  ) %>%
  ungroup()                                    #Restore data to original form

##View(LOQ_with_MDL)    # >>>CHANGE<<< commented out to show only final results view

#----------------------------------------------
#Summary Table from LOQ MDLs

LOQ_summary <- LOQ_with_MDL %>%                #Look at LOW_with_MDL data, and then
  group_by(xlodsetid, xanalyte, xunits) %>%    #Regroup by test set, xanaylte and unit of measurement.
  summarise(                                   #Create summary of
    n_reps = first(n_reps),                    #Replicate count (keep replicate count), pulls the first val
    sd_val = first(sd_val),                    #Standard deviation (keep standard deviation)
    MDL    = first(MDL),                       #MDL Value (keep the MDL value)
    .groups = "drop"                           #Return to regular data frame
  )

# View(LOQ_summary)     # >>>CHANGE<<< commented out to show only final results view

#----------------------------------------------
#Clean/Format data

#Fix spaces and parentheses, making sure that data formatting is not a future problem
mdlb <- janitor::clean_names(mdlb) %>%         # >>>CHANGE<<< assign back; previously not assigned
  #It is best to include all potential variables. There is no downside to having extra
  select(prep_id, run_batch, sample_id, matrix_gas, sample_type,
         date, time, analyte, units, concentration_blank_result)

#----------------------------------------------
#Define MDLb CALCULATION function for complete dataset (mdlbfunction)

#Create MDLb calc function, use EPA rules
calc_mdlb <- function(values, alpha = 0.01){                #Set alpha and define function so we can calculate input values
  values <- as.numeric(values)                              #Ensure everything is numeric
  values <- values[!is.na(values)]                          #Drop missing values (both these functions are similar to janitor)
  n <- length(values)                                       #Number of values
  if (n == 0) return(NA_real_)                              # >>>CHANGE<<< guard
  
  zeros <- sum(values == 0)                                 # >>>CHANGE<<< strict zero to mirror spreadsheets
  nonzeros <- n - zeros
  
  #Scenario 1: All zeros
  if (zeros == n) {                                         #IF all values are zero, then
    return(0)                                               #return output '0'
  }
  
  #Scenario 2A: Mix of zeros & non-zeros, n < 100
  if (zeros > 0 && nonzeros > 0 && n < 100) {               #IF mixed !0 and 0, and n is below 100, then
    return(max(values, na.rm = TRUE))                       #Return only the highest number
  }
  
  #Scenario 2B: Mix of zeros & non-zeros, n >= 100
  if (zeros > 0 && nonzeros > 0 && n >= 100) {              #If mixed !0 and 0, n is above 100, then
    return(as.numeric(quantile(values, probs = 0.99, type = 7, na.rm = TRUE)))  # >>>CHANGE<<< EPA 99th pct
  }
  
  #Scenario 3: No zeros
  if (nonzeros == n) {                                      #IF every single value is not 0, then
    xbar <- mean(values)                                    #Set xbar using the mean function on values
    if (xbar < 0) xbar <- 0                                 #IF xbar is > 0, Fail Safe
    s <- sd(values)                                         #Calculate standard deviation
    t_val <- qt(0.99, df = n - 1)                           #Student T value
    return(xbar + t_val * s)
  }
  
  return(NA_real_)}     #Fail safe: Ensure NAs dont cause collapse

#NA_Real = Return answer, unless no answer then return NA. Fail safe, ensures NA's don't cause collapse

#----------------------------------------------
#Call MDLb CALCULATION function

mdlb_results <- mdlb %>%                                 #Defines mdlb_results using mdlb and then (%>% means and then)
  group_by(analyte, units) %>%                           #Group using two things - anaytle and units
  summarise(
    n_reps = n(),                                        #Set n reps using standard r function, n()
    MDLb = calc_mdlb(concentration_blank_result),        #Set MDLb using the mdlb function on the conc_blank results column
    .groups = "drop"                                     #Drop our groupings we did, Fail safe (same function as ungroup)
  )
# View(mdlb_results)       # >>>CHANGE<<< commented out to show only final results view

#----------------------------------------------
#Calculate MDLs from Current Spiked Data

# >>>CHANGE<<< Robust column detection (handles case/spaces; keeps your downstream names)
sc <- janitor::clean_names(Spikecurrent)
compound_col  <- names(sc)[grepl("^(compound|analyte)$|compound|analyte", names(sc), ignore.case = TRUE)][1]
units_col_sc  <- names(sc)[grepl("unit", names(sc), ignore.case = TRUE)][1]
actual_col    <- names(sc)[grepl("actual", names(sc), ignore.case = TRUE)][1]
sample_id_col <- names(sc)[grepl("^sample[_ ]?id$|sample[_ ]?id", names(sc), ignore.case = TRUE)][1]
if (is.na(compound_col) || is.na(units_col_sc) || is.na(actual_col)) {
  stop("Spikecurrent is missing required columns (compound/analyte, units, actual).")
}

CurrentMDLs <- if (!is.na(sample_id_col)) {
  sc %>%
    filter(!is.na(.data[[actual_col]]), !is.na(.data[[compound_col]]), !is.na(.data[[units_col_sc]])) %>%
    group_by(.data[[compound_col]], .data[[units_col_sc]], .data[[sample_id_col]]) %>%
    summarise(actual = sum(.data[[actual_col]], na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[compound_col]], .data[[units_col_sc]]) %>%
    summarise(
      n_reps = n(),                           #Number of replicate samples
      sd_val = sd(actual, na.rm = TRUE),      #Define Standard Deviation
      MDL    = mdl_general(actual),           #Call EPA MDL Formula
      .groups = "drop"
    )
} else {
  sc %>%
    filter(!is.na(.data[[actual_col]]), !is.na(.data[[compound_col]]), !is.na(.data[[units_col_sc]])) %>%
    group_by(.data[[compound_col]], .data[[units_col_sc]]) %>%
    summarise(
      n_reps = n(),                           #Number of replicate samples
      sd_val = sd(.data[[actual_col]], na.rm = TRUE),  #Define Standard Deviation
      MDL    = mdl_general(.data[[actual_col]]),       #Call EPA MDL Formula
      .groups = "drop"
    )
}

# >>>CHANGE<<< Standardize to your expected names for the join
names(CurrentMDLs)[names(CurrentMDLs) == compound_col] <- "compound"
names(CurrentMDLs)[names(CurrentMDLs) == units_col_sc] <- "units"

# View(CurrentMDLs)        # >>>CHANGE<<< commented out to show only final results view

#----------------------------------------------
#Function for outputting what rule the blanks were calculated using

classify_epa_rule <- function(values, tol = .Machine$double.eps^0.5) {
  values <- as.numeric(values)
  values <- values[!is.na(values)]
  n <- length(values)
  if (n == 0) return(NA_character_)
  zeros <- sum(values == 0)                                 # >>>CHANGE<<< strict zero so labels mirror spreadsheets
  nonzeros <- n - zeros
  
  if (zeros == n) {
    return("1: all 0 values (MDL = 0)")
  } else if (zeros > 0 && nonzeros > 0 && n < 100) {
    return("2A: 0 and non-0 values, < 100 data pts (MDL = max)")
  } else if (zeros > 0 && nonzeros > 0 && n >= 100) {
    return("2B: 0 and non-0 values, >= 100 data pts (MDL = 99th percentile MDL)")
  } else if (nonzeros == n) {
    return("3: no 0 values (MDL = mean + t*sd)")
  } else {
    return(NA_character_)
  }
}

#----------------------------------------------
#Create final table

#Step 1: Call MDLb Calc Function, define values to be used in final table sheet
mdlb_results_full <- mdlb %>%             #Set mdlb_results by using mdlb data
  group_by(analyte, units) %>%            #Grouping by xanalyte
  summarise(
    Number_of_Data_Points_MDLb = n(),                                                    #Set n_reps = n(mdlb)
    MDLb = as.numeric(calc_mdlb(concentration_blank_result)),                             #Run MDLb Calc Function on the blank results
    
    Highest_Value = max(concentration_blank_result, na.rm = TRUE),                       #Defines highest value column
    Average_Value = mean(concentration_blank_result, na.rm = TRUE),                      #Defines average value column
    Std_Dev_MDLb = sd(concentration_blank_result, na.rm = TRUE),                         #Defines st dev value column
    
    Top_1_Percent = max(concentration_blank_result, na.rm = TRUE),                       # >>>CHANGE<<< DISPLAY MAX to match verified export
    Student_T_Used = ifelse(                                                             
      all(suppressWarnings(as.numeric(concentration_blank_result)) != 0) & n() >= 2,      # >>>CHANGE<<< only for Rule 3 groups
      qt(0.99, df = n() - 1),
      NA_real_
    ),
    .groups = "drop"
  )

#Step 2: Clean and Combine MDLb and MDLs
mdlb_clean <- mdlb_results_full %>%                 
  rename(                                      
    Compound = analyte,                        #Rename for consistency in nomenclature
    Units = units                              #Rename for consistency in nomenclature
  )

mdls_clean <- CurrentMDLs %>%                  
  janitor::clean_names() %>%                   
  rename(                                      
    Compound = compound,                       #Rename for consistency in nomenclature
    Units = units,                             #Rename for consistency in nomenclature
    MDLs = mdl,                                #Rename for consistency in nomenclature
    Number_of_Data_Points_MDLS = n_reps,       #Rename for consistency in nomenclature
    Std_Dev = sd_val                           #Rename for consistency in nomenclature
  )

#Calls Classify_EPA_Rule function
epa_rules <- mdlb %>%                                              
  group_by(analyte, units) %>%                                     
  summarise(                                                       
    EPA_Rule = classify_epa_rule(concentration_blank_result),      #This is a custom function
    .groups = "drop"                                               
  ) %>%                                                            
  rename(                                                          
    Compound = analyte,                                            
    Units = units                                                  
  )

# Step 3: Combine MDLb and MDLs, compute final values
# >>>CHANGE<<< inner join to match verified row set (analytes present in both)
final_mdls <- mdlb_clean %>% 
  inner_join(mdls_clean, by = c("Compound", "Units")) %>%           
  left_join(epa_rules,  by = c("Compound", "Units")) %>%           
  mutate(
    MDL_Final = pmax(MDLb, MDLs, na.rm = TRUE),                    #evaluate mdlb and mdls per anaylte (greater of the two)
    Comparison = case_when(                                        
      is.na(MDLb) | is.na(MDLs) ~ "insufficient data",
      MDLb > MDLs ~ "MDLb > MDLs",                                 
      MDLb < MDLs ~ "MDLs > MDLb",                                 
      TRUE ~ "MDLb = MDLs"                                         
    )
  ) %>%
  select(                               # << keep your original final-table format >>
    Compound,                          
    Number_of_Data_Points_MDLS,
    Number_of_Data_Points_MDLb,        
    Std_Dev,
    MDLs,
    MDLb,
    MDL_Final,                         
    Comparison,                        
    Highest_Value,                     
    Average_Value,                     
    Top_1_Percent,                     # DISPLAY only; MDLb calc uses 99th percentile when n>=100 per calc_mdlb()
    Student_T_Used,
    EPA_Rule                           
  )

# >>>CHANGE<<< Apply global rounding at the end (reporting only)
final_mdls <- round_df_mode(final_mdls, digits = ROUND_DIGITS, mode = ROUND_MODE)

# ----- ONLY show final results view -----
View(final_mdls)

# If you also want a file export while keeping the same view behavior:
# write.csv(final_mdls, "TO-15 MDL Results (R-only Minimal).csv", row.names = FALSE)

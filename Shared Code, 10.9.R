setwd('C:/Users/cleg1/OneDrive/Desktop/Setwd') #Designate a folder to search through, this will be consistent throughout
LOQ <- read.csv("LOQs 2023-2024.csv") #By using read.csv, we can import any data we can get in our folder into R.
Spikeold <- read.csv("Initial Spike Collection TO-15.csv")
Spikecurrent <- read.csv("MDLsongoing2.csv")
mdlb <- read.csv("MdlbDataAll.csv")
#We import the data set from excel, names like "MDLbDataAll.csv" will not be used again.
#However, these data sets: LOQ(old blanks), Spikeold, Spikecurrent, MDLb(ongoing blanks) are the names we will use for now
#We have things broken into parts, with some setup before spikeold will transform into Initial Spikes (with mdls), followed 
#by more of the same for the other datasets, as well as a finalized analysis.


#Return(NA_Real) what does it do?, define what leftjoin/rightjoin diff is. FIx classify epa rule and find spot
#install.packages(c("EnvStats", "readxl", "dplyr","janitor")) #First part of a two part sequence, we are installing these 
#packages onto R from a place in our computers files. Next we will do step 2
library(EnvStats) #We will be uploading these 3 from our library of downloaded stuff into our active working directory.
library(readxl) #There is also an extra first step, which is downloading these 3 packages from the internet to a computer
library(dplyr) #A lot of the unique functions of R don't come default. In order to complete 
#this project, we will specialize our data arrangement by downloading an official Environmental Stats package
#then run it alongside DPlyr and readxl.

#This is the second part of setup. Using the envstats package, we will create a mdl formula.
mdl <- function(y, alpha = 0.01) { #create function using alpha 0.01 value.
  n <- length(y) #It parses through different xanayltes row by row, setting length(y) equal to n.
  s <- sd(y, na.rm = TRUE) #s is set as standard dev, y will be inputted by the data we analyze.
  t_val <- qt(1 - alpha, df = n - 1) #alpha is pre-defined as 0.01, also set up degree of freedoms function.
  mdl_val <- t_val * s #Outputted student t * outputted standard deviation = MDL value
  return(mdl_val) #This will show our final answer for the custom made function
}
#This a designed to be a replicate of the official 2019 MDL formula
#Here is the MDl formula for reference, directly sourced from the official 2019 document.

#MDLs = tt(n−1, 1−∝=0.99)Ss
#  MDLs = the method detection limit based on spiked samples
#t(n-1, 1-α = 0.99) = this represents calculations consistent with a
#t-value appropriate for a single-tailed 99th percentile formula. This ensures
#our mdl is calculated with less than a 1% chance of being inaccurate.
#t statistic and a standard deviation estimate with n-1 degrees of freedom
#Ss = sample standard deviation of the replicate spiked sample analyses.
#make them look the same, translate using mdl sheet as basis.

freon113_spikes <- c(2.23, 2.05, 2.09, 2.20, 2.13, 2.18, 2.27, 2.09) #test case for single xanaytle
mdl_freon113 <- mdl(y = freon113_spikes) #We defined Freon spikes quantitative values, now we can use our mdl function
print(mdl_freon113) #mdl(y=) will run the mdl function, with y being = Freon spikes in this case. We can then print it.
#We want to check if the mdl formula works, in the future we will build upon our functions and look to do all xanalytes
#at once. Upon testing, we can see that we got exactly the same results as previously calculated on the spreadsheet.

#Doing everything this way wouldn't work well for a database. This upcoming formula is created in order to accomplish the
#task. Alike the previous MDL formula, but with modified code in order to be able to test the entire data set simultaneously.
mdl_general <- function(values) {   #Our new function is going to be called mdl_general
  values <- as.numeric(values) #We make sure all inputs are treated as numbers
  values <- values[!is.na(values)] #Remove missing NA values, if needed
  n <- length(values) #We need to assign n to the value length in order to keep track of replicate values. As you know
  if (n < 7) return(NA) #you need 7 replicates in order to do this testing properly & ensure we only test proper data
  s <- sd(values) #Standard values need to be calculated before we can get our final number
  t_val <- qt(0.99, df = n - 1) #df = n-1 is a typical t configuration, and with the alpha set & t_val defined
  return(t_val * s) #We can finally do Test_statistic * standard deviation, which will give us our EPA mdl value.
}

SpikesInitial <- Spikeold %>%  #%>% is computer talk for "and then"
  rowwise() %>% #We are telling r to process one row at a time, this is how it parses through xanayltes
  mutate( #We are mutating the data set, we will see an extra column in our future dataset.
    MDL = mdl_general(c_across("X.12":"X.19"))  #Apply our Mdl_General function to rows 12-19
  ) %>%
  ungroup() #ungroup undoes the "rowwise()" modification we no longer need.
#We take spikecrt and ask r to create a new column based on values S1-S8, aka x.12-x.19, being plugged into the
#general mdl formula as outlined by the EPA. We are able to apply this to the entire SPIKECRT and create a new
#datasheet almost instantly, SpikeCRT_with_MDL
#Lets use a function to take a look at our new datasheet.
View(Spikeold) #Wow! All MDL's are consistent with those already calculated


LOQ_with_MDL <- LOQ %>% #Look at LOQ data, and then
  group_by(xlodsetid, xanalyte, xunits) %>% #group columns xlodsetid, xanaylte, and xunits
  mutate(
    n_reps = n(),#count how many replicate rows per group
    sd_val = sd(xqcspmes, na.rm = TRUE), #Calculate SD of measured values
    MDL = mdl(xqcspmes) #Apply our MDL function to replicate measurements
  ) %>%
  ungroup() #restore data to original form
View(LOQ_with_MDL) #view data

LOQ_summary <- LOQ_with_MDL %>%
  group_by(xlodsetid, xanalyte, xunits) %>% #Regroup by test set, xanaylte and unit of measurement.
  summarise(
    n_reps = first(n_reps),   #keep replicate count
    sd_val = first(sd_val),   #keep standard deviation
    MDL    = first(MDL),      #Keep the MDL value
    .groups = "drop"       #return to regular data frame
  )
View(LOQ_summary) #We now have an effective summary table



library(dplyr)
library(janitor)#Time to start sweeping things clean!

mdlb <- janitor::clean_names(mdlb)
janitor::clean_names(mdlb) %>%   # janitor+clean_names fix spaces & parentheses, making sure data formatting
  select(prep_id, run_batch, sample_id, matrix_gas, sample_type, #isn't a future problem. It's best to include all
         date, time, analyte, units, concentration_blank_result) #potential variables, no real downside to having extra either
#to make sure the length of our script isn't a problem, we should be very consistent about cleaning our data and formatting it
#MDLb calc function, using EPA rules
calc_mdlb <- function(values, alpha = 0.01){   #set alpha and define function so we can calculate input values
  values <- as.numeric(values)          # ensure everything is numeric
  values <- values[!is.na(values)]      # drop missing values (both these functions are similar to janitor)
  n <- length(values) #now that all our data has been clean and alpha value is set, we include n in our calculation
  
  #Scenario 1: all zeros
  if (all(values == 0)) { #IF all values are zero, then
    return(0) #return output '0'
  }
  
  #Scenario 2A: mix of zeros & non-zeros, n < 100
  if (any(values == 0) && any(values != 0) && n < 100) { #if mixed !0 and 0, and n is below 100, then
    return(max(values)) #return only the highest number
  }
  
  #Scenario 2B: mix, n >100
  if (any(values == 0) && any(values != 0) && n >= 100) { #if mixed !0 and 0, n is above 100 then we need something
    #more complicated. Since we have over 100 values, doing it by EPA requirements requires some sorting
    sorted_vals <- sort(values) #clearly label our sorted values, and use sort, a function default in r
    rank_index <- ceiling(n * 0.99) #Record by rank, and also we have to set a ceiling using 0.99, based of the
    return(sorted_vals[rank_index]) #epa requirements for going over 100, we want to record and return the highest
    #few values
  }
  
  #scenario 3: no zeros
  if (all(values != 0)) { #IF every single value is not 0
    xbar <- mean(values) #set xbar using the mean function on values
    if (xbar < 0) xbar <- 0 #IF xbar is >0, Fail Safe
    s <- sd(values) #Calculate standard deviation
    t_val <- qt(0.99, df = n - 1)  #set t value formula with n-1 and a 0.99 alpha
    return(xbar + t_val * s) #Now that we have everything defined and set, execute the function with xbar values,
  }#added to the t-val we have set with n and our alpha of 0.99, along with var standard deviation, are all being used.
  return(NA_real_)}
#NA_Real = Return answer, unless no answer then return NA. Fail safe, ensures NA's don't cause collapse

#By dividing everything into buckets, we are able to perfectly follow the mdl blank rules of the EPA.
mdlb_results <- mdlb %>% #defines mdlb_results using mdlb and then (%>% means and then)
  group_by(analyte, units) %>% #group using two things - anaytle and units
  summarise(
    n_reps = n(),#set n reps using standard r function, n()
    MDLb = calc_mdlb(concentration_blank_result),#Set MDLb using the mdlb function on the conc_blank results column
    .groups = "drop" #Drop our groupings we did, Fail safe    
  )
View(mdlb_results)#Wonderful! We can create and anaylze large amounts of mdlb data


SpikeOngoing <- Spikecurrent %>% #Time to get back to spikes, with Spikes ongoing!
  filter(!is.na(Actual),!is.na(Compound), !is.na(Units)) %>% #Filter to make sure the data exists as it should
  group_by(Compound, Units, Sample.ID) %>% #Group using Compound and sample id (units should generally be the same)
  summarise(
    actual = sum(Actual, na.rm = TRUE),  # add replicates together
    .groups = "drop"
  ) %>%
  
  # Calculate MDLs across all replicates of each analyte/quarter
  group_by(Compound, Units) %>%
  summarise(
    n_reps = n(),                           # Number of replicate samples equal to n()
    sd_val = sd(actual, na.rm = TRUE),      # Set standard deviation
    MDL    = mdl_general(actual),           # use our EPA MDL formula on the actual(actual is where the final conc. is)
    .groups = "drop"                        # For the initial data sheet.
  )

View(SpikeOngoing)

#Now that we have SpikesOngoing and SpikesInitial along with blanks ongoing and initial, we can create our
#FINAL TABLE using comparisons between SpikesOngoing and BlanksOngoing to anaylze large chunks of analyte data.
mdlb_results <- mdlb %>% #set mdlb_results by using mdlb data, use janitor to clean the names
  janitor::clean_names() %>%
  group_by(analyte, units) %>% #before grouping by anaylte
  summarise(
    n_reps = n(),#Set n_reps = n(mdlb), the %>% from before allow us to use n()
    MDLb = as.numeric(calc_mdlb(concentration_blank_result)),# <- run calculation using our mdlb function on
    Highest_Value = max(concentration_blank_result, na.rm = TRUE), #the blank result. Define columns using functions
    Average_Value = mean(concentration_blank_result, na.rm = TRUE),#like mean, sd, and max.
    Std_Dev_MDLb = sd(concentration_blank_result, na.rm = TRUE), #Right now we are defining and naming these values
    Top_1_Percent = if(length(concentration_blank_result) >= 100) { #and later on we will execute the code necessary
      #to show these outputs in our final table sheet. Next we will sort the results and set a ceiling using 0.01
      sort(concentration_blank_result, decreasing = TRUE)[ceiling(0.01 * length(concentration_blank_result))]
    } else {
      max(concentration_blank_result, na.rm = TRUE)
    },
    Student_T_Used = qt(0.99, df = n() - 1),#add student_t value to be shown in the next table.
    .groups = "drop"
  )

mdlb_clean <- mdlb_results %>%
  rename(
    Compound = analyte, #using rename we will rename things so that they match between data sets,
    Units = units, #This is very important when doing comparisons!
    Number_of_Data_Points_MDLb = n_reps
  )

mdls_clean <-  SpikeOngoing%>% #Updating SpikesInitial to a "cleaned" version
  janitor::clean_names() %>% #We will clean everything and change it to a
  rename( #uppercased, "Student_T_Used" esque format
    Compound = compound,#start renaming
    Units = units,
    MDLs = mdl,
    Number_of_Data_Points_MDLS = n_reps,
    Std_Dev = sd_val
  )#Now everything in mdlb and mdls are "Units" & "Compound"
#along with a few specfications for each kind of mdl.

classify_epa_rule <- function(values) {#We will add an epa rule to make it akin to google doc
  v <- as.numeric(values)
  v <- v[!is.na(v)] #Make sure values are valid
  n <- length(v) #Set n and create different "if" scenarios to classify EPA outcomes
  
  
  if (n == 0) return(NA_character_)
  if (all(v == 0)) return("1: All zeros (Pass)")
  if (any(v == 0) && any(v != 0)) {
    if (n < 100) {
      return("2A: Mixed 0/non-0, n < 100 (MDL = max value)")
    } else {
      return("2B: Mixed 0/non-0, n >= 100 (MDL = 99th percentile)")
    }
  }
  if (all(v != 0)) return("3: No zeros (MDL = mean + t*SD)")
  return(NA_character_)
}


calc_mdlb <- function(values, alpha = 0.01) {
  values <- as.numeric(values) 
  values <- values[!is.na(values)]
  n <- length(values)
  
  # No values at all
  if (n == 0) return(NA_real_)
  
  # Scenario 1: all zeros
  if (all(values == 0)) return(0)
  
  # Scenario 2: mixture of zeros & non-zeros
  if (any(values == 0) && any(values != 0)) {
    if (n < 100) {
      return(max(values, na.rm = TRUE))
    } else {
      sv <- sort(values, na.last = NA)
      idx <- ceiling(n * 0.99)
      idx <- min(idx, length(sv))  # guard index
      return(sv[idx])
    }
  }
  
  # Scenario 3: all non-zeros
  if (all(values != 0)) {
    xbar <- mean(values, na.rm = TRUE)
    if (is.na(xbar)) return(NA_real_)
    if (xbar < 0) xbar <- 0
    s <- sd(values, na.rm = TRUE)
    t_val <- qt(1 - alpha, df = n - 1)
    return(xbar + t_val * s)
  }
  
  # default fallback
  return(NA_real_)
}#262-280 is VERY MESSY
epa_rules <- mdlb %>% #Lets add epa_rules into mdlb so we can get in our final summary!
  janitor::clean_names() %>% #Janitor appreciation month!
  group_by(analyte, units) %>% 
  summarise( #We can apply an epa_rule test
    EPA_Rule = classify_epa_rule(concentration_blank_result),#Function from earlier
    .groups = "drop"
  ) %>%
  rename(#Rename to "Units & "Compound"
    Compound = analyte,
    Units = units)

#Epa rules is converted to a table so that it can easily be added to final_mdls
epa_rules_summary <- mdlb %>%
  janitor::clean_names() %>%
  group_by(analyte, units) %>%
  summarise(EPA_Rule = classify_epa_rule(concentration_blank_result), .groups = "drop") %>%
  rename(Compound = analyte, Units = units)

#Combine MDLb and MDLs, compute final values
final_mdls <- mdlb_clean %>% #We will create final mdl by using mdlb_clean
  full_join(mdls_clean, by = c("Compound", "Units")) %>% #we clean mdls as well so that we can easily join them
  left_join(epa_rules_summary,  by = c("Compound", "Units")) %>% #to make it like the google sheet we will have to add columns  
  mutate(
    MDL_Final = pmax(MDLb, MDLs, na.rm = TRUE),#evaluate mdlb and mdls per anaylte
    Comparison = case_when(
      MDLb > MDLs ~ "MDLb > MDLs", #and compare them, making sure to define what happens,
      MDLb < MDLs ~ "MDLs > MDLb", #I.E. if MDL s is higher, then display MDLb < MDLs
      TRUE ~ "MDLb = MDLs"
    )
  ) %>%
  select(#select each of these things
    Compound,#things like compound get carried over from previous data
    Number_of_Data_Points_MDLS,        
    Number_of_Data_Points_MDLb,#but much of these names (which all have to be predefined before this to function)
    Std_Dev,                          
    MDLs,
    MDLb,
    MDL_Final,#are specially defined and created within the code.
    Comparison,#Which one is higher?
    Highest_Value,#Whats the highest value within the data set
    Average_Value,#Whats the average value within the dataset
    Top_1_Percent,#Show top 1% of data (usually will = mdl)
    Student_T_Used,
    EPA_Rule #custom!
  )

View(final_mdls)


#Enhance and compute final_mdls with flags
# Add a column that quantifies which side dominates and by how much
explain_final_mdl <- final_mdls %>%
  mutate(
    MDLb = as.numeric(MDLb),
    MDLs = as.numeric(MDLs),
    MDL_Final = pmax(MDLb, MDLs, na.rm = TRUE),
    MDL_Ratio = case_when(
      is.na(MDLb) & is.na(MDLs) ~ NA_real_,
      is.na(MDLb) ~ Inf,  # treat as MDLs only
      is.na(MDLs) ~ -Inf, # treat as MDLb only
      TRUE ~ MDLb / MDLs
    ),
    MDL_Flag = case_when(
      MDL_Ratio >= 2 ~ "MDLb >> MDLs (x2+)",
      MDL_Ratio >= 1.2 & MDL_Ratio < 2 ~ "MDLb moderately higher",
      MDL_Ratio > 0 & MDL_Ratio < 0.8333 ~ "MDLs moderately higher",
      MDL_Ratio <= 0.5 ~ "MDLs >> MDLb (x2+)",
      TRUE ~ "Similar"
    ),
    Likely_Cause = case_when(
      MDL_Flag == "MDLb >> MDLs (x2+)" ~ "Investigate blanks",
      MDL_Flag == "MDLs >> MDLb (x2+)" ~ "Investigate spikes: unstable spike prep or inconsistent spike recoveries",
      MDL_Flag == "Similar" ~ "No obvious problems",
      TRUE ~ NA_character_
    )
  )

View(explain_final_mdl)

#Analyte diagnostics (Mission: Find outliers,  check ratios)


# We want to anaylze why mdls are higher or mdlb are higher. We can correlate outside factors that commonly lead
#to mldb or mdls being high, such as  Trend anaylsis of mdlb if its rising over time, we can use time based anaylsis 
#and check if instruments need to be cleaned. Think of specific indications like this. A lot of the factors depend
#on cleaning. We can also look into which xanyalytes cleaning are most essential or has disruptions of influence
#If we can prove instrument needed cleaning. Trend analysis over time for different xanalytes.
#mark out ones where highest value is significantly higher than mdlb(aka second highest value)
#Way to analyze perhaps top 10%? Or just see all the points, ranked in order for specific analytes , try to anaylze how big of
#an outlier be. Look at top 10 and analyze trends.
#Next week questions about the stuff i sent




#looking at target vs actual, computer may not notice but sometimes target is different than other times for an
#analtyte. We need to filter by  the "target value" Seeing large variance in the same spike targets, 
#perhaps sort by actual vs theoretical. Could identify changes, but I don't know which data sets have a target value.

# We want to create a new distinction in flags, we want to use ongoing trend anaylsis for looking at the data
#Setup difference between trends and fixes, first task : identify where a fix needs to be made. This is where the theoretical
#outcome for a single analyte doesn't much its peers. FLAG. "Data != Homogenous" 2 output columns, trend and error analysis.
#Second error we could have, lets focus on mdlb >> mdls. This will trigger an "investigate" flag. Try 2 finish by tomorrow
#If possible, indicate both errors if multiple. 
#Extra flagging

#Two columns at the end of exported table, with one being logic message and one being error message.
#Focus on mdlb 2x >mdls, "Attention, blanks could be vetted for optimized result" identify where there are outlier spike concentrations, 
#one of spikes differing from peers, breaks mdl procedure rules. "Error, spikes inconsistent across analyte"
#Polish results well, continue updating final mdls while also creating functions that do extra anaylsis on spike / blanks
#and demonstrating them on MDLb and SpikesOngoing. 

#run through code, demo after deleting all
#Focus on first two flag objectives, consider going back and creating more diverse flag types

# Detect outlier spikes — spikes differing too much from others
detect_spike_outlier <- function(values, threshold = 2) {
  values <- as.numeric(values)
  values <- values[!is.na(values)]
  if (length(values) < 3) return(FALSE)  # not enough data
  
  mean_val <- mean(values)
  sd_val <- sd(values)
  any(abs(values - mean_val) > threshold * sd_val)
}
detect_spike_outlier(SpikesInitial)

# Classify logical and error messages for final table
logic_and_error_flags <- function(mdlb_results, SpikesInitial) {
  logic_msg <- "Pass"
  error_msg <- NA_character_
  
  # Logic flag: MDLb should not exceed 2x MDLs
  if (!is.na(mdlb) && !is.na(mdls) && mdlb > 2 * mdls) {
    logic_msg <- "Attention, blanks could be vetted for optimized result"
  }
  
  # Spike inconsistency flag
  spike_vals <- spikes_data %>% 
    filter(Compound == compound) %>%
    pull(Actual)
  
  if (length(spike_vals) >= 3 && detect_spike_outlier(spike_vals)) {
    error_msg <- "Error, spikes inconsistent across analyte"
  }
  
  return(list(logic = logic_msg, error = error_msg))
}


# Add logic and error message columns
final_mdls <- final_mdls %>%
  rowwise() %>%
  mutate(
    flag = list(logic_and_error_flags(mdlb_results, SpikesInitial)),
    Logic_Message = flag$logic,
    Error_Message = flag$error
  ) %>%
  ungroup() %>%
  select(Compound, Units, MDLb, MDLs, MDL_Final, Comparison, EPA_Rule,
         Logic_Message, Error_Message)

# Summary of spike precision across analytes
analyze_spike_precision <- function(spike_data) {
  spike_data %>%
    group_by(Compound) %>%
    summarise(
      Mean = mean(Actual, na.rm = TRUE),
      SD = sd(Actual, na.rm = TRUE),
      RSD = (SD / Mean) * 100,
      Outlier_Flag = detect_spike_outlier(Actual)
    ) %>%
    arrange(desc(RSD))
}

# Analyze blank distributions for MDL variability
analyze_blank_variation <- function(blank_data) {
  blank_data %>%
    group_by(analyte) %>%
    summarise(
      Mean_Conc = mean(concentration_blank_result, na.rm = TRUE),
      SD_Conc = sd(concentration_blank_result, na.rm = TRUE),
      Replicates = n()
    ) %>%
    arrange(desc(SD_Conc))
}

# Example demonstration:
spike_variation_report <- analyze_spike_precision(SpikesInitial)
blank_variation_report <- analyze_blank_variation(mdlb_results)


logic_and_error_flags <- function(mdlb_results, SpikesInitial) {
  logic_msg <- "Pass"
  error_msg <- NA_character_
  
  # Coerce to single numeric values
  MDLb <- as.numeric(MDLb[1])
  MDLs <- as.numeric(MDLs[1])
  
  # Check MDL ratio
  if (!is.na(MDLb) && !is.na(MDLs) && MDLb > 2 * MDLs) {
    logic_msg <- "Attention, blanks could be vetted for optimized result"
  }
  
  # Check spike inconsistency
  spike_vals <- spikes_data %>%
    filter(Compound == compound) %>%
    pull(Actual)
  
  if (length(spike_vals) >= 3 && any(abs(spike_vals - mean(spike_vals)) > 2 * sd(spike_vals))) {
    error_msg <- "Error, spikes inconsistent across analyte"
  }
  
  return(list(logic = logic_msg, error = error_msg))
}


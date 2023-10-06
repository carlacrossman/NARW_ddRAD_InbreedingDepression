##############################################################
#            Measuring Fitness via Fecundity                 #
#         in North Atlantic Right Whale Females              #
#             to assess inbreeding depression                #
##############################################################

#----------------------------------------------------#
#              Load Necessary Packages               #
#----------------------------------------------------#

library(ggplot2)
library(patchwork)
library(tidyverse)
library(stringr)

#----------------------------------------------------#
#                  Load Data Files                   #
#----------------------------------------------------#

# The code below imports the relevant data files. They were extracted tables provided by Phil Hamilton on February 15, 2023.


###--- Individuals ---###

# Load data file with a list of sample names
# Separate into column with species prefix and ID number
# Remove duplicates (ending in 'D' and remove EAU sample numbers)
# Remove the males for this analysis 1037 and 1131
# Keep only the nea number column

ind_samples <- read_csv("Fitness/ind_ids.csv", col_names = TRUE, show_col_types = FALSE)

ind_females <- ind_samples %>%
  separate(ID, into = c("species", "nea"), sep = 3) %>%
  filter(!(str_detect(nea, "D$") | species == "eau")) %>%
  filter(!(nea == 1037 | nea == 1131)) %>%
  mutate(nea = as.numeric(nea)) %>%
  select(nea)


###--- Calving Data---###

# Load data file containing all calving events including the mom, year and calf

calvingevent_temp <- read_csv("Fitness/data_files_feb2023/TblCalvingEvent-2023-02-15.csv", col_names = TRUE, show_col_types = FALSE) 


calvingevent <- calvingevent_temp %>%
  rename(nea = EGNo, year = CalvingYear) %>%
  select(nea, year, CalfNo)


###--- Whale Age ---###

# Load data file containing age class data for all years each nea was sighted

whale_age_temp <- read_csv("Fitness/data_files_feb2023/TblWhaleAge-2023-02-15.csv", col_names = TRUE, show_col_types = FALSE)

whale_age <- whale_age_temp %>%
  rename(nea = EGNo, year = AgeYear, age = Age, age_class = AgeClassCode, alive_status = AliveStatusCode)


###--- Whale ---###

# Load data file containing a list of all individuals with sex, mother, birth and death year if known, and first/last year sighted

whale_temp <- read_csv("Fitness/data_files_feb2023/TblWhale-2023-02-15.csv", col_names = TRUE, show_col_types = FALSE)

whale <- whale_temp %>%
  select(EGNo, GenderCode, FirstYearSighted, LastYearSighted, BirthYear, DeathYear, Mother, FirstCalvingYear) %>%
  rename(nea = EGNo, sex = GenderCode)


###--- Number of Adult Females Alive ---###
# Load data file with population size of NARW in different age cohortsfrom Pace model

total_alive_year <- read.csv("Fitness/PaceAFW.csv", header = TRUE)

total_females_alive_year <- total_alive_year %>%
  filter(grepl('NumAF[', Parameter, fixed = TRUE)) %>%
  filter(Year <= 2020)


# Key Data Files:

# - **ind_females** : list of individual females used in ddRAD study  
# - **calvingevent** : calving events for all females in the population
# - **whale_age** : age class data for all individuals in the population
# - **whale** : information on individual whale life history including nea, sex, birth and death years and first and last sighting years
# - **total_females_alive_year** : Number of adult females alive in each study year as determined by the Pace model



#----------------------------------------------------#
#                Load Helper Functions               #
#----------------------------------------------------#



###-------------- ALIVE --------------###

# Makes a list of the whales alive in a specified year.


alive = function(year) {
  
  #--- Initialize a vector for holding results ---#
  whales_year <-  NULL
  
  whales_year <-  whale_age[whale_age$year == year & (whale_age$alive_status == "A" | whale_age$alive_status == "R" | whale_age$alive_status == "OP" | whale_age$alive_status == "OSP"), 1]
  
  return(whales_year[["nea"]])
}


###-------------- ADULT --------------###

# Identify which individuals are alive and adult in a given year.


adult = function(year) {
  
  #--- Initialize a vector for holding results ---#
  age_year = NULL
  
  age_year = whale_age[whale_age$year == year & whale_age$age_class == "A" & (whale_age$alive_status == "A" | whale_age$alive_status == "R" | whale_age$alive_status == "OP" | whale_age$alive_status == "OSP"), "nea"]
  
  return(age_year[["nea"]])
}


###-------------- BIRTH --------------###

# Makes a list of the whales who gave birth in a specific year.


birth = function(time) {
  
  #--- Initialize a vector for holding results ---#
  birth_year = NULL
  
  birth_year <- calvingevent %>%
    filter(year == time) %>%
    select(nea)
  
  return(birth_year[["nea"]])
}



#----------------------------------------------------#
#        Prepare Summary Matrices and Lists          #
#----------------------------------------------------#



###--- Study Years ---###

study_years <- c(1990:2020) # Create a new vector that lists the years in our study period


###--- List all whales ---###

whale_list <- whale %>%
  select(nea) # Create a vector for all whales in the population so calculate survival for the species/population


###--- List of females ---###

female_whales <- whale[(whale[, "sex"] == "F"), 1] # Create a tibble of all whales that are females


###--- Matrix of which whales are alive in which years ---###

alive_by_year <- matrix(0, nrow = length(whale_list[["nea"]]), ncol = (length(study_years)))
colnames(alive_by_year) <- c(study_years)
rownames(alive_by_year) <- c(whale_list[["nea"]])

counter = 1
for (i in study_years) {
  alive_by_year[ , counter] <- as.numeric(whale_list[["nea"]] %in% alive(i))
  counter <- counter + 1
}


###--- Plot Number of Adult Females Over Time ---###

total_females_alive_year %>%
  ggplot(aes(Year, Mean)) +
  geom_point() +
  geom_line() +
  ylab("Mean # of Adult Females Alive")



#----------------------------------------------------#
#     Identify births for each female over time      #
#----------------------------------------------------#



# Create matrix of the female whales and study years
# Loop through years
#     - Report 0 or 1 for whether or not a female calved
#     - Record a list of the adults alive in a given year
# Loop through females
#     - If the female was not an adult in the study year, change the birth value to NA

births_by_year <- matrix(0, nrow = length(female_whales[["nea"]]), ncol = (length(study_years)))
colnames(births_by_year) <- c(study_years)
rownames(births_by_year) <- c(female_whales[["nea"]])

counter = 1
for (i in study_years) {
  births_by_year[ , counter] <- as.numeric(female_whales[["nea"]] %in% birth(i))
  adult_female_in_year <- female_whales[["nea"]] %in% adult(i)
  for (j in 1:length(female_whales[["nea"]])) {
    if(adult_female_in_year[j] == FALSE ) {
      births_by_year[j, counter] <- NA
    }
  }
  counter <- counter + 1
}


###--- Incorporate last year sighted ---###

# If the last year sighted was during the study period,
# For all years, after the last year sighted, replace values in the birth matrix with NA.
# This implies that the last year sighted could have been the last year alive.

for(i in row.names(births_by_year)){
  if(whale$LastYearSighted[whale$nea == i] < 2020 & whale$LastYearSighted[whale$nea == i] > 1990){
    births_by_year[i,as.character((whale$LastYearSighted[whale$nea == i]+1):2020)] <- NA
  }
}


###--- Calculate the number of births each year ---###

total_births_year <- data.frame(study_years, colSums(births_by_year, na.rm = TRUE))
colnames(total_births_year) <- c("year", "num_births")


###--- Plot number of births over time ---###

total_births_year %>%
  ggplot(aes(year, num_births)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(num_births), colour = "red"), linetype = "dashed" )



#----------------------------------------------------#
#             Annual Population Fecundity            #
#----------------------------------------------------#



###--- Calculate annual population fecundity ---###

annual_fecundity <- data.frame(study_years, total_births_year[["num_births"]] / total_females_alive_year[["Mean"]])
colnames(annual_fecundity) <- c("year", "pop_mean_fecundity")


###--- Plot annual fecundity ---###

annual_fecundity %>%
  ggplot(aes(year, pop_mean_fecundity)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(pop_mean_fecundity), colour = "red"), linetype = "dashed" )


###--- Calculate individual fecundity component for each year ---###

# Initialize an empty dataframe
# Calculate yearly contributions:
#   - For each whale, in each year starting with year 2 (1991)
#   - If the value for births for year-1, year, and year+1 is NA (the whale was not alive/adult), report NA
#   - If the sum of births for year-1, year, and year+1 is equal to or greater than 1, the female gave birth in that window. 
#     Calculate the fecundity contribution as 1 - the annuanl population fecundityof that year, divided by the number of whales alive - 1
#   - If the sum of births for year-1, year, and year+1 is equal to 0, the female did not give birth in that window.
#     Calculate the fecundity contribution as 0 - the annuanl population fecundityof that year, divided by the number of whales alive - 1
#   - Results in a data frame containing the fecundity contribution of each female in each year

fecundity_by_year <- data.frame(matrix(NA, nrow = length(female_whales[["nea"]]), ncol = (length(study_years))))
colnames(fecundity_by_year) <- c(study_years)

for (i in 1:length(female_whales[["nea"]])) {
  for (j in 2:(length(study_years)-1)) {
    if (sum(is.na(births_by_year[i, c(j-1, j, j+1)])) == 3) {
      fecundity_by_year[i, j] <- NA
    } else if (sum(births_by_year[i, c(j-1, j, j+1)], na.rm = TRUE) >= 1) {
        fecundity_by_year[i, j] <- (1 - (annual_fecundity[j, "pop_mean_fecundity"])) / (total_females_alive_year[j,"Mean"] - 1)
      } else if (sum(births_by_year[i, c(j-1, j, j+1)], na.rm = TRUE) == 0){
          fecundity_by_year[i, j] <- (0 - (annual_fecundity[j, "pop_mean_fecundity"])) / (total_females_alive_year[j,"Mean"] - 1)
        }
  }
}

row.names(fecundity_by_year) <- female_whales$nea

# This generates a data frame with each female represented by a row, each year by a column and 
# is populated with her individual fecundity contribution each year she was an adult and alive


###--- Drop females with less than 6 years of fecundity data ---###

# This requires females to have fecundity values for at least 6 years (therefore a minimum of 4 years alive & adult)
# Here we also calculate the mean fecundity contribution of each female

fecundity_by_year_6yr_min <- fecundity_by_year %>% 
  filter(rowSums(!is.na(.)) >= 6) %>%
  mutate(mean_fecundity = rowMeans(., na.rm = TRUE))

fecundity_by_year_6yr_min$nea <- row.names(fecundity_by_year_6yr_min)

###--- Calculate Total Contributions ---###

# total_fecundity_contribution <- data.frame(female_whales[["nea"]], rowSums(fecundity_by_year, na.rm = TRUE))
# colnames(total_fecundity_contribution) <- c("nea", "mean_ind_fecundity")

total_fecundity_contribution <- fecundity_by_year_6yr_min %>%
  select(nea, mean_fecundity)



#----------------------------------------------------#
#        Final Individual Fecundity Dataset          #
#----------------------------------------------------#



###--- Subset females from study list ---###

# Remove seven individual females from the list where ddRADseq data was poor
study_females <- ind_females %>%
  filter(!(nea %in% c(3320, 3294, 1620, 1140, 1201, 2042, 3292)))


###--- Subset fecundity data for study samples only ---###

ddRAD_fitness <- total_fecundity_contribution %>%
  filter(nea %in% study_females$nea)
write.table(ddRAD_fitness, file="Fitness/results/ddRAD_delifing_fitness.txt", sep="\t")


###--- Plot mean individual fecundity values to visualize distribution across samples ---###

ggplot(ddRAD_fitness, aes(as.factor(nea), mean_fecundity)) +
  geom_col() +
  ylab("Mean Individual Fecundity Contribution") +
  xlab("Individual") +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle=90,hjust=1))



#----------------------------------------------------#
#         SAMPLE PLOTS FOR A FEW INDIVIDUALS         #
#----------------------------------------------------#

# Here we generate sample plots to showing individual fecundity contributions for females with different reproductive histories

birth_plot <- total_births_year %>%
  ggplot(aes(year, num_births)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(num_births), colour = "red"), linetype = "dashed") +
  theme(legend.position = "none")

pop_fecundity <- annual_fecundity %>%
  ggplot(aes(year, pop_mean_fecundity)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(pop_mean_fecundity), colour = "red"), linetype = "dashed" ) +
  theme(legend.position = "none")

fecundity_by_year_long <- as.data.frame(fecundity_by_year_6yr_min) %>%
  mutate(nea = row.names(fecundity_by_year_6yr_min)) %>%
  gather(c(as.character(study_years)), key = "year", value="ind_fecund")



###--- nea 1407 ---###

ind_fecund_1407 <- fecundity_by_year_long %>%
  filter(nea == "1407") %>%
  ggplot(aes(as.numeric(year), ind_fecund)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(ind_fecund, na.rm = TRUE), colour = "red"), linetype = "dashed") +
  theme(legend.position = "none")

# births_by_year["1407",]
# alive_by_year["1407",]

p1_1407 <- birth_plot + 
  geom_vline(xintercept = c(1993, 2001), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1407], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ggtitle("nea 1407")

p2_1407 <- pop_fecundity + 
  geom_vline(xintercept = c(1993, 2001), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1407], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01)

p3_1407 <- ind_fecund_1407  + 
  geom_vline(xintercept = c(1993, 2001), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1407], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ylim(-0.015, 0.015)

plot_1407 <- p1_1407/p2_1407/p3_1407



###--- nea 1204 ---###

ind_fecund_1204 <- fecundity_by_year_long %>%
  filter(nea == "1204") %>%
  ggplot(aes(as.numeric(year), ind_fecund)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(ind_fecund, na.rm = TRUE), colour = "red"), linetype = "dashed") +
  theme(legend.position = "none")

# births_by_year["1204",]
# alive_by_year["1204",]

p1_1204 <- birth_plot + 
  geom_vline(xintercept = c(1991, 1995, 1999, 2002, 2005, 2009, 2013, 2019), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1204], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ggtitle("nea 1204")

p2_1204 <- pop_fecundity + 
  geom_vline(xintercept = c(1991, 1995, 1999, 2002, 2005, 2009, 2013, 2019), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1204], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01)

p3_1204 <- ind_fecund_1204  + 
  geom_vline(xintercept = c(1991, 1995, 1999, 2002, 2005, 2009, 2013, 2019), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1204], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ylim(-0.015, 0.015)

plot_1204 <- p1_1204/p2_1204/p3_1204



###--- nea 1719 ---###

ind_fecund_1719 <- fecundity_by_year_long %>%
  filter(nea == "1719") %>%
  ggplot(aes(as.numeric(year), ind_fecund)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(ind_fecund, na.rm = TRUE), colour = "red"), linetype = "dashed") +
  theme(legend.position = "none")

# births_by_year["1719",]
# alive_by_year["1719",]

p1_1719 <- birth_plot + 
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1719], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ggtitle("nea 1719")

p2_1719 <- pop_fecundity + 
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1719], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01)

p3_1719 <- ind_fecund_1719  + 
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1719], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ylim(-0.015, 0.015)

plot_1719 <- p1_1719/p2_1719/p3_1719



###--- nea 1620---###

ind_fecund_1620 <- fecundity_by_year_long %>%
  filter(nea == "1620") %>%
  ggplot(aes(as.numeric(year), ind_fecund)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(ind_fecund, na.rm = TRUE), colour = "red"), linetype = "dashed") +
  theme(legend.position = "none")

# births_by_year["1620",]
# alive_by_year["1620",]

p1_1620 <- birth_plot + 
  geom_vline(xintercept = c(1996, 2001, 2004, 2007, 2010, 2015), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1620], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ggtitle("nea 1620")

p2_1620 <- pop_fecundity + 
  geom_vline(xintercept = c(1996, 2001, 2004, 2007, 2010, 2015), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1620], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01)

p3_1620 <- ind_fecund_1620  + 
  geom_vline(xintercept = c(1996, 2001, 2004, 2007, 2010, 2015), colour = "blue") +
  geom_rect(aes(xmin = 1990, xmax = whale$LastYearSighted[whale$nea == 1620], ymin=-Inf, ymax = Inf), fill = "aquamarine", alpha = 0.01) +
  ylim(-0.015, 0.015)

plot_1620 <- p1_1620/p2_1620/p3_1620


###--- Summary Plots for Four Individuals ---###

(plot_1407 | plot_1719 | plot_1620 | plot_1204)



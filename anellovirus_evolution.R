library(dplyr)
library(data.table)
library(tidyverse)
library(paletteer)
library(reshape2)
library(lubridate)
library(gridExtra)
library(readxl)
library(broom)
library(ggpubr)
library(emmeans)

### Parsing the tables function - useful for other projects.
parse_my_table <- function(pattern) {
  # read file path
  all_paths <-
    list.files(pattern = pattern,
               full.names = TRUE)
  
  # read file content
  all_content <-
    all_paths %>%
    lapply(read.table,
           header = TRUE,
           sep = "\t",
           encoding = "UTF-8")
  
  # read file name
  all_filenames <- all_paths %>%
    basename() %>%
    as.list()
  
  # combine file content list and file name list
  all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)
  
  # unlist all lists
  all_result <- rbindlist(all_lists, fill = T)
  
  return(all_result)
}

all_result <- parse_my_table(".table")

### changing the names of variables (of the TPs)
all_result$V1 <- all_result$V1 %>% str_remove_all(fixed(".update.table")) %>% 
  str_remove_all(fixed("S3-")) %>% 
  str_remove_all(fixed("S4-"))

### renaming columns
all_result <- all_result %>% rename(TP = V1, Lineage = CHROM)
all_result <- all_result %>% mutate(TP = as.numeric(TP)) %>% select(-VAR)

### adding collection dates to the all_result dataframe
### I have the collection dates stored in another dataframe

years_TP <- read_csv("dates_col.csv")
years_TP$Date_of_collection <- mdy(years_TP$Date_of_collection)
years_TP$Date_of_collection <- as.Date(years_TP$Date_of_collection)

all_results_dates <- all_result %>% inner_join(years_TP, by = "TP")
all_results_dates$Date_of_collection <- all_results_dates$Date_of_collection %>% as.Date("%y-%m-%d")

### Adding the "change" column - which is the type of mutation and applying filters.
all_results_dates_mut <- all_results_dates %>%
  filter(TYPE == "SNP", AF > 0.05) %>%
  mutate(change = paste(REF, ALT))

### Selection of TPs based on the coverage - at least 90% coverage
### Loading a new dataset that contains samtools idxstats and mpileup results.

res_sam <- read.csv("results_samtools90.csv", stringsAsFactors = FALSE, header = TRUE)
reshape2::melt(res_sam) -> res_sam_melt

### cleanup
res_sam_melt$variable <- res_sam_melt$variable %>% str_remove(pattern = "TP")
res_sam_melt$variable <- res_sam_melt$variable %>% as.numeric()

### Selecting the time points per lineage based on coverage, copy number and presence of APOBEC
RPM <- read_excel("RPM.xlsx")
RPM_filtered <- RPM %>% filter(selected_TP == "yes") %>% select(-selected_TP, -note)

### Final selection of time points
selection_TP <- function(sel.lineage) {
  selection <- RPM_filtered %>% filter(Lineage == sel.lineage) %>% pull(TP)
  return(selection)
}

S1_01_selected_TP <- selection_TP("TTV-AMS-S1-01")
S1_23_selected_TP <- selection_TP("TTV-AMS-S1-23")
S1_41_selected_TP <- selection_TP("TTMV-AMS-S1-41")
S2_01_selected_TP <- selection_TP("TTV-AMS-S2-01")
S2_03_selected_TP <- selection_TP("TTV-AMS-S2-03")
S2_04_selected_TP <- selection_TP("TTV-AMS-S2-04")

### Making filtered data frames for each lineage, together with the follow-up period between the points.
## Dates of the first (reference) sample of all lineages
global_S1 <- as.Date("1985/11/15")
global_S2 <- as.Date("1987/09/07")
ref_S1_01 <- as.Date("1985/11/15")
ref_S1_23 <- as.Date("1985/11/15")
ref_S1_41 <- as.Date("1986/11/17")
ref_S2_01 <- as.Date("1990/06/15")
ref_S2_03 <- as.Date("2004/10/11")
ref_S2_04 <- as.Date("2005/10/17")

S1_01 <- all_results_dates_mut %>% 
  filter(Lineage == "TTV-AMS-S1-01", TP %in% S1_01_selected_TP) %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_01, 
         year_of_follow = as.numeric(day_of_follow / 365)) %>% 
  select(-year)
S1_01$year_of_follow <- format(round(S1_01$year_of_follow, 1), nsmall = 1)

S1_23 <- all_results_dates_mut %>% 
  filter(Lineage == "TTV-AMS-S1-23",TP %in% S1_23_selected_TP) %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_23, 
         year_of_follow = as.numeric(day_of_follow / 365)) %>% 
  select(-year)
S1_23$year_of_follow <- format(round(S1_23$year_of_follow, 1), nsmall = 1)

S1_41 <- all_results_dates_mut %>% 
  filter(Lineage == "TTMV-AMS-S1-41",TP %in% S1_41_selected_TP) %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_41, 
         year_of_follow = as.numeric(day_of_follow / 365)) %>% 
  select(-year)
S1_41$year_of_follow <- format(round(S1_41$year_of_follow, 1), nsmall = 1)

S2_01 <- all_results_dates_mut %>% 
  filter(Lineage == "TTV-AMS-S2-01",TP %in% S2_01_selected_TP) %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_01, 
         year_of_follow = as.numeric(day_of_follow / 365)) %>% 
  select(-year)
S2_01$year_of_follow <- format(round(S2_01$year_of_follow, 1), nsmall = 1)

S2_03 <- all_results_dates_mut %>% 
  filter(Lineage == "TTV-AMS-S2-03",TP %in% S2_03_selected_TP) %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_03, 
         year_of_follow = as.numeric(day_of_follow / 365)) %>% 
  select(-year)
S2_03$year_of_follow <- format(round(S2_03$year_of_follow, 1), nsmall = 1)

S2_04 <- all_results_dates_mut %>% 
  filter(Lineage == "TTV-AMS-S2-04",TP %in% S2_04_selected_TP) %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_04, 
         year_of_follow = as.numeric(day_of_follow / 365)) %>% 
  select(-year)
S2_04$year_of_follow <- format(round(S2_04$year_of_follow, 1), nsmall = 1)

### Putting these datasets together
all_sel_TPs <- rbind(S1_01, S1_23, S1_41, S2_01, S2_03, S2_04)
all_sel_TPs$year_of_follow <- as.numeric(all_sel_TPs$year_of_follow)

all_sel_TPs_S1 <- all_sel_TPs %>% 
  filter(TP <= 55) %>% 
  mutate(global_year = ymd(Date_of_collection), global_day_of_follow = global_year - global_S1, 
         global_year_of_follow = as.numeric(global_day_of_follow / 365)) %>%
  mutate(age = global_year_of_follow + 41) %>% 
  select(-global_year, -global_day_of_follow)

all_sel_TPs_S2 <- all_sel_TPs %>% 
  filter(TP > 55) %>% 
  mutate(global_year = ymd(Date_of_collection), global_day_of_follow = global_year - global_S2, 
         global_year_of_follow = as.numeric(global_day_of_follow / 365)) %>%
  mutate(age = global_year_of_follow + 35) %>% 
  select(-global_year, -global_day_of_follow)

all_sel_TPs <- rbind(all_sel_TPs_S1, all_sel_TPs_S2)
all_sel_TPs$global_year_of_follow <- format(round(all_sel_TPs$global_year_of_follow, 1), nsmall = 1)
all_sel_TPs$global_year_of_follow <- as.numeric(all_sel_TPs$global_year_of_follow)
all_sel_TPs <- all_sel_TPs %>% rename(lineage_year_of_follow = year_of_follow)
all_sel_TPs$age <- as.numeric(all_sel_TPs$age)
all_sel_TPs$age <- format(round(all_sel_TPs$age, 1), nsmall = 1)

write.csv(all_sel_TPs, "all_sel_TPs.csv")

### Graph displaying the selected lineages - how many years of follow

## add the first time points (just for the graphs)
missing_first_TPs <- tibble(
  Lineage = c("TTV-AMS-S1-01", "TTMV-AMS-S1-41", "TTV-AMS-S2-03", "TTV-AMS-S2-04"),
  TP = c(1, 4, 78, 80),
  age = c(41, 42.5, 52.1, 53.1)
)

all_sel_TPs %>% filter(lineage_year_of_follow > 0) %>% 
  ggplot(aes(lineage_year_of_follow, Lineage, fill = Lineage)) + 
  geom_point(size = 6, shape = 21, alpha = 0.9) +
  theme_minimal() +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35), limits = c(0, 35), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_blank(), 
        axis.text.x = element_text(vjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid"),
        plot.title = element_text(size = 20)) +
  labs(x = "Years of follow up", y = "Lineage", title = "A")

ggsave("Figure_years_of_follow_lin.png")
ggsave("Figure_years_of_follow_lin.svg")

all_sel_TPs %>% 
  ggplot(aes(as.numeric(age), Lineage, fill = Lineage)) + 
  geom_point(size = 6, shape = 21, alpha = 0.9) +
  geom_point(data = missing_first_TPs, aes(as.numeric(age), Lineage, fill = Lineage), size = 6, shape = 21, alpha = 0.9) +
  theme_minimal() +
  scale_x_continuous(breaks = c(35, 40, 45, 50, 55, 60, 65, 70, 75), limits = c(35, 75), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(vjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid"),
        plot.title = element_text(size = 20)) +
  labs(x = "Age of the subject", y = "")

ggsave("Figure_years_of_follow_lin_age.png")
ggsave("Figure_years_of_follow_lin_age.svg")

all_sel_TPs %>%
  group_by(Lineage) %>% 
  summarise(number_of_total_variations = n())

### Graph displaying the relationship between the years of follow up and number of total variations
all_sel_TPs %>%
  group_by(Lineage) %>% 
  summarise(max.age = max(age), number_of_total_variations = n()) %>% 
  ggplot(aes(as.numeric(max.age), number_of_total_variations, fill = Lineage)) +
  geom_point(size = 6, shape = 21, alpha = 0.9) +
  theme_minimal() +
  scale_x_continuous(limits = c(55, 75), breaks = seq(55, 75, 5), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid"),
        plot.title = element_text(size = 20)) +
  labs(x = "Age of the subject", y = "Total number of variations")

ggsave("total_number_of_var.svg")

all_sel_TPs %>%
  group_by(Lineage) %>% 
  distinct(POS, .keep_all = TRUE) %>% 
  summarise(max.age = max(age), number_of_variable_sites = n()) %>% 
  ggplot(aes(as.numeric(max.age), number_of_variable_sites, fill = Lineage)) +
  geom_point(size = 6, shape = 21, alpha = 0.9) +
  theme_minimal() +
  scale_x_continuous(limits = c(55, 75), breaks = seq(55, 75, 5), expand = c(0,0)) +
  theme(legend.position = "bottom", 
        legend.text = element_text(size =12), legend.title = element_blank(),
        axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid"),
                                        plot.title = element_text(size = 20)) +
  labs(x = "Years of follow up", y = "Number of variable sites")

ggsave("total_number_of_var_pos.svg")

### Making the linear regression of number of variable sites compared to reference in time
all_sel_TPs %>% 
  group_by(global_year_of_follow, Lineage) %>% 
  count(TYPE) %>% 
  ggplot(aes(global_year_of_follow, n)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ Lineage) +
  labs(y = "Number of substitutions", x = "Years of follow up") +
  theme_bw() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14), 
        panel.background = element_rect(fill = "white"))

ggsave("Figure_linreg.png")

### Making the same graph, but with age
all_sel_TPs$age <- as.numeric(all_sel_TPs$age)

all_sel_TPs %>% 
  group_by(age, Lineage) %>% 
  count(TYPE) %>% 
  ggplot(aes(age, n)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(method = "spearman", cor.coef.name = "rho", label.sep = "\n") +
  facet_wrap(~ Lineage) +
  labs(y = "Number of substitutions", x = "Age of the subject") +
  scale_x_continuous(breaks = c(35, 40, 45, 50, 55, 60, 65, 70, 75), limits = c(35, 75)) +
  theme_bw() +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 14), 
        panel.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12)) 

ggsave("Figure_linreg_age_R2.png")
ggsave("Figure_linreg_age_R2.svg")

### Correlation test- linear model using lm
all_sel_TPs_count <- all_sel_TPs %>% 
  group_by(Lineage, age) %>% 
  count(TYPE) %>% 
  rename(number_of_mutations = n)

mdl_TP_mut <- lm(data = all_sel_TPs_count, number_of_mutations ~ age * Lineage)
anova(mdl_TP_mut) ### generally there's a significant difference between the Lineages mutation rate.

pairwise.t.test(all_sel_TPs_count$number_of_mutations, all_sel_TPs_count$Lineage, p.adjust.method = "bonferroni")

## comparing the correlations (additionally)
pairs(emtrends(mdl_TP_mut, ~Lineage, var = "age"))

### Correlation test - Spearman
spearman_correlation_TP <- function(sel_lineage, cor.method) {
  selection <- all_sel_TPs %>% 
    filter(Lineage == sel_lineage) %>% 
    group_by(age) %>% 
    count(TYPE)
  
  cor.test(selection$n, selection$age, method = cor.method) %>% 
  tidy() %>% 
  cbind(tibble(Lineage = sel_lineage))
  
}

correlation_variable_sites_in_time <- rbind(
spearman_correlation_TP(sel_lineage = "TTV-AMS-S1-01", cor.method = "spearman"),
spearman_correlation_TP(sel_lineage = "TTV-AMS-S1-23", cor.method = "spearman"),
spearman_correlation_TP(sel_lineage = "TTMV-AMS-S1-41", cor.method = "spearman"),
spearman_correlation_TP(sel_lineage = "TTV-AMS-S2-01", cor.method = "spearman"),
spearman_correlation_TP(sel_lineage = "TTV-AMS-S2-03", cor.method = "spearman"),
spearman_correlation_TP(sel_lineage = "TTV-AMS-S2-04", cor.method = "spearman")
) %>% 
  select(Lineage, estimate, p.value) %>% 
  mutate(p.value.adjusted = p.adjust(p.value, method = "bonferroni"))

write.csv(correlation_variable_sites_in_time, "correlation_variable_sites_in_time.csv")

## Calculation of number of substitutions per each of the time point
mut_per_TP <- all_sel_TPs %>% 
  group_by(Lineage, TP) %>% 
  count(TP) %>% 
  ungroup() %>% 
  rename(number_of_mutations = n)

## Here I am joining the RPM table with the mutation per time point table
RPM_mut <- RPM_filtered %>%
  group_by(Lineage, TP) %>% 
  inner_join(mut_per_TP, by = c("TP", "Lineage")) %>% 
  ungroup()

suppl_figure1A <- RPM_mut %>% ggplot(aes(log(RPM), number_of_mutations)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(y = "Number of substitutions", x = "Log RPM", title = "A") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18))

suppl_figure1B <- RPM_mut %>% ggplot(aes(log(copies_per_mL), number_of_mutations)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(y = "Number of substitutions", x = "Log copies per mL", title = "B") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18))

grid.arrange(ncol = 1, nrow = 2, suppl_figure1A, suppl_figure1B)
suppl_figure1 <- arrangeGrob(suppl_figure1A, suppl_figure1B, nrow = 2)
ggsave(file = "suppl_figure1.svg", suppl_figure1, width = 2200, height = 2400, units = "px")
ggsave(file = "suppl_figure1.png", suppl_figure1, width = 2200, height = 2400, units = "px")

cor.test(RPM_mut$copies_per_mL, RPM_mut$number_of_mutations, method = "spearman")
cor.test(RPM_mut$RPM, RPM_mut$number_of_mutations, method = "spearman")


### DOTGRAPHS
### Graphs mutations per position in time
selected_TPs <- all_sel_TPs %>% select(TP, Lineage) %>% group_by(Lineage) %>%
  distinct(TP, .keep_all=TRUE) %>% ungroup() %>% arrange(TP)

### TTV-AMS-S1-01
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-01") %>% 
  ggplot() +
  geom_point(aes(x = as.factor(POS), y = as.factor(global_year_of_follow), 
                 colour = change, alpha = AF), size = 3) +
  theme(axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 8),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Year of follow-up", title = "TTV-AMS-S1-01")

all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-01") %>% 
  ggplot(aes(AF)) +
  geom_histogram() +
  labs(title = "TTV-AMS-S1-01")
### Many mutations of lower than 0.25 frequency

### TTV-AMS-S1-23
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-23") %>% 
  ggplot() +
  geom_point(aes(x = as.factor(POS), y = as.factor(global_year_of_follow), 
                 colour = change, alpha = AF), size = 3) +
  theme(axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 8),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Year of follow-up", title = "TTV-AMS-S1-23")

all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-23") %>% 
  ggplot(aes(AF)) +
  geom_histogram() +
  labs(title = "TTMV-AMS-S1-23")
### Also many mutations of lower than 0.25 frequency

### TTMV-AMS-S1-41
all_sel_TPs %>% 
  filter(Lineage == "TTMV-AMS-S1-41") %>% 
  ggplot() +
  geom_point(aes(x = as.factor(POS), y = as.factor(global_year_of_follow), 
                 colour = change, alpha = AF), size = 4) +
  theme(axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 10),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Year of follow-up", title = "TTMV-AMS-S1-41")

all_sel_TPs %>% 
  filter(Lineage == "TTMV-AMS-S1-41") %>% 
  ggplot(aes(AF)) +
  geom_histogram() +
  labs(title = "TTMV-AMS-S1-41")
### Here I see mainly 100% AF

### TTV-AMS-S2-01
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-01") %>% 
  ggplot() +
  geom_point(aes(x = as.factor(POS), y = as.factor(global_year_of_follow), 
                 colour = change, alpha = AF), size = 3) +
  theme(axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 8),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Year of follow-up", title = "TTV-AMS-S2-01")

all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-01") %>% 
  ggplot(aes(AF)) +
  geom_histogram() +
  labs(title = "TTV-AMS-S2-01")
### Also many mutations of lower than 0.25 frequency, but the 0.9 - 1 AF mutations dominate

### TTV-AMS-S2-03
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-03") %>% 
  ggplot() +
  geom_point(aes(x = as.factor(POS), y = as.factor(global_year_of_follow), 
                 colour = change, alpha = AF), size = 3) +
  theme(axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 8),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Year of follow-up", title = "TTV-AMS-S2-03")

all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-03") %>% 
  ggplot(aes(AF)) +
  geom_histogram() +
  labs(title = "TTV-AMS-S2-03")
### Also many mutations of lower than 0.25 frequency, but the 0.9 - 1 AF mutations dominate

### TTV-AMS-S2-04
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-04") %>% 
  ggplot() +
  geom_point(aes(x = as.factor(POS), y = as.factor(global_year_of_follow), 
                 colour = change), size = 6) +
  theme(axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 10),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Year of follow-up", title = "TTV-AMS-S2-04")

all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-04") %>% 
  ggplot(aes(AF)) +
  geom_histogram() +
  labs(title = "TTV-AMS-S2-04")
### All mutations are high frequencies

###HVR
all_sel_TPs %>% 
  ggplot(aes(POS)) +
  geom_histogram(bins = 30, color = "grey") +
  scale_x_continuous(breaks = seq(0, 2500, 100), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 900), expand = c(0,0)) +
  theme(axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 10),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "ORF1 position (nt)", y = "Total number of mutations in all lineages")

ggsave("HVR.png")
ggsave("HVR.svg")

### DOT-GRAPH with no factor for the ORF1 position 

### Assigning colors to values
change_colors <- c("A C" = "#F8766D", "A G" = "#DB6D00", "A T" = "#CD9600", 
                   "C A" = "#7CAE00", "C G" = "#00BA38", "C T" = "#00BFC4",
                   "G A" = "#619CFF", "G C" = "#00B0F6", "G T" = "#006DDB",
                   "T A" = "#C77CFF", "T C" = "#F564E3", "T G" = "#FF61CC")

legend_colors <- tibble(
  change = c("A C", "A G", "A T", "C A", "C G", "C T", 
             "G A", "G C", "G T", "T A", "T C", "T G"),
  values = seq(1, 1, 1)
)

ggplot(legend_colors, aes(change, values, fill = change)) +
  geom_tile(color = "white") + coord_equal() +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.title = element_blank(), axis.text.x = element_text(size = 16)) +
  labs(title = "Nucleotide change")  +
  scale_fill_manual(values = c("#F8766D", "#DB6D00", "#CD9600", 
                                "#7CAE00", "#00BA38", "#00BFC4",
                                "#619CFF", "#00B0F6", "#006DDB",
                                "#C77CFF", "#F564E3", "#FF61CC"))

ggsave("legend.svg")


### Graph type of mutation
all_sel_TPs %>% 
  group_by(Lineage, change) %>% 
  summarise(number_of_mutations = n()) %>%
  arrange(desc(number_of_mutations)) %>% 
  ggplot() +
  geom_bar(aes(x = Lineage, y = number_of_mutations, fill = change), position = "fill", stat = "identity", color = "grey") +
  theme_minimal() +
  scale_fill_manual(values = change_colors)+
  labs(x = "Lineage", y = "Proportion") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))

ggsave("suppl_figure2_types_mut.png")
ggsave("suppl_figure2_types_mut.svg")

### TTV-AMS-S1-01
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-01") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(100, 2300), breaks = seq(100, 2300, 100)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  geom_vline(xintercept = 800) +
  geom_vline(xintercept = 1300) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject", title = "TTV-AMS-S1-01") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S1-01_ORF1.png")
ggsave("TTV-AMS-S1-01_ORF1.svg")

### ZOOMIN
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-01") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(800, 1300), breaks = seq(800, 1300, 25), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S1-01_ORF1_zoom.png")
ggsave("TTV-AMS-S1-01_ORF1_zoom.svg")

### TTV-AMS-S1-23
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-23") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(100, 2300), breaks = seq(100, 2300, 100)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  geom_vline(xintercept = 800) +
  geom_vline(xintercept = 1300) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject", title = "TTV-AMS-S1-23") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S1-23_ORF1.png")
ggsave("TTV-AMS-S1-23_ORF1.svg")

### ZOOMIN
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S1-23") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(800, 1300), breaks = seq(800, 1300, 25), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S1-23_ORF1_zoom.png")
ggsave("TTV-AMS-S1-23_ORF1_zoom.svg")

### TTMV-AMS-S1-41
all_sel_TPs %>% 
  filter(Lineage == "TTMV-AMS-S1-41") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(49, 1600), breaks = seq(50, 1600, 50)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  geom_vline(xintercept = 800) +
  geom_vline(xintercept = 1300) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject", title = "TTMV-AMS-S1-41") +
  scale_color_manual(values = change_colors)

ggsave("TTMV-AMS-S1-41_ORF1.png")
ggsave("TTMV-AMS-S1-41_ORF1.svg")

### ZOOMIN (not necessary)
all_sel_TPs %>% 
  filter(Lineage == "TTMV-AMS-S1-41") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(800, 1300), breaks = seq(800, 1300, 25), expand = c(0,0)) +
  theme(legend.position = "none", 
        axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 12),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Age of the subject") +
  scale_color_manual(values = change_colors)

ggsave("TTMV-AMS-S1-41_ORF1_zoom.png")
ggsave("TTMV-AMS-S1-41_ORF1_zoom.svg")

### TTV-AMS-S2-01
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-01") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(50, 2300), breaks = seq(100, 2300, 100)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  geom_vline(xintercept = 800) +
  geom_vline(xintercept = 1300) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject", title = "TTV-AMS-S2-01") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S2-01_ORF1.png")
ggsave("TTV-AMS-S2-01_ORF1.svg")

### ZOOMIN
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-01") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(800, 1300), breaks = seq(800, 1300, 25), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S2-01_ORF1_zoom.png")
ggsave("TTV-AMS-S2-01_ORF1_zoom.svg")

### TTV-AMS-S2-03
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-03") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(50, 2300), breaks = seq(100, 2300, 100)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  geom_vline(xintercept = 800) +
  geom_vline(xintercept = 1300) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject", title = "TTV-AMS-S2-03") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S2-03_ORF1.png")
ggsave("TTV-AMS-S2-03_ORF1.svg")

### ZOOMIN
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-03") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change, alpha = AF), size = 3) +
  scale_x_continuous(limits = c(800, 1300), breaks = seq(800, 1300, 25), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S2-03_ORF1_zoom.png")
ggsave("TTV-AMS-S2-03_ORF1_zoom.svg")

### TTV-AMS-S2-04
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-04") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change), size = 3) +
  scale_x_continuous(limits = c(50, 2300), breaks = seq(100, 2300, 100)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 15), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 16),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  geom_vline(xintercept = 800) +
  geom_vline(xintercept = 1300) +
  labs(x = "ORF1 position (nt)", y = "Age of the subject", title = "TTV-AMS-S2-04") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S2-04_ORF1.png")
ggsave("TTV-AMS-S2-04_ORF1.svg")

### ZOOMIN
all_sel_TPs %>% 
  filter(Lineage == "TTV-AMS-S2-04") %>% 
  ggplot() +
  geom_point(aes(x = POS, y = as.factor(age), 
                 colour = change), size = 3) +
  scale_x_continuous(limits = c(800, 1300), breaks = seq(800, 1300, 25), expand = c(0,0)) +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0, size = 12), 
        axis.text.x = element_text(vjust = 0.5, angle = 90, size = 12),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill="transparent", colour = "grey50"),
        panel.grid.major = element_line(colour = "lightgray", linetype = "solid")) +
  labs(x = "Position (nt)", y = "Age of the subject") +
  scale_color_manual(values = change_colors)

ggsave("TTV-AMS-S2-04_ORF1_zoom.png")
ggsave("TTV-AMS-S2-04_ORF1_zoom.svg")

### Figures HVR no. mutations vs rest of ORF1

RPM_filtered %>%
  filter(TP < 56) %>% 
  ggplot(aes(TP, log(copies_per_mL), color = Lineage)) +
  scale_x_continuous(breaks = seq(0, 55, 5)) +
  geom_point(size = 5) +
  geom_label(aes(label = TP))

RPM_filtered %>%
  filter(TP < 56) %>% 
  ggplot(aes(TP, log(RPM), color = Lineage)) +
  scale_x_continuous(breaks = seq(0, 55, 5)) +
  geom_point(size = 5) +
  geom_label(aes(label = TP))
### Take TP52 for TTMV-S1-41 and TTV-S1-23, and 54 for TTV-AMS-S1-01

RPM_filtered %>%
  filter(TP >= 56) %>% 
  ggplot(aes(TP, log(copies_per_mL), color = Lineage)) +
  scale_x_continuous(breaks = seq(55, 110, 5)) +
  geom_point(size = 5) +
  geom_label(aes(label = TP))

RPM_filtered %>%
  filter(TP >= 56) %>% 
  ggplot(aes(TP, log(RPM), color = Lineage)) +
  scale_x_continuous(breaks = seq(55, 110, 5)) +
  geom_point(size = 5) +
  geom_label(aes(label = TP))
### Take 107 for TTV-AMS-S2-03, 106 for TTV-AMS-S2-03 and 92 for TTV-AMS-S2-04

### OR... do them all, it's not so much
HVR_part <- all_sel_TPs %>% 
  filter(POS >= 800 | POS <= 1300) %>% 
  group_by(Lineage, lineage_year_of_follow) %>%
  summarize(number_of_mutations = n())

non_HVR_part <- all_sel_TPs %>% 
  filter(POS < 800 | POS > 1300) %>% 
  group_by(Lineage, lineage_year_of_follow) %>%
  summarize(number_of_mutations = n())

lengths_of_ORF1 <- data.frame(
  Lineage = c("TTV-AMS-S1-01", "TTV-AMS-S1-23", "TTMV-AMS-S1-41",
              "TTV-AMS-S2-01", "TTV-AMS-S2-03", "TTV-AMS-S2-04"),
  ORF1_length = c(2367, 2319, 1551,
                  2210, 2226, 2289),
  length_HVR = c(500, 500, 500,
                 500, 500, 500)

  )

lengths_of_ORF1 <- lengths_of_ORF1 %>%
  mutate(length_no_HVR = ORF1_length - length_HVR)

HVR_part <- HVR_part %>% inner_join(lengths_of_ORF1, by = "Lineage") 

HVR_part <- HVR_part %>%
  filter(lineage_year_of_follow > 0) %>% 
  mutate(mutation_per_site_per_year = number_of_mutations / length_HVR / lineage_year_of_follow,
         site = "EHVR")
  
non_HVR_part <- non_HVR_part %>% inner_join(lengths_of_ORF1, by = "Lineage") 

non_HVR_part <- non_HVR_part %>%
  filter(lineage_year_of_follow > 0) %>% 
  mutate(mutation_per_site_per_year = number_of_mutations / length_no_HVR / lineage_year_of_follow, 
         site = "non_EHVR")

HVR_and_non_HVR <- rbind(HVR_part, non_HVR_part)

HVR_and_non_HVR %>% 
  ggplot(aes(site, log(mutation_per_site_per_year), fill = site)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.title = element_text(size = 14)) +
  labs(x = "Site", y = "log of number of mutations per site per year")


final_S1_01_selected_TP <- all_sel_TPs %>% filter(Lineage == "TTV-AMS-S1-01") %>% select(TP) %>% unique()
final_S1_23_selected_TP <- all_sel_TPs %>% filter(Lineage == "TTV-AMS-S1-23") %>% select(TP) %>% unique()
final_S1_41_selected_TP <- all_sel_TPs %>% filter(Lineage == "TTMV-AMS-S1-41") %>% select(TP) %>% unique()
final_S2_01_selected_TP <- all_sel_TPs %>% filter(Lineage == "TTV-AMS-S2-01") %>% select(TP) %>% unique()
final_S2_03_selected_TP <- all_sel_TPs %>% filter(Lineage == "TTV-AMS-S2-03") %>% select(TP) %>% unique()
final_S2_04_selected_TP <- all_sel_TPs %>% filter(Lineage == "TTV-AMS-S2-04") %>% select(TP) %>% unique()

### Which frequencies are the most common among the mutations at the same position?

more_than_1_mut <- all_sel_TPs %>% 
  group_by(Lineage, POS, TP) %>% 
  summarise(total_mutations_at_the_pos = n()) %>% 
  filter(total_mutations_at_the_pos > 1) %>% 
  ungroup()

more_than_1_mut_selection <- all_sel_TPs %>% 
  semi_join(more_than_1_mut)

more_than_1_mut_selection %>% 
  ggplot(aes(AF)) + 
  geom_histogram(bins = 10) +
  geom_vline(xintercept = 0.33)

all_sel_TPs %>% 
  filter(AF > 0.33)

### Boxplot synonymous and non-synonymous substitutions in HVR and outside HVR

changes_in_genomes <- read_excel("changes_in_genomes.xlsx")
changes_in_genomes <- changes_in_genomes %>% inner_join(lengths_of_ORF1, by = "Lineage")

### TTV-AMS-S1-01

changes_in_genomes_HVR_S1_01 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_EHVR, s_EHVR, length_HVR) %>% 
  filter(Lineage == "TTV-AMS-S1-01") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_01, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_EHVR / length_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_EHVR / length_HVR / lineage_year_of_follow,
         site = "EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

changes_in_genomes_no_HVR_S1_01 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_rest, s_rest, length_no_HVR) %>% 
  filter(Lineage == "TTV-AMS-S1-01") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_01, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_rest / length_no_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_rest / length_no_HVR / lineage_year_of_follow,
         site = "non_EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

### TTV-AMS-S1-23

changes_in_genomes_HVR_S1_23 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_EHVR, s_EHVR, length_HVR) %>% 
  filter(Lineage == "TTV-AMS-S1-23") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_23, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_EHVR / length_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_EHVR / length_HVR / lineage_year_of_follow,
         site = "EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

changes_in_genomes_no_HVR_S1_23 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_rest, s_rest, length_no_HVR) %>% 
  filter(Lineage == "TTV-AMS-S1-23") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_23, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_rest / length_no_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_rest / length_no_HVR / lineage_year_of_follow,
         site = "non_EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

### TTMV-AMS-S1-41

changes_in_genomes_HVR_S1_41 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_EHVR, s_EHVR, length_HVR) %>% 
  filter(Lineage == "TTMV-AMS-S1-41") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_41, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_EHVR / length_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_EHVR / length_HVR / lineage_year_of_follow,
         site = "EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

changes_in_genomes_no_HVR_S1_41 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_rest, s_rest, length_no_HVR) %>% 
  filter(Lineage == "TTMV-AMS-S1-41") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S1_41, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_rest / length_no_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_rest / length_no_HVR / lineage_year_of_follow,
         site = "non_EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)


### TTV-AMS-S2-01

changes_in_genomes_HVR_S2_01 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_EHVR, s_EHVR, length_HVR) %>% 
  filter(Lineage == "TTV-AMS-S2-01") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_01, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_EHVR / length_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_EHVR / length_HVR / lineage_year_of_follow,
         site = "EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

changes_in_genomes_no_HVR_S2_01 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_rest, s_rest, length_no_HVR) %>% 
  filter(Lineage == "TTV-AMS-S2-01") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_01, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_rest / length_no_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_rest / length_no_HVR / lineage_year_of_follow,
         site = "non_EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

### TTV-AMS-S2-03

changes_in_genomes_HVR_S2_03 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_EHVR, s_EHVR, length_HVR) %>% 
  filter(Lineage == "TTV-AMS-S2-03") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_03, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_EHVR / length_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_EHVR / length_HVR / lineage_year_of_follow,
         site = "EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

changes_in_genomes_no_HVR_S2_03 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_rest, s_rest, length_no_HVR) %>% 
  filter(Lineage == "TTV-AMS-S2-03") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_03, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_rest / length_no_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_rest / length_no_HVR / lineage_year_of_follow,
         site = "non_EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

### TTV-AMS-S2-04

changes_in_genomes_HVR_S2_04 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_EHVR, s_EHVR, length_HVR) %>% 
  filter(Lineage == "TTV-AMS-S2-04") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_04, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_EHVR / length_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_EHVR / length_HVR / lineage_year_of_follow,
         site = "EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)

changes_in_genomes_no_HVR_S2_04 <- changes_in_genomes %>% 
  select(Lineage, TP, ns_rest, s_rest, length_no_HVR) %>% 
  filter(Lineage == "TTV-AMS-S2-04") %>% 
  inner_join(years_TP, by = "TP") %>% 
  mutate(year = ymd(Date_of_collection), day_of_follow = year - ref_S2_04, 
         lineage_year_of_follow = as.numeric(day_of_follow / 365), 
         s_per_site_per_year = s_rest / length_no_HVR / lineage_year_of_follow,
         ns_per_site_per_year = ns_rest / length_no_HVR / lineage_year_of_follow,
         site = "non_EHVR") %>% 
  select(Lineage, TP, s_per_site_per_year, ns_per_site_per_year, site)


changes_in_genomes_HVR_and_non_HVR <- rbind(changes_in_genomes_HVR_S1_01, changes_in_genomes_no_HVR_S1_01,
                                            changes_in_genomes_HVR_S1_23, changes_in_genomes_no_HVR_S1_23,
                                            changes_in_genomes_HVR_S1_41, changes_in_genomes_no_HVR_S1_41,
                                            changes_in_genomes_HVR_S2_01, changes_in_genomes_no_HVR_S2_01,
                                            changes_in_genomes_HVR_S2_03, changes_in_genomes_no_HVR_S2_03,
                                            changes_in_genomes_HVR_S2_04, changes_in_genomes_no_HVR_S2_04)


HVR_and_non_HVR %>% 
  ggplot(aes(site, log(mutation_per_site_per_year), fill = site)) +
  geom_boxplot() +
  facet_wrap(~Lineage) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12)) +
  labs(x = "Site", y = "log of number of mutations per site per year") +
  scale_fill_manual(values = c("#8C8C8C", "#BFBFBF"))+
  stat_compare_means(aes(site, log(mutation_per_site_per_year)))
ggsave("mut_per_site_HVR_and_non_HVR.png")
ggsave("mut_per_site_HVR_and_non_HVR.svg")

ggplot(changes_in_genomes_HVR_and_non_HVR) +
  geom_boxplot(aes(site, log(s_per_site_per_year), fill = site)) +
  facet_wrap(~ Lineage) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12)) +
  labs(x = "Site", y = "log of number of synonymous mutations per site per year")+
  scale_fill_manual(values = c("#8C8C8C", "#BFBFBF")) +
  stat_compare_means(aes(site, log(s_per_site_per_year)))
ggsave("s_per_site_year_HVR_and_non_HVR.png")
ggsave("s_per_site_year_HVR_and_non_HVR.svg")


ggplot(changes_in_genomes_HVR_and_non_HVR) +
  geom_boxplot(aes(site, log(ns_per_site_per_year), fill = site))  +
  facet_wrap(~ Lineage) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12)) +
  labs(x = "Site", y = "log of number of non-synonymous mutations per site per year")+
  scale_fill_manual(values = c("#8C8C8C", "#BFBFBF")) +
  stat_compare_means(aes(site, log(ns_per_site_per_year)))
ggsave("ns_per_site_year_HVR_and_non_HVR.png")
ggsave("ns_per_site_year_HVR_and_non_HVR.svg")

all_sel_TPs %>% filter(Lineage == "TTV-AMS-S1-01") %>% count(TP, sort = TRUE)

### GRAPH WITH THE ORF1 DOMAINS
selection_domains <- read_excel("orf1_picture.xlsx")
domains <- read_excel("orf1_picture.xlsx", sheet = 2)

domains %>% ggplot() +
  geom_col(aes(length, Lineage, fill = region), color = "black") +
  geom_point(data = selection_domains, mapping = aes(Codon, Lineage, color = Selection), shape = "|", size = 20) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 800, 100)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14)) +
  labs(x = "ORF1 position (amino acids)", y = "") +
  scale_fill_manual(values = c("#AEEEEE","#FF8C00","#3FA0FF","#458B00","#3FA0FF","#FF8C00","#B66DFF")) +
  scale_color_manual(values = c("#7CCD7C", "#CD3700")) +
  geom_vline(xintercept = 265, size = 1, linetype = "dashed") +
  geom_vline(xintercept = 433, size = 1, linetype = "dashed")

ggsave("domains.svg")

### Legend for this graph

### Assigning colors to values
domain_colors <- c("C-Term" = "#AEEEEE", "Jelly-Roll" = "#FF8C00",
                  "P1" = "#3FA0FF", "P2" = "#458B00",
                  "Arginine Rich" = "#B66DFF")

legend_domain_colors <- tibble(
  domain = c("C-Term", "Jelly-Roll", "P1", "P2", "Arginine Rich"),
  values = seq(1, 1, 1)
)

ggplot(legend_domain_colors, aes(values, domain, fill = domain)) +
  geom_tile(color = "white") + coord_equal() +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.title = element_blank(), axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 20)) +
  labs(title = "ORF1 domains")  +
  scale_fill_manual(values = domain_colors)

ggsave("legend_domains.svg")


# Title: Reproducible Analysis of the Scientific Article "Compensatory Proximal Adjustments Characterize Effective Reaching Movements After Stroke"
# Author: Silke Wolf
# Contact: si.wolf@uke.de
# Date: June 17, 2025
# Repository: https://github.com/xeni-lab/reaching_movements_after_stroke

#
# Description:
# This R script contains the data analysis, visualizations, and results related  
# to the scientific article "Compensatory Proximal Adjustments Characterize Effective Reaching Movements After Stroke". 
# It aims to make the analyses and plots reproducible and transparent, allowing others to verify and build upon 
# the work.


rm(list=ls())

# load packages ####

library(readr)
library(tidyr)
library(dplyr)
library(lme4)
library(lmtest)
library(lmerTest)
library(stats)
library(ggplot2)
library(ggpubr)
library(psych)
library(compareGroups)
library(gridExtra)


# Load data (insert path, where rds files are stored) ####

# * Main data ####
data <- read_rds(".../1_data.rds")


# * Additional data for plots 1a, 2a and 2d ####

# path length two example trials
plot1a_dat <- readRDS(".../2_plot1a.rds")

# peak velocity of hand for two example trials
plot1a_pv <- readRDS(".../3_plot1a_pv.rds")

# for labeling
subjects <- data.frame(subject = c("control18", "patient13"), label = c("Control", "Stroke"))

# joint angle data for two example trials
plot2a_dat <- read_rds(".../4_plot2a.rds")

# peak angular velocity elbow flexion for two example trials
plot2a_pv_el <- read_rds(".../5_plot2a_pv_el.rds")

# peak angular velocity shoulder rotation for two example trials
plot2a_pv_shz <- read_rds(".../6_plot2a_pv_shz.rds")

# # joint angle data for two example trials
plot2d_dat <- read_rds(".../7_plot2d.rds")

# peak angular velocity shoulder flexion for two example trials
plot2d_pv_shx <- read_rds(".../8_plot2d_shx.rds")

# peak angular velocity shoulder abduction for two example trials
plot2d_pv_shy <- read_rds(".../9_plot2d_shy.rds")

# A) Description of variables in 'data' ####

# subject = id
# condition = control or patient
# age = age in years
# sex = 0 = male, 1 = female

# nihss_total = sum score NIHSS at inclusion
# ratio_mean_FA = left FA mean / right FA mean
# arat_total_affected = ARAT score affected side
# ashw_mean_affected = Ashworth spasticity - averaged for the affected limb
# bbt_te_rel = BBT - ratio of blocks moved by the more affected hand to the less affected hand
# fma_total = Fugl-Meyer score affected side 
# rel_nhpt_affected = Ratio of Nine Hole Peg Test measured in pegs/second

# block = block 1 - 8 of experiment
# trial = trial 1 - 20 per block
# side = the side which is moved (right or left)
# trial_number = consecutive trial number, entire experiment (max. 160)

# move_time = movement time in seconds 
# time_peak_vel_hand_in_sec = time to peak velocity (TTPVHD) for each individual trial of the hand, in seconds
# rel_time_peak_vel_hand = time to peak velocity (TTPVHD) for each individual trial of the hand, expressed as a percentage from 0 to 100%
# time_peak_acc_hand_in_sec = time to peak acceleration (TTPAHD) for each individual trial of the hand, in seconds
# rel_time_peak_acc_hand = time to peak acceleration (TTPAHD) for each individual trial of the hand, expressed as a percentage from 0 to 100%
# number_vel_peaks = number of velocity peaks per trial = smoothness

# rel_time_pv_sh_x = time to peak angular velocity shoulder (TTPVSH) flexion/extension as a percentage from 0% to 100%
# rel_time_pv_sh_y = time to peak angular velocity shoulder (TTPVSH) abd-/adduction as a percentage from 0% to 100%
# rel_time_pv_sh_z_neg = time to peak angular velocity shoulder (TTPVSH) internal/external rotation as a percentage from 0% to 100%
# rel_time_pv_el_x = time to peak angular velocity elbow (TTPVEL) flexion/extension as a percentage from 0% to 100%

# sh_x_el_x = relative time interval TTPVEL and TTPVSH flexion 
# sh_y_el_x = relative time interval TTPVEL and TTPVSH abduction 
# sh_z_el_x = relative time interval TTPVEL and TTPVSH rotation
# sh_x_sh_y = relative time interval TTPVSH flexion and TTPVSH abduction
# sh_x_sh_z = relative time interval TTPVSH flexion and TTPVSH rotation
# sh_y_sh_z = relative time interval TTPVSH abduction and TTPVSH rotation

# TAS = time after stroke (months)
# lesion_location = lesion location


# *** Table 1 #####

table1 <- data %>% 
  group_by(subject) %>% 
  mutate(sex = as.factor(sex)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(age, sex, TAS, lesion_location, ratio_mean_FA, nihss_total, fma_total, arat_total_affected, rel_nhpt_affected, bbt_te_rel, ashw_mean_affected, condition)

table1 %>% 
  print(n=26)

summary <- compareGroups::compareGroups(condition ~ ., data = table1, 
                                    method = c(fma_total = 2, arat_total_affected = 2), p.corrected = T, Q1 = 0, Q3 = 1)

summary_tab <- compareGroups::createTable(summary, type = 1, sd.type = 2, show.n = T, show.all =  T, digits.p = 3)
summary_tab

# *** Table 2 ####

# skim relevant variables 

data %>% 
  filter(side == "right") %>% 
  group_by(condition) %>% 
  select(move_time, rel_time_peak_vel_hand, rel_time_peak_acc_hand, number_vel_peaks) %>% 
  skimr::skim()

data %>% 
  filter(side == "right") %>% 
  group_by(condition) %>% 
  select(rel_time_pv_el_x, rel_time_pv_sh_z_neg, rel_time_pv_sh_x, rel_time_pv_sh_y) %>% 
  skimr::skim()

data %>% 
  filter(side == "right") %>% 
  group_by(condition) %>% 
  select(sh_x_el_x, sh_y_el_x, sh_z_el_x) %>% 
  skimr::skim()

data %>% 
  filter(side == "right") %>% 
  group_by(condition) %>% 
  select(sh_x_sh_y, sh_x_sh_z, sh_y_sh_z) %>% 
  skimr::skim()


# B) Models for Table 3 ####

# *** Endpoint variables ####

# 1. MOVEMENT TIME ####

# 1.1 simple linear regression ####

#1
dur_m <- lm(move_time ~ condition + side, data)
summary(dur_m)

# 1.2 full mixed-effect model ####

#2
dur_m2 <- lmer(move_time ~ side + trial_number + condition + (1 | subject) , data)
summary(dur_m2)

dur_m2_no_stroke <- lmer(move_time ~ side + trial_number + (1 | subject) , data)
summary(dur_m2_no_stroke)

dur_comp <- lmtest::lrtest(dur_m2_no_stroke, dur_m2)

# 2. TTP VELOCITY HAND ####

# 2.1 simple linear regression ####

#1
ttpv_m <- lm(rel_time_peak_vel_hand ~ condition + side, data)
summary(ttpv_m)

# 2.2 full mixed-effect model ####

#2
ttpv_m2 <- lmer(rel_time_peak_vel_hand ~ side + trial_number + condition + (1 | subject) , data)
summary(ttpv_m2)

ttpv_m2_no_stroke <- lmer(rel_time_peak_vel_hand ~ side + trial_number + (1 | subject) , data)
summary(ttpv_m2_no_stroke)

ttpv_comp <- lmtest::lrtest(ttpv_m2_no_stroke, ttpv_m2)


# 3. TTP ACCELERATION HAND ####

# 3.1 simple linear regression ####

#1
ttpa_m <- lm(rel_time_peak_acc_hand ~ condition + side, data)
summary(ttpa_m)

# 3.2 full mixed-effect model ####

#2
ttpa_m2 <- lmer(rel_time_peak_acc_hand ~ side + trial_number + condition + (1 | subject) , data)
summary(ttpa_m2)

ttpa_m2_no_stroke <- lmer(rel_time_peak_acc_hand ~ side + trial_number + (1 | subject) , data)
summary(ttpa_m2_no_stroke)

ttpa_comp <- lmtest::lrtest(ttpa_m2_no_stroke, ttpa_m2)

# 4. SMOOTHNESS ####

# 4.1 simple linear regression ####

#1
peaks_m <- lm(number_vel_peaks ~ condition + side, data)
summary(peaks_m)

# 4.2 full mixed-effect model ####

#2
peaks_m2 <- lmer(number_vel_peaks ~ side + trial_number + condition + (1 | subject) , data)
summary(peaks_m2)

peaks_m2_no_stroke <- lmer(number_vel_peaks ~ side + trial_number + (1 | subject) , data)
summary(peaks_m2_no_stroke)

peaks_comp <- lmtest::lrtest(peaks_m2_no_stroke, peaks_m2)


# Compile results full model ####

# Initialize a list to store results
results_list <- list()

# Define a function to extract and round the needed values
extract_and_round_values <- function(model, effect_name, digits=3) {
  summary_model <- summary(model)
  coef_values <- summary_model$coefficients
  effect_values <- coef_values[effect_name, ]
  rounded_values <- round(effect_values, digits)
  return(rounded_values)
}

# Extract and round values for each model and store them in the list
results_list[["dur_m2"]] <- extract_and_round_values(dur_m2, "conditionpatient")
results_list[["ttpv_m2"]] <- extract_and_round_values(ttpv_m2, "conditionpatient")
results_list[["ttpa_m2"]] <- extract_and_round_values(ttpa_m2, "conditionpatient")
results_list[["peaks_m2"]] <- extract_and_round_values(peaks_m2, "conditionpatient")

# Convert list to a data frame
results_df <- do.call(rbind, results_list)

# Rename the rows for clarity
rownames(results_df) <- c("dur_m2", "ttpv_m2", "ttpa_m2", "peaks_m2")

# Remove the "df" column if it exists
if("df" %in% colnames(results_df)) {
  results_df <- results_df[, !colnames(results_df) %in% "df"]
}

# Display the compiled results
print(results_df)


# Compile results model comparison ####


# Extract and compile the results of each comparison
extract_comparison <- function(comp) {
  # Extract chi-square statistic, degrees of freedom, and p-value
  chi_sq <- round(comp$`Chisq`[2], 3)
  df <- comp$Df[2]
  p_value <- round(comp$`Pr(>Chisq)`[2], 3)
  
  results <- c(Chi_Square = chi_sq, DF = df, P_Value = p_value
  )
  
  return(results)
}

# A list to store comparison results
comp_results <- list()

# Extract and store the results for each comparison
comp_results[["dur_comp"]] <- extract_comparison(dur_comp)
comp_results[["ttpv_comp"]] <- extract_comparison(ttpv_comp)
comp_results[["ttpa_comp"]] <- extract_comparison(ttpa_comp)
comp_results[["peaks_comp"]] <- extract_comparison(peaks_comp)


# Convert the list to a data frame
comp_results_df <- do.call(rbind, comp_results)

# Rename the rows for clarity
rownames(comp_results_df) <- c("dur_comp", "ttpv_comp", "ttpa_comp", "peaks_comp")

# Display the compiled results
print(comp_results_df)

# *** Inter-joint variables ####


# 1. SHOULDER FLEXION - ELBOW FLEXION ####

# 1.1 simple linear regression ####

#1
sx_ex <- lm(sh_x_el_x ~ condition + side, data)
summary(sx_ex)

# 1.2 full mixed-effect model ####

#2
sx_ex_m2 <- lmer(sh_x_el_x ~ condition + side + trial_number + (1 | subject) , data)
summary(sx_ex_m2)

sx_ex_m2_no_stroke <- lmer(sh_x_el_x ~  side + trial_number + (1 | subject) , data)
summary(sx_ex_m2_no_stroke)

sxex_com <- lmtest::lrtest(sx_ex_m2_no_stroke, sx_ex_m2)

# 2. SHOULDER ABDUCTION - ELBOW FLEXION ####

# 2.1 simple linear regression ####

#1
sy_ex <- lm(sh_y_el_x ~ condition + side, data)
summary(sy_ex)

# 2.2 full mixed-effect model ####

#2
sy_ex_m2 <- lmer(sh_y_el_x ~ condition + side + trial_number + (1 | subject) , data)
summary(sy_ex_m2)

sy_ex_m2_no_stroke <- lmer(sh_y_el_x ~  side + trial_number + (1 | subject) , data)
summary(sy_ex_m2_no_stroke)

syex_com <- lmtest::lrtest(sy_ex_m2_no_stroke, sy_ex_m2)

# 3. SHOULDER ROTATION - ELBOW FLEXION ####

# 3.1 simple linear regression ####

#1
sz_ex <- lm(sh_z_el_x ~ condition + side, data)
summary(sz_ex)

# 3.2 full mixed-effect model ####

#2
sz_ex_m2 <- lmer(sh_z_el_x ~ condition + side + trial_number + (1 | subject) , data)
summary(sz_ex_m2)

sz_ex_m2_no_stroke <- lmer(sh_z_el_x ~  side + trial_number + (1 | subject) , data)
summary(sz_ex_m2_no_stroke)

szex_com <- lmtest::lrtest(sz_ex_m2_no_stroke, sz_ex_m2)

# *** Intra-joint variables ####

# 1. SHOULDER FLEXION - SHOULDER ABDUCTION ####

# 1.1 simple linear regression ####

#1
sx_sy <- lm(sh_x_sh_y ~ condition + side, data)
summary(sx_sy)

# 1.2 full mixed-effect model ####

#2
sx_sy_m2 <- lmer(sh_x_sh_y ~ condition + side + trial_number + (1 | subject) , data)
summary(sx_sy_m2)

sx_sy_m2_no_stroke <- lmer(sh_x_sh_y ~  side + trial_number + (1 | subject) , data)
summary(sx_sy_m2_no_stroke)

sxsy_com <- lmtest::lrtest(sx_sy_m2_no_stroke, sx_sy_m2)


# 2. SHOULDER FLEXION - SHOULDER ROTATION ####

# 2.1 simple linear regression ####

#1
sx_sz <- lm(sh_x_sh_z ~ condition + side, data)
summary(sx_sz)

# 2.2 full mixed-effect model ####

#2
sx_sz_m2 <- lmer(sh_x_sh_z ~ condition + side + trial_number + (1 | subject) , data)
summary(sx_sz_m2)

sx_sz_m2_no_stroke <- lmer(sh_x_sh_z ~  side + trial_number + (1 | subject) , data)
summary(sx_sz_m2_no_stroke)

sxsz_com <- lmtest::lrtest(sx_sz_m2_no_stroke, sx_sz_m2)


# 3. SHOULDER ABDUCTION - SHOULDER ROTATION ####

# 3.1 simple linear regression ####

#1
sy_sz <- lm(sh_y_sh_z ~ condition + side, data)
summary(sy_sz)

# 3.2 full mixed-effect model ####

#2
sy_sz_m2 <- lmer(sh_y_sh_z ~ condition + side + trial_number + (1 | subject) , data)
summary(sy_sz_m2)

sy_sz_m2_no_stroke <- lmer(sh_y_sh_z ~  side + trial_number + (1 | subject) , data)
summary(sy_sz_m2_no_stroke)

sysz_com <- lmtest::lrtest(sy_sz_m2_no_stroke, sy_sz_m2)



# Compile results full model ####

# Initialize a list to store results
results_list <- list()

# Extract and round values for each model and store them in the list
results_list[["sx_ex_m2"]] <- extract_and_round_values(sx_ex_m2, "conditionpatient")
results_list[["sy_ex_m2"]] <- extract_and_round_values(sy_ex_m2, "conditionpatient")
results_list[["sz_ex_m2"]] <- extract_and_round_values(sz_ex_m2, "conditionpatient")
results_list[["sx_sy_m2"]] <- extract_and_round_values(sx_sy_m2, "conditionpatient")
results_list[["sx_sz_m2"]] <- extract_and_round_values(sx_sz_m2, "conditionpatient")
results_list[["sy_sz_m2"]] <- extract_and_round_values(sy_sz_m2, "conditionpatient")


# Convert list to a data frame
results_df <- do.call(rbind, results_list)

# Rename the rows for clarity
rownames(results_df) <- c("sx_ex_m2", "sy_ex_m2", "sz_ex_m2",
                          "sx_sy_m2", "sx_sz_m2", "sy_sz_m2")

# Remove the "df" column if it exists
if("df" %in% colnames(results_df)) {
  results_df <- results_df[, !colnames(results_df) %in% "df"]
}

# Display the compiled results
print(results_df)


# Compile results model comparison ####

# A list to store comparison results
comp_results <- list()

# Extract and store the results for each comparison
comp_results[["sxex_com"]] <- extract_comparison(sxex_com)
comp_results[["syex_com"]] <- extract_comparison(syex_com)
comp_results[["szex_com"]] <- extract_comparison(szex_com)
comp_results[["sxsy_com"]] <- extract_comparison(sxsy_com)
comp_results[["sxsz_com"]] <- extract_comparison(sxsz_com)
comp_results[["sysz_com"]] <- extract_comparison(sysz_com)


# Convert the list to a data frame
comp_results_df <- do.call(rbind, comp_results)

# Rename the rows for clarity
rownames(comp_results_df) <- c("sxex_com", "syex_com", "szex_com", 
                               "sxsy_com", "sxsz_com", "sysz_com")

# Display the compiled results
print(comp_results_df)


# C) Plots ####

# plot 1 ####

# 1a - examples ####

plot1a <- plot1a_dat %>% 
  ggplot(aes(x = mov_time, y = dist_sum)) +
  geom_line(size=1) +
  geom_point(data = plot1a_pv, aes(x = mov_time, y = dist_sum, shape = "Peak velocity"), size = 4, color = "black") +
  theme_minimal() +
  xlim(0, 1.2) +
  scale_x_continuous(breaks = c(0, .5, 1)) + 
  xlab("Time (s)") + 
  ylab("3D hand path length (mm)") +
  scale_shape_manual(values = c("Peak velocity" = 17)) +
  scale_linetype_manual(values = c("solid", "solid")) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), size=16),
    axis.title.y = element_text(margin = margin(r = 10), size=16),
    axis.line = element_line(color = "black"),  # Set axes lines color
    axis.ticks = element_line(color = "black"),  # Set ticks color
    axis.text = element_text(color = "black", size=14),  # Set tick labels color
    panel.grid = element_blank(),  # Remove grid
    legend.text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  ) +
  theme(axis.title.x = element_text(margin = margin(t = 10)),  
        axis.title.y = element_text(margin = margin(r = 10))) +
  facet_grid(rows = vars(subject)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()
  ) +
  # Add annotations
  geom_text(data = subjects, mapping = aes(x = 0, y = 680, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6) +
  # Add custom legend
  guides(shape = guide_legend(title = "", override.aes = list(size = 4, color = "black"))) +
  theme(legend.position = c(0.9, 0.85))

#plot1a

# 1b - frequencies peak velocity ####

# Annotation data frames
control_annotation <- data.frame(condition = "control", label = "Control")
stroke_annotation <- data.frame(condition = "patient", label = "Stroke")

summary(data)

plot1b <- data %>% 
  filter(side == "right") %>% 
  ggplot(aes(x =  rel_time_peak_vel_hand)) +
  geom_density(aes(y = ..count..),
               fill = NA, 
               color = "black", 
               size = 1) +
  xlim(0, 100) +
  scale_y_continuous(breaks = c(0, 50, 100)) +  
  theme_minimal() +
  xlab("Relative time (%)") + 
  ylab("Number of trials at peak velocity time points") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), size=16),
    axis.title.y = element_text(margin = margin(r = 10), size=16),
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black", size=14),  
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")  
  ) +
  facet_grid(rows = vars(condition)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  # Add annotations for Control
  geom_text(data = control_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6) +
  # Add annotations for Stroke
  geom_text(data = stroke_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6)
#plot1b

# 1c - box plots ####

plot1c <- data %>% 
  mutate(condition = factor(condition, levels = c("patient", "control"))) %>% 
  filter(side == "right") %>% 
  ggplot(aes(y=condition, x=rel_time_peak_vel_hand, group=condition)) +
  geom_boxplot(width=.3, size=1, position = position_dodge(width = 2)) +
  xlim(0,100) +
  theme_minimal() +
  ylab("") + 
  xlab("Relative time (%)") +
  theme(axis.title.x = element_text(margin = margin(t = 10), size=16),  
        axis.title.y = element_text(margin = margin(r = 10), size=16),
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),  
        axis.text = element_text(color = "black", size=14),  
        axis.text.y = element_text(face = "bold", size = 18),
        panel.grid = element_blank()  
  ) +
  theme(legend.position = "none") +
  scale_y_discrete(labels = c("control" = "Control", "patient" = "Stroke"))

#plot1c

ggpubr::ggarrange(plot1a, NULL, plot1b, NULL, plot1c, ncol = 5, widths = c(1, 0.3, 1, 0, 1))

# plot 2 ####

# 2a - examples ####

plot2a <- plot2a_dat %>%
  ggplot(aes(x = mov_time)) +
  geom_line(aes(y = ElbowAngles_x, color = "Elbow flexion/extension"), size = 1) +
  geom_line(aes(y = ShoulderAngles_z, color = "Shoulder rotation"), size = 1) +
  geom_point(data = plot2a_pv_el, aes(x = mov_time, y = ElbowAngles_x, color = "Elbow flexion/extension", shape = "Peak Velocity"), size = 4) +
  geom_point(data = plot2a_pv_shz, aes(x = mov_time, y = ShoulderAngles_z, color = "Shoulder rotation", shape = "Peak Velocity"), size = 4) +
  theme_minimal() +
  xlim(0, 1.2) +
  scale_x_continuous(breaks = c(0, .5, 1)) + 
  xlab("Time (s)") +
  ylab(expression(paste("Joint angle (", degree, ")"))) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), size = 16),
    axis.title.y = element_text(margin = margin(r = 10), size = 16),
    axis.line = element_line(color = "black"),  # Set axes lines color
    axis.ticks = element_line(color = "black"),  # Set ticks color
    axis.text = element_text(color = "black", size = 14),  # Set tick labels color
    panel.grid = element_blank(),  # Remove grid
    legend.text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  ) +
  facet_grid(rows = vars(subject)) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()) +
  # Add annotations
  geom_text(data = subjects, mapping = aes(x = 0, y = 120, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6) +
  # Custom legend
  scale_color_manual(
    values = c("Elbow flexion/extension" = "#f57c00", "Shoulder rotation" = "#c0ca33"),
    name = ""
  ) +
  scale_shape_manual(
    values = c("Peak Velocity" = 17),
    name = ""
  ) +
  theme(legend.position = "none")


#plot2a

# 2b - frequencies peak velocity ####

# Annotation data frames
control_annotation <- data.frame(condition = "control", label = "Control")
stroke_annotation <- data.frame(condition = "patient", label = "Stroke")

el_x_sh_z <- data %>% 
  filter(side == "right") %>% 
  select(condition, rel_time_pv_sh_z_neg, rel_time_pv_el_x) %>%
  pivot_longer(
    cols = starts_with("rel_time_pv"),
    names_to = "joint",
    values_to = "value"
  ) %>%
  mutate(joint = case_when(
    joint == "rel_time_pv_sh_z_neg" ~ "sh_z",
    joint == "rel_time_pv_el_x" ~ "el_x"
  ))



plot2b <- el_x_sh_z %>% 
  ggplot(aes(x = value, color = joint, fill = joint)) +
  geom_density(aes(y = ..count..),
               alpha=.3,
               size = 1) +
  scale_color_manual(values = c("#f57c00", "#c0ca33")) +
  scale_fill_manual(values = c("#f57c00", "#c0ca33")) +
  xlim(0, 100) +
  scale_y_continuous(breaks = c(0, 50, 100)) + 
  theme_minimal() +
  xlab("Relative time (%)") + 
  ylab("Number of trials at peak velocity time points") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), size=16),
    axis.title.y = element_text(margin = margin(r = 10), size=16),
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black", size=14),  
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")   
  ) +
  facet_grid(rows = vars(condition)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  theme(legend.position = "none") +
  # Add annotations for Control
  geom_text(data = control_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6) +
  # Add annotations for Stroke
  geom_text(data = stroke_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6)


#plot2b

# 2c - box plots ####

plot2c <- el_x_sh_z %>% 
  mutate(condition = factor(condition, levels = c("patient", "control"))) %>% 
  ggplot(aes(y=condition, x=value, color=joint)) +
  geom_boxplot(width=.5, size=1, position=position_dodge(width=0.7)) +
  scale_color_manual(values = c("#f57c00", "#c0ca33")) +
  xlim(0,100) +
  theme_minimal() +
  ylab("") + 
  xlab("Relative time (%)") +
  theme(axis.title.x = element_text(margin = margin(t = 10), size=16),  
        axis.title.y = element_text(margin = margin(r = 10), size=16),
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),  
        axis.text = element_text(color = "black", size=14),  
        axis.text.y = element_text(face = "bold", size = 18),
        panel.grid = element_blank()  
  ) +
  theme(legend.position = "none") +
  scale_y_discrete(labels = c("control" = "Control", "patient" = "Stroke"))

#plot2c

ggpubr::ggarrange(plot2a, NULL, plot2b, NULL, plot2c, ncol = 5, widths = c(1, 0.3, 1, 0, 1))

# 2d - examples ####

plot2d <- plot2d_dat %>%
  ggplot(aes(x = mov_time)) +
  geom_line(aes(y = ShoulderAngles_x, color = "Shoulder flexion/extension"), size = 1) +
  geom_line(aes(y = ShoulderAngles_y, color = "Shoulder ab-/adduction"), size = 1) +
  geom_point(data = plot2d_pv_shx, aes(x = mov_time, y = ShoulderAngles_x, color = "Shoulder flexion/extension", shape = "Peak Velocity"), size = 4) +
  geom_point(data = plot2d_pv_shy, aes(x = mov_time, y = ShoulderAngles_y, color = "Shoulder ab-/adduction", shape = "Peak Velocity"), size = 4) +
  theme_minimal() +
  xlim(0, 1.2) +
  scale_x_continuous(breaks = c(0, .5, 1)) + 
  xlab("Time (s)") +
  ylab(expression(paste("Joint angle (", degree, ")"))) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), size = 16),
    axis.title.y = element_text(margin = margin(r = 10), size = 16),
    axis.line = element_line(color = "black"),  # Set axes lines color
    axis.ticks = element_line(color = "black"),  # Set ticks color
    axis.text = element_text(color = "black", size = 14),  # Set tick labels color
    panel.grid = element_blank(),  # Remove grid
    legend.text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  ) +
  facet_grid(rows = vars(subject)) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()) +
  # Add annotations
  geom_text(data = subjects, mapping = aes(x = 0, y = 130, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6) +
  # Custom legend
  scale_color_manual(
    values = c("Shoulder flexion/extension" = "#5e35b1", "Shoulder ab-/adduction" = "#00acc1"),
    name = ""
  ) +
  scale_shape_manual(
    values = c("Peak Velocity" = 17),
    name = ""
  ) +
  theme(legend.position = "none")

#plot2d

# 2e - frequencies peak velocity ####

# Annotation data frames
control_annotation <- data.frame(condition = "control", label = "Control")
stroke_annotation <- data.frame(condition = "patient", label = "Stroke")

sh_x_sh_y <- data %>% 
  filter(side == "right") %>% 
  select(condition, rel_time_pv_sh_x, rel_time_pv_sh_y) %>%
  pivot_longer(
    cols = starts_with("rel_time_pv"),
    names_to = "joint",
    values_to = "value"
  ) %>%
  mutate(joint = case_when(
    joint == "rel_time_pv_sh_x" ~ "sh_x",
    joint == "rel_time_pv_sh_y" ~ "sh_y"
  ))

plot2e <- sh_x_sh_y %>% 
  ggplot(aes(x = value, color = joint, fill = joint)) +
  #geom_histogram(alpha=.3,
  #               binwidth = 2,
  #               size = 1,
  #               position = "identity") +
  geom_density(aes(y = ..count..),
               alpha=.3,
               size = 1) +
  scale_color_manual(values = c("#5e35b1", "#00acc1")) +
  scale_fill_manual(values = c("#5e35b1", "#00acc1")) +
  xlim(0, 100) +
  scale_y_continuous(breaks = c(0, 50, 100)) + 
  theme_minimal() +
  xlab("Relative time (%)") + 
  ylab("Number of trials at peak velocity time points") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), size=16),
    axis.title.y = element_text(margin = margin(r = 10), size=16),
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black", size=14),  
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")
  ) +
  facet_grid(rows = vars(condition)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  theme(legend.position = "none") +
  # Add annotations for Control
  geom_text(data = control_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6) +
  # Add annotations for Stroke
  geom_text(data = stroke_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6)

#plot2e

# 2f - box plots ####

plot2f <- sh_x_sh_y %>% 
  mutate(condition = factor(condition, levels = c("patient", "control"))) %>% 
  ggplot(aes(y=condition, x=value, color=joint)) +
  geom_boxplot(width=.5, size=1, position=position_dodge(width=0.7)) +
  scale_color_manual(values = c("#5e35b1", "#00acc1")) +
  xlim(0,100) +
  theme_minimal() +
  ylab("") + 
  xlab("Relative time (%)") +
  theme(axis.title.x = element_text(margin = margin(t = 10), size=16),  
        axis.title.y = element_text(margin = margin(r = 10), size=16),
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),  
        axis.text = element_text(color = "black", size=14),  
        axis.text.y = element_text(face = "bold", size = 18),
        panel.grid = element_blank()  
  ) +
  theme(legend.position = "none") +
  scale_y_discrete(labels = c("control" = "Control", "patient" = "Stroke"))

#plot2f

ggpubr::ggarrange(plot2d, NULL, plot2e, NULL, plot2f, ncol = 5, widths = c(1, 0.3, 1, 0, 1))

# plot S1 ####

custom_labels <- c("control" = "Healthy", "patient" = "Stroke")

data %>% 
  filter(side == "right") %>% 
  ggplot(aes()) +
  geom_density(aes(x=rel_time_pv_sh_x, y = ..count..), color="#5e35b1", fill="#5e35b1", alpha=.3, size = 1) +
  geom_density(aes(x=rel_time_pv_sh_y, y = ..count..), color="#00acc1", fill="#00acc1", alpha=.3, size = 1) +
  geom_density(aes(x=rel_time_pv_sh_z_neg, y = ..count..), color="#c0ca33", fill="#c0ca33", alpha=.3, size = 1, bw = "nrd0") +
  geom_density(aes(x=rel_time_pv_el_x,  y = ..count..), color="#f57c00", fill="#f57c00", alpha=.3, size = 1) +
  xlim(0, 100) +
  scale_y_continuous(breaks = c(0, 50, 100)) + 
  theme_minimal() +
  xlab("Relative time (%)") + 
  ylab("Number of trials at peak velocity time points") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), size=16),
    axis.title.y = element_text(margin = margin(r = 10), size=16),
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black", size=14),  
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")   
  ) +
  facet_grid(rows = vars(condition)) + 
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  theme(legend.position = "none") +
  # Add annotations for Control
  geom_text(data = control_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6) +
  # Add annotations for Stroke
  geom_text(data = stroke_annotation, mapping = aes(x = 0, y = 100, label = label),
            hjust = 0, vjust = 1, fontface = "bold", inherit.aes = FALSE, size = 6)

# plot S2 ####

# * only coefficients ####

# prepare data 

cor_data <- data %>%
  filter(condition == "patient") %>% 
  filter(side == "right") %>% 
  select(subject, arat_total_affected, ashw_mean_affected, bbt_te_rel, rel_nhpt_affected, fma_total, ratio_mean_FA,
         move_time, time_peak_vel_hand_in_sec, time_peak_acc_hand_in_sec, number_vel_peaks,
         rel_time_peak_vel_hand, rel_time_peak_acc_hand,
         sh_x_el_x, sh_y_el_x, sh_z_el_x, sh_x_sh_y, sh_x_sh_z, sh_y_sh_z) %>% 
  group_by(subject) %>% 
  mutate(mean_dur = mean(move_time), 
         mean_peaks = mean(number_vel_peaks),
         mean_ttpa = mean(rel_time_peak_acc_hand),
         mean_ttpv = mean(rel_time_peak_vel_hand),
         mean_sh_x_el_x = mean(sh_x_el_x),
         mean_sh_y_el_x = mean(sh_y_el_x),
         mean_sh_z_el_x = mean(sh_z_el_x),
         mean_sh_x_sh_y = mean(sh_x_sh_y),
         mean_sh_x_sh_z = mean(sh_x_sh_z),
         mean_sh_y_sh_z = mean(sh_y_sh_z)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(-c(subject, trial_number, move_time, 
            rel_time_peak_vel_hand, rel_time_peak_acc_hand,
            time_peak_vel_hand_in_sec, time_peak_acc_hand_in_sec, number_vel_peaks,
            sh_x_el_x, sh_y_el_x, sh_z_el_x, sh_x_sh_y, sh_x_sh_z, sh_y_sh_z))


cor_results <- psych::corr.test(cor_data, method = "spearman")
cor_matrix <- cor_results$r
p_matrix <- cor_results$p

new_names <- c("ARAT", "Ashworth", "BBT", "NHPT", 
               #"SIS", 
               "UFMA", 
               "Ratio mean FA",
               "Movement time", "Smoothness", "TTPa", "TTPv", 
               "TI EL f SH f", "TI EL f SH a", "TI EL f SH r", 
               "TI SH f SH a", "TI SH f SH r", "TI SH a SH r")


colnames(cor_matrix) <- new_names
rownames(cor_matrix) <- new_names


corr_df <- as.data.frame(as.table(cor_matrix))
colnames(corr_df) <- c("Var1", "Var2", "value")


upper <- corr_df %>% 
  filter(match(Var1, unique(corr_df$Var1)) > match(Var2, unique(corr_df$Var2))) %>% 
  ggplot(aes(Var1, Var2, fill=value)) +
  geom_tile(color = "#f2f4f4", height = 1, width = 1) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = .01, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)) +  
  geom_text(aes(label = sprintf("%.2f", value), 
                fontface = "plain"),
            color = "black",  
            size = 3) + 
  coord_fixed() +
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank()) +
  scale_y_discrete(limits = rev(unique(corr_df$Var2)), position = "right") +
  scale_fill_gradient2(low = "#ef5350", high = "#5c6bc0", mid = "white", 
                       limit = c(-1, 1), name = "Correlation\ncoefficient") +
  guides(fill = FALSE) +
  scale_x_discrete(position = "top")

upper

# * only p-values ####

colnames(p_matrix) <- new_names
rownames(p_matrix) <- new_names

p_df <- as.data.frame(as.table(p_matrix))
colnames(p_df) <- c("Var1", "Var2", "p_value")

p_df <- p_df %>%
  filter(Var1 != Var2) 


lower <- p_df %>% 
  filter(match(Var1, unique(corr_df$Var1)) > match(Var2, unique(corr_df$Var2))) %>% 
  ggplot(aes(Var2, Var1)) +
  geom_tile(color = "#d6dbdf", fill = "white", height = 1, width = 1) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color="black"),   
        axis.text.y = element_text(size = 9, , color="black")  
  ) + 
  geom_text(aes(label = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)), 
                fontface = ifelse(p_value < 0.05, "bold", "plain")),  
            color = "black",  
            size = 2.7) + 
  
  
  coord_fixed() +
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank()) +
  scale_y_discrete(limits = rev(unique(p_df$Var2)), position = "left") +
  scale_x_discrete(position = "bottom") 


lower

# plot S3 ####

new_names <- c("ARAT", "Ashworth", "BBT", "NHPT", "UFMA", 
                "Ratio_mean_FA",
                "Movement_time", "Smoothness", "TTPa", "TTPv", 
                "TI_EL_f_SH_f", "TI_EL_f_SH_a", "TI_EL_f_SH_r", 
                "TI_SH_f_SH_a", "TI_SH_f_SH_r", "TI_SH_a_SH_r")

plot_data_new <- cor_data
colnames(plot_data_new) <- new_names

summary(plot_data_new)

# define clinical variables to plot

variables_to_plot <- c("ARAT", "BBT", "NHPT", "UFMA")

# select kinematic variables for correlation - select one of the three lines

kinematic_variables <- c(
  #"Movement_time", "Smoothness", "TTPa", "TTPv"  # endpoint
  #"TI_EL_f_SH_f", "TI_EL_f_SH_a", "TI_EL_f_SH_r" # inter-joint
  "TI_SH_f_SH_a", "TI_SH_f_SH_r", "TI_SH_a_SH_r"  # intra-joint
)


# create the scatter plots for the defined combinations

plot_list <- list() # empty list

for (var1 in variables_to_plot) {
  for (var2 in kinematic_variables) {
    # calculate correlation value and p-value with psych::corr.test
    corr_result <- psych::corr.test(plot_data_new[[var1]], plot_data_new[[var2]], method = "spearman")
    
    # access to the correlation value and p-value
    if (length(corr_result$r) == 1) {
      # if only one correlation was calculated (between two variables)
      cor_value <- corr_result$r
      p_value <- corr_result$p
    } else {
      # if several values are expected in a matrix
      cor_value <- corr_result$r[1, 2]  
      p_value <- corr_result$p[1, 2]     
    }
    
    if (!is.na(cor_value) && !is.na(p_value)) { 
      # format p-value
      p_value_formatted <- ifelse(p_value < 0.001, "<0.001", formatC(p_value, format = "f", digits = 3))
      
      # replace underscores with spaces for the graph labels
      var1_name_plot <- gsub("_", " ", var1)
      var2_name_plot <- gsub("_", " ", var2)
      
      # create the scatter plot with ggplot2
      p <- ggplot(plot_data_new, aes_string(x = var1, y = var2)) +
        geom_point(color = "#34495e", size = 2, alpha = .6) +
        geom_smooth(method = "lm", size = 1, linetype = "solid", color = "blue", alpha = 0.2, se = TRUE) +
        
        # using r and p-values, create the title
        labs(title = paste(var1_name_plot, "/", var2_name_plot, 
                           "\n(r =", formatC(cor_value, format = "f", digits = 2), 
                           ", p =", p_value_formatted, ")"),
             x = var1_name_plot,
             y = var2_name_plot) +
        
        theme_minimal() +
        theme(
          plot.title = element_text(size = 9),
          axis.title.x = element_text(size = 8),  
          axis.title.y = element_text(size = 8),  
          axis.text = element_text(size = 8)
        ) 
      
      # add the plot to the plot list
      plot_list[[length(plot_list) + 1]] <- p
    }
  }
}

for (i in 1:length(plot_list)) {
  print(paste("Plot Nummer:", i))  # output of the plot number
  print(plot_list[[i]])              # display plot
}


#selected_indices <- setdiff(1:16, c()) # endpoint
selected_indices <- setdiff(1:12, c()) # inter-joint | intra-joint


selected_plots <- plot_list[selected_indices]
grid.arrange(grobs = selected_plots, ncol = 4)

# plot S4 ####

# insert

# variable                  name

#move_time                  Movement time
#rel_time_peak_vel_hand     TTPv
#rel_time_peak_acc_hand     TTPa
#number_vel_peaks           Smoothness

#sh_x_el_x                  TI EL f SH f
#sh_y_el_x                  TI EL f SH a
#sh_z_el_x                  TI EL f SH r

#sh_x_sh_y                  TI SH f SH a
#sh_x_sh_z                  TI SH f SH r
#sh_y_sh_z                  TI SH a SH r



s4 <- data %>% 
  filter(side == "right") %>% 
  group_by(subject) %>% 
  mutate(mean = mean(move_time)) %>% #specify here which variable to plot against NHPT, select from the column "variable" above 
  slice(1)

# Custom colors for conditions
custom_colors <- c("#1a237e", "#b71c1c") 


cor_results <- s4 %>%
  group_by(condition) %>%
  do(cor_test = psych::corr.test(select(.,  rel_nhpt_affected, mean), method = "spearman"))

cor_results <- cor_results %>%
  rowwise() %>%
  mutate(
    cor_value = cor_test$r[1, 2],    
    p_value = cor_test$p[1, 2],      
    y_position = if_else(condition == unique(s4$condition)[1], 1.8, 1.65) # give meaningful values here for the point at which the r and p values are to be shown in relation to the y-axis
  )


ggplot(s4, aes(x =  rel_nhpt_affected, y = mean, color = condition, shape = rev(condition))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = T, alpha = 0.1, aes(fill = condition, color = condition), linetype = "longdash") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(x = "NHPT", y = "Movement time") + #specify here which variable to plot against NHPT, select from the column "name" above 
  geom_text(data = cor_results, aes(
    label = paste0("r = ", round(cor_value, 2), ", p = ", sprintf("%.3f", p_value)), 
    x = 0.4, 
    y = y_position, 
    color = condition), 
    hjust = 0,
    size=5) +
  theme(
    text = element_text(size = 14),          
    axis.title = element_text(size = 16),      
    axis.text = element_text(size = 14)
  ) +
  theme(legend.position = "none")

# ~ ####



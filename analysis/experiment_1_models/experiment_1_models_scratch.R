
#'
#' Scratch file for experiment 1 modeling
#' TODO move the final versions of these analyses to a clean script
#'


# INIT ----
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(tidyverse)


# GLOBALS ----

# Directories
DATA_DIR = '../../data' # folder containing data files
E1_DATA_FILE = 'experiment_1/physics_continual_learning_trialdata_e2_full.csv' # name of file containing full dataset
NONSOCIAL_DATA_FILE = 'experiment_1_nonsocial/nonsocial_trialdata.csv' # name of file containing trial dataset for nonsocial control


# Experiment parameters
TRIAL_BLOCKS = 8
CONDITION_LOOKUP = c(
  'bad' = 'unreliable',
  'good' = 'reliable',
  'improve' = 'improving',
  'worsen' = 'worsening',
  'nonsocial' = 'nonsocial'
)
CONDITION_ORDER = c('unreliable', 'reliable', 'improving', 'worsening', 'nonsocial')

# Figure variables
COLORS = c(
  'unreliable' = '#E74C3C', # unreliable
  'reliable' = '#2980B9', # reliable
  'improving' = '#7D3C98', # improving
  'worsening' = '#AF7AC5', # worsening
  'nonsocial' = 'gray' # nonsocial
)

PLOT_THEME = theme(
  # titles
  plot.title = element_text(face = 'bold', size = 32, family = 'Charter', margin = margin(b = 0.5, unit = 'line')),
  axis.title.y = element_text(face = 'bold', size = 32, family = 'Charter', margin = margin(r = 0.5, unit = 'line')),
  axis.title.x = element_text(face = 'bold', size = 32, family = 'Charter', margin = margin(t = 0.5, unit = 'line')),
  legend.title = element_text(face = 'bold', size = 24, family = 'Charter'),
  # axis text
  axis.text.x = element_text(size = 20, face = 'bold', angle = 0, vjust = 1, family = 'Charter', margin = margin(t = 0.5, unit = 'line'), color = 'black'),
  axis.text.y = element_text(size = 20, face = 'bold', family = 'Charter', margin = margin(r = 0.5, unit = 'line'), color = 'black'),
  # legend text
  legend.text = element_text(size = 24, face = 'bold', family = 'Charter', margin = margin(b = 0.5, unit = 'line')),
  # facet text
  strip.text = element_text(size = 12, family = 'Charter'),
  # backgrounds, lines
  panel.background = element_blank(),
  strip.background = element_blank(),
  panel.grid = element_line(color = 'gray'),
  axis.line = element_line(color = 'black'),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  # positioning
  legend.position = 'bottom',
  legend.key = element_rect(colour = 'transparent', fill = 'transparent')
)


# ANALYSIS FUNCTIONS ----
# NB: these are same as in `experiment_1_analysis_scratch.R`


# Determine the participant's final paddle angle on each trial
# Uses either the agent's suggestion or the final angle in the sequence of
# captured paddle angles
get_final_angle = function(paddle_rho, paddle_tr) {
  paddle_tr = strsplit(paddle_tr, '\\[')[[1]][2]
  paddle_tr = strsplit(paddle_tr, '\\]')[[1]][1]
  paddle_tr = str_replace_all(paddle_tr, ' ', '')
  if(paddle_tr == '') {
    return(as.numeric(paddle_rho))
  } else {
    paddle_vec = strsplit(paddle_tr, ',')[[1]]
    return(as.numeric(paddle_vec[length(paddle_vec)]))
  }
}

# Function for getting unsigned difference (radians) between `start_angle` and `end_angle`
get_angle_difference = function(start_angle, end_angle) {
  return(abs(((start_angle - end_angle + pi) %% (2 * pi)) - pi))
}

# Difference between `participant_angle` and `true_angle` towards / away from `agent_angle`
# If distance of `participant_angle` away from `true_angle` is *towards* `agent_angle`, this value is positive.
# If distance of `participant_angle` away from `true_angle` is *away from* `agent_angle`, this value is negative
get_signed_participant_error = function(agent_angle, participant_angle, true_angle) {
  agent_error_signed = ((agent_angle - true_angle + pi) %% (2 * pi)) - pi
  participant_error_signed = ((participant_angle - true_angle + pi) %% (2 * pi)) - pi
  if(agent_error_signed * participant_error_signed > 0) {
    return(abs(participant_error_signed))
  } else {
    return(-abs(participant_error_signed))
  }
}

# Function for transforming radian measures to degrees
rad_to_deg = function(rad) {
  return(rad * (180 / pi))
}



# DATA PROCESSING ----

# Read in trial data
trial_data = read_csv(paste(DATA_DIR, E1_DATA_FILE, sep = '/'))
glimpse(trial_data)

# Add derived columns
# NB: this takes 5-10s to run
# TODO remove any columns below that aren't used in the final analyses
trial_data = trial_data |>
  rowwise() |>
  mutate(
    trial_block = ceiling((trialInd + 1) / ((max(trial_data$trialInd) + 1) / TRIAL_BLOCKS)),
    condition_str = factor(CONDITION_LOOKUP[agent_cond], levels = CONDITION_ORDER),
    paddleIntervention = !trustedAgent,
    final_par = get_final_angle(paddle_rho, paddle_tr),
    intervene_dist = get_angle_difference(paddle_rho, final_par),
    intervene_dist_degrees = rad_to_deg(intervene_dist),
    bot_error = get_angle_difference(groundtruthAngle, paddle_rho),
    bot_error_degrees = rad_to_deg(bot_error),
    human_error = get_angle_difference(groundtruthAngle, final_par),
    human_error_degrees = rad_to_deg(human_error),
    signed_human_error = get_signed_participant_error(paddle_rho, final_par, groundtruthAngle),
    signed_human_error_degrees = rad_to_deg(signed_human_error)
  )
# sanity check
glimpse(trial_data)


# Read in nonsocial data
nonsocial_data = read_csv(paste(DATA_DIR, NONSOCIAL_DATA_FILE, sep = '/'))
glimpse(nonsocial_data)

# Add derived columns
# NB: similar columns here but slightly different from above, not copy-paste
nonsocial_data = nonsocial_data |>
  rowwise() |>
  mutate(
    trial_block = ceiling((trialInd + 1) / ((max(trial_data$trialInd) + 1) / TRIAL_BLOCKS)),
    condition_str = factor(CONDITION_LOOKUP[condition], levels = CONDITION_ORDER),
    paddle_intervention = updatedPaddle,
    final_participant_angle = get_final_angle(paddle_rho, paddle_tr),
    intervene_dist = get_angle_difference(paddle_rho, final_participant_angle),
    intervene_dist_degrees = rad_to_deg(intervene_dist),
    human_error = get_angle_difference(groundtruthAngle, final_participant_angle),
    human_error_degrees = rad_to_deg(human_error)
  )
# sanity check
glimpse(nonsocial_data)









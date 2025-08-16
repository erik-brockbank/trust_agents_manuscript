
#'
#' Scratch file for experiment 1 figures
#' TODO move the final versions of these figures to a clean script
#'



# INIT ----
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggExtra) # for ggMarginal plot below
library(tidyverse)


# GLOBALS ----

# Directories
DATA_DIR = '../../data/experiment_1' # folder containing data files
DATA_FILE = 'physics_continual_learning_trialdata_e2_full.csv' # name of file containing full dataset
SURVEY_FILE = 'physics_continual_learning_surveydata_e2_full.csv' # name of file containing survey responses
FIGURE_PATH = '../../results/experiment_1' # output directory

# Experiment parameters
TRIAL_BLOCKS = 8
CONDITION_LOOKUP = c(
  'bad' = 'unreliable',
  'good' = 'reliable',
  'improve' = 'improving',
  'worsen' = 'worsening'
)
CONDITION_ORDER = c('unreliable', 'reliable', 'improving', 'worsening')

# Figure variables
COLORS = c(
  'unreliable' = '#E74C3C', # unreliable
  'reliable' = '#2980B9', # reliable
  'improving' = '#7D3C98', # improving
  'worsening' = '#AF7AC5' # worsening
)

DEFAULT_PLOT_THEME = theme(
  # titles
  plot.title = element_text(size = 32, family = 'Avenir', margin = margin(b = 0.5, unit = 'line')),
  axis.title.y = element_text(size = 24, family = 'Avenir', margin = margin(r = 0.5, unit = 'line')),
  axis.title.x = element_text(size = 24, family = 'Avenir', margin = margin(t = 0.5, unit = 'line')),
  legend.title = element_text(size = 24, family = 'Avenir'),
  # axis text
  axis.text.x = element_text(size = 20, angle = 0, vjust = 1, family = 'Avenir', margin = margin(t = 0.5, unit = 'line'), color = 'black'),
  axis.text.y = element_text(size = 20, family = 'Avenir', margin = margin(r = 0.5, unit = 'line'), color = 'black'),
  # legend text
  legend.text = element_text(size = 24, family = 'Avenir', margin = margin(b = 0.5, unit = 'line')),
  # facet text
  strip.text = element_text(size = 12, family = 'Avenir'),
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

# Read in data
trial_data = read_csv(paste(DATA_DIR, DATA_FILE, sep = '/'))
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


# FIGURE: Intervention rate by trial block, condition ----
# Calculate 95% CIs beforehand (customizable Y axis ranges)
fig_intervention = trial_data |>
  group_by(condition_str, gameID, trial_block) |>
  # Calculate subject-level intervention rate in each block
  summarize(
    subject_intervention_rate = mean(as.numeric(paddleIntervention)),
    .groups = 'drop'
  ) |>
  group_by(condition_str, trial_block) |>
  # Calculate 95% CI across subjects in each block
  summarize(
    ci_stats = list(Hmisc::smean.cl.boot(subject_intervention_rate)),
    .groups = 'drop'
  ) |>
  tidyr::unnest_wider(ci_stats) |>
  ggplot(aes(x = trial_block, y = Mean, color = condition_str)) +
  geom_point(
    size = 6
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0,
    linewidth = 1
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  scale_x_continuous(
    name = 'block number',
    breaks = seq(TRIAL_BLOCKS),
    labels = as.character(seq(TRIAL_BLOCKS))
  ) +
  scale_y_continuous(
    name = 'intervention rate',
    breaks = seq(0.5, 1, by = 0.1),
    labels = as.character(seq(0.5, 1, by = 0.1)),
    limits = c(0.495, 1.0)
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    legend.position = 'none'
  )

fig_intervention
ggsave(
  fig_intervention,
  filename = 'raw_intervention_rate.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 8,
  height = 8,
)



# FIGURE: Signed error by condition and trial block ----
# Calculate 95% CIs beforehand (allows customizable Y axis ranges)
fig_signed_error = trial_data |>
  group_by(condition_str, gameID, trial_block) |>
  # Calculate subject-level mean signed error in each block
  summarize(
    mean_signed_error = mean(signed_human_error_degrees),
    .groups = 'drop'
  ) |>
  group_by(condition_str, trial_block) |>
  # Calculate 95% CI across subjects in each block
  summarize(
    ci_stats = list(Hmisc::smean.cl.boot(mean_signed_error)),
    .groups = 'drop'
  ) |>
  tidyr::unnest_wider(ci_stats) |>
  ggplot(aes(x = trial_block, y = Mean, color = condition_str)) +
  geom_point(
    size = 4
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0,
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.5) +
  annotate(
    geom = 'polygon', x = c(-Inf, -Inf, Inf, Inf), y = c(0, Inf, Inf, 0),
    fill = 'green', # these colors get modified "post-production"
    alpha = 0.1
  ) +
  annotate(
    geom = 'polygon', x = c(-Inf, -Inf, Inf, Inf), y = c(0, -Inf, -Inf, 0),
    fill = 'red', # these colors get modified "post-production"
    alpha = 0.1
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  scale_x_continuous(
    name = 'block number',
    breaks = seq(TRIAL_BLOCKS),
    labels = as.character(seq(TRIAL_BLOCKS))
  ) +
  scale_y_continuous(
    name = 'signed error magnitude (deg.)',
    breaks = seq(-5, 10, by = 5),
    labels = as.character(seq(-5, 10, by = 5)),
    limits = c(-5, 12.9)
  ) +
  facet_wrap(~condition_str, ncol = 1) +
  DEFAULT_PLOT_THEME +
  theme(
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.line.x = element_blank(),
    strip.text.x = element_text(size = 20, family = 'Avenir'),
    legend.position = 'none',
  )

fig_signed_error
ggsave(
  fig_signed_error,
  filename = 'raw_signed_error.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 8,
  height = 8,
)



# FIGURE: Critical trial intervention rate by condition ----
# Calculate 95% CIs beforehand (allows customizable Y axis ranges)
fig_critical_trial_intervention_rate = trial_data |>
  filter(criticalTrial) |>
  group_by(condition_str, gameID) |>
  # Subject-level intervention rate across all critical trials
  summarize(
    subject_intervention_rate = mean(as.numeric(paddleIntervention)),
    .groups = 'drop'
  ) |>
  group_by(condition_str) |>
  # Calculate 95% CI across subjects
  summarize(
    ci_stats = list(Hmisc::smean.cl.boot(subject_intervention_rate)),
    .groups = 'drop'
  ) |>
  tidyr::unnest_wider(ci_stats) |>
  ggplot(
    aes(
      x = condition_str,
      y = Mean,
      fill = condition_str,
    )
  ) +
  geom_bar(
    aes(y = Mean),
    stat = 'identity',
    width = 0.75,
    alpha = 0.75,
    color = 'black'
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0,
    linewidth = .75,
  ) +
  geom_hline(yintercept = 1.0, linewidth = 0.5, linetype = 'dashed') +
  scale_x_discrete(
    name = element_blank(),
    breaks = CONDITION_ORDER,
    labels = CONDITION_ORDER
  ) +
  scale_y_continuous(
    name = 'critical trial intervention rate',
    breaks = seq(0, 1.0, by = 0.25),
    labels = as.character(seq(0, 1.0, by = 0.25)),
  ) +
  scale_fill_manual(
    name = element_blank(),
    values = COLORS
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    legend.position = 'none'
  )

fig_critical_trial_intervention_rate
ggsave(
  fig_critical_trial_intervention_rate,
  filename = 'raw_critical_trial_intervention_rate.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 8,
  height = 7,
)



# FIGURE: Bot error and intervention magnitude by condition ----

distance_figure = trial_data |>
  # filter(!criticalTrial) |>
  ggplot(aes(x = bot_error_degrees, y = intervene_dist_degrees, color = condition_str, fill = condition_str)) +
  # NB: we add scatterplot points here for density estimation, but don't actually show them
  geom_point(
    size = 0.75,
    alpha = 0 # hide the points
  ) +
  geom_smooth(
    method = 'lm',
    linewidth = 1.5,
  ) +
  geom_abline(
    linetype = 'dashed',
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = mean(trial_data$bot_error_degrees[trial_data$criticalTrial]),
    linetype = 'solid',
    linewidth = 0.25,
    # color = 'darkgray'
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  scale_fill_manual(
    name = element_blank(),
    values = COLORS
  ) +
  scale_x_continuous(
    name = 'bot error (deg.)',
    breaks = seq(0, 180, by = 60),
    labels = seq(0, 180, by = 60),
    limits = c(-1, 180)
  ) +
  scale_y_continuous(
    name = 'intervention magnitude (deg.)',
    breaks = seq(0, 180, by = 60),
    labels = seq(0, 180, by = 60),
    limits = c(-1, 180)
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    legend.position = 'none'
  )

distance_figure
# Add marginal density plots
# NB: we pipe directly to the save call below so the margin
ggMarginal(
  distance_figure,
  groupColour = T,
  type = 'density',
  margins = 'both',
  size = 5,
  linewidth = 1,
  # kernel density estimation bandwidth (higher == smooth)
  bw = 10
) |>
ggsave(
  filename = 'raw_bot_error_intervention_magnitude.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 8,
  height = 7,
)






#'
#' Scratch file for experiment 2 analyses
#' TODO move the final versions of these figures and analyses to a clean script
#'



# INIT ----
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(lme4)
library(tidyverse)


# GLOBALS ----

# Directories
DATA_DIR = '../../data/experiment_2' # folder containing data files
DATA_FILE = 'trust_agents_trialdata_bias_complete.csv' # name of file containing full dataset
SURVEY_FILE = 'trust_agents_surveydata_bias_complete.csv' # name of file containing survey responses
FIGURE_PATH = '../../results/experiment_2' # output directory

# Figure variables
TRIAL_BLOCKS = 8
CONDITION_LOOKUP = c(
  'accurate' = 'no bias',
  'heavy' = 'heavy',
  'light' = 'light'
)
CONDITION_ORDER = c('heavy', 'light', 'no bias')

COLORS = c(
  'no bias' = '#2980B9',
  'heavy' = '#7D3C98',
  'light' = '#AF7AC5'
)

# Same theme as E1
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

# Difference between two angles towards / away from "anchor angle"
# If distance of `end_angle` away from `start_angle` is *towards* `anchor_angle`, this value is positive.
# If distance of `end_angle` away from `start_angle` is *away from* `anchor_angle`, this value is negative
get_signed_error = function(anchor_angle, start_angle, end_angle) {
  anchor_distance_signed = ((anchor_angle - end_angle + pi) %% (2 * pi)) - pi
  angle_distance_signed = ((start_angle - end_angle + pi) %% (2 * pi)) - pi
  if(anchor_distance_signed * angle_distance_signed > 0) {
    return(abs(angle_distance_signed))
  } else {
    return(-abs(angle_distance_signed))
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
    trial_block = ceiling((trialInd + 1) / ((max(trial_data$trialInd) + 1) / TRIAL_BLOCKS)), # NB: trialInd is 0-indexed
    condition_str = factor(CONDITION_LOOKUP[agent_cond], levels = CONDITION_ORDER),
    paddleIntervention = !trustedAgent,
    final_par = get_final_angle(paddle_rho, paddle_tr),
    intervene_dist = get_angle_difference(paddle_rho, final_par),
    intervene_dist_degrees = rad_to_deg(intervene_dist),
    bot_error = get_angle_difference(groundtruthAngle, paddle_rho), # orig
    bot_error_degrees = rad_to_deg(bot_error),
    human_error = get_angle_difference(groundtruthAngle, final_par),
    human_error_degrees = rad_to_deg(human_error),
    human_error_reduction_degrees = bot_error_degrees - human_error_degrees,
    improved = as.numeric(human_error_reduction_degrees > 0),
    # human error signed relative to bot suggestion (away from GT towards bot -> positive)
    signed_human_error = get_signed_error(paddle_rho, final_par, groundtruthAngle),
    signed_human_error_degrees = rad_to_deg(signed_human_error),
    # bot error signed relative to *launch* (away from GT towards launch -> positive)
    signed_bot_error = get_signed_error(launching_rho, paddle_rho, groundtruthAngle),
    signed_bot_error_degrees = rad_to_deg(signed_bot_error),
    # human correction signed relative to ground truth (away from suggestion towards GT -> positive)
    signed_human_correction = get_signed_error(groundtruthAngle, final_par, paddle_rho),
    signed_human_correction_degrees = rad_to_deg(signed_human_correction),
    # proportion of error resolved
    bot_error_proportion_resolved = (bot_error_degrees - human_error_degrees) / bot_error_degrees
  )

# sanity check
glimpse(trial_data)



# FIGURE: Intervention rate by trial block, condition ----
# Calculate 95% CIs beforehand (customizable Y axis ranges)
fig_intervention = trial_data |>
  # NB: filtering out mystery round trials
  filter(!criticalTrial) |>
  group_by(condition_str, gameID, trial_block) |>
  # Calculate subject-level intervention rate in each block
  summarize(
    subject_intervention_rate = mean(as.numeric(paddleIntervention)),
    .groups = 'drop'
  ) |>
  group_by(condition_str, trial_block) |>
  # Calculate 95% CI across subjects
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
    breaks = seq(0.4, 1, by = 0.1),
    labels = as.character(seq(0.4, 1, by = 0.1)),
    limits = c(0.37, 1)
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    legend.position = 'none'
  )

# Save figure
fig_intervention
ggsave(
  fig_intervention,
  filename = 'raw_intervention_rate.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 8,
  height = 7,
)


# FIGURE: Intervention magnitude by bot error decile, condition ----
# Calculate 95% CIs beforehand (customizable Y axis ranges)
fig_intervention_magnitude = trial_data |>
  filter(!criticalTrial) |>
  group_by(condition_str) |>
  mutate(
    error_range = cut_interval(bot_error_degrees, n = 10, labels = seq(1, 10, by = 1)),
    .groups = 'drop'
  ) |>
  group_by(condition_str, gameID, error_range) |>
  # Calculate subject-level intervention rate in each block
  summarize(
    subject_intervention_magnitude = mean(intervene_dist_degrees),
    .groups = 'drop'
  ) |>
  group_by(condition_str, error_range) |>
  # Calculate 95% CI across subjects
  summarize(
    ci_stats = list(Hmisc::smean.cl.boot(subject_intervention_magnitude)),
    .groups = 'drop'
  ) |>
  tidyr::unnest_wider(ci_stats) |>
  ggplot(aes(x = error_range, y = Mean, color = condition_str)) +
  geom_point(
    size = 6
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0,
    linewidth = 1
  ) +
  scale_x_discrete(
    name = 'bot error decile',
    breaks = as.character(seq(1, 10, by = 1)),
    labels = seq(1, 10, by = 1)
  ) +
  scale_y_continuous(
    name = 'intervention magnitude (deg.)',
    breaks = seq(0, 40, by = 10),
    labels = as.character(seq(0, 40, by = 10)),
    limits = c(0, 45)
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  facet_wrap(
    ~ condition_str,
    # scales = 'free'
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    # strip.background = element_blank(),
    # strip.text.x = element_blank(),
    strip.text.x = element_text(size = 24, family = 'Avenir'),
    legend.position = 'none'
  )

# Save figure
fig_intervention_magnitude
ggsave(
  fig_intervention_magnitude,
  filename = 'raw_intervention_magnitude.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 14,
  height = 7,
)



# Sanity check: individual points that the 95% CIs above are based on:
# average subject's intervention magnitude by bot error decile
trial_data |>
  filter(!criticalTrial) |>
  group_by(condition_str) |>
  mutate(
    error_range = cut_interval(bot_error_degrees, n = 10, labels = seq(1, 10, by = 1)),
    .groups = 'drop'
  ) |>
  group_by(condition_str, gameID, error_range) |>
  # Calculate subject-level intervention rate in each block
  summarize(
    subject_intervention_magnitude = mean(intervene_dist_degrees),
    .groups = 'drop'
  ) |>
  ggplot(aes(x = error_range, y = subject_intervention_magnitude, color = condition_str)) +
  geom_point(size = 1) +
  scale_x_discrete(
    name = 'bot error decile',
    labels = seq(1, 10, by = 1)
  ) +
  scale_y_continuous(
    name = 'intervention magnitude (deg.)',
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  facet_wrap(
    ~ condition_str,
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    # strip.background = element_blank(),
    # strip.text.x = element_blank(),
    strip.text.x = element_text(size = 24, family = 'Avenir'),
    legend.position = 'none'
  )



# FIGURE: Mystery round intervention rate by trial block, condition ----
# Calculate 95% CIs beforehand (customizable Y axis ranges)
fig_mystery_intervention = trial_data |>
  # NB: mystery round trials only
  filter(criticalTrial) |>
  group_by(condition_str, gameID, trial_block) |>
  # Calculate subject-level intervention rate in each block
  summarize(
    # NB: this is only ever 0 or 1 bc there's one trial per subject per trial block
    subject_intervention_rate = mean(as.numeric(paddleIntervention)),
    .groups = 'drop'
  ) |>
  group_by(condition_str, trial_block) |>
  # Calculate 95% CI across subjects
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
    name = 'intervention rate \n(mystery round trials only)',
    breaks = seq(0.25, 1, by = 0.25),
    labels = as.character(seq(0.25, 1, by = 0.25)),
    limits = c(0.15, 1)
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    legend.position = 'none'
  )

# Save figure
fig_mystery_intervention
ggsave(
  fig_mystery_intervention,
  filename = 'raw_mystery_round_intervention_rate.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 8,
  height = 7,
)






# FIGURE: Signed intervention distance ----

# As with the signed error magnitude in E1, what was the distribution of mystery round interventions
# in E2 in the correct or incorrect direction? Which side of 0 is the distribution on?


# If distance of `participant_angle` away from `agent_angle` is *towards* `true_angle`, this value is positive.
# If distance of `participant_angle` away from `agent_angle` is *away from* `true_angle`, this value is negative.
# TODO: compare the output of this to just checking whether participant error is < original bot error,
# returning positive intervention distance if so, negative if not
get_signed_intervention_distance = function(agent_angle, participant_angle, true_angle) {

  participant_intervention_signed = ((participant_angle - agent_angle + pi) %% (2 * pi)) - pi
  ball_offset_signed = ((true_angle - agent_angle + pi) %% (2 * pi)) - pi

  if(participant_intervention_signed * ball_offset_signed > 0) {
    return(abs(participant_intervention_signed))
  } else {
    return(-abs(participant_intervention_signed))
  }
}

# Add signed intervention distance
trial_data = trial_data |>
  rowwise() |>
  mutate(
    signed_intervention_distance = get_signed_intervention_distance(
      paddle_rho, final_par, groundtruthAngle
    ),
    signed_intervention_distance_degrees = rad_to_deg(signed_intervention_distance)
  )


# Calculate 95% CIs beforehand (customizable Y axis ranges)
fig_signed_intervention_dist = trial_data |>
  filter(criticalTrial) |>
  group_by(condition_str, gameID, trial_block) |>
  # Calculate subject-level mean intervention magnitude in each block
  summarize(
    subject_signed_intervention_dist = mean(signed_intervention_distance_degrees),
    .groups = 'drop'
  ) |>
  group_by(condition_str, trial_block) |>
  # Calculate 95% CI across subjects
  summarize(
    ci_stats = list(Hmisc::smean.cl.boot(subject_signed_intervention_dist)),
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
  scale_x_continuous(
    name = 'block number',
    breaks = seq(TRIAL_BLOCKS),
    labels = as.character(seq(TRIAL_BLOCKS))
  ) +
  scale_y_continuous(
    name = 'signed intervention magnitude (deg.)',
    breaks = seq(-5, 20, by = 5),
    labels = as.character(seq(-5, 20, by = 5)),
    limits = c(-5, 21)
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  facet_wrap(~condition_str, ncol = 1) +
  DEFAULT_PLOT_THEME +
  theme(
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),
    strip.text.x = element_text(size = 20, family = 'Avenir'),
    legend.position = 'none',
  )


fig_signed_intervention_dist
ggsave(
  fig_signed_intervention_dist,
  filename = 'raw_mystery_round_signed_intervention_dist.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 8,
  height = 7,
)



# FIGURE: Intervention distance, bot error deciles ----
# Aggregate: mystery round intervention magnitude by bot error decile, condition
# Calculate 95% CIs beforehand (customizable Y axis ranges)
fig_mystery_round_intervention_magnitude = trial_data |>
  filter(criticalTrial) |>
  group_by(condition_str) |>
  mutate(
    error_range = cut_interval(bot_error_degrees, n = 10, labels = seq(1, 10, by = 1)),
    .groups = 'drop'
  ) |>
  group_by(condition_str, gameID, error_range) |>
  # Calculate subject-level intervention rate in each block
  summarize(
    subject_intervention_magnitude = mean(intervene_dist_degrees),
    .groups = 'drop'
  ) |>
  group_by(condition_str, error_range) |>
  # Calculate 95% CI across subjects
  summarize(
    ci_stats = list(Hmisc::smean.cl.boot(subject_intervention_magnitude)),
    .groups = 'drop'
  ) |>
  tidyr::unnest_wider(ci_stats) |>
  ggplot(aes(x = error_range, y = Mean, color = condition_str)) +
  geom_point(
    size = 6
  ) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0,
    linewidth = 1
  ) +
  scale_x_discrete(
    name = 'bot error decile',
    breaks = as.character(seq(1, 10, by = 1)),
    labels = seq(1, 10, by = 1)
  ) +
  scale_y_continuous(
    name = 'intervention magnitude (deg.) \n(mystery round trials only)',
    breaks = seq(0, 40, by = 10),
    labels = as.character(seq(0, 40, by = 10)),
    limits = c(0, 41)
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  facet_wrap(
    ~ condition_str,
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    # strip.background = element_blank(),
    # strip.text.x = element_blank(),
    strip.text.x = element_text(size = 24, family = 'Avenir'),
    legend.position = 'none'
  )

# Save figure
fig_mystery_round_intervention_magnitude
ggsave(
  fig_mystery_round_intervention_magnitude,
  filename = 'raw_mystery_round_intervention_magnitude.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 14,
  height = 7,
)


# Sanity check: Intervention distance, individual points
# Individual points: mystery round intervention magnitude by bot error, condition
trial_data |>
  filter(criticalTrial) |>
  ggplot(
    aes(
      x = bot_error_degrees,
      y = intervene_dist_degrees,
      color = condition_str
    )
  ) +
  geom_point(
    alpha = 0.5,
    size = 1.5
  ) +
  geom_smooth(
    method = 'lm',
    linewidth = 1,
    color = 'black'
  ) +
  scale_x_continuous(
    name = 'bot error (deg.)',
  ) +
  scale_y_continuous(
    name = 'intervention magnitude (deg.) \n(mystery round trials only)',
  ) +
  scale_color_manual(
    name = element_blank(),
    values = COLORS
  ) +
  facet_wrap(
    ~ condition_str,
    scales = 'free'
  ) +
  DEFAULT_PLOT_THEME +
  theme(
    strip.text.x = element_text(size = 24, family = 'Avenir'),
    legend.position = 'none'
  )


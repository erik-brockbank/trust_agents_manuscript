
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
FIGURE_PATH = '../../results/model_figures' # output directory

# Experiment parameters
N_TRIALS = 96
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
# TODO remove any columns below that aren't used in the final analyses
trial_data = trial_data |>
  rowwise() |>
  mutate( # NB: this takes 5-10s to run
    trial_block = ceiling((trialInd + 1) / (N_TRIALS / TRIAL_BLOCKS)),
    condition_str = factor(CONDITION_LOOKUP[agent_cond], levels = CONDITION_ORDER),
    paddleIntervention = !trustedAgent,
    final_participant_angle = get_final_angle(paddle_rho, paddle_tr),
    intervene_dist = get_angle_difference(paddle_rho, final_participant_angle),
    intervene_dist_degrees = rad_to_deg(intervene_dist),
    bot_error = get_angle_difference(groundtruthAngle, paddle_rho),
    bot_error_degrees = rad_to_deg(bot_error),
    human_error = get_angle_difference(groundtruthAngle, final_participant_angle),
    human_error_degrees = rad_to_deg(human_error),
    signed_human_error = get_signed_participant_error(paddle_rho, final_participant_angle, groundtruthAngle),
    signed_human_error_degrees = rad_to_deg(signed_human_error)
  )
# sanity check
glimpse(trial_data)


# Read in nonsocial data
nonsocial_data = read_csv(paste(DATA_DIR, NONSOCIAL_DATA_FILE, sep = '/'))
glimpse(nonsocial_data)

# Add derived columns
nonsocial_data = nonsocial_data |>
  rowwise() |>
  mutate( # NB: similar columns here but slightly different from above, not copy-paste
    trial_block = ceiling((trialInd + 1) / (N_TRIALS / TRIAL_BLOCKS)),
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


# DATA PROCESSING ----

# Remove anybody that didn't complete 96 trials

# Summarize: incomplete subjects in trial_data
incomplete_ids = trial_data |>
  group_by(gameID) |>
  summarize(
    complete_trials = n()
  ) |>
  ungroup() |>
  filter(
    complete_trials < N_TRIALS
  ) |>
  arrange(complete_trials)
incomplete_ids

# Remove incomplete subjects
trial_data = trial_data |>
  filter(
    !gameID %in% unique(incomplete_ids$gameID)
  )

# Summarize: incomplete subjects in nonsocial_data
incomplete_ids_nonsocial = nonsocial_data |>
  group_by(gameID) |>
  summarize(
    complete_trials = n()
  ) |>
  ungroup() |>
  filter(
    complete_trials < N_TRIALS
  ) |>
  arrange(complete_trials)
incomplete_ids_nonsocial

# Remove incomplete subjects
nonsocial_data = nonsocial_data |>
  filter(
    !gameID %in% unique(incomplete_ids_nonsocial$gameID)
  )


# Remove all trials in which people did not move the paddle

# Summarize: how many trials did subjects skip?
nonsocial_data |>
  group_by(gameID) |>
  summarize(
    num_skipped = sum(!paddle_intervention),
    num_trials = n(),
    prop_skipped = num_skipped / num_trials
  ) |>
  ungroup() |>
  arrange(desc(num_skipped)) |>
  print(n=nrow(nonsocial_data))

# Filter out trials with no paddle intervention
nonsocial_data = nonsocial_data |>
  filter(updatedPaddle)


# MODEL ----

DECAY_SLOPE = -2.0

# > Calculate variance of bot error ----

# Add squared bot error columns and initialize additional columns needed for calculating variance of bot error
trial_data = trial_data |>
  rowwise() |>
  mutate(
    prev_trials = trialInd, # since trialInd is 0-indexed, it also encodes the number of previous trials
    # bot error colums: radians
    bot_error_sq = bot_error^2, # squared bot error in each trial
    prev_bot_error_sq = 0, # initialize supporting column for calculating sum of squared bot error in previous trials
    sum_prev_bot_error_sq = 0, # initialize sum of squared bot error in previous trials
    var_prev_bot_error = 0, # initialize variance of bot error in previous trials
    # bot error colums: degrees
    bot_error_degrees_sq = bot_error_degrees^2, # squared bot error (degrees) in each trial
    prev_bot_error_degrees_sq = 0, # initialize supporting column for calculating sum of squared bot error (degrees) in previous trials
    sum_prev_bot_error_degrees_sq = 0, # initialize sum of squared bot error (degrees) in previous trials
    var_prev_bot_error_degrees = 0, # initialize variance of bot error (degrees) in previous trials
    # weighted error: initial trials
    early_sum_prev_bot_error_sq = 0, # *weighted* sum of squared bot error across previous (early) trials
    early_var_prev_bot_error = 0, # variance of *weighted* bot error across previous (early) trials
    early_sum_prev_bot_error_degrees_sq = 0, # *weighted* sum of squared bot error (degrees) across previous (early) trials
    early_var_prev_bot_error_degrees = 0, # variance of *weighted* bot error (degrees) across previous (early) trials
    # weighted error: recent trials
    late_sum_prev_bot_error_sq = 0, # *weighted* sum of squared bot error across previous (recent) trials
    late_var_prev_bot_error = 0, # variance of *weighted* bot error across previous (recent) trials
    late_sum_prev_bot_error_degrees_sq = 0, # *weighted* sum of squared bot error (degrees) across previous (recent) trials
    late_var_prev_bot_error_degrees = 0 # variance of *weighted* bot error (degrees) across previous (recent) trials
  )

# Fill in supporting columns for calculating sum of squared bot error across previous trials
trial_data = trial_data |>
  group_by(gameID) |>
  arrange(gameID, trialInd) |>
  mutate(
    prev_bot_error_sq = ifelse(is.na(lag(bot_error_sq, 1)), 0, lag(bot_error_sq, 1)),
    prev_bot_error_degrees_sq = ifelse(is.na(lag(bot_error_degrees_sq, 1)), 0, lag(bot_error_degrees_sq, 1))
  )

# Fill in sum of squared bot error across previous trials
trial_data = trial_data |>  # NB: ensure that data is arranged by subject, increasing trial index
  arrange(gameID, sort(trialInd))

for (subj in unique(trial_data$gameID)) { # NB: this takes 60-120s to run
  for (idx in sort(unique(trial_data$trialInd))) {
    trial_data$sum_prev_bot_error_sq[trial_data$gameID == subj & trial_data$trialInd == idx] = sum(trial_data[trial_data$gameID == subj,]$prev_bot_error_sq[1:(idx+1)])
    trial_data$sum_prev_bot_error_degrees_sq[trial_data$gameID == subj & trial_data$trialInd == idx] = sum(trial_data[trial_data$gameID == subj,]$prev_bot_error_degrees_sq[1:(idx+1)])
    # Fill in *weighted* sum of squared bot error across previous (early) trials
    trial_data$early_sum_prev_bot_error_sq[trial_data$gameID == subj & trial_data$trialInd == idx] = sum(trial_data[trial_data$gameID == subj,]$prev_bot_error_sq[1:(idx+1)]*(seq(from=0, to=idx)^DECAY_SLOPE), na.rm=T)
    trial_data$early_sum_prev_bot_error_degrees_sq[trial_data$gameID == subj & trial_data$trialInd == idx] = sum(trial_data[trial_data$gameID == subj,]$prev_bot_error_degrees_sq[1:(idx+1)]*(seq(from=0, to=idx)^DECAY_SLOPE), na.rm=T)
    # Fill in *weighted* sum of squared bot error across previous (recent) trials
    trial_data$late_sum_prev_bot_error_sq[trial_data$gameID == subj & trial_data$trialInd == idx] = sum(trial_data[trial_data$gameID == subj,]$prev_bot_error_sq[1:(idx+1)]*(rev(seq(from=1, to=idx+1)^DECAY_SLOPE)), na.rm=T)
    trial_data$late_sum_prev_bot_error_degrees_sq[trial_data$gameID == subj & trial_data$trialInd == idx] = sum(trial_data[trial_data$gameID == subj,]$prev_bot_error_degrees_sq[1:(idx+1)]*(rev(seq(from=1, to=idx+1)^DECAY_SLOPE)), na.rm=T)
  }
}

# Fill in variance of bot error across previous trials
trial_data = trial_data |>
  rowwise() |>
  mutate(
    var_prev_bot_error = ifelse(prev_trials < 2,
                                NA, # NB: prev error variance is undefined over first 2 trials
                                sum_prev_bot_error_sq / (prev_trials-1)),
    var_prev_bot_error_degrees = ifelse(prev_trials < 2,
                                        NA, # # NB: prev error variance is undefined over first 2 trials
                                        sum_prev_bot_error_degrees_sq / (prev_trials-1)),
    # Fill in *weighted* variance of bot error across previous (early) trials
    early_var_prev_bot_error = ifelse(prev_trials < 2,
                                      NA,
                                      early_sum_prev_bot_error_sq / sum(seq(from=1, to=prev_trials-1)^DECAY_SLOPE)),
    early_var_prev_bot_error_degrees = ifelse(prev_trials < 2,
                                              NA,
                                              early_sum_prev_bot_error_degrees_sq / sum(seq(from=1, to=prev_trials-1)^DECAY_SLOPE)),
    # Fill in *weighted* variance of bot error across previous (early) trials
    late_var_prev_bot_error = ifelse(prev_trials < 2,
                                     NA,
                                     late_sum_prev_bot_error_sq / sum(seq(from=1, to=prev_trials-1)^DECAY_SLOPE)),
    late_var_prev_bot_error_degrees = ifelse(prev_trials < 2,
                                             NA,
                                             late_sum_prev_bot_error_degrees_sq / sum(seq(from=1, to=prev_trials-1)^DECAY_SLOPE))
  )


# sanity checks
trial_data |>
  filter(gameID == subj) |>
  select(trialInd,
         prev_bot_error_degrees_sq, sum_prev_bot_error_degrees_sq, var_prev_bot_error_degrees,
         early_sum_prev_bot_error_degrees_sq, early_var_prev_bot_error_degrees,
         late_sum_prev_bot_error_degrees_sq, late_var_prev_bot_error_degrees
  ) |>
  print(n=nrow(trial_data))


# > Calculate variance of human error ----

# Calculate squared error on each trial
nonsocial_data = nonsocial_data |>
  rowwise() |>
  mutate(
    human_error_sq = human_error^2, # squared human error in each trial
    human_error_degrees_sq = human_error_degrees^2, # squared human error (degrees) in each trial
  )

# Calculate variance of human error on this trial index: sum of human squared error on this trial index / (number of participants-1)
nonsocial_data_trial_ind_summary = nonsocial_data |>
  group_by(trialInd) |>
  summarize(
    var_human_error_trial_idx = sum(human_error_sq) / (n()-1),
    var_human_error_degrees_trial_idx = sum(human_error_degrees_sq) / (n()-1)
  ) |>
  ungroup()

# Calculate average nonsocial response for each trial angle (used to generate 'internal estimate' below)
nonsocial_data_trial_angle_summary = nonsocial_data |>
  group_by(launching_rho, paddle_rho_original) |> # NB: paddle_rho_original is ground truth est. from simulation
  summarize(
    avg_nonsocial_trial_estimate = mean(final_participant_angle),
    avg_nonsocial_trial_estimate_degrees = rad_to_deg(avg_nonsocial_trial_estimate),
    sd_nonsocial_trial_estimate = sd(final_participant_angle), # NB: probably not necessary
    sd_nonsocial_trial_estimate_degrees = sd(rad_to_deg(final_participant_angle))
  ) |>
  ungroup()


# > Calculate variance weights ----

# Combine nonsocial trial index summary above with trial data
trial_data = trial_data |>
  inner_join(
    nonsocial_data_trial_ind_summary,
    by = c('trialInd')
  )

# Combine nonsocial trial angle summary above with trial data
trial_data = trial_data |>
  inner_join(
    nonsocial_data_trial_angle_summary,
    by = c('launching_rho', 'paddle_rho_original')
  )

# Calculate weights using the columns added above
# NB: below is redundant (human and bot weights sum to 1)
trial_data = trial_data |>
  rowwise() |>
  mutate(
    inv_var_prev_bot_error = 1 / var_prev_bot_error,
    inv_var_prev_bot_error_degrees = 1 / var_prev_bot_error_degrees,
    inv_var_human_error_trial_idx = 1 / var_human_error_trial_idx,
    inv_var_human_error_degrees_trial_idx = 1 / var_human_error_degrees_trial_idx,
    weight_participant_estimate = inv_var_human_error_trial_idx / sum(inv_var_human_error_trial_idx, inv_var_prev_bot_error),
    weight_participant_estimate_degrees = inv_var_human_error_degrees_trial_idx / sum(inv_var_human_error_degrees_trial_idx, inv_var_prev_bot_error_degrees),
    weight_bot_estimate = inv_var_prev_bot_error / sum(inv_var_human_error_trial_idx, inv_var_prev_bot_error),
    weight_bot_estimate_degrees = inv_var_prev_bot_error_degrees / sum(inv_var_human_error_degrees_trial_idx, inv_var_prev_bot_error_degrees),
    # calculate weights for *weighted* error across previous (early) trials
    early_inv_var_prev_bot_error = 1 / early_var_prev_bot_error,
    early_inv_var_prev_bot_error_degrees = 1 / early_var_prev_bot_error_degrees,
    early_weight_participant_estimate = inv_var_human_error_trial_idx / sum(inv_var_human_error_trial_idx, early_inv_var_prev_bot_error),
    early_weight_participant_estimate_degrees = inv_var_human_error_degrees_trial_idx / sum(inv_var_human_error_degrees_trial_idx, early_inv_var_prev_bot_error_degrees),
    early_weight_bot_estimate = early_inv_var_prev_bot_error / sum(inv_var_human_error_trial_idx, early_inv_var_prev_bot_error),
    early_weight_bot_estimate_degrees = early_inv_var_prev_bot_error_degrees / sum(inv_var_human_error_degrees_trial_idx, early_inv_var_prev_bot_error_degrees),
    # calculate weights for *weighted* error across previous (recent) trials
    late_inv_var_prev_bot_error = 1 / late_var_prev_bot_error,
    late_inv_var_prev_bot_error_degrees = 1 / late_var_prev_bot_error_degrees,
    late_weight_participant_estimate = inv_var_human_error_trial_idx / sum(inv_var_human_error_trial_idx, late_inv_var_prev_bot_error),
    late_weight_participant_estimate_degrees = inv_var_human_error_degrees_trial_idx / sum(inv_var_human_error_degrees_trial_idx, late_inv_var_prev_bot_error_degrees),
    late_weight_bot_estimate = late_inv_var_prev_bot_error / sum(inv_var_human_error_trial_idx, late_inv_var_prev_bot_error),
    late_weight_bot_estimate_degrees = late_inv_var_prev_bot_error_degrees / sum(inv_var_human_error_degrees_trial_idx, late_inv_var_prev_bot_error_degrees)
  )


# > Calculate weighted paddle placement estimates ----

# Combine agent partner estimate, participant estimate
trial_data = trial_data |>
  rowwise() |>
  mutate(
    # baseline model (weight based on variance)
    weighted_paddle_estimate_baseline = sum((weight_participant_estimate*avg_nonsocial_trial_estimate), (weight_bot_estimate*paddle_rho)),
    weighted_paddle_estimate_degrees_baseline = sum((weight_participant_estimate_degrees*avg_nonsocial_trial_estimate_degrees), (weight_bot_estimate_degrees*rad_to_deg(paddle_rho))),
    # participant estimate only
    weighted_paddle_estimate_participant = avg_nonsocial_trial_estimate,
    weighted_paddle_estimate_degrees_participant = avg_nonsocial_trial_estimate_degrees,
    # agent partner estimate only
    weighted_paddle_estimate_bot = paddle_rho,
    weighted_paddle_estimate_degrees_bot = rad_to_deg(paddle_rho),
    # weighted paddle estimates using *weighted* error across previous (early) trials
    weighted_paddle_estimate_early = sum((early_weight_participant_estimate*avg_nonsocial_trial_estimate), (early_weight_bot_estimate*paddle_rho)),
    weighted_paddle_estimate_degrees_early = sum((early_weight_participant_estimate_degrees*avg_nonsocial_trial_estimate_degrees), (early_weight_bot_estimate_degrees*rad_to_deg(paddle_rho))),
    # weighted paddle estimates using *weighted* error across previous (late) trials
    weighted_paddle_estimate_late = sum((late_weight_participant_estimate*avg_nonsocial_trial_estimate), (late_weight_bot_estimate*paddle_rho)),
    weighted_paddle_estimate_degrees_late = sum((late_weight_participant_estimate_degrees*avg_nonsocial_trial_estimate_degrees), (late_weight_bot_estimate_degrees*rad_to_deg(paddle_rho)))
  )


# > Calculate error of weighted paddle placement estimates ----

trial_data = trial_data |>
  rowwise() |>
  mutate(
    # baseline model (weight based on variance)
    weighted_paddle_estimate_error_baseline = get_angle_difference(final_participant_angle, weighted_paddle_estimate_baseline),
    weighted_paddle_estimate_error_degrees_baseline = rad_to_deg(weighted_paddle_estimate_error_baseline),
    # participant estimate only
    weighted_paddle_estimate_error_participant = get_angle_difference(final_participant_angle, weighted_paddle_estimate_participant),
    weighted_paddle_estimate_error_degrees_participant = rad_to_deg(weighted_paddle_estimate_error_participant),
    # agent partner estimate only
    weighted_paddle_estimate_error_bot = get_angle_difference(final_participant_angle, weighted_paddle_estimate_bot),
    weighted_paddle_estimate_error_degrees_bot = rad_to_deg(weighted_paddle_estimate_error_bot),
    # weighted paddle estimates using *weighted* error across previous (early) trials
    weighted_paddle_estimate_error_early = get_angle_difference(final_participant_angle, weighted_paddle_estimate_early),
    weighted_paddle_estimate_error_degrees_early = rad_to_deg(weighted_paddle_estimate_error_early),
    # weighted paddle estimates using *weighted* error across previous (late) trials
    weighted_paddle_estimate_error_late = get_angle_difference(final_participant_angle, weighted_paddle_estimate_late),
    weighted_paddle_estimate_error_degrees_late = rad_to_deg(weighted_paddle_estimate_error_late)
  )


# FIGURE: model comparison ----

# Figure: compare error across models
model_comparison_fig = trial_data |>
  select(
    gameID, trialInd, trial_block, condition_str,
    weighted_paddle_estimate_error_degrees_bot,
    weighted_paddle_estimate_error_degrees_participant,
    weighted_paddle_estimate_error_degrees_baseline,
    weighted_paddle_estimate_error_degrees_early,
    weighted_paddle_estimate_error_degrees_late
  ) |>
  filter( # NB: these are all NA in same trials (first 2 rows for each subject)
    !is.na(weighted_paddle_estimate_error_degrees_baseline),
    !is.na(weighted_paddle_estimate_error_degrees_early),
    !is.na(weighted_paddle_estimate_error_degrees_late)
  ) |>
  pivot_longer(
    cols = starts_with('weighted_paddle_estimate_error'),
    names_to = c('.value', 'model'),
    names_pattern = 'weighted_paddle_estimate_error_(.*)_(.*)'
  ) |>
  mutate(
    model = factor(
      model,
      levels = c(
        'bot',
        'participant',
        'baseline',
        'early',
        'late'
      )
    )
  ) |>
  ggplot(
    aes(x = model, y = degrees)
  ) +
  stat_summary(
    fun.data = 'mean_cl_boot',
    geom = 'pointrange'
  ) +
  scale_x_discrete(
    name = element_blank(),
    labels = c(
      'bot' = 'partner \nonly',
      'participant' = 'participant \nonly',
      'baseline' = 'cue \nintegration',
      'early' = 'cue \nintegration \n(early trials)',
      'late' = 'cue \nintegration \n(late trials)'
    )
  ) +
  scale_y_continuous(
    name = 'error (degrees)',
    breaks = seq(10, 20, by = 2),
    labels = seq(10, 20, by = 2),
    # limits = c(10, 20)
  ) +
  PLOT_THEME
# Save figure
model_comparison_fig
# ggsave(
#   model_comparison_fig,
#   filename = 'raw_model_comparison_summary.pdf',
#   path = FIGURE_PATH,
#   device = cairo_pdf,
#   width = 12,
#   height = 7,
# )


# Figure: compare error across models (separate by condition)
model_comparison_condition_fig = trial_data |>
  select(
    gameID, trialInd, trial_block, condition_str,
    weighted_paddle_estimate_error_degrees_bot,
    weighted_paddle_estimate_error_degrees_participant,
    weighted_paddle_estimate_error_degrees_baseline,
    weighted_paddle_estimate_error_degrees_early,
    weighted_paddle_estimate_error_degrees_late
  ) |>
  filter( # NB: these are all NA in same trials (first 2 rows for each subject)
    !is.na(weighted_paddle_estimate_error_degrees_baseline),
    !is.na(weighted_paddle_estimate_error_degrees_early),
    !is.na(weighted_paddle_estimate_error_degrees_late)
  ) |>
  pivot_longer(
    cols = starts_with('weighted_paddle_estimate_error'),
    names_to = c('.value', 'model'),
    names_pattern = 'weighted_paddle_estimate_error_(.*)_(.*)'
  ) |>
  mutate(
    model = factor(
      model,
      levels = c(
        'bot',
        'participant',
        'baseline',
        'early',
        'late'
      )
    )
  ) |>
  ggplot(
    aes(x = model, y = degrees)
  ) +
  stat_summary(
    aes(fill = condition_str),
    color = 'black',
    fun.data = 'mean_cl_boot',
    geom = 'bar',
    position = position_dodge(width = 0.9),
    alpha = 0.75,
    width = 0.75,
    # show.legend = F
  ) +
  stat_summary(
    # NB: aes throws a warning but we need this to align the errorbars without changing color
    aes(fill = condition_str),
    fun.data = 'mean_cl_boot',
    geom = 'errorbar',
    position = position_dodge(width = 0.9),
    width = 0,
    show.legend = F
  ) +
  stat_summary(
    aes(x = model, y = degrees, color = 'average'),
    fun.data = 'mean_cl_boot',
    geom = 'pointrange',
    size = 2,
    # linewidth = 1,
    show.legend = F
  ) +
  scale_x_discrete(
    name = element_blank(),
    labels = c(
      'bot' = 'partner \nonly',
      'participant' = 'participant \nonly',
      'baseline' = 'cue \nintegration',
      'early' = 'cue \nintegration \n(early trials)',
      'late' = 'cue \nintegration \n(late trials)'
    )
  ) +
  scale_y_continuous(
    name = 'error (degrees)',
    breaks = seq(0, 40, by = 10),
    labels = seq(0, 40, by = 10),
  ) +
  scale_color_manual(
    name = element_blank(),
    values = append(COLORS, c('average' = 'black'))
  ) +
  scale_fill_manual(
    name = element_blank(),
    values = COLORS
  ) +
  PLOT_THEME

# Save figure
model_comparison_condition_fig
ggsave(
  # NB: save call removes legend
  model_comparison_condition_fig + theme(legend.position = 'none'),
  filename = 'raw_model_comparison_by_condition.pdf',
  path = FIGURE_PATH,
  device = cairo_pdf,
  width = 12,
  height = 8,
)








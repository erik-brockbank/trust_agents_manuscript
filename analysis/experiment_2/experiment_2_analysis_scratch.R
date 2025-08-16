
#'
#' Scratch file for experiment 2 analyses
#' TODO move the final versions of these analyses to a clean script
#'


# INIT ----
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(emmeans)
library(lme4)
library(LaplacesDemon) # invlogit function for interpreting lmer model coeffs
library(pwr)
library(tidyverse)


# GLOBALS ----

# Directories
DATA_DIR = '../../data/experiment_2' # folder containing data files
DATA_FILE = 'trust_agents_trialdata_bias_complete.csv' # name of file containing full dataset
SURVEY_FILE = 'trust_agents_surveydata_bias_complete.csv' # name of file containing survey responses

# Experiment parameters
TRIAL_BLOCKS = 8
CONDITION_LOOKUP = c(
  'accurate' = 'no bias',
  'heavy' = 'heavy',
  'light' = 'light'
)
CONDITION_ORDER = c('heavy', 'light', 'no bias')


# ANALYSIS FUNCTIONS ----

# Determine the participant's final paddle angle on each trial
# Uses either the agent's suggestion or the final angle in the sequence of
# captured paddle angles
get_final_angle = function(paddle_rho, paddle_tr) {
  paddle_tr = strsplit(paddle_tr, "\\[")[[1]][2]
  paddle_tr = strsplit(paddle_tr, "\\]")[[1]][1]
  paddle_tr = str_replace_all(paddle_tr, " ", "")
  if(paddle_tr == "") {
    return(as.numeric(paddle_rho))
  } else {
    paddle_vec = strsplit(paddle_tr, ",")[[1]]
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
survey_data = read_csv(paste(DATA_DIR, SURVEY_FILE, sep = '/'))
glimpse(trial_data)
glimpse(survey_data)

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


# ANALYSIS: Demographics ----

# N participants
length(unique(survey_data$gameID))
length(unique(trial_data$gameID))
setdiff(unique(survey_data$gameID), unique(trial_data$gameID)) # should be empty

# Age
summary(survey_data$age)
sd(survey_data$age)

# Gender
table(survey_data$gender)

# Education
table(survey_data$edu)

# Completion time
survey_data |>
  group_by(gameID) |>
  mutate(time_elapsed_min = (time_elapsed / 1000) / 60) |>
  ungroup() |>
  summarize(
    mean_time = mean(time_elapsed_min),
    sd_time = sd(time_elapsed_min)
  )

# N participants per condition
trial_data |>
  group_by(condition_str) |>
  summarize(
    subjects = n_distinct(gameID)
  )


# ANALYSIS: Power ----

# > Functions ----

# Simulate a person responding on 8 Mystery Round trials with correct probability `prob`
person_sim = function(prob) {
  rbinom(n = 8, size = 1, prob = prob)
}

# Simulate `n` participants responding on 8 Mystery Round trials with correct probability `prob`
exp_sim = function(n, prob) {
  replicate(n, person_sim(prob))
}

# Return the proportion of correct responses for simulated subjects in `responses`
# NB: `responses` is output of `exp_sim` above
get_exp_probs = function(responses) {
  prob_fn = function(x) {
    sum(x) / length(x)
  }
  apply(responses, 2, prob_fn)
}

# Run t-test for a set of correct response proportions from a simulated set of subjects
# Return p-value for easy significance check
# NB: `response_proportions` is output of `get_exp_probs` above
run_exp_stats = function(response_proportions) {
  t.test(
    response_proportions,
    mu = 0.5,
    alternative = 'greater' # NB: can toggle comment this line
  )$p.value
}

# Wrapper for getting p value from t-test of simulated experiment responses
analysis_sim = function(n, prob) {
  run_exp_stats(get_exp_probs(exp_sim(n, prob)))
}

# Wrapper for running repeated experiment simulations + t-tests
# Run `k_sims` experiment simulations with `n_subj` participants,
# Mystery Round correct probability `prob`
# Returns the simulated power for those repeated experiment simulations
# This is the number of simulations in which p < .05 / total simulations
power_sim = function(n_subj, prob, k_sims) {
  p_vals = replicate(k_sims, analysis_sim(n_subj, prob))
  # How often do we reject the null out of total number of tests?
  sum(p_vals < 0.05) / length(p_vals)
}



# > Results ----

N = 50
P = 0.6
SIMS = 10000 # 10000 => 5-10s


# TESTING
# Test individual experiment functions
simulated_responses = exp_sim(N, P) # simulated set of trial responses
response_probabilities = get_exp_probs(simulated_responses) # proportions for those responses
run_exp_stats(response_probabilities) # p-value for t-test based on these proportions
# Test wrapper for individual experiment functions
analysis_sim(N, P) # p-value for t-test based on simulated results


# SIMULATION
# Below we run `SIMS` iterations of the experiment simulation,
# with `N` subjects per experiment,
# and `P` probability of each subject choosing correctly on each mystery round
# Returns the power to detect true probability `P`
power_sim(N, P, SIMS)

# Summary of results:
# N = 50, P = 0.6 => power = approx. 0.99


# Quantifying effect sizes for P values based on above
# Reference: https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html
ES.h(p1 = P, p2 = 0.5)
# Summary: P = 0.6 has a Cohen's h of 0.2, a "small" effect size


# ANALYSIS: Intervention rate ----

# Role of condition and condition * block interaction on intervention rates
m_rand = glmer(
  as.numeric(paddleIntervention) ~ (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  family = 'binomial',
)

m_condition = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  family = 'binomial',
)

m_block = glmer(
  as.numeric(paddleIntervention) ~ sessionBlock + (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  family = 'binomial',
)

m_block_condition = glmer(
  as.numeric(paddleIntervention) ~ sessionBlock + as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  family = 'binomial',
)

m_block_condition_int = glmer(
  as.numeric(paddleIntervention) ~ sessionBlock * as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  family = 'binomial',
)


# Main effect of condition
anova(m_condition, m_rand, test = 'LRT')
emmeans(m_condition, 'condition_str')
# convert emmeans to probability space: mean + 95% CI
invlogit(2.4792); invlogit(2.058); invlogit(2.900) # heavy
invlogit(2.1924); invlogit(1.816); invlogit(2.569)  # light
invlogit(-0.0759); invlogit(-0.439); invlogit(0.288) # no bias
# Compare between conditions
emmeans(m_condition, specs = pairwise ~ condition_str)

# Effect of session block, session block * condition interaction
anova(m_block_condition_int, m_block_condition, m_block, m_rand, test = 'LRT')
# Comparing slopes on session block for each condition (slope estimate for each condition at top, contrast below)
emtrends(m_block_condition_int, specs = pairwise ~ condition_str, var = 'sessionBlock')


# ANALYSIS: Relationship between bot error and intervention magnitude ----

# Motivation: in this experiment, people need to make continuous corrections to counteract bias,
# not just fixed adjustment. Did they do this?

m_rand = lmer(
  intervene_dist_degrees ~ (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  REML = F,
  control = lmerControl(optimizer = 'Nelder_Mead')
)

m_error = lmer(
  intervene_dist_degrees ~ bot_error_degrees + (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  REML = F,
)

m_condition_error = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) + bot_error_degrees + (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  REML = F,
  control = lmerControl(optimizer = 'Nelder_Mead')
)

m_condition_error_int = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) * bot_error_degrees + (1+sessionBlock|gameID),
  data = trial_data |> filter(!criticalTrial), # excluse mystery round trials
  REML = F,
  control = lmerControl(optimizer = 'Nelder_Mead')
)

# Did condition * bot error interaction improve model fit?
anova(m_condition_error_int, m_condition_error, m_error, m_rand, test = 'LRT')
# Compare bot error slopes across conditions
# (Top: slopes > 0 in each condition? Bottom: difference in slopes across conditions)
# NB: the slope estimates vary wildly across conditions...
emtrends(m_condition_error_int, specs = pairwise ~ condition_str, var = 'bot_error_degrees')


# ANALYSIS: Mystery round intervention rate ----

m_rand = glmer(
  as.numeric(paddleIntervention) ~ (1|gameID), # incl. launching rho as "item" intercept
  data = trial_data |> filter(criticalTrial),
  family = 'binomial'
)

m_condition = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  family = 'binomial'
)


# Do conditions differ in their intervention rates on critical trials?
anova(m_condition, m_rand, test = 'LRT')
# What is the estimated intervention rate in each condition?
emmeans(m_condition, 'condition_str')
# Convert from logit space to probability space
invlogit(2.41); invlogit(1.73); invlogit(3.096) # heavy
invlogit(2.02); invlogit(1.41); invlogit(2.634) # light
invlogit(-1.33); invlogit(-1.95); invlogit(-0.704) # no bias
# Which contrasts are significant?
emmeans(m_condition, specs = pairwise ~ condition_str)


# Do mystery round intervention rates vary over the experiment?
m_block = glmer(
  as.numeric(paddleIntervention) ~ sessionBlock + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  family = 'binomial',
)

m_condition_block = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) + sessionBlock + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  family = 'binomial'
)

m_condition_block_int = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) * sessionBlock + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  family = 'binomial',
)


# Interaction between condition and session block
anova(m_condition_block_int, m_condition_block, m_block, m_rand, test = 'LRT')
# Comparing slopes on session block for each condition (slope estimate for each condition at top, contrast below)
emtrends(m_condition_block_int, specs = pairwise ~ condition_str, var = 'sessionBlock')


# ANALYSIS: Mystery round signed intervention magnitude ----

# TODO move the function here up top, add to trial_data with the other new columns
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


# Models
m_rand = lmer(
  signed_intervention_distance_degrees ~ (1|gameID),
  data = trial_data |> filter(criticalTrial),
  REML = F,
)

m_condition = lmer(
  signed_intervention_distance_degrees ~ as.factor(condition_str) + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  REML = F,
)

m_block = lmer(
  signed_intervention_distance_degrees ~ sessionBlock + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  REML = F,
)

m_block_condition = lmer(
  signed_intervention_distance_degrees ~ sessionBlock + as.factor(condition_str) + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  REML = F,
)

m_block_condition_int = lmer(
  signed_intervention_distance_degrees ~ sessionBlock * as.factor(condition_str) + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  REML = F,
)

# Main effect of condition
anova(m_condition, m_rand)
emmeans(m_condition, ~ condition_str)
# Compare condition values
emmeans(m_condition, specs = pairwise ~ condition_str)
# Interaction between condition and session block
anova(m_block_condition_int, m_block_condition, m_block, m_rand, test = 'LRT')
# Comparing slopes on session block for each condition (slope estimate for each condition at top, contrast below)
emtrends(m_block_condition_int, specs = pairwise ~ condition_str, var = 'sessionBlock')



# ANALYSIS: Mystery round relationship between bot error and intervention magnitude ----

# Motivation: can people make continuous corrections to counteract bias, *even in mystery round trials*??

m_rand = lmer(
  intervene_dist_degrees ~ (1|gameID),
  data = trial_data |> filter(criticalTrial), # mystery round trials only
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
)

m_error = lmer(
  intervene_dist_degrees ~ bot_error_degrees + (1|gameID),
  data = trial_data |> filter(criticalTrial), # mystery round trials only
  REML = F,
)

m_condition_error = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) + bot_error_degrees + (1|gameID),
  data = trial_data |> filter(criticalTrial), # mystery round trials only
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
)

m_condition_error_int = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) * bot_error_degrees + (1|gameID),
  data = trial_data |> filter(criticalTrial), # mystery round trials only
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
)

# Did condition * bot error interaction improve model fit?
anova(m_condition_error_int, m_condition_error, m_error, m_rand, test = 'LRT')
# Compare bot error slopes across conditions
# (Top: slopes > 0 in each condition? Bottom: difference in slopes across conditions)
# NB: the slope estimates vary wildly across conditions...
emtrends(m_condition_error_int, specs = pairwise ~ condition_str, var = 'bot_error_degrees')




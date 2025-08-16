
#'
#' Scratch file for experiment 1 analyses
#' TODO move the final versions of these analyses to a clean script
#'


# INIT ----
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(emmeans)
library(lme4)
library(LaplacesDemon) # invlogit function for interpreting lmer model coeffs
library(tidyverse)


# GLOBALS ----

# Directories
DATA_DIR = '../../data/experiment_1' # folder containing data files
DATA_FILE = 'physics_continual_learning_trialdata_e2_full.csv' # name of file containing full dataset
SURVEY_FILE = 'physics_continual_learning_surveydata_e2_full.csv' # name of file containing survey responses

# Experiment parameters
TRIAL_BLOCKS = 8
CONDITION_LOOKUP = c(
  'bad' = 'unreliable',
  'good' = 'reliable',
  'improve' = 'improving',
  'worsen' = 'worsening'
)
CONDITION_ORDER = c('unreliable', 'reliable', 'improving', 'worsening')


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
    signed_human_error = get_signed_participant_error(paddle_rho, final_par, groundtruthAngle),
    signed_human_error_degrees = rad_to_deg(signed_human_error)
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


# ANALYSIS: Trial summary ----

# Percent correct overall
trial_data |>
  group_by(gameID) |>
  summarize(subj_acc = sum(correct) / n()) |>
  ungroup() |>
  summarize(
    mean(subj_acc),
    sd(subj_acc)
  )

# Percent correct across trial blocks
trial_data |>
  group_by(gameID, trial_block) |>
  summarize(subj_block_acc = sum(correct) / n()) |>
  ungroup() |>
  group_by(trial_block) |>
  summarize(
    mean(subj_block_acc),
    sd(subj_block_acc)
  )

# RMSE overall
trial_data |>
  group_by(gameID) |>
  summarize(subj_rmse = sqrt(mean(human_error_degrees^2))) |>
  ungroup() |>
  summarize(
    mean(subj_rmse),
    sd(subj_rmse)
  )

# TODO delete if not reporting condition-level differences
# Accuracy differs by condition
# summary(aov(
#   subj_acc ~ factor(condition_str),
#   data = trial_data |> group_by(gameID, condition_str) |> summarize(subj_acc = sum(correct) / n()) |> ungroup()
# ))
# trial_data |> group_by(gameID, condition_str) |> summarize(subj_acc = sum(correct) / n()) |> ungroup() |>
#   ggplot(aes(x = condition_str, y = subj_acc)) +
#   geom_jitter() +
#   DEFAULT_PLOT_THEME
# # RMSE differs by condition
# summary(aov(
#   subj_rmse ~ factor(condition_str),
#   data = trial_data |> group_by(gameID, condition_str) |> summarize(subj_rmse = sqrt(mean(human_error_degrees^2))) |> ungroup()
# ))
# trial_data |> group_by(gameID, condition_str) |> summarize(subj_rmse = sqrt(mean(human_error_degrees^2))) |> ungroup() |>
#   ggplot(aes(x = condition_str, y = subj_rmse)) +
#   geom_jitter() +
#   DEFAULT_PLOT_THEME


# ANALYSIS: Intervention rates ----

# Role of condition and condition * block interaction on intervention rates
# NB: models take 5-10s to fit
m_rand = glmer(
  as.numeric(paddleIntervention) ~ (1+sessionBlock|gameID),
  data = trial_data,
  family = 'binomial',
  control = glmerControl(optimizer = 'Nelder_Mead')
)

m_condition = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data,
  family = 'binomial',
  control = glmerControl(optimizer = 'Nelder_Mead')
)

m_block = glmer(
  as.numeric(paddleIntervention) ~ sessionBlock + (1+sessionBlock|gameID),
  data = trial_data,
  family = 'binomial',
  control = glmerControl(optimizer = 'Nelder_Mead')
)

m_block_condition = glmer(
  as.numeric(paddleIntervention) ~ sessionBlock + as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data,
  family = 'binomial',
  control = glmerControl(optimizer = 'Nelder_Mead')
)

m_block_condition_int = glmer(
  as.numeric(paddleIntervention) ~ sessionBlock * as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data,
  family = 'binomial',
  # control = glmerControl(optimizer = 'Nelder_Mead')
  control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
)


# Main effect of condition
anova(m_condition, m_rand, test = 'LRT')
emmeans(m_condition, 'condition_str')
# convert emmeans to probability space: mean + 95% CI
invlogit(2.789); invlogit(2.488); invlogit(3.09) # unreliable
invlogit(0.982); invlogit(0.692); invlogit(1.27)  # reliable
invlogit(1.948); invlogit(1.665); invlogit(2.23) # improving
invlogit(1.670); invlogit(1.288); invlogit(2.05) # worsening
# Compare between conditions
emmeans(m_condition, specs = pairwise ~ condition_str)

# Effect of session block, session block * condition interaction
anova(m_block_condition_int, m_block_condition, m_block, m_rand, test = 'LRT')
# Comparing slopes on session block for each condition (slope estimate for each condition at top, contrast below)
emtrends(m_block_condition_int, specs = pairwise ~ condition_str, var = 'sessionBlock')


# ANALYSIS: Signed error magnitude ----

m_rand = lmer(
  signed_human_error_degrees ~ (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead', optCtrl = list(maxfun = 1000000))
  # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  control = lmerControl(optimizer = 'bobyqa')
  )

m_condition = lmer(
  signed_human_error_degrees ~ as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  control = lmerControl(optimizer = 'Nelder_Mead')
  # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
)

m_block = lmer(
  signed_human_error_degrees ~ sessionBlock + (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
  control = lmerControl(optimizer = 'bobyqa')
  # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
)

m_block_condition = lmer(
  signed_human_error_degrees ~ sessionBlock + as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
  control = lmerControl(optimizer = 'bobyqa')
  # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
)

m_block_condition_int = lmer(
  signed_human_error_degrees ~ sessionBlock * as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
  control = lmerControl(optimizer = 'bobyqa')
  # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
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



# ANALYSIS: Critical trial interventions ----

# Intervention rate
m_rand = glmer(
  as.numeric(paddleIntervention) ~ (1|gameID), # TODO try including slope here?? Or (1|launching_rho) as "item" intercept?
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
invlogit(2.12); invlogit(1.747); invlogit(2.48) # unreliable
invlogit(1.35); invlogit(0.995); invlogit(1.70) # reliable
invlogit(2.24); invlogit(1.857); invlogit(2.63) # improving
invlogit(1.49); invlogit(1.128); invlogit(1.85) # worsening
# Which contrasts are significant?
emmeans(m_condition, specs = pairwise ~ condition_str)


# TODO delete if not reporting 1st and 2nd half critical trial interventions
# Do critical trial intervention rates vary over the experiment?
# m_block = glmer(
#   as.numeric(paddleIntervention) ~ sessionBlock + (1|gameID) + (1|launching_rho),
#   data = trial_data |> filter(criticalTrial),
#   family = 'binomial'
# )
# m_condition_block = glmer(
#   as.numeric(paddleIntervention) ~ as.factor(condition_str) + sessionBlock + (1|gameID) + (1|launching_rho),
#   data = trial_data |> filter(criticalTrial),
#   family = 'binomial'
# )
# m_condition_block_int = glmer(
#   as.numeric(paddleIntervention) ~ as.factor(condition_str) * sessionBlock + (1|gameID) + (1|launching_rho),
#   data = trial_data |> filter(criticalTrial),
#   family = 'binomial',
#   # control = glmerControl(optimizer = 'Nelder_Mead', optCtrl = list(maxfun = 1000000))
#   control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
# )

# Change over time
# anova(m_block, m_rand, test = 'LRT')
# summary(m_block)
# Interaction between condition and session block
# anova(m_condition_block_int, m_condition_block, m_condition, m_block, m_rand, test = 'LRT')
# Comparing slopes on session block for each condition (slope estimate for each condition at top, contrast below)
# emtrends(m_condition_block_int, specs = pairwise ~ condition_str, var = 'sessionBlock')

# Do the changes in intervention rate over time remove the differences between conditions?
# Compare critical trial intervention rates across conditions for first and second half of trials
# m_rand_first = glmer(
#   as.numeric(paddleIntervention) ~ (1|gameID),
#   data = trial_data |> filter(criticalTrial) |> filter(sessionBlock <= 4),
#   family = 'binomial'
# )
# m_rand_second = glmer(
#   as.numeric(paddleIntervention) ~ (1|gameID),
#   data = trial_data |> filter(criticalTrial) |> filter(sessionBlock > 4),
#   family = 'binomial',
# )
#
# m_condition_first = glmer(
#   as.numeric(paddleIntervention) ~ as.factor(condition_str) + (1|gameID),
#   data = trial_data |> filter(criticalTrial) |> filter(sessionBlock <= 4),
#   family = 'binomial'
# )
# m_condition_second = glmer(
#   as.numeric(paddleIntervention) ~ as.factor(condition_str) + (1|gameID),
#   data = trial_data |> filter(criticalTrial) |> filter(sessionBlock > 4),
#   family = 'binomial'
# )
#
# anova(m_condition_first, m_rand_first, test = 'LRT')
# anova(m_condition_second, m_rand_second, test = 'LRT')
# emmeans(m_condition_first, specs = pairwise ~ condition_str)$contrasts
# emmeans(m_condition_second, specs = pairwise ~ condition_str)$contrasts
#


# Intervention magnitude
m_rand = lmer(
  intervene_dist_degrees ~ (1|gameID), # TODO try including slope here?? Or (1|launching_rho) as "item" intercept?
  data = trial_data |> filter(criticalTrial),
  REML = F
)

m_condition = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) + (1|gameID),
  data = trial_data |> filter(criticalTrial),
  REML = F
)

# Do conditions differ in their intervention magnitudes on critical trials?
anova(m_condition, m_rand, test = 'LRT')
# What is the estimated intervention magnitude in each condition?
emmeans(m_condition, 'condition_str')
# Which contrasts are significant?
emmeans(m_condition, specs = pairwise ~ condition_str)$contrasts



# ANALYSIS: Relationship between bot error and intervention magnitude across conditions ----

# NB: this analysis *seeds* the critical trial analyses below
# First: interaction between condition and bot error in predicting intervention magnitude
# (i.e., did the relationship between bot error and intervention magnitude differ across conditions?)
m_rand = lmer(
  intervene_dist_degrees ~ (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
  control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  # control = lmerControl(optimizer = 'bobyqa')
)

m_condition = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) + (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
  control = lmerControl(optimizer = 'bobyqa')
)

m_condition_error = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) + bot_error_degrees + (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  # control = lmerControl(optimizer = 'Nelder_Mead')
  # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  control = lmerControl(optimizer = 'bobyqa')
)

m_condition_error_int = lmer(
  intervene_dist_degrees ~ as.factor(condition_str) * bot_error_degrees + (1+sessionBlock|gameID),
  data = trial_data,
  REML = F,
  control = lmerControl(optimizer = 'Nelder_Mead')
)

# Did condition * bot error interaction improve model fit?
anova(m_condition_error_int, m_condition_error, m_condition, m_rand, test = 'LRT')
# Compare bot error slopes across conditions
# (Top: slopes > 0 in each condition? Bottom: difference in slopes across conditions)
# NB: the slope estimates vary wildly across conditions...
emtrends(m_condition_error_int, specs = pairwise ~ condition_str, var = 'bot_error_degrees')




# Do we obtain similar results for intervention RATES?
# NB: these models take 5-10s to fit
m_rand = glmer(
  as.numeric(paddleIntervention) ~ (1+sessionBlock|gameID),
  # data = trial_data |> filter(paddleIntervention),
  data = trial_data,
  family = 'binomial',
  # control = glmerControl(optimizer = 'Nelder_Mead')
  # control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  # control = glmerControl(optimizer = 'bobyqa')
)

m_condition = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) + (1+sessionBlock|gameID),
  # data = trial_data |> filter(paddleIntervention),
  data = trial_data,
  family = 'binomial',
  # control = glmerControl(optimizer = 'Nelder_Mead')
  # control = glmerControl(optimizer = 'bobyqa')
)

m_condition_error = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) + bot_error_degrees + (1+sessionBlock|gameID),
  # data = trial_data |> filter(paddleIntervention),
  data = trial_data,
  family = 'binomial',
  # control = glmerControl(optimizer = 'Nelder_Mead')
  # control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  # control = glmerControl(optimizer = 'bobyqa')
)

m_condition_error_int = glmer(
  as.numeric(paddleIntervention) ~ as.factor(condition_str) * bot_error_degrees + (1+sessionBlock|gameID),
  # data = trial_data |> filter(paddleIntervention),
  data = trial_data,
  family = 'binomial',
  # control = glmerControl(optimizer = 'Nelder_Mead')
  # control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))
  control = glmerControl(optimizer = 'bobyqa')
)

# Did condition * bot error interaction improve model fit?
anova(m_condition_error_int, m_condition_error, m_condition, m_rand, test = 'LRT')
# Compare bot error slopes across conditions
# (Top: slopes > 0 in each condition? Bottom: difference in slopes across conditions)
# NB: the slope estimates vary wildly across conditions...
emtrends(m_condition_error_int, specs = pairwise ~ condition_str, var = 'bot_error_degrees')





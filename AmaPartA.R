source("https://edin.ac/4boOrMs")
get_my_data("B293230")

# Packages ────────────────────────────────────────────────────────────────
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)

# Data Wrangling ──────────────────────────────────────────────────────────
cat("N participants before cleaning:", n_distinct(gardenpath$participant), "\n")

# Check for invalid entries
n_invalid_age <- gardenpath |>
  filter(age > 25) |>
  distinct(participant) |>
  nrow()
cat("Number of participants with age > 25:", n_invalid_age, "\n")

# Identify participants with valid number of timepoints
valid_participants <- gardenpath |>
  group_by(participant) |>
  summarise(n_trials = n(), .groups = "drop") |>
  filter(n_trials <= 7) |>
  pull(participant) |>
  as.character()

# Filter dataset
gardenpath <- gardenpath |>
  filter(
    as.character(participant) %in% valid_participants,
    age <= 25,
    !is.na(RSES)
  )
gardenpath <- gardenpath |>
  mutate(
    gtype = factor(gtype, levels = 1:3, labels = c("Paired", "Community", "Solitary")),
    prevexp = factor(prevexp, levels = 0:1, labels = c("No", "Yes")),
    allotment = factor(allotment),
    participant = factor(participant)
  )

cat("N participants after cleaning:", n_distinct(gardenpath$participant), "\n")
cat("Total RSES observations:", sum(!is.na(gardenpath$RSES)), "\n")

# Summary of timepoints per participant
final_trials_per_person <- gardenpath |>
  group_by(participant) |>
  summarise(n_trials = n(), .groups = "drop") |>
  count(n_trials) |>
  rename("Number of Timepoints" = n_trials,
         "Number of Participants" = n)
print(final_trials_per_person)

# Descriptives ────────────────────────────────────────────────────────────
desc <- gardenpath |>
  group_by(gtype, months) |>
  summarise(
    mean_RSES = mean(RSES, na.rm = TRUE),
    sd_RSES   = sd(RSES, na.rm = TRUE),
    n         = n(),
    se        = sd_RSES / sqrt(n),
    .groups = "drop"
  )
print(desc, n = 21)

# Visualisation ───────────────────────────────────────────────────────────
# Plot 1: Individual trajectories
ggplot(gardenpath, aes(x = months, y = RSES, group = participant, colour = gtype)) +
  geom_line(alpha = 0.15) +
  stat_summary(aes(group = gtype), fun = mean, geom = "line",
               linewidth = 1.5, colour = "black", linetype = "dashed") +
  scale_colour_manual(values = c("Paired" = "purple", 
                                 "Community" = "darkorange", 
                                 "Solitary" = "green")) +
  facet_wrap(~gtype) +
  scale_x_continuous(breaks = c(0,3,6,9,12,15,18)) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(
    title = "Individual RSES Trajectories by Gardening Session Type",
    subtitle = "Dashed line = group mean trajectory",
    x = "Months Since Baseline",
    y = "RSES Score",
    colour = "Session Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # legend redundant with facet labels

# Plot 2: Mean trajectories with SE (cleaner than ribbons)
ggplot(gardenpath, aes(x = months, y = RSES, colour = gtype)) +
  scale_colour_manual(values = c("Paired" = "purple", 
                                 "Community" = "darkorange", 
                                 "Solitary" = "green")) +
  stat_summary(geom = "pointrange") +
  stat_summary(geom = "line") +
  scale_x_continuous(breaks = c(0,3,6,9,12,15,18)) +
  labs(
    title = "Mean RSES Over Time by Gardening Session Type",
    x = "Months Since Baseline",
    y = "Mean RSES Score",
    colour = "Session Type"
  ) +
  theme_minimal()

# Plot 3: By allotment (from your classmate - useful to check clustering)
ggplot(gardenpath, aes(x = months, y = RSES, colour = gtype)) +
  scale_colour_manual(values = c("Paired" = "purple", 
                                 "Community" = "darkorange", 
                                 "Solitary" = "green")) +
  stat_summary(geom = "pointrange") +
  stat_summary(geom = "line") +
  facet_wrap(~allotment) +
  scale_x_continuous(breaks = c(0,3,6,9,12,15,18)) +
  labs(
    title = "Mean RSES Over Time by Allotment",
    x = "Months Since Baseline",
    y = "Mean RSES Score",
    colour = "Session Type"
  ) +
  theme_minimal()

# Model Building ──────────────────────────────────────────────────────────
m1 <- lmer(RSES ~ months + (1 + months | participant) + (1 | allotment),
           data = gardenpath, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa"))

m2 <- lmer(RSES ~ months + gtype + (1 + months | participant) + (1 | allotment),
           data = gardenpath, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa"))

m3 <- lmer(RSES ~ months * gtype + (1 + months | participant) + (1 | allotment),
           data = gardenpath, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa"))

m4 <- lmer(RSES ~ months * gtype + age + prevexp + (1 + months | participant) + (1 | allotment),
           data = gardenpath, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa"))

# Model comparison
anova(m1, m2, m3, m4)


# Final Model ─────────────────────────────────────────────────────────────
m_final <- lmer(RSES ~ months * gtype + age + prevexp +
                  (1 + months | participant) + (1 | allotment),
                data = gardenpath, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa"))

summary(m_final)
confint(m_final, method = "Wald")


# Post-hoc ────────────────────────────────────────────────────────────────
# Compare rates of change (slopes) across gardening types
emtrends(m_final, pairwise ~ gtype, var = "months")

# Pairwise differences at key timepoints
pairs(emmeans(m_final, ~ gtype | months, at = list(months = c(0, 6, 12, 18))))

# Assumptions ─────────────────────────────────────────────────────────────
# Residuals vs fitted
plot(m_final)

# Scale-location (homoscedasticity)
plot(m_final, form = sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))

# Normality of residuals
qqnorm(resid(m_final)); qqline(resid(m_final))

# Normality of random effects
lattice::qqmath(ranef(m_final, condVar = TRUE))

# Model fit check
check_predictions(m_final, iterations = 1000)
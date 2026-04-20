# ── 0. Packages ────────────────────────────────────────────────────────────────
library(lavaan)
library(tidyverse)
library(psych)
library(semPlot)

# ── 1. Data Preparation ────────────────────────────────────────────────────────
# Initial data check
summary(locus)
describe(locus)

# Check distributions before reverse coding
multi.hist(locus[, -1], global = FALSE)  # exclude ID column

# Reverse score negatively worded items:
# loc4: "It is just luck that decides what happens in my life"
# cope3: "When I have a problem, my feelings tend to get in the way of me solving it"
locus <- locus |>
  mutate(
    loc4 = 6 - loc4,
    cope3 = 6 - cope3
  )

# ── 2. CFA: Locus of Control ───────────────────────────────────────────────────
# Initial model
loc_mod <- "
loc =~ loc1 + loc2 + loc3 + loc4 + loc5 + loc6"

loc_mod_est <- cfa(loc_mod, missing = "ML", data = locus)
fitmeasures(loc_mod_est)[c("rmsea", "srmr", "tli", "cfi")]
modindices(loc_mod_est, sort = TRUE, maximum.number = 5)

# Modified model - adding loc3 ~~ loc4 based on MI (61.457)
# Justified: both items relate to future outcomes being controllable
loc_mod2 <- "
loc =~ loc1 + loc2 + loc3 + loc4 + loc5 + loc6
loc3 ~~ loc4"

loc_mod_est2 <- cfa(loc_mod2, missing = "ML", data = locus)
fitmeasures(loc_mod_est2)[c("rmsea", "srmr", "tli", "cfi")]
modindices(loc_mod_est2, sort = TRUE, maximum.number = 5)
parameterestimates(loc_mod_est2, standardized = TRUE)
# Final fit: RMSEA = 0.000, SRMR = 0.010, TLI = 1.005, CFI = 1.000 ✅

# ── 3. CFA: Coping Strategies ──────────────────────────────────────────────────
# Initial model
cope_mod <- "
cope =~ cope1 + cope2 + cope3 + cope4 + cope5 + cope6 + cope7 + cope8"

cope_mod_est <- cfa(cope_mod, missing = "ML", data = locus)
fitmeasures(cope_mod_est)[c("rmsea", "srmr", "tli", "cfi")]
modindices(cope_mod_est, sort = TRUE, maximum.number = 5)

# Modified model - adding cope4 ~~ cope7 based on MI (296.535)
# Justified: both items describe structured step-by-step problem solving
cope_mod2 <- "
cope =~ cope1 + cope2 + cope3 + cope4 + cope5 + cope6 + cope7 + cope8
cope4 ~~ cope7"

cope_mod_est2 <- cfa(cope_mod2, missing = "ML", data = locus)
fitmeasures(cope_mod_est2)[c("rmsea", "srmr", "tli", "cfi")]
modindices(cope_mod_est2, sort = TRUE, maximum.number = 5)
parameterestimates(cope_mod_est2, standardized = TRUE)
# Final fit: RMSEA = 0.035, SRMR = 0.019, TLI = 0.990, CFI = 0.993 ✅

# ── 4. CFA: Anxiety and Depression ────────────────────────────────────────────
# Initial model - two correlated factors as specified in the brief
anx_dep_mod <- "
anx =~ anx1 + anx2 + anx3 + anx4
dep =~ dep1 + dep2 + dep3 + dep4
anx ~~ dep"

anx_dep_est <- cfa(anx_dep_mod, missing = "ML", data = locus)
fitmeasures(anx_dep_est)[c("rmsea", "srmr", "tli", "cfi")]
modindices(anx_dep_est, sort = TRUE, maximum.number = 5)

# Modified model - adding dep1 ~~ dep2 based on MI (16.963)
# Justified: both items describe emotional/somatic expressions of sadness
anx_dep_mod2 <- "
anx =~ anx1 + anx2 + anx3 + anx4
dep =~ dep1 + dep2 + dep3 + dep4
anx ~~ dep
dep1 ~~ dep2"

anx_dep_est2 <- cfa(anx_dep_mod2, missing = "ML", data = locus)
fitmeasures(anx_dep_est2)[c("rmsea", "srmr", "tli", "cfi")]
modindices(anx_dep_est2, sort = TRUE, maximum.number = 5)
parameterestimates(anx_dep_est2, standardized = TRUE)
# Final fit: RMSEA = 0.029, SRMR = 0.019, TLI = 0.991, CFI = 0.994 ✅

# ── 5. Full SEM with Mediation ─────────────────────────────────────────────────
full_mod <- "
  # Measurement models (with modifications from CFAs)
  loc =~ loc1 + loc2 + loc3 + loc4 + loc5 + loc6
  loc3 ~~ loc4

  cope =~ cope1 + cope2 + cope3 + cope4 + cope5 + cope6 + cope7 + cope8
  cope4 ~~ cope7

  anx =~ anx1 + anx2 + anx3 + anx4
  dep =~ dep1 + dep2 + dep3 + dep4
  dep1 ~~ dep2

  # Structural paths (mediation)
  cope ~ a*loc
  anx ~ b1*cope + c1*loc
  dep ~ b2*cope + c2*loc

  # Correlated outcomes
  anx ~~ dep

  # Indirect effects
  ind_anx := a*b1
  ind_dep := a*b2

  # Total effects
  total_anx := ind_anx + c1
  total_dep := ind_dep + c2"

full_mod_sem <- sem(full_mod, data = locus, missing = "ML",
                    std.lv = TRUE, se = "bootstrap", bootstrap = 1000)

summary(full_mod_sem, ci = TRUE)
fitmeasures(full_mod_sem)[c("rmsea", "srmr", "tli", "cfi")]
parameterestimates(full_mod_sem, boot.ci.type = "perc", ci = TRUE)

# ── 6. Path Diagram ────────────────────────────────────────────────────────────
semPaths(full_mod_sem,
         what = "std",
         layout = "tree2",
         rotation = 2,
         edge.label.cex = 0.7,
         node.label.cex = 0.8,
         label.cex = 0.8,
         sizeMan = 5,
         sizeLat = 10,
         color = list(
           lat = "lightblue",
           man = "white"
         ),
         edge.color = "black",
         whatLabels = "std",
         residuals = FALSE,
         intercepts = FALSE,
         nCharNodes = 0,
         exoCov = FALSE,
         mar = c(6, 6, 6, 6),
         label.prop = 0.5,
         edge.label.position = 0.5,
         curvePivot = TRUE,
         pastel = TRUE)

title("SEM: Locus of Control → Coping → Internalising Problems")
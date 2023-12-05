rm(list=ls())

library(babynames)
library(dplyr)
library(tidyr)
library(ggplot2)

# I'd like to falsify a neutral model by testing whether it generates the same
# rate of decay in name survivorship as we see in the population.

# Let's do it for only female names (not entirely sure why this matters).
b <- babynames[babynames$sex == "F", ]

# Set parameters for the simulation
# For now, these are low to test.  In reality, we'll probably want to start in
# like 1920, do 100 reps and 80 years or so.
base.year <- 1880
n.repss <- 5
n.years <- 5

# Filter for just names in base year.
b.base <- b[b$year == base.year, ]

# Define a function to calculate survivorship over time vs base year.
CalcSurv <- function(x) {sum(b.base$name %in% x) / nrow(b.base)}

# Calculate survivorship over time for real data, fit a linear model of the
#relationship, plot and extract the slope.
surv.r <- b[b$year > base.year, c(1, 3)] %>%
  group_by(year) %>%
  summarize(survivorship = CalcSurv(name))
surv.m <- lm(data = surv.r, survivorship ~ year)
plot(surv.r$year, surv.r$survivorship, xlab = "Year", ylab = "Survivorship")
abline(surv.m)
slope.r <- surv.m$coefficients[2]

# Capture population by year, to use in sampling.
pop <- summarize(group_by(b, year), pop_size = sum(n))

# To sample in our neutral model, we need to start with vector containing all of
# the names found in the base year, each appearing the number of times it
# appeared in the population that year
names.s <- c()
for (i in 1:nrow(b.base)) {
  temp <- rep(b.base[[i, 3]], each = b.base[[i, 4]])
  names.s <- append(names.s, temp)
}

# Create 3 objects to use in simulation. 
# 1) A DF to contain the results for each year. This will contain the year, 
# samples, population per the sampling (to confirm behavior), and survivorship.
# 2) An empty list to store each simulation run.
# 3) An empty DF to store the slope of each model
seed <- data.frame(
  year = base.year,
  data = I(list(names.s)),
  pop = length(names.s),
  surv = 1)
complete.runs = list()
slopes <- data.frame(run = integer(), slope = integer())

# Starting one after the base year, take a number of samples from the previous
# year equal to the population in that year.  Complete this exercise for the
# preset number of years and reps, saving each rep to the complete.runs object.
for (i in 1:n.repss) {
  run.hist <- seed
  for (j in 1:n.years) {
    smp <- sample(run.hist[[j, 2]], pop[[j + 1, 2]], replace = TRUE)
    run.hist.row <- data.frame(
      year = j + base.year,
      data = I(list(smp)),
      pop = length(smp),
      surv = CalcSurv(smp))
    run.hist <- rbind(run.hist, run.hist.row)
    if(j == n.years) {complete.runs[[i]] <- run.hist[-1, ]}
  }
}

# For each simulation, fit a linear model, extract the slope, and save it to
# the slope vector.
for (i in 1:n.repss) {
  temp <- data.frame(
    run = i,
    slope = lm(
      complete.runs[[i]]$surv ~ complete.runs[[i]]$year,
      data = complete.runs[[i]])$coefficients[[2]])
  slopes <- rbind(slopes, temp)
}

# Plot a histogram of the slopes for all runs of the simulation vs the slope
# of the actual population.
ggplot(slopes, aes(x = slope)) +
  geom_histogram() +
  geom_point(aes(x = slope.r, y = 0), colour = "blue", size = 6, shape = 8)

# We can see that clearly the rate of decline of name survivorship in the
# neutral model is much faster than in the actual data. In fact, we can  reject
# the neutral model as being outside the 95% confidence intervals. Why might
# that be? Well, perhaps, in reality, parents don't actually sample names
# randomly but instead have a preference for rare names and against common
# names. Some names have more fitness than others!

# To test this, instead of sampling randomly, we can weight our sampling so that
# rare names get a boost (are fitter) and common names get downgraded (are less)
# fit.  So, let's write a function to derive weights for a vector of names.
# This is a bit arbitrary, but let's say that the top 10% of names by commonness
# get a 10% downgrade and the bottom 10% get a 10% boost.

# Takes in a vector of names and returns a vector of weights.
CreateWeightsVector <- function(x) {
  weights.by.name <- data.frame(names = x) %>%
    group_by(names) %>%
    summarize(count = n()) %>%
    arrange(desc(count)) %>%
    mutate(
      weight = ifelse(
        row_number() <= round(length(unique(x)) / 10),
        1 / length(unique(x)) * .9,
        ifelse(
          row_number() >= length(unique(x)) - round(length(unique(x)) / 10) + 1,
          1 / length(unique(x)) * 1.1,
          1 / length(unique(x)))
      ))
  weight.vector <- left_join(
    data.frame(names = x),
    weights.by.name[, c(1, 3)],
    by = "names")[, 2]
  return(weight.vector)
  }

# Again, we need to start with 3 global objects for use in the simulation.
seed.w <- data.frame(
  year = base.year,
  data = I(list(names.s)),
  weights = I(list(CreateWeightsVector(names.s))),
  pop = length(names.s),
  surv = 1)
complete.runs.w = list()
slopes.w <- data.frame(run = integer(), slope = integer())

# Then, we rerun our simulations exactly the same, except now with weighting.
for (i in 1:n.repss) {
  run.hist <- seed.w
  for (j in 1:n.years) {
    smp <- sample(
      run.hist[[j, 2]],
      pop[[j + 1, 2]],
      prob = run.hist[[j, 3]],
      replace = TRUE)
    run.hist.row <- data.frame(
      year = j + base.year,
      data = I(list(smp)),
      weights = I(list(CreateWeightsVector(smp))),
      pop = length(smp),
      surv = CalcSurv(smp))
    run.hist <- rbind(run.hist, run.hist.row)
    if(j == n.years) {complete.runs.w[[i]] <- run.hist[-1, ]}
  }
}

# Everything else is exactly the same.  We fit linear models and plot a
# histogram of their slopes.
for (i in 1:n.repss) {
  temp <- data.frame(
    run = i,
    slope = lm(
      complete.runs.w[[i]]$surv ~ complete.runs.w[[i]]$year,
      data = complete.runs.w[[i]])$coefficients[[2]])
  slopes.w <- rbind(slopes.w, temp)
}
ggplot(slopes.w, aes(x = slope)) +
  geom_histogram() +
  geom_point(aes(x = slope.r, y = 0), colour = "blue", size = 6, shape = 8)

# Now population result lies within the simulation 95% confidence intervals!
# TODO: Rewrite this text after a run with full parameters set.
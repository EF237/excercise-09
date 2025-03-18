library(tidyverse)

f <- "https://raw.githubusercontent.com/difiore/ada-datasets/main/Street_et_al_2017.csv"
  
d <- read_csv(f, col_names = TRUE)

install.packages("skimr")
library(skimr)

# Step 1 
skim(d)

# Step 2 

par(mfrow = c(2,2))
  
 # Group_size
plot(d$Group_size, d$ECV, main = "ECV vs Group Size",
     xlab = "Group Size", ylab = "ECV")

 # Longevity
plot(d$Longevity, d$ECV, main = "ECV vs Longevity",
     xlab = "Longevity", ylab = "ECV")

  # Weaning
plot(d$Weaning, d$ECV, main = "ECV vs Weaning Period",
     xlab = "Weaning Period", ylab = "ECV")

  # Reproductive Lifespan
plot(d$Repro_lifespan, d$ECV, main = "ECV vs Reproductive Lifespan",
     xlab = "Reproductive Lifespan", ylab = "ECV")

# Step 3 - Group Size by Hand

filtered_Group_size <- na.omit(d[, c("ECV", "Group_size")])

  # beta 1 by hand
(beta1 <- sum((filtered_Group_size$ECV - mean(filtered_Group_size$ECV)) * (filtered_Group_size$Group_size - mean(filtered_Group_size$Group_size)))/sum((filtered_Group_size$Group_size - mean(filtered_Group_size$Group_size))^2))
(beta1 <- cov(filtered_Group_size$Group_size, filtered_Group_size$ECV)/var(filtered_Group_size$Group_size))
(beta1 <- cor(filtered_Group_size$Group_size, filtered_Group_size$ECV) * (sd(filtered_Group_size$ECV)/sd(filtered_Group_size$Group_size)))
  # output - 2.463071

  # beta0 by hand
(beta0 <- mean(filtered_Group_size$ECV) - beta1 * mean(filtered_Group_size$Group_size))
 #ouput - 30.35652

EG <- na.omit(d[, c("ECV", "Group_size")])

# beta 1 by hand (3 options)
(beta1 <- sum((EG$ECV - mean(EG$ECV)) * (EG$Group_size - mean(EG$Group_size)))/sum((EG$Group_size - mean(EG$Group_size))^2))
#or
(beta1 <- cov(EG$Group_size, EG$ECV)/var(EG$Group_size))
#or
(beta1 <- cor(EG$Group_size, EG$ECV) * (sd(EG$ECV)/sd(EG$Group_size)))
# output - 2.463071

# beta0 by hand
(beta0 <- mean(EG$ECV) - beta1 * mean(EG$Group_size))
#ouput - 30.35652





# Step 4 - Group Size by lm()

(m <- lm(ECV ~ Group_size, data = EG))

summary(m)

par(mfrow = c(2,2))

plot(m)

# Step 5

  # Catarrhini
c = d |>
  filter (Taxonomic_group == "Catarrhini")

(c.m <- lm(ECV ~ Group_size, data = c))

summary(c.m)

par(mfrow = c(2,2))

plot(c.m)

  # Platyrrhini
p = d |>
  filter (Taxonomic_group == "Platyrrhini")

(p.m <- lm(ECV ~ Group_size, data = p))

summary(p.m)

par(mfrow = c(2,2))

plot(p.m)

  # Strepsirhini
s = d |>
  filter (Taxonomic_group == "Strepsirhini")

(s.m <- lm(ECV ~ Group_size, data = s))

summary(s.m)

par(mfrow = c(2,2))

plot(s.m)

##################### Step 6 - linear regression of ECV on social group size

# Remove NAs
filtered_Group_size <- d[complete.cases(d$ECV, d$Group_size), ]

# Assign variables
X <- filtered_Group_size$Group_size
Y <- filtered_Group_size$ECV
n <- length(X)

# Compute means
X_bar <- mean(X)
Y_bar <- mean(Y)

# Compute slope (beta1)
beta1 <- sum((X - X_bar) * (Y - Y_bar)) / sum((X - X_bar)^2)

# Compute intercept (beta0)
beta0 <- Y_bar - beta1 * X_bar

# Compute predicted values
Y_hat <- beta0 + beta1 * X

# Compute residuals
residuals <- Y - Y_hat

# Compute residual sum of squares (RSS)
rss <- sum(residuals^2)

# Compute standard error of beta1
stderror <- sqrt(rss / (n - 2)) / sqrt(sum((X - X_bar)^2))

# Print correct standard error
stderror

m <- lm(ECV ~ Group_size, data = filtered_Group_size)
summary(m)$coefficients["Group_size", "Std. Error"]





















  # 95% Cl

t_critical <- qt(0.975, (n-2)) 
(lower <- beta1 - t_critical * stderror)
(upper <- beta1 + t_critical * stderror)

  #Verify
confint(m)

  # P value associated with coefficent
t_stat <- beta1 / stderror
(p <- 2 * (1 - pt(abs(t_stat), (n-2))))

  #Verify
summary(m)$coefficients["Group_size", "Pr(>|t|)"]


# Step 7 - permutation



n <- 1000

perm <- vector()

for (i in 1:n) {
  perm_ECV <- sample(EG$ECV)
  perm_model <- lm(perm_ECV ~ EG$Group_size)
  perm[i] <- coef(perm_model)["EG$Group_size"]
}

(perm_p <- mean(abs(perm) >= abs(beta1)))









n <- 10000

m_orig <- lm(ECV ~ Group_size, data = filtered_Group_size)
beta1_obs <- coef(m_orig)["Group_size"]  

perm_slopes <- numeric(n_permutations)

for (i in 1:n_permutations) {
  permuted_ECV <- sample(filtered_Group_size$ECV)  # Shuffle ECV
  perm_model <- lm(permuted_ECV ~ filtered_Group_size$Group_size)  # Fit model
  perm_slopes[i] <- coef(perm_model)["filtered_Group_size$Group_size"]  # Store slope
}







perm <- vector(length = n)  

for (i in 1:n) {
  perm_ECV <- sample(EG$ECV)  
  
  perm_model <- lm(perm_ECV ~ EG$Group_size, data = EG) 
  
  perm[i] <- coef(perm_model)[2]
}

# Compute p-value (two-tailed)
perm_p <- mean(abs(perm) >= abs(beta1))

print(perm_p)













# Compute p-value (two-tailed test)
p_value <- mean(abs(perm_slopes) >= abs(beta1_obs))

# Print results
cat("Observed Beta1:", beta1, "\n")
cat("Permutation-based p-value:", p, "\n")

cat("Permutation-based p-value:", perm_p, "\n")

p
# Plot the null distribution
hist(perm_slopes, breaks = 30, col = "lightblue", main = "Permutation-Based Null Distribution of Beta1",
     xlab = "Slope Coefficient (Beta1)")
abline(v = beta1_obs, col = "red", lwd = 2, lty = 2)  








# Step 8 - bootstrapping

n_boot <- 1000  

boot_slopes <- numeric(n_boot)  

n <- nrow(EG)  

for (i in 1:n_boot) {
  boot <- EG[sample(1:n, size = n, replace = TRUE), ] 
  boot_model <- lm(ECV ~ Group_size, data = boot)  
  boot_slopes[i] <- coef(boot_model)[2]  
}

# Quantile Method

ci_quantile <- quantile(boot_slopes, probs = c(0.025, 0.975))

# Theory Method
boot_se <- sd(boot_slopes) 

(ci_theory <- c(beta1 - 1.96 * boot_se, beta1 + 1.96 * boot_se))




zero_in_ci_quantile <- ci_quantile[1] <= 0 & ci_quantile[2] >= 0
zero_in_ci_theory <- ci_theory[1] <= 0 & ci_theory[2] >= 0
















se_beta1 <- summary(m)$coefficients["Group_size", "Std. Error"]

 # p-value for beta1
p_value_beta1 <- summary(m)$coefficients["Group_size", "Pr(>|t|)"]

 # 95% confidence interval
ci_beta1 <- confint(m)["Group_size", ]  

 # Print results
se_beta1
p_value_beta1
ci_beta1

n <- nrow(d)  # Sample size
residuals_squared <- sum(residuals(m)^2)  
x_var <- sum((d$Group_size - mean(d$Group_size))^2)  

se_beta1_manual <- sqrt(residuals_squared / (n - 2) / x_var)
se_beta1_manual


t_critical <- qt(0.975, df = n - 2)
beta1 <- coef(m)["Group_size"]

ci_lower <- beta1 - t_critical * se_beta1_manual
ci_upper <- beta1 + t_critical * se_beta1_manual

c(ci_lower, ci_upper)



t_value <- beta1 / se_beta1_manual
p_value_manual <- 2 * (1 - pt(abs(t_value), df = n - 2))
p_value_manual





# Step 7 - permutation



# Step 8 - bootstrapping



library(tidyr)

summary(m)
broom::tidy(m)
confint(m)
broom::glance(m)










abline(ceof(m1), col = "blue")



install.packages("lmodel2")
library(lmodel2)

m2 = lmodel2(height~weight, data = d,
             range.y = "relative", range.x = "relative",
             nperm = 1000)

m2

plot(x=d$weight, y=d$height)

betas <- broom::tidy(m2) |>
  filter(method == "OLS") |>
  pull(estimate)
abline(betas, col = "blue")

betas <- broom::tidy(m2)|>
  filter(method == "RMA") |>
  pull(estimate)
abline(betas, col = "red")

betas <- broom::tidy(m2)|>
  filter(method == "MA") |>
  pull(estimate)
abline(betas, col = "green")

or

plot(m2, "OLS")
plot(m2, "MA")
plot(m2, "RMA")

f <- "https://raw.githubusercontent.com/difiore/ada-datasets/main/Street_et_al_2017.csv"
d <- read_csv(f, col_names = TRUE)

m = lm(formula = ECV ~ Group_size, data = d)

m = lm(ECV ~ Group_size, d)

par(mfrow = c(2,2))

plot(d$Group_size, d$ECV) 
plot(d$Longevity, d$ECV) 
plot(d$Weaning, d$ECV) 
plot(d$Repro_lifespan, d$ECV) 

m1 = lm(formula = ECV~ Group_size, data = d)
m2 = lm(formula = ECV~ Longevity, data = d)
m3 = lm(formula = ECV~ Weaning, data = d)
m4 = lm(formula = ECV~ Repro_lifespan, data = d)

broom::tidy(m1)
confint(m1)
broom::tidy(m2)
confint(m2)
broom::tidy(m3)
confint(m3)
broom::tidy(m4)
confint(m4)


# Step 5

c = d |>
  filter (Taxonomic_group == "Catarrhini")

p = d |>
  filter (Taxonomic_group == "Platyrrhini")

s = d |>
  filter (Taxonomic_group == "Strepsirhini")

# Step 6



# Step 7
install.packages("infer")


library(broom)
library(infer)

alpha <- 0.05
confidence_level <- 1 - alpha
p_lower <- alpha/2
p_upper <- 1 - (alpha/2)
degrees_of_freedom <- nrow(EG) - 2
critical_value <- qt(p_upper, df = degrees_of_freedom)

original.slope <- lm(data = EG, ECV ~ Group_size) |>
  tidy(conf.int = TRUE, conf.level = confidence_level) |>
  mutate(lower = estimate - std.error * critical_value, upper = estimate + std.error *
           critical_value) |>
  filter(term == "Group_size")

permuted.slope <- EG |>
  specify(ECV ~ Group_size) |>
  hypothesize(null = "independence") |>
  generate(reps = 1000, type = "permute") |>
  calculate(stat = "slope")

permuted.slope.summary <- permuted.slope |>
  summarize(estimate = mean(stat), std.error = sd(stat), lower = estimate - std.error *
              critical_value, upper = estimate + std.error * critical_value, perm.lower = quantile(stat,
                                                                                                   p_lower), perm.upper = quantile(stat, p_upper))

get_ci(permuted.slope, level = 1 - alpha, type = "percentile")

get_ci(permuted.slope, level = 1 - alpha, type = "se", point_estimate = pull(permuted.slope.summary,
                                                                             estimate))

p.value <- permuted.slope |>
  mutate(abs_stat = abs(stat)) |>
  summarize(estimate = mean(abs_stat >= abs(pull(original.slope, estimate))))

p.value

(p.value <- permuted.slope |>
    get_p_value(obs_stat = original.slope$estimate, direction = "two_sided"))







boot.slope <- EG |>
  specify(ECV ~ Group_size) |>
  generate(reps = 1000, type = "bootstrap") |>
  calculate(stat = "slope")

boot.slope.summary <- boot.slope |>
  summarize(estimate = mean(stat), std.error = sd(stat), lower = estimate - std.error *
              critical_value, upper = estimate + std.error * critical_value, boot.lower = quantile(stat,
                                                                                                   p_lower), boot.upper = quantile(stat, p_upper))

boot.slope.summary


for (i in 1:n_boot) {
  boot[[i]] <- mean(samples(EG, n, replace = TRUE))
  # boot <- EG[sample(1:n, size = n, replace = TRUE), ]
  # boot_model <- lm(ECV ~ Group_size, data = boot)
  # boot_slopes[i] <- coef(boot_model)[2]
}

# Quantile Method
(CI.percentile <- get_ci(boot.slope, level = 1 - alpha, type = "percentile"))

# Theory Method
(CI.theory <- get_ci(boot.slope, level = 1 - alpha, type = "se", point_estimate = pull(boot.slope.summary,
                                                                                       estimate)))






library(dplyr)

set.seed(210) 

inv.logit<-function(l) exp(l)/(1+exp(l))

# population of 500,000 subjects
n = 500000 

beta = c(log(1.221403), log(1.1)) # P(Y|X,Z)
gamma = c(log(1.349859), log(2), log(2), 0.3) #P(D|X,Y,Z)

# generate the genetic variant X with MAF = 0.3
x = rbinom(n, size=2, prob=0.3)

# generate the standardized continuous covariate Z correlated with X
z = rnorm(n, mean=0.5*x - 0.3, sd=1)

# generate the binary secondary trait Y
y = rbinom(n, size=1, prob=inv.logit(-1 + beta[1]*x + beta[2]*z))

# Generate disease D. Adjust alpha to get prevalence of interest (0.1). 
alpha= -2.96
d = rbinom(n, size=1, prob=inv.logit(alpha + gamma[1]*x + gamma[2]*y + gamma[3]*z + gamma[4]*x*y))

df <- data.frame(x, y, z, d)

model <- glm(y ~ x + z, data=df %>% sample_n(10000), family=binomial)
exp(cbind(OR = coef(model), confint(model)))

# Sample case control
# case-control sampling 
result = NULL
n.sample = 10000
for(i in 1:100) {
  df_cases <- df %>% filter(d == TRUE) %>% sample_n(n.sample)
  df_cases$weight <- sum(df$d)/n.sample
  df_controls <- df %>% filter(d == FALSE) %>% sample_n(nrow(df_cases))
  df_controls$weight <- (n - sum(df$d))/n.sample
  
  df.c <- rbind(df_cases, df_controls)

  unweighted <- glm(y ~ x + z, data=df.c, family=binomial)
  weighted <- glm(y ~ x + z, data=df.c, family=quasibinomial, weights=weight)
  result = rbind(result, cbind(exp(coef(unweighted))[3], exp(coef(weighted))[3]))
}
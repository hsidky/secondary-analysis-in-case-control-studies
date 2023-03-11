library(dplyr)

set.seed(210) 

inv.logit<-function(l) exp(l)/(1+exp(l))

# population of 500,000 subjects
n = 500000 

beta = c(log(1.221403), log(1.1)) # P(Y|X,Z)
gamma = c(log(1.349859), log(2), log(2), 0.3) #P(D|X,Y,Z)

generate.data <- function(n, alpha) {
  # generate the genetic variant X with MAF = 0.3
  x = rbinom(n, size=2, prob=0.3)
  
  # generate the standardized continuous covariate Z correlated with X
  z = rnorm(n, mean=0.5*x - 0.3, sd=1)
  
  # generate the binary secondary trait Y
  y = rbinom(n, size=1, prob=inv.logit(-1 + beta[1]*x + beta[2]*z))
  
  # Generate disease D. Adjust alpha to get prevalence of interest (0.1). 
  d = rbinom(n, size=1, prob=inv.logit(alpha + gamma[1]*x + gamma[2]*y + gamma[3]*z + gamma[4]*x*y))
  
  return(data.frame(x, y, z, d))
}

df1 <- generate.data(n, -2.96) # Prevalence 0.1.
df2 <- generate.data(n, -2.00) # Prevalence 0.2.
df3 <- generate.data(n, -1.50) # Prevalence 0.3

df <- rbind(df1, df2, df3)

model <- glm(y ~ x + z, data=df %>% sample_n(100000), family=binomial)
cbind(OR = coef(model), confint(model))

model <- glm(d ~ x + y + z + x:y, data=df %>% sample_n(100000), family=binomial)
cbind(OR = coef(model), confint(model))

# Sample case control
# case-control sampling 
result = NULL
n.sample = 5000
for(i in 1:100) {
  
  df.c = NULL
  for(df.s in list(df1, df2, df3)) {
    df_cases <- df.s %>% filter(d == TRUE) %>% sample_n(n.sample)
    df_cases$weight <- sum(df.s$d == TRUE)/n.sample*(sum(df1$d) + sum(df2$d) + sum(df3$d))/sum(df.s$d)
    
    df_controls <- df.s %>% filter(d == FALSE) %>% sample_n(n.sample)
    df_controls$weight <- sum(df.s$d == FALSE)/n.sample*(sum(df1$d) + sum(df2$d) + sum(df3$d))/sum(df.s$d)
    
    df.c <- rbind(df.c, df_cases, df_controls)
  }

  unweighted <- glm(y ~ x + z, data=df.c, family=binomial)
  weighted <- glm(y ~ x + z, data=df.c, family=quasibinomial, weights=weight)
  result = rbind(result, cbind(coef(unweighted)[3], coef(weighted)[3]))
}

est.bias <- (colMeans(result) - beta[2])/beta[2]*100
est.bias

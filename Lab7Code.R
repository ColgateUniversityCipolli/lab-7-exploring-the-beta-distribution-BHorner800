#############################################################################
#Lab 7
#Ben Horner
#Math240
#############################################################################
#############################################################################
#Libraries
#############################################################################
library(ggplot2)
library(tidyverse)
library(e1071)
library(cumstats)
library(patchwork)
#############################################################################
#Task 1: Describe The Population Distribution
#############################################################################
summarize.beta = function(alpha, beta){
  ab.values = paste(c(alpha,beta), collapse = ",")
  mean = (alpha)/(alpha + beta)
  var = (alpha*beta)/(((alpha + beta)^2)*(alpha + beta + 1))
  skew = (2*(beta - alpha)*sqrt(alpha + beta + 1))/((alpha + beta + 2)*sqrt(alpha*beta))
  kurt = 6*((alpha-beta)^2 * (alpha + beta + 1) - alpha*beta*(alpha + beta + 2))/(
    alpha*beta*(alpha + beta + 2)*(alpha + beta + 3))
  
  beta.summaries = data.frame(beta_distribution = ab.values, 
                              mean = mean,
                              varience = var,
                              skewness = skew,
                              excess_kurtosis = kurt)
  return(beta.summaries)
}

plot.beta = function(alpha, beta){
  fig.dat <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   # generate a grid of points
    mutate(beta.pdf = dbeta(x, alpha, beta),                      # compute the beta PDF
           norm.pdf = dnorm(x,                                    # Gaussian distribution with
                            mean = alpha/(alpha+beta),            # same mean and variance
                            sd = sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))))
  
  ggplot(data= fig.dat)+                                              # specify data
    geom_line(aes(x=x, y=beta.pdf, color="Beta")) +                 # plot beta dist
    geom_line(aes(x=x, y=norm.pdf, color="Gaussian(0.2857, 0.0255)")) +  # plot guassian dist
    geom_hline(yintercept=0)+                                            # plot x axis
    theme_bw()+                                                          # change theme
    xlab("x")+                                                           # label x axis
    ylab("Density")+                                                     # label y axis
    scale_color_manual("", values = c("black", "grey"))+                 # change colors
    theme(legend.position = "bottom")+
    ggtitle(paste("Beta(alpha=", alpha, ", beta=", beta, ")"))
}

#List of alpha beta pairs to loop through
alpha.list = c(2, 5, 5, 0.5)
beta.list = c(5, 5, 2, 0.5)
beta.summaries = data.frame()
beta.plots = list()

#Loop must be of same length as list of alpha and beta
for (i in 1:4){ 
  alpha = alpha.list[i]
  beta = beta.list[i]
  new.row = summarize.beta(alpha, beta)
  beta.summaries = rbind(beta.summaries, new.row) #saves the summaries of each beta
  new.plot = plot.beta(alpha, beta)
  beta.plots[[i]] = new.plot #saves each beta plot as an element in this list
}
p1 = plot.beta(2, 5)
p2 = plot.beta(5, 5)
p3 = plot.beta(5, 2)
p4 = plot.beta(0.5, 0.5)
(p1+p2)/(p3+p4)

#############################################################################
#Task 2: Compute the moments
#############################################################################
#creating the integrans for centered and uncentered functions
uncent.integrand = function(x, alpha, beta, k){
  (x^k)*dbeta(x, shape1 = alpha, shape2 = beta)
}
cent.integrand = function(x, alpha, beta, k){
  ux = integrate(uncent.integrand, lower = 0, upper = 1, alpha = alpha, beta = beta, k = 1)$value
  ((x-ux)^k) * dbeta(x, shape1 = alpha, shape2 = beta)
}

#calculating the centered and uncentered beta moment
beta.moment = function(alpha, beta, k, centered = T){
  #uncentered moment
  if (centered == F){
    uncent.moment = integrate(uncent.integrand, lower = 0, upper = 1, alpha = alpha, beta = beta, k = k)
    return(uncent.moment$value)
  }
  #centered moment
  if (centered == T){
    cent.moment = integrate(cent.integrand, lower = 0, upper = 1, alpha = alpha, beta = beta, k = k)
    return(cent.moment$value)
  }
}

#calculating population-level characteristics using the moment
pop.chars = function(alpha, beta){ 
  mean = beta.moment(alpha, beta, 1, F)
  var = beta.moment(alpha, beta, 2)
  skew = (beta.moment(alpha, beta, 3))/((beta.moment(alpha, beta, 2))^1.5)
  kurt = (beta.moment(alpha, beta, 4)/(beta.moment(alpha, beta, 2)^2))-3
  results = c(mean, var, skew, kurt)
}


#############################################################################
#Task 3: Do Data Summaries Help
#############################################################################
##########################
# Generate a sample
##########################
beta.sample.summary = function(alpha, beta, seed){
  set.seed(seed) # Set seed so we all get the same results.
  sample.size <- 500 # Specify sample details
  beta.sample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter
  
  beta.sample.tibble = tibble(values = beta.sample) %>%
    summarize(mean = mean(beta.sample),
              variance = var(beta.sample),
              skewness = skewness(beta.sample),
              kurtosis = kurtosis(beta.sample)
              
    )
  df <- data.frame(beta_sample = beta.sample) #makes it into df instead of vector. 
  #Plotting
  ggplot(df, aes(x = beta_sample)) + #assigns the data as the data frame and aes = x asis
    geom_histogram(aes(y=after_stat(density)), binwidth = 0.05, fill = "lightblue", color = "black") + #makes histogram
    geom_density(color = "red", size = 1) +  #Density curve
    stat_function(fun = dbeta, args = list(shape1 = alpha, shape2 = beta), color = "blue", size = 1) +  #True prob density
    labs(title = "Histogram of Beta Sample", #titles and axis labels
         x = "Beta Distribution Values", 
         y = "Probability Density")
}

x = (beta.sample.summary(2, 5, 7272)) #testing it


#############################################################################
#Task 4: Is Sample Size Important
#############################################################################
set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha = 2
beta = 5
beta.sample <- rbeta(n = sample.size,  # sample size
                     shape1 = alpha,   # alpha parameter
                     shape2 = beta) 

c.mean = cummean(beta.sample)
c.var  = cumvar(beta.sample)
c.skew = cumskew(beta.sample)
c.kurt = cumkurt(beta.sample)

beta.pop.data = summarize.beta(alpha, beta)
plot.on.og.sample <- plot.on.og.sample +
  geom_line(data = new.data, aes(x=x, y=y), color = i)


#############################################################################
#Task 5: How can we model the variation
#############################################################################
beta.sample.summary = function(alpha, beta, seed, n){
  '
  input -> alpha, beta, seed, n (number of samples
  output -> beta sample summary (mean, var, skew, kurt)
  '
  set.seed(seed) # Set seed so we all get the same results.
  sample.size <- n # Specify sample details
  beta.sample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter
  
  beta.sample.tibble = tibble(values = beta.sample) %>%
    summarize(mean = mean(beta.sample),
              variance = var(beta.sample),
              skewness = skewness(beta.sample),
              excess_kurtosis = kurtosis(beta.sample))-3
  return(beta.sample.tibble)
}

#Looping Through sample to see variance in data
n = 500 #sample size
a = 2 #alpha
b = 5 #beta
seed = 7272
stats.summary = data.frame() #where we are going to be storing the data

for (i in 1:1000){
  seed = 7272 + i
  new.row = beta.sample.summary(a, b, seed, n)
  stats.summary = rbind(stats.summary, new.row) #saves the summaries of each beta
}


plot.mean = ggplot(stats.summary, aes(x = mean)) + #assigns the data as the data frame and aes = x asis
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.005, fill = "lightblue", color = "black") + #makes histogram
  geom_density(color = "red", size = 1) +  #Density curve
  labs(title = "Histogram of mean", #titles and axis labels
       x = "Value of Statistic", 
       y = "Count")

plot.var = ggplot(stats.summary, aes(x = variance)) + #assigns the data as the data frame and aes = x asis
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.001, fill = "lightblue", color = "black") + #makes histogram
  geom_density(color = "red", size = 1) +  #Density curve
  labs(title = "Histogram of Variance", #titles and axis labels
       x = "Value of Statistic", 
       y = "Count")

plot.skew = ggplot(stats.summary, aes(x = skewness)) + #assigns the data as the data frame and aes = x asis
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.05, fill = "lightblue", color = "black") + #makes histogram
  geom_density(color = "red", size = 1) +  #Density curve
  labs(title = "Histogram of Skewness", #titles and axis labels
       x = "Value of Statistic", 
       y = "Count")

plot.kurt = ggplot(stats.summary, aes(x = excess_kurtosis)) + #assigns the data as the data frame and aes = x asis
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.05, fill = "lightblue", color = "black") + #makes histogram
  geom_density(color = "red", size = 1) +  #Density curve
  labs(title = "Histogram of Excess Kurtosis", #titles and axis labels
       x = "Value of Statistic", 
       y = "Count")

stat.plot.summary = (plot.mean + plot.var) / (plot.skew + plot.kurt)
stat.plot.summary

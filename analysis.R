

# Bayesian reanalysis of the article: 

# Effectiveness of therapeutic heparin versus prophylactic
# heparin on death, mechanical ventilation, or intensive care unit
# admission in moderately ill patients with covid-19 admitted to
# hospital: RAPID randomised clinical trial

# http://dx.doi.org/10.1136/bmj.n2400

# Install/load packages
renv::use("tidyverse",
          "ggdist")

library(tidyverse)

# Primary composite outcome
# Transform it to the log scale and calculate the standard error (SE)
# https://training.cochrane.org/handbook/current/chapter-06#section-6-3-2

UCL = 1.1
LCL = 0.43

SE = (log(UCL) - log(LCL))/3.92
mean = log(0.69)

#### Bayesian normal conjugate analysis ----

# Spiegelhalter DJ, Abrams KR, Myles JP. Bayesian Approaches to Clinical Trials
# and Health Care Evaluation. Wiley; 2004.

post.normal.mean <- function(prior.mean, prior.var, data.mean, data.var)
{
  post.mean.numerator <- prior.mean/prior.var + data.mean/data.var
  post.mean.denominator <- 1/prior.var + 1/data.var
  post.mean <-  post.mean.numerator/post.mean.denominator
  post.var <- (1/(1/prior.var + 1/data.var))
  newlist <- dplyr::tibble(post.mean, post.var)
  return(newlist)
}

# Calculate posterior distribution assuming a vague prior (Normal(0, 10^2)),
# log scale
vague = post.normal.mean(prior.mean = 0, prior.var = 100,
                         data.mean = mean, data.var = SE^2)




set.seed(123) # set seed for reproducibility
N = 10e4

# Generate 10,000 samples
dist = dplyr::tibble(posterior = rnorm(N,
                                mean = vague$post.mean,
                                sd = sqrt(vague$post.var)))

# Calculate the median and 95% credible interval
width = dist %>% 
  ggdist::median_qi(exp(posterior)) # exp() to transform into the linear scale

width

# Calculate posterior probabilities

probs = 
  dist %>% 
  summarise(less1 = mean(posterior < log(1)),
            less09 = mean(posterior < log(0.9))) %>% 
  mutate(across(everything(), ~round(., 2)))

probs

# Visualize the posterior distribution

dist %>% 
  ggplot(aes(x = exp(posterior))) +
  ggdist::stat_halfeye(fill = "#A263A9",
                       point_interval = ggdist::median_qi,
                       .width = 0.95) +
  annotate("rect", xmin = 0.9, xmax = 1/0.9, ymin = -Inf, ymax = Inf,
           alpha = .7, fill = "white") +
  annotate("text", x = 1, y = 0.09,
           label = "ROPE", size = 5) +
  geom_vline(xintercept = 0.9, linetype = 2 , color = "gray30") +
  geom_vline(xintercept = 1/0.9, linetype = 2 , color = "gray30") +
  scale_x_continuous(breaks = c(1,seq(0.3, 1.5, 0.2))) +
  coord_cartesian(x = c(0.2, 1.6)) +
  labs(x = "\nOdds Ratio",
       y = "Density",
       title = "Primary Composite Outcome: Posterior Distribution Assuming a Vague Prior",
       subtitle = paste0("\n \nPr(< 1.0)    = ",probs$less1,
                         "\n \nPr(< 0.9)    = ",probs$less09),
       caption = "ROPE: range of practical equivalence") +
  theme(
    plot.title.position = 'plot',
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.background = element_blank(),
    panel.grid.major.x = element_line(color = "gray80", size = 0.3),
    plot.margin = margin(20, 20, 20, 20)
  )
  
### Converting from odds ratio to risk difference ----

#Doi SA, Furuya-Kanamori L, Xu C, Lin L, Chivese T, Thalib L.
# Questionable utility of the relative risk in clinical research: a call for
# change to practice. J Clin Epidemiol. 2020;0(0).
# doi:10.1016/j.jclinepi.2020.08.019

# Let's first convert OR to absolute risk in the therapeutic heparin group

Risk_fun = function(OR){
  # Equation 8 in https://doi.org/10.1016/j.jclinepi.2020.08.019,
  # where Rc is the risk in control arm
  
  z = (Rc*OR)/(Rc*(OR - 1) + 1)
  
  z
}

Rc = 0.219

# Calculate the median and 95% credible interval

dist %>% 
  mutate(Risk = 1000*Risk_fun(exp(posterior)), # Exp to convert to Odds Ratio
         # Subtract by the risk in control arm to calculate the risk difference
         RD = -(Risk - 1000*Rc)) %>% 
  ggdist::median_qi(RD)

# Calculate posterior probabilities
  
dist %>% 
  mutate(Risk = Risk_fun(exp(posterior)),
         RD = -(Risk - Rc)) %>% 
  summarise(NNT50 = mean(RD > 0.02),
            NNT20 = mean(RD > 0.05))


# Exemplar runs

# Both functions contamination_perm_test() and skewness_comparison() must be installed, run respective scripts
source("Supplement1_test_perm_contamination.R")
source("Supplement2_skewness_comparison.R")

# You can generate data with a known proportion of false-negative individuals with this script and try permutation test for contaminated data on them.

count<-1000               # sets sample size
inf.prop<-0.5             # sets proportion of seropositive individuals
fixed.effect<-(-3)        # sets the effect of infection on simulated trait
healthy.average<-101.5    # sets the average trait value in noninfected group
sd<-12                    # sets the standard deviation of within group 

false.negatives<-5        # sets proportion of false negative individuals

# Computes counts in respektive groups using set properties
count.infected<-round(inf.prop*count)
count.healthy<-count-count.infected

# Creates a variables that indicates infection
infected<-c(rep(F,count.healthy),rep(T,count.infected))

# Calculates how many false negative individuals will be in the seronegative group
reassign<-round(count.healthy*(false.negatives/100))
# Creates a vector of actual infection
really.infected<-c(rep(F,count.healthy-reassign),rep(T,count.infected+reassign))

# Generates the population (all healthy individuals)
trait<-rnorm(count,healthy.average,sd)

# Modifies infected individuals
trait[(sum(!really.infected)+1):count]<-trait[(sum(!really.infected)+1):count]+fixed.effect

# Sorts infeceted individuals such that the most changed ones are grouped with healthy individulas according to the identification vector. Those extreme individuals will be in effect marked as seronegative 
trait<-c(trait[really.infected==F],sort(trait[really.infected==T],decreasing=(sign(fixed.effect)==1)))

# Scrambles both vectors
scramble<-sample(1:count)

infected<-infected[scramble]
trait<-trait[scramble]

# Trial data are now ready



# Executes permutation test for contaminated data with default argument values
contamination_perm_test(trait,infected)

# Executes permutation test for contaminated data with the skewness analysis
contamination_perm_test(trait,infected,skewness.analysis=T)

# Executes permutation test for contaminated data, the hypothesized difference is changed manually. This comes useful when you have a reason to suspect a paradoxical switch in group means.
contamination_perm_test(trait,infected,higher.healthy=F)

# Executes permutation test on only 100 permutation runs per test - it is quicker, but less accurate
contamination_perm_test(trait,infected,runs=100)

# Executes two tailed permutation test for contaminated data
contamination_perm_test(trait,infected,two.tailed=T)

# Executes skewness comparison to estimate proportion and distribution tail of possible contamination prior to the permutation test
skewness_comparison(trait,infected)

# Executes permutation test for contaminated data, levels of contamination are set manually
contamination_perm_test(trait,infected,percentages=seq(1,15,1))






# This script contains contamination_perm_test() function
# contamination_perm_test(trait, identification, percentages=c(0,5,10), higher.healthy=(mean(trait[identification==F])>mean(trait[identification==T])), runs=10000, skewness.analysis=F)

# Arguments are described below.
# It is necessary to define function calculating skewness index first:

# Function that returns Fisher-Pearson coefficient of skewness
# Input is a vector of real numbers.

FPskewness<-function(x){
  return((sum((x-mean(x))^3)/length(x))/((sqrt(sum((x-mean(x))^2)/length(x)))^3))
}

# contamination_perm_test
# Function that delegates the parameters to either one-tailed or two-tailed version of the test described below

contamination_perm_test<-function(trait,identification,percentages=c(0,5,10),higher.healthy=(mean(trait[identification==F])>mean(trait[identification==T])),runs=10000,two.tailed=F,skewness.analysis=F){
  if(two.tailed==F){
    contamination_perm_test_one(trait=trait,identification=identification,percentages=percentages,higher.healthy=higher.healthy,runs=runs,skewness.analysis=skewness.analysis)
  }else{
    contamination_perm_test_two(trait=trait,identification=identification,percentages=percentages,higher.healthy=higher.healthy,runs=runs,skewness.analysis=skewness.analysis)
  }
}

# The Function works with the following arguments:

# trait - Vector of trait values
# identification - Logical vector of assumed presence (T) or absence (F) of infection
# percentages - Vector of percentages of false negatives amongst negative subjects (contamination levels) for which the permutation test for contaminated data will be run.
# two.tailed - Specifies the version of the test, two.tailed=F is the default.
# higher healthy - Logical. This parameter indicates whether we assume healthy individuals to show higher (T) or lower (F) trait values. When not specified, the script attempts to assume this relationship from the original group means.
# runs - Number of resamplings used in the permutation test

# In this scenario the seropositive group mean is subtracted from the seronegative group mean and the one-tailed permutation test is conducted accordingly.


# One-tailed version of the test

contamination_perm_test_one<-function(trait,identification,percentages=c(0,5,10),higher.healthy=(mean(trait[identification==F])>mean(trait[identification==T])),runs=10000,skewness.analysis=F){
  
  if(length(trait)!=length(identification)){
    stop("The vectors of trait values and infection indication are of different lengths.")
  }
  
  higher<-(mean(trait[identification==F])>mean(trait[identification==T]))
  set.higher<-higher.healthy
  
  orig.means<-tapply(trait,identification,mean)
  orig.means<-data.frame(orig.means)
  
  names(orig.means)<-"Original mean values"
  rownames(orig.means)[which(rownames(orig.means)=="FALSE")]<-"Identified as healthy"
  rownames(orig.means)[which(rownames(orig.means)=="TRUE")]<-"Identified as infected"
  
  higher.report<-ifelse(higher==T,
                        "In the original sample, individuals identified as healthy showed higher \naverage trait value.",
                        "In the original sample, individuals identified as infected showed higher \naverage trait value."
  )
  
  concord<-ifelse(higher==set.higher,"Consequently,","Despite that,")
  
  set.report1<-paste(concord,ifelse(set.higher==T,
                                    "healthy individuals were hypothesised to have higher \naverage trait value in a contamination-free sample. \n",
                                    "infected individuals were hypothesised to have higher \naverage trait value in a contamination-free sample. \n"
  ))
  
  set.report2<-paste(ifelse(set.higher==T,
                            "\nFor each contamination level respective proportion of seronegative \nindividuals with the lowest trait values was relabeled as seropositive \nin the original sample as well as in each permutation test run.",
                            "\nFor each contamination level respective proportion of seronegative \nindividuals with the highest trait values was relabeled as seropositive \nin the original sample as well as in each permutation test run."
  ))
  
  trait<-c(trait[identification==F],trait[identification==T])
  infected<-sort(identification)
  
  count.healthy<-sum(!identification)
  count.infected<-sum(identification)
  
  Nperc<-length(percentages)
  
  vector.ident<-list() 
  
  for(i in 1:Nperc){
    reassign<-round(count.healthy*(percentages[i]/100))
    identification<-c(rep(F,count.healthy-reassign),rep(T,count.infected+reassign))
    vector.ident[[i]]<-identification
  }
  
  which.test<-paste("One-tailed permutation test for contaminated data was executed. \nProportion of differences (mean of non-infected - mean of infected)",
                    ifelse(set.higher==T,"\nhigher","\nlower"),
                    "than the observed difference is returned as an equivalent \nof p-value.\n",collapse=" ")
  
# Sort healthy individuals to indicate possible false-negatives
  trait2<-c(sort(trait[infected==F],decreasing=higher.healthy),trait[infected==T])
  
  dist.reals<-1:Nperc
  names(dist.reals)<-paste(as.character(percentages), "%")
  
  contamination<-paste(as.character(percentages), "%")
  names(contamination)<-paste(as.character(percentages), "%")
  
  mean.healthy<-dist.reals
  mean.infected<-dist.reals
  
  sd.healthy<-dist.reals
  sd.infected<-dist.reals
  
  N.healthy<-dist.reals
  N.infected<-dist.reals
  
  mean.dist.perm<-dist.reals
  p.vals.perm<-dist.reals
  
# Compute group means in non-permuted sample
  for(i in 1:Nperc){
    mean.healthy[i]<-mean(trait2[vector.ident[[i]]==F])
    mean.infected[i]<-mean(trait2[vector.ident[[i]]==T])
    sd.healthy[i]<-sd(trait2[vector.ident[[i]]==F])
    sd.infected[i]<-sd(trait2[vector.ident[[i]]==T])
    N.healthy[i]<-sum(vector.ident[[i]]==F)
    N.infected[i]<-sum(vector.ident[[i]]==T)
    dist.reals[i]<-(mean(trait2[vector.ident[[i]]==F])-mean(trait2[vector.ident[[i]]==T]))
  }
  
  cohen<-abs(dist.reals)/((sd.healthy*N.healthy+sd.infected*N.infected)/(N.healthy+N.infected))
  
# Skewness computation
  skewness.healthy<-FPskewness(trait[infected==F])
  skewness.infected<-FPskewness(trait[infected==T])
  
  skew.diff<-abs(skewness.healthy-skewness.infected)
  
  skewness<-c(skewness.healthy,skewness.infected)
  skewness<-data.frame(skewness)
  
  names(skewness)<-"Fisher-Pearson coefficient of skewness"
  rownames(skewness)<-c("Identified as healthy","Identified as infected")
  
  signs<-sign(dist.reals)
  
# Permutation test with skewness add-on
  perm.dist<-array(NA,dim=c(Nperc,runs))
  rand.skew<-NA
  
  for(run in 1:runs){
    trait2<-sample(trait)
    rand.skew[run]<-abs(FPskewness(trait2[infected==F])-FPskewness(trait2[infected==T]))
    
    trait2<-c(sort(trait2[infected==F],decreasing=higher.healthy),trait2[infected==T])
    
    for(i in 1:Nperc){
      perm.dist[i,run]<-mean(trait2[vector.ident[[i]]==F])-mean(trait2[vector.ident[[i]]==T])
    }
  }
  
  skew.p<-sum(rand.skew>skew.diff)/runs
  
  skew.higher<-ifelse(skewness.healthy>skewness.infected,"test-negative","test-positive")
  skew.guess.higher<-ifelse(skewness.healthy>skewness.infected,FALSE,TRUE)
  healthy.positive<-ifelse(skewness.healthy>0,TRUE,FALSE)
  infected.positive<-ifelse(skewness.infected>0,TRUE,FALSE)
  skew.sig<-ifelse(skew.p<0.05,TRUE,FALSE)
  
  skew.message<-paste(
    ifelse(healthy.positive==infected.positive,
           paste(
             "The distribution of trait value was",
             ifelse(healthy.positive,"positively","negatively"),
             "skewed \nin both groups.",
             "The Fisher-Pearson coefficient of skewness \nwas higher in",skew.higher,"group.")
           ,
           paste("The distribution of individuals identified as healthy \nwas skewed",
                 ifelse(healthy.positive,"positively,","negatively,"),
                 "the distribution of individuals \nidentified as infected",
                 ifelse(infected.positive,"positively.","negatively."))
    )
    ,
    paste("\n\nThe difference between the coefficients of skewness was",
          ifelse(skew.sig," \nstatistically significant.\n",", \nhowever, not statistically significant.\n"), 
          "(Two-tailed permutation test of skewness difference \non ",runs," runs was executed.)",sep="")
    ,
    ifelse(skew.sig==FALSE,
           "\n\nThis might question the assumption of data contamination \nsince we would expect a difference in skewness between \nthe groups in contaminated data. \nProceed with caution.",
           paste("\n\nThis supports the assumption of data contamination.",
                 "\nBased on the difference in skewness we would assume \ncontamitation of healthy group by false negative \nsubjects from the",
                 ifelse(skew.higher=="test-positive","lower","upper"),
                 "tail of the distribution \nof infected individuals, which would lead to overall",
                 ifelse(skew.higher=="test-positive","\ndecrease","\nincrease"),
                 "of test-negative group mean."))
    ,
    ifelse(skew.sig==FALSE,"",
           paste(ifelse(set.higher==skew.guess.higher,
                        paste("\n\nThe skewness analysis provides further support to the hypothesis \nof",
                              ifelse(set.higher,"higher","lower"),
                              "mean in non-contaminated group of healthy \nindividuals, which was used in current permuation test \nfor contaminated data.\n"),
                        paste("\n\nThe skewness analysis, however, does not support the hypothesis \nof",
                              ifelse(set.higher,"higher","lower"),
                              "mean in non-contaminated group of healthy \nindividuals, which was used in current permuation test \nfor contaminated data. Proceed with caution.\n")
           )))
  )
  
  skewness<-rbind(skewness[1],"",skew.p)
  
  rownames(skewness)[c(3,4)]<-c("","p-value")
  
  skewness[c(1,2,4),1]<-format(round(as.numeric(skewness[c(1,2,4),1]),3))
  
  mean.dist.perm<-rowMeans(perm.dist)
  
  if(higher.healthy==T){
    for(i in 1:Nperc){
      p.vals.perm[i]<-sum(dist.reals[i]<perm.dist[i,])/runs
    }
  }else{
    
    for(i in 1:Nperc){
      p.vals.perm[i]<-sum(dist.reals[i]>perm.dist[i,])/runs
    }
  }
  
  mean.healthy<-format(round(mean.healthy,2))
  mean.infected<-format(round(mean.infected,2))
  sd.healthy<-format(round(sd.healthy,2))
  sd.infected<-format(round(sd.infected,2))
  cohen<-format(round(cohen,2))
  dist.reals<-format(round(dist.reals,2))
  
  mean.dist.perm<-format(round(mean.dist.perm,2))
  p.vals.perm<-format(round(p.vals.perm,3))
  
  
  res.table<-rbind(contamination,mean.healthy,mean.infected,dist.reals,mean.dist.perm,sd.healthy,sd.infected,N.healthy,N.infected,cohen,p.vals.perm)
  res.table<-as.data.frame(res.table)
  
  colnames(res.table)<-NULL
  rownames(res.table)<-c("contamination","non-infeceted mean","infected mean","mean difference","expected difference","non-infeceted SD","infected SD","non-infeceted N","infected N","Cohen's d","p-value")
  
  final.message<-ifelse(all(signs>0),
                        "\nThe mean difference was positive in all \nhypothesised contamination levels.",
                        ifelse(all(signs<0),
                               "\nThe mean difference was negative in all \nhypothesised contamination levels.",
                               ifelse(signs[1]>0,
                                      "\nThe mean difference started as positive, but turned negative \nwith growing hypothesised contamination level. \nThe results should be interpreted with caution. \nRunning the test that assumess the opposite relationship \nbetween group means (higher.healthy=T) or a two tailed test \nis worth consideration.",
                                      "\nThe mean difference started as negative, but turned positive \nwith growing hypothesised contamination level. \nThe results should be interpreted with caution. \nRunning the test that assumess the opposite relationship \nbetween group means (higher.healthy=F) or a two tailed test \nis worth consideration."
                               )))
  
  cat("\nSample characteristics:\n")
  print(orig.means)
  cat(paste("\n",higher.report,"\n\n",sep=""))
  cat(set.report1)
  
  if(skewness.analysis==T){
    cat("\n\nSkewness report:\n")
    print(skewness)
    cat("\n")
    cat(skew.message)
  }
  
  cat("\n\nPermutation test for contaminated data:\n")
  cat(paste(set.report2,"\n\n",sep=""))
  cat(which.test)
  cat(paste("\n",runs,"sample permutations were performed.\n"))
  print(res.table)
  cat(paste(final.message,"\n\n",sep=""))
  
  results<-list(orig.means,higher.report,set.report1,skewness,skew.message,set.report2,which.test,res.table,final.message)
  
  return(invisible(results))
}



# Two-tailed version of the test:

contamination_perm_test_two<-function(trait,identification,percentages=c(0,5,10),higher.healthy=(mean(trait[identification==F])>mean(trait[identification==T])),runs=10000,skewness.analysis=F){
  
  if(length(trait)!=length(identification)){
    stop("The vectors of trait values and infection indication are of different lengths.")
  }
  
  higher<-(mean(trait[identification==F])>mean(trait[identification==T]))
  set.higher<-higher.healthy
  
  orig.means<-tapply(trait,identification,mean)
  orig.means<-data.frame(orig.means)
  
  names(orig.means)<-"Original mean values"
  rownames(orig.means)[which(rownames(orig.means)=="FALSE")]<-"Identified as healthy"
  rownames(orig.means)[which(rownames(orig.means)=="TRUE")]<-"Identified as infected"
  
  higher.report<-ifelse(higher==T,
                        "In the original sample, individuals identified as healthy showed higher \naverage trait value.",
                        "In the original sample, individuals identified as infected showed higher \naverage trait value."
  )
  
  concord<-ifelse(higher==set.higher,"Consequently,","Despite that,")
  
  set.report1<-paste(concord,ifelse(set.higher==T,
                                    "healthy individuals were hypothesised to have higher \naverage trait value in a contamination-free sample. \n",
                                    "infected individuals were hypothesised to have higher \naverage trait value in a contamination-free sample. \n"
  ))
  
  set.report2<-paste(ifelse(set.higher==T,
                            "\nFor each contamination level respective proportion of seronegative \nindividuals with the lowest trait values was relabeled as seropositive \nin the original sample as well as in each permutation test run.",
                            "\nFor each contamination level respective proportion of seronegative \nindividuals with the highest trait values was relabeled as seropositive \nin the original sample as well as in each permutation test run."
  ))
  
  trait<-c(trait[identification==F],trait[identification==T])
  infected<-sort(identification)
  
  count.healthy<-sum(!identification)
  count.infected<-sum(identification)
  
  Nperc<-length(percentages)
  
  vector.ident<-list() 
  
  for(i in 1:Nperc){
    reassign<-round(count.healthy*(percentages[i]/100))
    identification<-c(rep(F,count.healthy-reassign),rep(T,count.infected+reassign))
    vector.ident[[i]]<-identification
  }
  
  which.test<-paste("Two-tailed permutation test for contaminated data was executed. \nProportion of differences (in absolute value)",
                "higher than the \nobserved difference is returned as an equivalent of p-value.\n",collapse=" ")
  
# Sort healthy individuals to indicate possible false-negatives
  trait2<-c(sort(trait[infected==F],decreasing=higher.healthy),trait[infected==T])
  
  dist.reals<-1:Nperc
  names(dist.reals)<-paste(as.character(percentages), "%")
  
  contamination<-paste(as.character(percentages), "%")
  names(contamination)<-paste(as.character(percentages), "%")
  
  mean.healthy<-dist.reals
  mean.infected<-dist.reals
  
  sd.healthy<-dist.reals
  sd.infected<-dist.reals
  
  N.healthy<-dist.reals
  N.infected<-dist.reals
  
  mean.dist.perm<-dist.reals
  p.vals.perm<-dist.reals
  
# Compute group means in non-permuted sample
  for(i in 1:Nperc){
    mean.healthy[i]<-mean(trait2[vector.ident[[i]]==F])
    mean.infected[i]<-mean(trait2[vector.ident[[i]]==T])
    sd.healthy[i]<-sd(trait2[vector.ident[[i]]==F])
    sd.infected[i]<-sd(trait2[vector.ident[[i]]==T])
    N.healthy[i]<-sum(vector.ident[[i]]==F)
    N.infected[i]<-sum(vector.ident[[i]]==T)
    dist.reals[i]<-abs((mean(trait2[vector.ident[[i]]==F])-mean(trait2[vector.ident[[i]]==T])))
  }
  
  cohen<-abs(dist.reals)/((sd.healthy*N.healthy+sd.infected*N.infected)/(N.healthy+N.infected))
  
# Skewness computation
  skewness.healthy<-FPskewness(trait[infected==F])
  skewness.infected<-FPskewness(trait[infected==T])
  
  skew.diff<-abs(skewness.healthy-skewness.infected)
  
  skewness<-c(skewness.healthy,skewness.infected)
  skewness<-data.frame(skewness)
  
  names(skewness)<-"Fisher-Pearson coefficient of skewness"
  rownames(skewness)<-c("Identified as healthy","Identified as infected")
  
  signs<-sign(dist.reals)
  
# Permutation test with skewness add-on
  perm.dist<-array(NA,dim=c(Nperc,runs))
  rand.skew<-NA
  
  for(run in 1:runs){
    trait2<-sample(trait)
    rand.skew[run]<-abs(FPskewness(trait2[infected==F])-FPskewness(trait2[infected==T]))
    
    trait2<-c(sort(trait2[infected==F],decreasing=higher.healthy),trait2[infected==T])
    
    for(i in 1:Nperc){
      perm.dist[i,run]<-abs(mean(trait2[vector.ident[[i]]==F])-mean(trait2[vector.ident[[i]]==T]))
    }
  }
  
  skew.p<-sum(rand.skew>skew.diff)/runs
  
  skew.higher<-ifelse(skewness.healthy>skewness.infected,"test-negative","test-positive")
  skew.guess.higher<-ifelse(skewness.healthy>skewness.infected,FALSE,TRUE)
  healthy.positive<-ifelse(skewness.healthy>0,TRUE,FALSE)
  infected.positive<-ifelse(skewness.infected>0,TRUE,FALSE)
  skew.sig<-ifelse(skew.p<0.05,TRUE,FALSE)
  
  skew.message<-paste(
    ifelse(healthy.positive==infected.positive,
           paste(
             "The distribution of trait value was",
             ifelse(healthy.positive,"positively","negatively"),
             "skewed \nin both groups.",
             "The Fisher-Pearson coefficient of skewness \nwas higher in",skew.higher,"group.")
           ,
           paste("The distribution of individuals identified as healthy \nwas skewed",
                 ifelse(healthy.positive,"positively,","negatively,"),
                 "the distribution of individuals \nidentified as infected",
                 ifelse(infected.positive,"positively.","negatively."))
    )
    ,
    paste("\n\nThe difference between the coefficients of skewness was",
          ifelse(skew.sig," \nstatistically significant.\n",", \nhowever, not statistically significant.\n"), 
          "(Two-tailed permutation test of skewness difference \non ",runs," runs was executed.)",sep="")
    ,
    ifelse(skew.sig==FALSE,
           "\n\nThis might question the assumption of data contamination \nsince we would expect a difference in skewness between \nthe groups in contaminated data. \nProceed with caution.",
           paste("\n\nThis supports the assumption of data contamination.",
                 "\nBased on the difference in skewness we would assume \ncontamitation of healthy group by false negative \nsubjects from the",
                 ifelse(skew.higher=="test-positive","lower","upper"),
                 "tail of the distribution \nof infected individuals, which would lead to overall",
                 ifelse(skew.higher=="test-positive","\ndecrease","\nincrease"),
                 "of test-negative group mean."))
    ,
    ifelse(skew.sig==FALSE,"",
           paste(ifelse(set.higher==skew.guess.higher,
                        paste("\n\nThe skewness analysis brings further support to the hypothesis \nof",
                              ifelse(set.higher,"higher","lower"),
                              "mean in non-contaminated group of healthy \nindividuals, which was used in current permuation test \nfor contaminated data.\n"),
                        paste("\n\nThe skewness analysis, however, does not support the hypothesis \nof",
                              ifelse(set.higher,"higher","lower"),
                              "mean in non-contaminated group of healthy \nindividuals, which was used in current permuation test \nfor contaminated data. Proceed with caution.\n")
           )))
  )
  
  skewness<-rbind(skewness[1],"",skew.p)
  
  rownames(skewness)[c(3,4)]<-c("","p-value")
  
  skewness[c(1,2,4),1]<-format(round(as.numeric(skewness[c(1,2,4),1]),3))
  
  mean.dist.perm<-rowMeans(perm.dist)
  
  for(i in 1:Nperc){
    p.vals.perm[i]<-sum(dist.reals[i]<perm.dist[i,])/runs
  }
  
  mean.healthy<-format(round(mean.healthy,2))
  mean.infected<-format(round(mean.infected,2))
  sd.healthy<-format(round(sd.healthy,2))
  sd.infected<-format(round(sd.infected,2))
  cohen<-format(round(cohen,2))
  dist.reals<-format(round(dist.reals,2))
  
  mean.dist.perm<-format(round(mean.dist.perm,2))
  p.vals.perm<-format(round(p.vals.perm,3))
  
  res.table<-rbind(contamination,mean.healthy,mean.infected,dist.reals,mean.dist.perm,sd.healthy,sd.infected,N.healthy,N.infected,cohen,p.vals.perm)
  res.table<-as.data.frame(res.table)
  
  colnames(res.table)<-NULL
  rownames(res.table)<-c("contamination","non-infeceted mean","infected mean","mean difference","expected difference","non-infeceted SD","infected SD","non-infeceted N","infected N","Cohen's d","p-value")
  
  final.message<-paste("\nTwo-tailed permutation test for contaminated data was executed. \nIt was assumed that",
                       ifelse(set.higher==T,"healthy","infected"), 
                       "individuals have on average higher \ntrait value if all false negative individuals are relocated correctly.")
  
  cat("\nSample characteristics:\n")
  print(orig.means)
  cat(paste("\n",higher.report,"\n\n",sep=""))
  cat(set.report1)
  
  if(skewness.analysis==T){
    cat("\n\nSkewness report:\n")
    print(skewness)
    cat("\n")
    cat(skew.message)
  }
  
  cat("\n\nPermutation test for contaminated data:\n")
  cat(paste(set.report2,"\n\n",sep=""))
  cat(which.test)
  cat(paste("\n",runs,"sample permutations were performed.\n"))
  print(res.table)
  cat(paste(final.message,"\n\n",sep=""))
  
  results<-list(orig.means,higher.report,set.report1,skewness,skew.message,set.report2,which.test,res.table,final.message)
  
  return(invisible(results))
}




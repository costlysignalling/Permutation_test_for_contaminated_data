# This script contains skewness_comparison() function
# skewness_comparison(trait,identification,percentages=seq(0,50,1),higher.healthy=(mean(trait[identification==F])>mean(trait[identification==T])),runs=10000)

# Arguments are described below.
# It is necessary to define function calculating skewness index first:

# Function that returns Fisher-Pearson coefficient of skewness
# Input is a vector of real numbers.

FPskewness<-function(x){
  return((sum((x-mean(x))^3)/length(x))/((sqrt(sum((x-mean(x))^2)/length(x)))^3))
}

# Skewness comparison
# Function that reports how the relocation of seronegative individuals changes the skewness of the distribution in both seronegative and seropositive groups.
# trait - Numerical vector of trait values
# identification - Logical vector of assumed presence (T) or absence (F) of infection
# percentages - Numerical vector of percentages of false negative amongst negative subjects (contamination levels) for which the permutation test for contaminated data will be run
# higher healthy - Logical. This parameter indicates whether we assume healthy individuals to show higher (T) or lower (F) trait values. When not specified, the script assumes this relationship based on group means with no hypothesized contamination.
# runs - Number of resamplings used in the permutation test

skewness_comparison<-function(trait,identification,percentages=seq(0,50,1),higher.healthy=(mean(trait[identification==F])>mean(trait[identification==T])),runs=10000){
  
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
  
  
  run.report<-paste("Two-tailed permutation test of skewness difference \non ",runs," runs was executed on each contamination level.\n\n",sep="")
  
  set.report2<-paste("\nFor each contamination level respective proportion of seronegative \nindividuals with", 
                     ifelse(set.higher==T,"lowest","highest"),
                     "trait value was relabeled as seropositive \nand the difference between the group skewness was measured. \nReferential skewness differences from permutation runs were based \non random non-contaminated sample with group sizes corresponding \nto respective contamination levels.\n\n"
  )
  
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
  
  # Sort healthy individuals to indicate possible false-negatives
  trait2<-c(sort(trait[infected==F],decreasing=higher.healthy),trait[infected==T])
  
  skewness.healthy<-1:Nperc
  names(skewness.healthy)<-paste(as.character(percentages), "%")
  skewness.infected<-skewness.healthy
  
  contamination<-paste(as.character(percentages), "%")
  names(contamination)<-paste(as.character(percentages), "%")
  
  # Compute group means in non-permuted sample
  for(i in 1:Nperc){
    skewness.healthy[i]<-FPskewness(trait2[vector.ident[[i]]==F])
    skewness.infected[i]<-FPskewness(trait2[vector.ident[[i]]==T])
  }
  
  skew.diff<-abs(skewness.healthy-skewness.infected)
  
  # Permutation test of skewness difference 
  
  skew.diff.perm<-array(NA,dim=c(Nperc,runs))
  
  for(run in 1:runs){
    trait2<-sample(trait)
    
    for(i in 1:Nperc){
      skew.diff.perm[i,run]<-abs(FPskewness(trait2[vector.ident[[i]]==F])-FPskewness(trait2[vector.ident[[i]]==T]))
    }
  }
  
  p.vals.skew<-NA
  
  for(i in 1:Nperc){
    p.vals.skew[i]<-sum(skew.diff[i]<skew.diff.perm[i,])/runs
  }
  
  guess.perc<-percentages[which(skew.diff==min(skew.diff))]
  
  possible<-percentages[p.vals.skew>0.05]
  
  if(length(possible)==0){
    ok.report<-paste("\nThe difference in Fisher-Pearson coefficients of skewness between \ngroups was significant at all investigated contamination levels.\nThis may be caused by extreme proportion of false negative \nindividuals, insufficient number of relocated fractions, \ndifferent shapes of distributions of healthy and infeceted \nindividuals, or, most likely, by wrong setting of higher group \nmean in non-contaminated sample. \nTry running this comparison with parameter higher.healthy =",ifelse(higher.healthy==TRUE,"FALSE","TRUE"),"\n\n")
  }else{
    if(min(possible)==max(possible)){
      ok.report<-paste("\nThe difference in Fisher-Pearson coefficients of skewness between \ngroups was not significant at",
                       min(possible),
                       "% of relocated individuals.\n\n")
    }else{
      ok.report<-paste("\nThe difference in Fisher-Pearson coefficients of skewness between \ngroups was not significant between",
                       min(possible),
                       "% and",
                       max(possible),
                       "% of relocated individuals.\n\n")
    }
  }
  
  best.guess<-paste("The difference between group skewness was smallest \nwhen",
                    guess.perc,
                    "% of seronegative individuals with",
                    ifelse(higher.healthy,"lowest","highest"),
                    "trait value \nwas relocated to seropositive group.\n\n")
  
  
  skewness.healthy<-format(round(skewness.healthy,2))
  skewness.infected<-format(round(skewness.infected,2))
  p.vals.skew<-format(round(p.vals.skew,3))
  
  skewness.res<-rbind(contamination,skewness.healthy,skewness.infected,"",p.vals.skew)
  
  skewness.res<-as.data.frame(skewness.res)
  
  colnames(skewness.res)<-NULL
  rownames(skewness.res)<-c("contamination","skewness healthy","skewness infected","","p-value")
  
  
  cat("\nSample characteristics:\n")
  print(orig.means)
  cat(paste("\n",higher.report,"\n\n",sep=""))
  cat(set.report1)
  cat(set.report2)
  
  cat(run.report)
  
  cat(paste("Skewness comparison:","\n",sep=""))
  
  print(skewness.res)
  
  cat(ok.report)
  cat(best.guess)
  
  results<-list(orig.means,higher.report,set.report1,skewness.res)
  
  return(invisible(results))
}
source("http://bioconductor.org/biocLite.R")
biocLite()
setwd("C:/Users/elifdogandar/Desktop/data")
library(affy)
affy.data = ReadAffy()  #reading data into R
#image of the raw arrays
op = par(mfrow = c(2,4))
for(i in 1:8){
  image(affy.data[ ,i]) 
}      
#density plot for raw data
plotDensity.AffyBatch(affy.data[,1:8], col = c("black","black","red","red","blue","blue","yellow","yellow"), log = TRUE,
                      which="pm",
                      ylab = "density",
                      xlab = NULL, main="density plot of unnormalized data")   
#boxplot of the raw data
boxplot(affy.data,las=2, main="boxplot of unnormalized data",col = c("purple","purple","red","red","blue","blue","yellow","yellow"))  #boxplot of raw data

#mvaplot for some array pairs
mva.pairs(pm(affy.data)[,c(1,2)]) #for the pair brain 1 and brain 2 
mva.pairs(pm(affy.data)[,c(1,3)]) #for the pair brain 1 and fetal brain 1  
mva.pairs(pm(affy.data)[,c(1,4)]) #for the pair brain 1 and fetal brain 2 

library(simpleaffy)
#assessing quality measures for raw data
data.qc = qc(affy.data)
avbg(data.qc)  # the average intensities of the background probes on each array
sfs(data.qc)   # the scale factors: factors used to equalize the mean intensities of the arrays
percent.present(data.qc) #the percent present calls: the percentage of spots that generate a significant signal 
#(significantly higher than background) according to the Affymetrix detection algorithm.
ratios(data.qc)  #the 3'/5' ratios of the quality control probe sets  

pm(affy.data)[1:5,]   #perfect match intensities of the first 5 observations
mm(affy.data)[1:5,]   #mismatch intensities of the first 5 observations
pm(affy.data)[1:5,]-mm(affy.data)[1:5,]    #differences between perfect and mismatch intensities

#Normalization by using rma
data.rma = rma(affy.data)
data.matrix = exprs(data.rma)
colnames(data.matrix) = c("brain1", "brain2", 
                        "fetal.brain1", "fetal.brain2",
                        "fetal.liver1", "fetal.liver2", 
                        "liver1", "liver2")

#densities of the normalized data
d1 <- density(data.matrix[,1])
d2 <- density(data.matrix[,2])
d3 <- density(data.matrix[,3])
d4 <- density(data.matrix[,4])
d5 <- density(data.matrix[,5])
d6 <- density(data.matrix[,6])
d7 <- density(data.matrix[,7])
d8 <- density(data.matrix[,8])
plot(range(d1$x, d2$x,d3$x), range(d1$y, d2$y,d3$y), type = "n", xlab = "x",
     ylab = "Density",main="Density Plot of Normalized Data")
lines(d1, col = "red")
lines(d2, col = "blue")
lines(d3, col = "pink")
lines(d4, col = "black")
lines(d5, col = "yellow")
lines(d6, col = "brown")
lines(d7, col = "orange")
lines(d8, col = "purple")

#boxplot of the normalized data
boxplot(data.matrix,las=2,col = c("purple","purple","red","red","blue","blue","yellow","yellow")) 


install.packages("ggplot2")
library(ggplot2)
#comparison plot of background normalized and raw data
bgcorr = pm(bg.correct(affy.data,method="rma"))
pmexp = pm(affy.data)
sampleNames = vector()
logs = vector()
corrlogs = vector()
for (i in 1:8)
{
  logs = c(logs,log2(pmexp[,i]))
  corrlogs = c(corrlogs,log2(bgcorr[,i]))
}
sampleNames = c("brain.1", "brain.2", "fetal.brain.1", "fetal.brain.2",
                "fetal.liver.1", "fetal.liver.2", "liver.1", "liver.2")
corrData = data.frame(logInt=logs,bgcorr_logInt=corrlogs,sampleName=sampleNames)
dataScatter = ggplot(corrData, aes(logInt,bgcorr_logInt))
dataScatter + geom_point() + geom_abline(intercept=0,slope=1,colour='red') + facet_grid(.~sampleName)

#some maplots for normalized data
mva.pairs(data.matrix[,c(1,2)]) 
mva.pairs(data.matrix[,c(1,3)]) 
mva.pairs(data.matrix[,c(1,4)]) 

#Principal Component Analysis
color=c('green','green','red','red','black','black','blue','blue')
data.PC = prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x[1:8],col=color, ylim=c(-150,100),main="PCA Plot")
text(data.PC$x[1:8],pos=1,labels=c("brain.1", "brain.2", "fetal.brain.1", "fetal.brain.2",
              "fetal.liver.1", "fetal.liver.2", "liver.1", "liver.2"))

#Use the t test to test for a difference in means with all genes.

brain.p.value.all.genes = apply(data.matrix, 1, function(x) { t.test(x[1:2], x[3:4]) $p.value } )
liver.p.value.all.genes = apply(data.matrix, 1, function(x) { t.test(x[5:6], x[7:8]) $p.value } )

sorteddata1 <- brain.p.value.all.genes[order(brain.p.value.all.genes)]  #sorting the raw p-values 
brain.fdr.pvals = p.adjust(sorteddata1, method="fdr", n=length(sorteddata1))  #adjusted p-values according to fdr procedure
brain.bon.pvals = p.adjust(sorteddata1, method="bonferroni",n=length(sorteddata1))  #adjusted p-values according to bonferroni procedure
brain.hochberg.pvals = p.adjust(sorteddata1, method="hochberg",n=length(sorteddata1))  #adjusted p-values according to hochberg procedure
n <- 12626

#########Bonferrroni
activebon <- rep(0, n)  # it will be 1 if the gene difference is significant 0 otherwise
counterbon <- 0         # it will give the number of significantly different genes

## making comparison with Bonferroni
bon_comp <- function(x,y,n){
  for(i in 1:n){
    if (x[i]<y[i]){
      activebon[i]=1
      counterbon=counterbon+1 
    }
  }
  return(activebon)
}
activebon= bon_comp(sorteddata1,brain.bon.pvals ,n)
counterbon=sum(activebon)

#########Hochberg
activehoch <- rep(0, n)     
counterhoch <- 0              


## making comparison with Hochberg
hoch_comp <- function(x,y,n){
  for(i in 1:n){
    if (x[(n+1-i)]<y[(n+1-i)]){
      activehoch[1:(n+1-i)]=1
      counterhoch=n+1-i
      break
    }
  }
  return(activehoch)
}
activehoch= hoch_comp(sorteddata1,brain.hochberg.pvals,n)
counterhoch= sum(activehoch)

##########Benjamini-Hochberg
activebenj <- rep(0, n)      
counterbenj <- 0

## Making comparison with Benjamini-Hochberg
benj_comp <- function(x,y,n){
  for(i in 1:n){
    if (x[(n+1-i)]<y[(n+1-i)]){
      activebenj[1:(n+1-i)]=1
      counterbenj=n+1-i
      break
    }
  }
  return(activebenj)
}
activebenj= benj_comp(sorteddata1,brain.fdr.pvals ,n)
counterbenj= sum(activebenj)

counterbon
counterhoch
counterbenj

#hypothesis testing with treshold value 0.025

#comparison of adult brain and fetal brain tissues
brain.p.value.all.genes = apply(data.matrix, 1, function(x) {
  t.test(x[1:2], x[3:4]) $p.value } )

counter.brain <- 0    #number of significant genes
sorteddata1 <- brain.p.value.all.genes[order(brain.p.value.all.genes)]
active.brain <- rep(0, n)  #indices of significant genes in sorted data
for(i in 1:n){
  if (sorteddata1[i]<0.025){
    active.brain[1:i]=1
  }
}
counter.brain= sum(active.brain)

#comparison of adult liver and fetal liver tissues
liver.p.value.all.genes = apply(data.matrix, 1, function(x) {
  t.test(x[5:6], x[7:8]) $p.value } )
sorteddata2 <- liver.p.value.all.genes[order(liver.p.value.all.genes)]
counter.liver <- 0
active.liver <- rep(0, n)
for(i in 1:n){
  if (sorteddata2[i]<0.025){
    active.liver[1:i]=1
  }
}
counter.liver= sum(active.liver)

#comparison between adult and fetal tissues
adult.p.value.all.genes = apply(data.matrix, 1, function(x) {
  t.test(x[c(1,2,7,8)], x[c(3,4,5,6)]) $p.value } )
sorteddata3 <- adult.p.value.all.genes[order(adult.p.value.all.genes)]
counter.adult <- 0
active.adult <- rep(0, n)
for(i in 1:n){
  if (sorteddata3[i]<0.025){
    active.adult[1:i]=1
  }
}
counter.adult= sum(active.adult)

#comparison between brain and liver tissues
brli.p.value.all.genes = apply(data.matrix, 1, function(x) {
  t.test(x[c(1,2,3,4)], x[c(5,6,7,8)]) $p.value } )
sorteddata4 <- brli.p.value.all.genes[order(brli.p.value.all.genes)]
counter.brli <- 0
active.brli <- rep(0, n)
for(i in 1:n){
  if (sorteddata4[i]<0.025){
    active.brli[1:i]=1
  }
}
counter.brli= sum(active.brli)


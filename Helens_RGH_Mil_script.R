library(lattice)

setwd("/Users/helencockerton/Desktop/QTL/RgxHp_Mildew/")
locM<-read.csv("locM.csv")
V1 = "RGH_mil_12"

rm(QTL)
R_QTL_S_V1<-paste("R_QTL_S_",V1,".csv",sep="")
QTL<-read.csv(R_QTL_S_V1, header=TRUE)
head(QTL)

QTL$Rname<- gsub('-', '.', QTL$Rname)
x = paste("qtl_model= step(lm(",V1,"~")
for (i in 1:length(QTL[, 4]))
{
    n =i+1
    print(n)
    x<- paste(x, QTL[i,"Rname"], '+')
    print(x)
}

x<-paste(x, ", data=locM))")
x<-gsub("+ ,",",",fixed = TRUE, x)
print(x,quote = FALSE)

filename=paste0("qtl_model_",V1,".R")
write(x, file=filename)
# source(filename)
# summary(qtl_model)

locM$Affx.88876663 <- as.character(locM$Affx.88876663)
for (i in 1:length(locM$Affx.88876663))
{
    locM$Affx.88876663[i] <- if (locM$Affx.88876663[i] == "np") locM$Affx.88876663[i] <- "A" else "B"
}

locM$Affx.88876663 <- as.factor(locM$Affx.88876663)

locM$Affx.88864419 <- as.character(locM$Affx.88864419)
for (i in 1:length(locM$Affx.88864419)) {
  locM$Affx.88864419[i] <- if (locM$Affx.88864419[i] == "np") locM$Affx.88864419[i] <- "A" else "B"}
locM$Affx.88864419 <- as.factor(locM$Affx.88864419)

locM$Affx.88834576 <- as.character(locM$Affx.88834576)
for (i in 1:length(locM$Affx.88834576)) {
  locM$Affx.88834576[i] <- if (locM$Affx.88834576[i] == "np") locM$Affx.88834576[i] <- "A" else "B"}
locM$Affx.88834576 <- as.factor(locM$Affx.88834576)

#                 hh          hk            kh          kk
# Affx-88860589 58.51621622  50.28648649    55.67837838 62.25731707
TukeyHSD(aov(RGH_mil_12~Affx.88860589, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = RGH_mil_12 ~ Affx.88860589, data = locM)
#
# $Affx.88860589
# diff         lwr       upr     p adj
# hk-hh -8.229730 -15.1554750 -1.303984 0.0127450
# kh-hh -2.837838  -9.7635831  4.087907 0.7114867
# kk-hh  3.741101  -3.0136120 10.495814 0.4770859
# kh-hk  5.391892  -1.5338534 12.317637 0.1842265
# kk-hk 11.970831   5.2161177 18.725543 0.0000518
# kk-kh  6.578939  -0.1757741 13.333652 0.0593179

locM$Affx.88860589 <- as.character(locM$Affx.88860589)
for (i in 1:length(locM$Affx.88860589)) {
  locM$Affx.88860589[i] <- if (locM$Affx.88860589[i] == "hk") locM$Affx.88860589[i] <- "B" else "A"}
locM$Affx.88860589 <- as.factor(locM$Affx.88860589)

#                 hh          hk            kh          kk
# Affx-88850184 53.51111111  59.465     61.575  53.27264151
TukeyHSD(aov(RGH_mil_12~Affx.88850184, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
#Fit: aov(formula = RGH_mil_12 ~ Affx.88850184, data = locM)

# $Affx.88850184
# diff         lwr        upr     p adj
# hk-hh  5.9538889  -2.1139787 14.0217565 0.2251193
# kh-hh  8.0638889   0.5618059 15.5659718 0.0298242
# kk-hh -0.2384696  -7.4294652  6.9525260 0.9997697
# kh-hk  2.1100000  -5.1601742  9.3801742 0.8747675
# kk-hk -6.1923585 -13.1410707  0.7563538 0.0990995
# kk-kh -8.3023585 -14.5853007 -2.0194163 0.0042572

locM$Affx.88850184 <- as.character(locM$Affx.88850184)
for (i in 1:length(locM$Affx.88850184)) {
  locM$Affx.88850184[i] <- if (locM$Affx.88850184[i] == "hh") locM$Affx.88850184[i] <- "B" else if (locM$Affx.88850184[i] == "kk") locM$Affx.88850184[i] <- "B" else "A"}
locM$Affx.88850184 <- as.factor(locM$Affx.88850184)

#                 hh          hk          kh            kk
# Affx-88864387 62.825  55.74210526 60.2625 52.69534884
TukeyHSD(aov(RGH_mil_12~Affx.88864387, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = RGH_mil_12 ~ Affx.88864387, data = locM)
#
# $Affx.88864387
# diff        lwr        upr     p adj
# hk-hh  -7.082895 -14.506579  0.3407891 0.0673841
# kh-hh  -2.562500 -11.049165  5.9241646 0.8613551
# kk-hh -10.129651 -17.903168 -2.3561341 0.0049803
# kh-hk   4.520395  -2.520260 11.5610491 0.3441344
# kk-hk  -3.046756  -9.209135  3.1156218 0.5740990
# kk-kh  -7.567151 -14.975748 -0.1585548 0.0432978

locM$Affx.88864387 <- as.character(locM$Affx.88864387)
for (i in 1:length(locM$Affx.88864387)) {
  locM$Affx.88864387[i] <- if (locM$Affx.88864387[i] == "kk") locM$Affx.88864387[i] <- "B" else "A"}
locM$Affx.88864387 <- as.factor(locM$Affx.88864387)

source(filename)
summary(qtl_model)

# Model after step
# qtl_model<- lm(formula = RGH_mil_12 ~ Affx.88826196 + Affx.88860589 + Affx.88834745 +
# Affx.88850184 + Affx.88856302 + Affx.88848483 + Affx.88864419 +
#   Affx.88904022 + Affx.88876663 + Affx.88899370 + Affx.88892871,
# data = locM)

# Delete non sig
qtl_model<-lm(formula = RGH_mil_12 ~ Affx.88826196 + Affx.88860589 +
                Affx.88850184 + Affx.88856302 + Affx.88848483 + Affx.88864419 +
                Affx.88904022 + Affx.88876663 + Affx.88892871,
              data = locM)
summary(qtl_model)

outpdftitle<-paste("res",V1,".pdf", sep="")
outpdfgraph<-plot(qtl_model$residuals,qtl_model$fitted)
pdf(outpdftitle)
outpdfgraph
dev.off()

locM$V1<-paste("res",V1,".pdf", sep="")
name <-paste("RxH",V1)

outpdftitle<-paste("fitted",V1,".pdf", sep="")
locM_V<-which( colnames(locM)== V1 )
locM_V1<- paste("locM[as.numeric(rownames(data.frame(qtl_model$fitted))),",locM_V,"]",sep="")

outpdfgraph<- xyplot(as.vector(qtl_model$fitted) ~ locM[as.numeric(rownames(data.frame(qtl_model$fitted))),locM_V], type = c("p","r"), col= "black", col.line = "black", ylab="Predicted",xlab="Observed", main=name)
pdf(outpdftitle)
outpdfgraph
dev.off()

#Select only significant variable to work with

y <- variable.names(qtl_model)
y<- as.data.frame(y)

y <- lapply(y, function(x) {
  gsub("kk", "", x)
})
y <- lapply(y, function(x) {
  gsub("hk", "", x)
})
y <- lapply(y, function(x) {
  gsub("kh", "", x)
})
y <- lapply(y, function(x) {
  gsub("np", "", x)
})
y <- lapply(y, function(x) {
  gsub("lm", "", x)
})
y <- lapply(y, function(x) {
  gsub("B", "", x)
})
y<- as.data.frame(y)
y<- unique(y)
y <- y[-1,]
y<- as.data.frame(y)
y

# Creat models without each varriable in turn to calculate PRE score

for(d in 1:length(y[, 1]))
{
  qtl_modelc = paste("qtl_modelc",d,"= lm(",V1,"~",sep="")
  for(p in 1:length(y[, 1]))
  {
    if(p==d) next
    else
      qtl_modelc<- paste(qtl_modelc, y[p,1], '+')
    d+1
  }
  file_name1<-paste(V1,"qtl_modelc.R",sep="")
  qtl_modelc<-paste(qtl_modelc, ", data=locM)\n")
  qtl_modelc<-gsub("+ ,",",",fixed = TRUE, qtl_modelc)
  write(qtl_modelc, file=file_name1,append=TRUE)
}

source(file_name1)

# create formulas to calcualte rss
x<- paste("rss <- sum(residuals(qtl_model)^2)\n")
for (i in 1:length(y[, 1])) {
  file_name3<-paste(V1,"rss.R",sep="")
  x<-paste(x,'rss',i,"<- sum(residuals(qtl_modelc",i,")^2)\n", sep="")
  write(x, file=file_name3,append=TRUE)
}
source(file_name3)

# create formulas to calcualte PRE

for (i in 1:length(y[, 1])) {
  file_name2<-paste(V1,"X.R",sep="")
  x<-paste('X',i,"<- (rss",i,"-rss)/rss",i,"*100\n", sep="")
  write(x, file=file_name2,append=TRUE)
}

source(file_name2)

# create list of PRE
x<-"h=c("
for (i in 1:length(y[, 1])) {
  x<-paste(x,'X',i,",", sep="")
}
x<-paste(x,")", sep="")
x<-gsub(",)",")",fixed = TRUE, x)
write(x, file="h.r")
source("h.r")

# subset markers of interest
QTLs<- QTL[which(QTL$Rname %in% y[,1]),]
# Add PRE scores
QTLs$PRE=h
# Save file containing markers
Significant_QTL_V1.csv<-paste("Significant_QTL_",V1,".csv",sep="")
filename<-paste("Significant_QTL_",V1,".csv", sep="")
write.csv(QTLs, file=filename)

#######
rm(QTL)
# End for RGH_mil_12
V1 = "RGH_mil_13"

locM<-read.csv("locM.csv")

R_QTL_S_V1<-paste("R_QTL_S_",V1,".csv",sep="")
QTL<-read.csv(R_QTL_S_V1, header=TRUE)
head(QTL)

QTL$Rname<- gsub('-', '.', QTL$Rname)
x = paste("qtl_model= step(lm(",V1,"~")
for (i in 1:length(QTL[, 4])) {
  n =i+1
  print(n)
  x<- paste(x, QTL[i,"Rname"], '+')
  print(x)}
x<-paste(x, ", data=locM))")
x<-gsub("+ ,",",",fixed = TRUE, x)
print(x,quote = FALSE)

filename=paste0("qtl_model_",V1,".R")
write(x, file=filename)
# source(filename)
# summary(qtl_model)

locM$Affx.88810984 <- as.character(locM$Affx.88810984)
for (i in 1:length(locM$Affx.88810984)) {
  locM$Affx.88810984[i] <- if (locM$Affx.88810984[i] == "np") locM$Affx.88810984[i] <- "A" else "B"}
locM$Affx.88810984 <- as.factor(locM$Affx.88810984)

locM$Affx.88863329 <- as.character(locM$Affx.88863329)
for (i in 1:length(locM$Affx.88863329)) {
  locM$Affx.88863329[i] <- if (locM$Affx.88863329[i] == "lm") locM$Affx.88863329[i] <- "A" else "B"}
locM$Affx.88863329 <- as.factor(locM$Affx.88863329)

locM$Affx.88873309 <- as.character(locM$Affx.88873309)
for (i in 1:length(locM$Affx.88873309)) {
  locM$Affx.88873309[i] <- if (locM$Affx.88873309[i] == "np") locM$Affx.88873309[i] <- "A" else "B"}
locM$Affx.88873309 <- as.factor(locM$Affx.88873309)

locM$Affx.88879735 <- as.character(locM$Affx.88879735)
for (i in 1:length(locM$Affx.88879735)) {
  locM$Affx.88879735[i] <- if (locM$Affx.88879735[i] == "lm") locM$Affx.88879735[i] <- "A" else "B"}
locM$Affx.88879735 <- as.factor(locM$Affx.88879735)

#                 hh          hk            kh          kk
# Affx-88857591 29.02307692  35.11666667    33.04482759 35.62222222
TukeyHSD(aov(RGH_mil_13~Affx.88857591, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = RGH_mil_13 ~ Affx.88857591, data = locM)
#
# $Affx.88857591
# diff       lwr       upr     p adj
# hk-hh  6.0935897  1.162669 11.024510 0.0086954
# kh-hh  4.0217507 -1.486735  9.530237 0.2339342
# kk-hh  6.5991453  1.445924 11.752367 0.0060017
# kh-hk -2.0718391 -7.810307  3.666628 0.7846274
# kk-hk  0.5055556 -4.892804  5.903915 0.9949231
# kk-kh  2.5773946 -3.353181  8.507970 0.6724024

locM$Affx.88857591 <- as.character(locM$Affx.88857591)
for (i in 1:length(locM$Affx.88857591)) {
  locM$Affx.88857591[i] <- if (locM$Affx.88857591[i] == "hh") locM$Affx.88857591[i] <- "B" else "A"}
locM$Affx.88857591 <- as.factor(locM$Affx.88857591)

#                 hh          hk          kh          kk
# Affx-88878897  30.84827586  36.8754717    30.82692308 30.86862745
TukeyHSD(aov(RGH_mil_13~Affx.88878897, data=locM))
#   Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = RGH_mil_13 ~ Affx.88878897, data = locM)
#
# $Affx.88878897
# diff         lwr       upr     p adj
# hk-hh  6.02719584   0.5429644 11.511427 0.0250355
# kh-hh -0.02135279  -6.4340607  6.391355 0.9999998
# kk-hh  0.02035159  -5.5017792  5.542482 0.9999997
# kh-hk -6.04854862 -11.7336043 -0.363493 0.0322200
# kk-hk -6.00684425 -10.6642002 -1.349488 0.0055552
# kk-kh  0.04170437  -5.6799206  5.763329 0.9999976

locM$Affx.88878897 <- as.character(locM$Affx.88878897)
for (i in 1:length(locM$Affx.88878897)) {
  locM$Affx.88878897[i] <- if (locM$Affx.88878897[i] == "hk") locM$Affx.88878897[i] <- "A" else "B"}
locM$Affx.88878897 <- as.factor(locM$Affx.88878897)

#                 hh          hk          kh            kk
# Affx-88857591 29.02307692  35.11666667    33.04482759 35.62222222
TukeyHSD(aov(RGH_mil_13~Affx.88857591, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = RGH_mil_13 ~ Affx.88857591, data = locM)
#
# $Affx.88857591
# diff       lwr       upr     p adj
# hk-hh  6.0935897  1.162669 11.024510 0.0086954
# kh-hh  4.0217507 -1.486735  9.530237 0.2339342
# kk-hh  6.5991453  1.445924 11.752367 0.0060017
# kh-hk -2.0718391 -7.810307  3.666628 0.7846274
# kk-hk  0.5055556 -4.892804  5.903915 0.9949231
# kk-kh  2.5773946 -3.353181  8.507970 0.6724024

locM$Affx.88857591 <- as.character(locM$Affx.88857591)
for (i in 1:length(locM$Affx.88857591)) {
  locM$Affx.88857591[i] <- if (locM$Affx.88857591[i] == "hh") locM$Affx.88857591[i] <- "B" else "A"}
locM$Affx.88857591 <- as.factor(locM$Affx.88857591)
source(filename)
summary(qtl_model)

# qtl_model <-lm(formula = RGH_mil_13 ~ Affx.88810984 + Affx.88809245 + Affx.88826196 +
# Affx.88836959 + Affx.88857344 + Affx.88857559 + Affx.88863329 +
#   Affx.88873309 + Affx.88877358 + Affx.88878897 + Affx.88904022 +
#   Affx.88890499 + Affx.88898654, data = locM)

# Removing NS

qtl_model <-lm(formula = RGH_mil_13 ~ Affx.88810984 + Affx.88809245 + Affx.88826196 +
                 Affx.88836959 + Affx.88857344 + Affx.88857559 + Affx.88863329 +
                 Affx.88873309 + Affx.88877358 + Affx.88878897 + Affx.88904022 +
                  Affx.88898654, data = locM)
summary(qtl_model)

outpdftitle<-paste("res",V1,".pdf", sep="")
outpdfgraph<-plot(qtl_model$residuals,qtl_model$fitted)
pdf(outpdftitle)
outpdfgraph
dev.off()

locM$V1<-paste("res",V1,".pdf", sep="")
name <-paste("RxH",V1)

outpdftitle<-paste("fitted",V1,".pdf", sep="")
locM_V<-which( colnames(locM)== V1 )
locM_V1<- paste("locM[as.numeric(rownames(data.frame(qtl_model$fitted))),",locM_V,"]",sep="")

outpdfgraph<- xyplot(as.vector(qtl_model$fitted) ~ locM[as.numeric(rownames(data.frame(qtl_model$fitted))),locM_V], type = c("p","r"), col= "black", col.line = "black", ylab="Predicted",xlab="Observed", main=name)
pdf(outpdftitle)
outpdfgraph
dev.off()

#Select only significant variable to work with

y <- variable.names(qtl_model)
y<- as.data.frame(y)

y <- lapply(y, function(x) {
  gsub("kk", "", x)
})
y <- lapply(y, function(x) {
  gsub("hk", "", x)
})
y <- lapply(y, function(x) {
  gsub("kh", "", x)
})
y <- lapply(y, function(x) {
  gsub("np", "", x)
})
y <- lapply(y, function(x) {
  gsub("lm", "", x)
})
y <- lapply(y, function(x) {
  gsub("B", "", x)
})
y<- as.data.frame(y)
y<- unique(y)
y <- y[-1,]
y<- as.data.frame(y)
y

# Creat models without each varriable in turn to calculate PRE score

for(d in 1:length(y[, 1]))
{
  qtl_modelc = paste("qtl_modelc",d,"= lm(",V1,"~",sep="")
  for(p in 1:length(y[, 1]))
  {
    if(p==d) next
    else
      qtl_modelc<- paste(qtl_modelc, y[p,1], '+')
    d+1
  }
  file_name1<-paste(V1,"qtl_modelc.R",sep="")
  qtl_modelc<-paste(qtl_modelc, ", data=locM)\n")
  qtl_modelc<-gsub("+ ,",",",fixed = TRUE, qtl_modelc)
  write(qtl_modelc, file=file_name1,append=TRUE)
}

source(file_name1)

# create formulas to calcualte rss
x<- paste("rss <- sum(residuals(qtl_model)^2)\n")
for (i in 1:length(y[, 1])) {
  file_name3<-paste(V1,"rss.R",sep="")
  x<-paste(x,'rss',i,"<- sum(residuals(qtl_modelc",i,")^2)\n", sep="")
  write(x, file=file_name3,append=TRUE)
}
source(file_name3)

# create formulas to calcualte PRE

for (i in 1:length(y[, 1])) {
  file_name2<-paste(V1,"X.R",sep="")
  x<-paste('X',i,"<- (rss",i,"-rss)/rss",i,"*100\n", sep="")
  write(x, file=file_name2,append=TRUE)
}

source(file_name2)

# create list of PRE
x<-"h=c("
for (i in 1:length(y[, 1])) {
  x<-paste(x,'X',i,",", sep="")
}
x<-paste(x,")", sep="")
x<-gsub(",)",")",fixed = TRUE, x)
write(x, file="h.r")
source("h.r")

# subset markers of interest
QTLs<- QTL[which(QTL$Rname %in% y[,1]),]
# Add PRE scores
QTLs$PRE=h
# Save file containing markers
Significant_QTL_V1.csv<-paste("Significant_QTL_",V1,".csv",sep="")
filename<-paste("Significant_QTL_",V1,".csv", sep="")
write.csv(QTLs, file=filename)

#######
rm(QTL)
#End of RGH_mil_13

V1 = "RGH_mil_14"

locM<-read.csv("locM.csv")

R_QTL_S_V1<-paste("R_QTL_S_",V1,".csv",sep="")
QTL<-read.csv(R_QTL_S_V1, header=TRUE)
head(QTL)

QTL$Rname<- gsub('-', '.', QTL$Rname)
x = paste("qtl_model= step(lm(",V1,"~")
for (i in 1:length(QTL[, 4])) {
  n =i+1
  print(n)
  x<- paste(x, QTL[i,"Rname"], '+')
  print(x)}
x<-paste(x, ", data=locM))")
x<-gsub("+ ,",",",fixed = TRUE, x)
print(x,quote = FALSE)

filename=paste0("qtl_model_",V1,".R")
write(x, file=filename)
# source(filename)
# summary(qtl_model)

locM$Affx.88810984 <- as.character(locM$Affx.88810984)
for (i in 1:length(locM$Affx.88810984)) {
  locM$Affx.88810984[i] <- if (locM$Affx.88810984[i] == "np") locM$Affx.88810984[i] <- "A" else "B"}
locM$Affx.88810984 <- as.factor(locM$Affx.88810984)

locM$Affx.88836048 <- as.character(locM$Affx.88836048)
for (i in 1:length(locM$Affx.88836048)) {
  locM$Affx.88836048[i] <- if (locM$Affx.88836048[i] == "lm") locM$Affx.88836048[i] <- "A" else "B"}
locM$Affx.88836048 <- as.factor(locM$Affx.88836048)

locM$Affx.88832542 <- as.character(locM$Affx.88832542)
for (i in 1:length(locM$Affx.88832542)) {
  locM$Affx.88832542[i] <- if (locM$Affx.88832542[i] == "np") locM$Affx.88832542[i] <- "A" else "B"}
locM$Affx.88832542 <- as.factor(locM$Affx.88832542)

locM$Affx.88852176 <- as.character(locM$Affx.88852176)
for (i in 1:length(locM$Affx.88852176)) {
  locM$Affx.88852176[i] <- if (locM$Affx.88852176[i] == "np") locM$Affx.88852176[i] <- "A" else "B"}
locM$Affx.88852176 <- as.factor(locM$Affx.88852176)

locM$Affx.88893505 <- as.character(locM$Affx.88893505)
for (i in 1:length(locM$Affx.88893505)) {
  locM$Affx.88893505[i] <- if (locM$Affx.88893505[i] == "np") locM$Affx.88893505[i] <- "A" else "B"}
locM$Affx.88893505 <- as.factor(locM$Affx.88893505)

locM$Affx.88895830 <- as.character(locM$Affx.88895830)
for (i in 1:length(locM$Affx.88895830)) {
  locM$Affx.88895830[i] <- if (locM$Affx.88895830[i] == "lm") locM$Affx.88895830[i] <- "A" else "B"}
locM$Affx.88895830 <- as.factor(locM$Affx.88895830)


#                 hh          hk          kh          kk
# Affx-88826246  103.504032  82.95227307    105.077703  82.79513863
TukeyHSD(aov(RGH_mil_14~Affx.88826246, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# #Fit: aov(formula = RGH_mil_14 ~ Affx.88826246, data = locM)
#
# $Affx.88826246
# diff        lwr       upr     p adj
# hk-hh -20.5517589 -36.546897 -4.556621 0.0057889
# kh-hh   1.5736710 -15.767312 18.914654 0.9953741
# kk-hh -20.7088933 -38.159328 -3.258458 0.0128746
# kh-hk  22.1254299   6.982404 37.268456 0.0012014
# kk-hk  -0.1571344 -15.425377 15.111108 0.9999931
# kk-kh -22.2825643 -38.955431 -5.609698 0.0037189

locM$Affx.88826246 <- as.character(locM$Affx.88826246)
for (i in 1:length(locM$Affx.88826246)) {
  locM$Affx.88826246[i] <- if (locM$Affx.88826246[i] == "kk") locM$Affx.88826246[i] <- "B" else if (locM$Affx.88826246[i] == "hk") locM$Affx.88826246[i] <- "B" else "A"}
locM$Affx.88826246 <- as.factor(locM$Affx.88826246)

#                 hh          hk          kh          kk
# Affx-88838175  101.6726188  96.29545476   83.19285741 84.75
TukeyHSD(aov(RGH_mil_14~Affx.88838175, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = RGH_mil_14 ~ Affx.88838175, data = locM)
#
# $Affx.88838175
# diff       lwr        upr     p adj
# hk-hh  -5.377164 -21.29321 10.5388837 0.8165500
# kh-hh -18.479761 -35.36564 -1.5938819 0.0258921
# kk-hh -16.922619 -33.44093 -0.4043116 0.0423930
# kh-hk -13.102597 -29.81313  3.6079308 0.1792146
# kk-hk -11.545455 -27.88447  4.7935568 0.2608762
# kk-kh   1.557143 -15.72799 18.8422706 0.9954734

locM$Affx.88838175 <- as.character(locM$Affx.88838175)
for (i in 1:length(locM$Affx.88838175)) {
  locM$Affx.88838175[i] <- if (locM$Affx.88838175[i] == "kk") locM$Affx.88838175[i] <- "B" else if (locM$Affx.88838175[i] == "kh") locM$Affx.88838175[i] <- "B" else "A"}
locM$Affx.88838175 <- as.factor(locM$Affx.88838175)

#                 hh          hk          kh          kk
# Affx-88853896  84.45500019  88.74999975   85.72265625 110.2852567
TukeyHSD(aov(RGH_mil_14~Affx.88853896, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
#Fit: aov(formula = RGH_mil_14 ~ Affx.88853896, data = locM)
#
# $Affx.88853896
# diff        lwr      upr     p adj
# hk-hh  4.295000 -11.060860 19.65086 0.8864035
# kh-hh  1.267656 -14.885468 17.42078 0.9969924
# kk-hh 25.830256  10.586665 41.07385 0.0001168
# kh-hk -3.027344 -20.146873 14.09219 0.9677231
# kk-hk 21.535257   5.271126 37.79939 0.0041407
# kk-kh 24.562600   7.543701 41.58150 0.0014195

locM$Affx.88853896 <- as.character(locM$Affx.88853896)
for (i in 1:length(locM$Affx.88853896)) {
  locM$Affx.88853896[i] <- if (locM$Affx.88853896[i] == "hh") locM$Affx.88853896[i] <- "B" else if (locM$Affx.88853896[i] == "kh") locM$Affx.88853896[i] <- "B" else "A"}
locM$Affx.88853896 <- as.factor(locM$Affx.88853896)

#                 hh          hk          kh          kk
# Affx-88866836  104.0052091  94.51694931   94.42999962 82.47303903
TukeyHSD(aov(RGH_mil_14~Affx.88866836, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
#Fit: aov(formula = RGH_mil_14 ~ Affx.88866836, data = locM)
#
# $Affx.88866836
# diff       lwr       upr     p adj
# hk-hh  -9.48825981 -27.41102  8.434504 0.5169706
# kh-hh  -9.57520950 -30.73055 11.580129 0.6432043
# kk-hh -21.53217010 -39.85690 -3.207444 0.0141067
# kh-hk  -0.08694969 -17.75307 17.579171 0.9999992
# kk-hk -12.04391028 -26.19803  2.110208 0.1250744
# kk-kh -11.95696059 -30.03075  6.116831 0.3177701

locM$Affx.88866836 <- as.character(locM$Affx.88866836)
for (i in 1:length(locM$Affx.88866836)) {
  locM$Affx.88866836[i] <- if (locM$Affx.88866836[i] == "kk") locM$Affx.88866836[i] <- "C" else if (locM$Affx.88866836[i] == "hh") locM$Affx.88866836[i] <- "A" else "B"}
locM$Affx.88866836 <- as.factor(locM$Affx.88866836)

source(filename)
summary(qtl_model)


# # qtl_model <- lm(formula = RGH_mil_14 ~ Affx.88810984 + Affx.88809245 + Affx.88826246 +
# Affx.88838175 + Affx.88836048 + Affx.88832542 + Affx.88852176 +
#   Affx.88853896 + Affx.88857559 + Affx.88866836 + Affx.88874172 +
#   Affx.88904022 + Affx.88893505 + Affx.88895830 + Affx.88901393,
# data = locM)

qtl_model <- lm(formula = RGH_mil_14 ~ Affx.88810984 + Affx.88826246 +
                  Affx.88838175 + Affx.88836048 + Affx.88832542 + Affx.88852176 +
                  Affx.88853896 + Affx.88866836 +
                  Affx.88904022 + Affx.88893505 + Affx.88895830 + Affx.88901393,
                data = locM)

summary(qtl_model)

outpdftitle<-paste("res",V1,".pdf", sep="")
outpdfgraph<-plot(qtl_model$residuals,qtl_model$fitted)
pdf(outpdftitle)
outpdfgraph
dev.off()

locM$V1<-paste("res",V1,".pdf", sep="")
name <-paste("RxH",V1)

outpdftitle<-paste("fitted",V1,".pdf", sep="")
locM_V<-which( colnames(locM)== V1 )
locM_V1<- paste("locM[as.numeric(rownames(data.frame(qtl_model$fitted))),",locM_V,"]",sep="")

outpdfgraph<- xyplot(as.vector(qtl_model$fitted) ~ locM[as.numeric(rownames(data.frame(qtl_model$fitted))),locM_V], type = c("p","r"), col= "black", col.line = "black", ylab="Predicted",xlab="Observed", main=name)
pdf(outpdftitle)
outpdfgraph
dev.off()

#Select only significant variable to work with

y <- variable.names(qtl_model)
y<- as.data.frame(y)

y <- lapply(y, function(x) {
  gsub("kk", "", x)
})
y <- lapply(y, function(x) {
  gsub("hk", "", x)
})
y <- lapply(y, function(x) {
  gsub("kh", "", x)
})
y <- lapply(y, function(x) {
  gsub("np", "", x)
})
y <- lapply(y, function(x) {
  gsub("lm", "", x)
})
y <- lapply(y, function(x) {
  gsub("B", "", x)
})
y <- lapply(y, function(x) {
  gsub("C", "", x)
})
y<- as.data.frame(y)
y<- unique(y)
y <- y[-1,]
y<- as.data.frame(y)
y

# Creat models without each varriable in turn to calculate PRE score

for(d in 1:length(y[, 1]))
{
  qtl_modelc = paste("qtl_modelc",d,"= lm(",V1,"~",sep="")
  for(p in 1:length(y[, 1]))
  {
    if(p==d) next
    else
      qtl_modelc<- paste(qtl_modelc, y[p,1], '+')
    d+1
  }
  file_name1<-paste(V1,"qtl_modelc.R",sep="")
  qtl_modelc<-paste(qtl_modelc, ", data=locM)\n")
  qtl_modelc<-gsub("+ ,",",",fixed = TRUE, qtl_modelc)
  write(qtl_modelc, file=file_name1,append=TRUE)
}

source(file_name1)

# create formulas to calcualte rss
x<- paste("rss <- sum(residuals(qtl_model)^2)\n")
for (i in 1:length(y[, 1])) {
  file_name3<-paste(V1,"rss.R",sep="")
  x<-paste(x,'rss',i,"<- sum(residuals(qtl_modelc",i,")^2)\n", sep="")
  write(x, file=file_name3,append=TRUE)
}
source(file_name3)

# create formulas to calcualte PRE

for (i in 1:length(y[, 1])) {
  file_name2<-paste(V1,"X.R",sep="")
  x<-paste('X',i,"<- (rss",i,"-rss)/rss",i,"*100\n", sep="")
  write(x, file=file_name2,append=TRUE)
}

source(file_name2)

# create list of PRE
x<-"h=c("
for (i in 1:length(y[, 1])) {
  x<-paste(x,'X',i,",", sep="")
}
x<-paste(x,")", sep="")
x<-gsub(",)",")",fixed = TRUE, x)
write(x, file="h.r")
source("h.r")

# subset markers of interest
QTLs<- QTL[which(QTL$Rname %in% y[,1]),]
# Add PRE scores
QTLs$PRE=h
# Save file containing markers
Significant_QTL_V1.csv<-paste("Significant_QTL_",V1,".csv",sep="")
filename<-paste("Significant_QTL_",V1,".csv", sep="")
write.csv(QTLs, file=filename)

#######
#End of 14

locM<-read.csv("locM.csv")
V1 = "RGH_mil_16"

rm(QTL)
R_QTL_S_V1<-paste("R_QTL_S_",V1,".csv",sep="")
QTL<-read.csv(R_QTL_S_V1, header=TRUE)
head(QTL)

QTL$Rname<- gsub('-', '.', QTL$Rname)
x = paste("qtl_model= step(lm(",V1,"~")
for (i in 1:length(QTL[, 4])) {
  n =i+1
  print(n)
  x<- paste(x, QTL[i,"Rname"], '+')
  print(x)}
x<-paste(x, ", data=locM))")
x<-gsub("+ ,",",",fixed = TRUE, x)
print(x,quote = FALSE)

filename=paste0("qtl_model_",V1,".R")
write(x, file=filename)
# source(filename)
# summary(qtl_model)

locM$Affx.88837970 <- as.character(locM$Affx.88837970)
for (i in 1:length(locM$Affx.88837970)) {
  locM$Affx.88837970[i] <- if (locM$Affx.88837970[i] == "lm") locM$Affx.88837970[i] <- "A" else "B"}
locM$Affx.88837970 <- as.factor(locM$Affx.88837970)


#                 hh          hk            kh          kk
# Affx-88890287 52.43952381  53.56481482    49.51666667 45.42956989
TukeyHSD(aov(RGH_mil_16~Affx.88890287, data=locM))
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
#Fit: aov(formula = RGH_mil_16 ~ Affx.88890287, data = locM)
#
# $Affx.88890287
# diff        lwr        upr     p adj
# hk-hh  1.125291  -6.116159  8.3667409 0.9773114
# kh-hh -2.922857  -7.624157  1.7784425 0.3701108
# kk-hh -7.009954 -11.788674 -2.2312341 0.0012386
# kh-hk -4.048148 -11.334345  3.2380485 0.4709913
# kk-hk -8.135245 -15.471634 -0.7988557 0.0235303
# kk-kh -4.087097  -8.933356  0.7591624 0.1294589

locM$Affx.88890287 <- as.character(locM$Affx.88890287)
for (i in 1:length(locM$Affx.88890287)) {
  locM$Affx.88890287[i] <- if (locM$Affx.88890287[i] == "kk") locM$Affx.88890287[i] <- "B" else "A"}
locM$Affx.88890287 <- as.factor(locM$Affx.88890287)

source(filename)
summary(qtl_model)

# Model after step
# qtl_model<- lm(formula = RGH_mil_16 ~ Affx.88816105 + Affx.88826833 + Affx.88851610 +
# Affx.88846808 + Affx.88837970 + Affx.88890287, data = locM)

# Delete non sig
qtl_model<-lm(formula = RGH_mil_16 ~ Affx.88826833 + Affx.88851610 +
                Affx.88846808 + Affx.88837970 + Affx.88890287, data = locM)
summary(qtl_model)

outpdftitle<-paste("res",V1,".pdf", sep="")
outpdfgraph<-plot(qtl_model$residuals,qtl_model$fitted)
pdf(outpdftitle)
outpdfgraph
dev.off()

locM$V1<-paste("res",V1,".pdf", sep="")
name <-paste("RxH",V1)

outpdftitle<-paste("fitted",V1,".pdf", sep="")
locM_V<-which( colnames(locM)== V1 )
locM_V1<- paste("locM[as.numeric(rownames(data.frame(qtl_model$fitted))),",locM_V,"]",sep="")

outpdfgraph<- xyplot(as.vector(qtl_model$fitted) ~ locM[as.numeric(rownames(data.frame(qtl_model$fitted))),locM_V], type = c("p","r"), col= "black", col.line = "black", ylab="Predicted",xlab="Observed", main=name)
pdf(outpdftitle)
outpdfgraph
dev.off()

#Select only significant variable to work with

y <- variable.names(qtl_model)
y<- as.data.frame(y)

y <- lapply(y, function(x) {
  gsub("kk", "", x)
})
y <- lapply(y, function(x) {
  gsub("hk", "", x)
})
y <- lapply(y, function(x) {
  gsub("kh", "", x)
})
y <- lapply(y, function(x) {
  gsub("np", "", x)
})
y <- lapply(y, function(x) {
  gsub("lm", "", x)
})
y <- lapply(y, function(x) {
  gsub("B", "", x)
})
y<- as.data.frame(y)
y<- unique(y)
y <- y[-1,]
y<- as.data.frame(y)
y

# Creat models without each varriable in turn to calculate PRE score

for(d in 1:length(y[, 1]))
{
  qtl_modelc = paste("qtl_modelc",d,"= lm(",V1,"~",sep="")
  for(p in 1:length(y[, 1]))
  {
    if(p==d) next
    else
      qtl_modelc<- paste(qtl_modelc, y[p,1], '+')
    d+1
  }
  file_name1<-paste(V1,"qtl_modelc.R",sep="")
  qtl_modelc<-paste(qtl_modelc, ", data=locM)\n")
  qtl_modelc<-gsub("+ ,",",",fixed = TRUE, qtl_modelc)
  write(qtl_modelc, file=file_name1,append=TRUE)
}

source(file_name1)

# create formulas to calcualte rss
x<- paste("rss <- sum(residuals(qtl_model)^2)\n")
for (i in 1:length(y[, 1])) {
  file_name3<-paste(V1,"rss.R",sep="")
  x<-paste(x,'rss',i,"<- sum(residuals(qtl_modelc",i,")^2)\n", sep="")
  write(x, file=file_name3,append=TRUE)
}
source(file_name3)

# create formulas to calcualte PRE

for (i in 1:length(y[, 1])) {
  file_name2<-paste(V1,"X.R",sep="")
  x<-paste('X',i,"<- (rss",i,"-rss)/rss",i,"*100\n", sep="")
  write(x, file=file_name2,append=TRUE)
}

source(file_name2)

# create list of PRE
x<-"h=c("
for (i in 1:length(y[, 1])) {
  x<-paste(x,'X',i,",", sep="")
}
x<-paste(x,")", sep="")
x<-gsub(",)",")",fixed = TRUE, x)
write(x, file="h.r")
source("h.r")

# subset markers of interest
QTLs<- QTL[which(QTL$Rname %in% y[,1]),]
# Add PRE scores
QTLs$PRE=h
# Save file containing markers
Significant_QTL_V1.csv<-paste("Significant_QTL_",V1,".csv",sep="")
filename<-paste("Significant_QTL_",V1,".csv", sep="")
write.csv(QTLs, file=filename)

#######
rm(QTL)
# End for RGH_mil_16

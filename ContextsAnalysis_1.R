library(faraway)
library(foreign)
library(vioplot)
library(lme4)
library(latticeExtra)


### Define special functions

coeftestS4 <- function (x, level = 0.95,full=FALSE)
  {
    if(full){
      ll.model <- -logLik(x)
      aic.model <- summary(x)@AICtab[1]
    }
    level2 <- 1-((1-level)/2)
    coefs <- summary(x)@coefs[,1]
    serrors <- summary(x)@coefs[,2]
    result <- cbind(coefs, coefs+t(outer(c(-1,1),qnorm(level2)*serrors)))
    if(full){
      ll.null <- -logLik(update(x, ~ +(1|CountryDistrict)+(1|country)))
      cat("LogLik: ", 2*(ll.null-ll.model),"\nAIC: ",as.numeric(aic.model),"\n")
    }
    colnames(result) <- c("Coefficients",paste(level*100,"% Conf.",sep=""),"Interval")
    return(result)
  }


stdize  <- function(x){
  (x-mean(x,na.rm=TRUE))/(sd(x,na.rm=TRUE))
                                        #x
}

## Read data in

wasdb.full <- read.csv("ContextsData.csv")


## Estimate Multilevel Models
cat("\n##############################################\n")
cat("\nHopeless Votes\n")
was.model.nb <- lmer(logitHopeless~stdize(laggedvol)
                     +stdize(npshareofm)
                     +stdize(log(magnitude))
                     +I(experience<=1)
                     +I(experience==2)
                     +I(experience==3)
                     +I(experience==4)
                     +as.factor(elec_sys)
                     +compensatory
                     +federal
                     +presidential
                     +stdize(eth_ling85)
                     +(1|CountryDistrict)
                     + (1|country),
                     data=wasdb.full,control=list(maxIter=1000))

print(coeftestS4(was.model.nb, full=TRUE))

cat("\nSfRatio \n")
sf.model.nb <- lmer(logitSF ~stdize(laggedvol)
                    +stdize(npshareofm)
                    +stdize(log(magnitude))
                    +I(experience<=1)
                    +I(experience==2)
                    +I(experience==3)
                    +I(experience==4)
                    +as.factor(elec_sys)
                    +compensatory
                    +federal
                    +presidential
                    +stdize(eth_ling85)
                    +(1|CountryDistrict)
                    +(1|country),
                    data=wasdb.full,control=list(maxIter=1000))
print(coeftestS4(sf.model.nb, full=TRUE))

cat("\nENPlosing*Hopeless \n")
enpwas.model.nb <- lmer(I(((wastage))*enplosing) ~ stdize(laggedvol)
                        +stdize(npshareofm)
                        +stdize(log(magnitude))
                        +I(experience<=1)
                        +I(experience==2)
                        +I(experience==3)
                        +I(experience==4)
                        +as.factor(elec_sys)
                        +compensatory
                        +federal
                        +presidential
                        +stdize(eth_ling85)
                        +(1|CountryDistrict)
                        +(1|country),
                        data=wasdb.full,control=list(maxIter=1000))
print(coeftestS4(enpwas.model.nb,full=TRUE))


cat("\nWasted Votes\n")
hop.model.nb <- lmer(logitWastage~stdize(laggedvol)
                     +stdize(npshareofm)
                     +stdize(log(magnitude))
                     +I(experience<=1)
                     +I(experience==2)
                     +I(experience==3)
                     +I(experience==4)
                     +as.factor(elec_sys)
                     +compensatory
                     +federal
                     +presidential
                     +stdize(eth_ling85)
                     +(1|CountryDistrict)
                     + (1|country),
                     data=wasdb.full,control=list(maxIter=1000))

print(coeftestS4(hop.model.nb, full=TRUE))

cat("\nHopeless Votes (double rate)\n")
was.model.nb2 <- lmer(logitHopeless~stdize(laggedvol)
                                          +stdize(npshareofm)
                                          +stdize(log(magnitude))
                                          +I(experience2<=1)
                                          +I(experience2==2)
                                          +I(experience2==3)
                                          +I(experience2==4)
                                          +as.factor(elec_sys)
                                          +compensatory
                                          +federal
                                          +presidential
                                          +stdize(eth_ling85)
                                          +(1|CountryDistrict)
                                          + (1|country),
                                          data=wasdb.full,control=list(maxIter=1000))

print(coeftestS4(was.model.nb2, full=TRUE))

cat("\nSfRatio (double rate) \n")
sf.model.nb2 <- lmer(logitSF ~stdize(laggedvol)
                                        +stdize(npshareofm)
                                        +stdize(log(magnitude))
                                        +I(experience2<=1)
                                        +I(experience2==2)
                                        +I(experience2==3)
                                        +I(experience2==4)
                                        +as.factor(elec_sys)
                                        +compensatory
                                        +federal
                                        +presidential
                                        +stdize(eth_ling85)
                                        +(1|CountryDistrict)
                                        +(1|country),
                                        data=wasdb.full,control=list(maxIter=1000))
print(coeftestS4(sf.model.nb2, full=TRUE))

cat("\nENPlosing*Hopeless (double rate)\n")
enpwas.model.nb2 <- lmer(I(((wastage))*enplosing) ~ stdize(laggedvol)
                        +stdize(npshareofm)
                        +stdize(log(magnitude))
                        +I(experience2<=1)
                        +I(experience2==2)
                        +I(experience2==3)
                        +I(experience2==4)
                        +as.factor(elec_sys)
                        +compensatory
                        +federal
                        +presidential
                        +stdize(eth_ling85)
                        +(1|CountryDistrict)
                        +(1|country),
                        data=wasdb.full,control=list(maxIter=1000))
print(coeftestS4(enpwas.model.nb2,full=TRUE))


cat("\nNP as the DV\n")
np.model.nb <- lmer(npshareofm ~ stdize(laggedvol)
                    +stdize(log(magnitude))
                    +I(experience<=1)
                    +I(experience==2)
                    +I(experience==3)
                    +I(experience==4)
                    +as.factor(elec_sys)
                    +compensatory
                    +federal
                    +presidential
                    +stdize(eth_ling85)
                    +(1|CountryDistrict)
                    +(1|country),
                    data=wasdb.full,control=list(maxIter=1000))
print(coeftestS4(np.model.nb))


cat("\nPlotting...")


#### PLOTS ####

###Regression Outcome Plot

test <- coeftestS4(was.model.nb)[c(2:8),]

all.results <- array(NA,c(dim(test),3))
all.results[,,1] <- coeftestS4(was.model.nb)[c(2:8),]
all.results[,,2] <- coeftestS4(sf.model.nb)[c(2:8),]
all.results[,,3] <- coeftestS4(enpwas.model.nb)[c(2:8),]


dimnames(all.results) <- list(c("Lagged Volatility","New Parties / M", "Logged Magnitude","One or Fewer Prior Elections", "Two Prior Elections","Three Prior Elections", "Four Prior Elections"),NULL,c("Hopeless Votes","SF Ratio","Coordination Product"))


y.vars <- rep(dimnames(all.results)[[3]],each=dim(all.results)[1])
y.vars <- factor(y.vars,levels=c("Coordination Product", "SF Ratio", "Hopeless Votes"))
x.vals <- c(all.results[,1,])
l.bounds <- c(all.results[,2,])
u.bounds <- c(all.results[,3,])
cond.vars <- dimnames(all.results)[[1]]


trellis.device(device="pdf",file="~/MyCloud/HierViability/RegResults.pdf",width=5,height=12,col=TRUE)
figure <-dotplot(y.vars~x.vals | cond.vars,
        layout = c(1,dim(all.results)[1],1),
        panel = function(x,y,subscripts){
          panel.abline(v=0,
                       col="gray70",
                       lty=2)
          panel.abline(h=as.numeric(y),
                       col="gray70",
                       lwd=0.9)
          panel.segments(l.bounds[subscripts],
                         as.numeric(y),
                         u.bounds[subscripts],
                         as.numeric(y),
                         lwd=1.9,
                         col="black")
          panel.xyplot(x,
                       y,
                       col="black",
                       pch=19,
                       cex=0.5
                       )
        },
        subscripts=TRUE,
        strip=strip.custom(bg="black",par.strip.text=list(col="white")),
        xlab="Coefficient Value",
        ylab="Response Variable",
        index.cond=list(c(3,2,4,1,6,7,5))
        )
print(figure)
graphics.off()

####Violin Plots

x.vals <- c(wasdb.full$experience,wasdb.full$npshareofm,wasdb.full$laggedvol,wasdb.full$magnitude)
y.vals <- wasdb.full$country
cond.vals <- rep(c("Experience","New Parties","Volatility","District Magnitude"),each=dim(wasdb.full)[1])

trellis.device(device="pdf",file="~/MyCloud/HierViability/AllViolin.pdf",width=10,height=10,col=TRUE)

fig1 <- bwplot(wasdb.full$country ~ wasdb.full$experience, box.ratio=15,
                panel = function(..., box.ratio) {
                  panel.violin(..., col = "gray70",
                               varwidth = FALSE, box.ratio = box.ratio)
                  panel.bwplot(...,
                               col="black", cex=0.5,coef=0,
                               fill = "white", box.ratio = .1)
                } )
fig2 <- bwplot(wasdb.full$country ~ wasdb.full$npshareofm, box.ratio=15,
                panel = function(..., box.ratio) {
                  panel.violin(..., col = "gray70",
                               varwidth = FALSE, box.ratio = box.ratio)
                  panel.bwplot(...,
                               col="black", cex=0.5,coef=0,
                               fill = "white", box.ratio = .1)
                } )
fig3 <- bwplot(wasdb.full$country ~ wasdb.full$laggedvol, box.ratio=15,
                panel = function(..., box.ratio) {
                  panel.violin(..., col = "gray70",
                               varwidth = FALSE, box.ratio = box.ratio)
                  panel.bwplot(...,
                               col="black", cex=0.5,coef=0,
                               fill = "white", box.ratio = .1)
                } )
fig4 <- bwplot(wasdb.full$country ~ log(wasdb.full$magnitude), box.ratio=15,
                panel = function(..., box.ratio) {
                  panel.violin(..., col = "gray70",
                               varwidth = FALSE, box.ratio = box.ratio)
                  panel.bwplot(...,
                               col="black",cex=0.5,coef=0,
                               fill = "white", box.ratio = .1)
                } )
figure2 <- c("Prior Electoral Experience"=fig1,
            "New Parties / M"=fig2,
            "Lagged Volatility"=fig3,
            "Logged District Magnitude"=fig4,
            layout=c(2,2,1),
            x.same=FALSE,
            y.same=TRUE
            )
figure2$strip <- figure$strip
figure2$par.settings <- list(box.rectangle=list(col="white"),
          box.umbrella=list(lty=1, col="white"))
figure2$xlab <- ""
figure2$index.cond <- list(c(3,4,1,2))

print(figure2)
graphics.off()

###Violin of response 
trellis.device(device="pdf",file="~/MyCloud/HierViability/ViolinResp.pdf",width=15,height=7,col=TRUE)
fig1 <- bwplot(wasdb.full$country ~ wasdb.full$logitWas, box.ratio=15,
                panel = function(..., box.ratio) {
                  panel.violin(..., col = "gray70",
                               varwidth = FALSE, box.ratio = box.ratio)
                  panel.bwplot(...,
                               col="black", cex=0.5,coef=0,
                               fill = "white", box.ratio = .1)
                } )
fig2 <- bwplot(wasdb.full$country ~ wasdb.full$logitSF, box.ratio=15,
                panel = function(..., box.ratio) {
                  panel.violin(..., col = "gray70",
                               varwidth = FALSE, box.ratio = box.ratio)
                  panel.bwplot(...,
                               col="black", cex=0.5,coef=0,
                               fill = "white", box.ratio = .1)
                } )
fig3 <- bwplot(wasdb.full$country ~ I(wasdb.full$enplosing*wasdb.full$wastage), box.ratio=15,
                panel = function(..., box.ratio) {
                  panel.violin(..., col = "gray70",
                               varwidth = FALSE, box.ratio = box.ratio)
                  panel.bwplot(...,
                               col="black", cex=0.5,coef=0,
                               fill = "white", box.ratio = .1)
                } )
figure3 <- c("Hopeless Votes"=fig1,
            "SF Ratio"=fig2,
            "Coordination Product"=fig3,
            layout=c(3,1,1),
            x.same=FALSE,
            y.same=TRUE
            )
figure3$strip <- figure$strip
figure3$par.settings <- list(box.rectangle=list(col="white"),
          box.umbrella=list(lty=1, col="white"))
figure3$xlab <- ""

print(figure3)
graphics.off()



#### Predicted Values

plot.preds <- function(model,vals,level=0.85,new.data,x.lab,y.lab,y.lim,title,min.num,max.num,preds=FALSE, print.graph=FALSE){
  indices <- c(1,25,50,75,100)
  pred.vals <- new.data %*% model@fixef
  bm <- diag(sqrt(new.data %*% summary(model)@vcov %*% t(new.data)))
  bm.QNorm <- qnorm(level)*bm
  pred <- cbind(ilogit(pred.vals),ilogit(pred.vals-bm.QNorm),ilogit(pred.vals+bm.QNorm))
  if(preds){return(pred)}
  graph.result <- xyplot(pred[,1]~vals,
                         ylim=y.lim,
                         ylab=y.lab,
                         xlab=x.lab,
                         panel=function(...){
                           panel.polygon(x=c(vals,rev(vals)),
                                        y=c(pred[,3],rev(pred[,2])),
                                         ylim=y.lim,
                                         col='gray60',
                                         border=NA)
                           panel.xyplot(vals,pred[,1],
                                        type='l',
                                       ylim=y.lim,
                                       lwd=2.3,
                                       col='white')
                         },
                         strip=strip.custom(var.name=title,
                           bg="black",
                           par.strip.text=list(col="white")),
                         scale=list(relation="free",
                           x=list(at=vals[indices],
                             labels=round(seq(min.num,max.num,length.out=100),2)[indices])
                           )
                         )
  if(print.graph){
    print(graph.result)
  }else{
    return(graph.result)
  }
}


minVol <-  min(wasdb.full$laggedvol)#
minVolstd <-  min(stdize(wasdb.full$laggedvol))
maxVol <-  max(wasdb.full$laggedvol)#
maxVolstd <-  max(stdize(wasdb.full$laggedvol))

minNP <-  min(wasdb.full$npshareofm)#
maxNP <-  max(wasdb.full$npshareofm)#
minNPstd <-  min(stdize(wasdb.full$npshareofm))#
maxNPstd <-  max(stdize(wasdb.full$npshareofm))#


minMag <-  min(wasdb.full$magnitude)#
maxMag <-  max(wasdb.full$magnitude)#
minMagstd <-  min(stdize(log(wasdb.full$magnitude)))#
maxMagstd <-  max(stdize(log(wasdb.full$magnitude)))



##Volatility#

see.vals <- seq(minVolstd,maxVolstd,length.out=100)#
new.data <- cbind(1,see.vals,0,0,0,0,0,0,1,0,0,0,1,0)#

#pdf("PredVolHop.pdf")
volhop <- plot.preds(model=was.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed Lagged Volatility",
           y.lab="Share of Hopeless Votes",
           y.lim=c(0,0.14),
           title="Hopeless Votes",
           min.num=minVol,max.num=maxVol)
#dev.off()
#pdf("PredVolSF.pdf")
volsf <- plot.preds(model=sf.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed Lagged Volatility",
           y.lab="SF Ratio",
           y.lim=c(0,1),
           title="SF Ratio",
           min.num=minVol,max.num=maxVol)
#dev.off()
#pdf("PredVolWENLP.pdf")
volcp <- plot.preds(model=enpwas.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed Lagged Volatility",
           y.lab="Coordination product",
           y.lim=c(0.55,0.75),
           title="Coordination Product",
           min.num=minVol,max.num=maxVol)
#dev.off()




##New Parties#
see.vals <- seq(minNPstd,maxNPstd,length.out=100)#
new.data <- cbind(1,0,see.vals,0,0,0,0,0,1,0,0,0,1,0)#

#pdf("PredNPHop.pdf")
nphop <- plot.preds(model=was.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed New Parties (as share of M)",
           y.lab="Share of Hopeless Votes",
           y.lim=c(0,0.4),
           title="Hopeless Votes",
           min.num=minNP,max.num=maxNP)
#dev.off()
#pdf("PredNPSF.pdf")
npsf <- plot.preds(model=sf.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed New Parties (as share of M)",
           y.lab="SF Ratio",
           y.lim=c(0,1),
           title="SF Ratio",
           min.num=minNP,max.num=maxNP)
#dev.off()
#pdf("PredNPWENLP.pdf")
npcp <- plot.preds(model=enpwas.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed New Parties (as share of M)",
           y.lab="Coordination product",
           y.lim=c(0.5,1),
           title="Coordination Product",
           min.num=minNP,max.num=maxNP)
#dev.off()

##Magnitude#
see.vals <- seq(minMagstd,maxMagstd,length.out=100)#
new.data <- cbind(1,0,0,see.vals,0,0,0,0,1,0,0,0,1,0)#

#pdf("PredMagHop.pdf")
maghop <- plot.preds(model=was.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed District Magnitude",
           y.lab="Share of Hopeless Votes",
           y.lim=c(0,0.1),
           title="Hopeless Votes",
           min.num=minMag,max.num=maxMag)
#dev.off()

#pdf("PredMagSF.pdf")
magsf <- plot.preds(model=sf.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed District Magnitude",
           y.lab="SF Ratio",
           y.lim=c(0,1),
           title="SF Ratio",
min.num=minMag,max.num=maxMag)
#dev.off()

#pdf("PredMagWENLP.pdf")
magcp <- plot.preds(model=enpwas.model.nb,
           vals=see.vals,
           new.data=new.data,
           x.lab="Observed District Magnitude",
           y.lab="Coordination product",
           y.lim=c(0.4,0.75),
           title="Coordination Product",
           min.num=minMag,max.num=maxMag)
#dev.off()



##Experience#

plot.exp <- function(ylab,ylim,pred,title,print.graph=FALSE){
  graph.result <- xyplot(pred[,1]~(1:5),#
                         xlab="Number of Prior Elections Under the Same Rules",#
                         ylab=ylab,#
                         ylim=ylim,#
                         strip=strip.custom(var.name=title,
                           bg="black",
                           par.strip.text=list(col="white")),
                         scale=list(relation="free",
                           x=list(at=1:5,#
                             labels=c("One or Fewer", "Two","Three","Four","More than Four"))),                    
                         panel=function(...){
                           panel.polygon(c(1:5,rev(1:5)),#
                                         c(pred[,3],rev(pred[,2])),#
                                         ylim=ylim,
                                         col='gray60',#
                                         border=NA)
                           panel.lines(1:5,pred[,1],
                                       type="l",
                                       lwd=2,
                                       col='white')#
                         }
                         )
  if(print.graph){
    print(graph.result)
  }else{
    return(graph.result)
  }
}

new.data <- cbind(1,0,0,0,1,0,0,0,1,0,0,0,1,0)#
pred0 <- plot.preds(model=was.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,1,0,0,1,0,0,0,1,0)#
pred1 <- plot.preds(model=was.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,1,0,1,0,0,0,1,0)#
pred2 <- plot.preds(model=was.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,0,1,1,0,0,0,1,0)#
pred3 <- plot.preds(model=was.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,0,0,1,0,0,0,1,0)#
pred4 <- plot.preds(model=was.model.nb,
                    new.data=new.data,
                    preds=TRUE)

pred <- rbind(pred0,pred1,pred2,pred3,pred4)#

#pdf("PredExpHop.pdf")
exphop <- plot.exp(ylab="Share of Hopeless Votes",
         ylim=c(0,0.3),
         pred=pred,
         title="Hopeless Votes")
#dev.off()#

new.data <- cbind(1,0,0,0,1,0,0,0,1,0,0,0,1,0)#
pred0 <- plot.preds(model=enpwas.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,1,0,0,1,0,0,0,1,0)#
pred1 <- plot.preds(model=enpwas.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,1,0,1,0,0,0,1,0)#
pred2 <- plot.preds(model=enpwas.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,0,1,1,0,0,0,1,0)#
pred3 <- plot.preds(model=enpwas.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,0,0,1,0,0,0,1,0)#
pred4 <- plot.preds(model=enpwas.model.nb,
                    new.data=new.data,
                    preds=TRUE)

pred <- rbind(pred0,pred1,pred2,pred3,pred4)#

#pdf("PredExpWENLP.pdf")
expcp <- plot.exp(ylab="Coordination product",
         ylim=c(0.5,0.85),
         pred=pred,
         title="Coordination Product")
#dev.off()#

new.data <- cbind(1,0,0,0,1,0,0,0,1,0,0,0,1,0)#
pred0 <- plot.preds(model=sf.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,1,0,0,1,0,0,0,1,0)#
pred1 <- plot.preds(model=sf.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,1,0,1,0,0,0,1,0)#
pred2 <- plot.preds(model=sf.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,0,1,1,0,0,0,1,0)#
pred3 <- plot.preds(model=sf.model.nb,
                    new.data=new.data,
                    preds=TRUE)
new.data <- cbind(1,0,0,0,0,0,0,0,1,0,0,0,1,0)#
pred4 <- plot.preds(model=sf.model.nb,
                    new.data=new.data,
                    preds=TRUE)

pred <- rbind(pred0,pred1,pred2,pred3,pred4)#

#pdf("PredExpSF.pdf")
expsf <- plot.exp(ylab="SF Ratio",
         ylim=c(0,1),
         pred=pred,
         title="SF Ratio")
#dev.off()#
trellis.device(device="pdf",file="~/MyCloud/HierViability/PredMag.pdf",width=13,height=6,col=TRUE)
pred.graph1 <- c("Hopeless Votes"=maghop,
                  "SF Ratio"=magsf,
                 "Coordination Product"=magcp,
                  layout=c(3,1,1),
                  x.same=TRUE,
                  y.same=NA)
pred.graph1$xlab <- "Observed District Magnitude"
pred.graph1$strip <- figure$strip
pred.graph1$y.limits <- list(unlist(maghop$y.limits),unlist(magsf$y.limits),unlist(magcp$y.limits))
pred.graph1$ylab=""
print(pred.graph1)
dev.off()

trellis.device(device="pdf",file="~/MyCloud/HierViability/PredVol.pdf",width=13,height=6,col=TRUE)
pred.graph2 <- c("Hopeless Votes"=volhop,
                 "SF Ratio"=volsf,
                 "Coordination Product"=volcp,
                 layout=c(3,1,1),
                 x.same=TRUE,
                 y.same=NA)
pred.graph2$xlab <- "Observed Lagged Volatility"
pred.graph2$strip <- figure$strip
pred.graph2$y.limits <- list(unlist(volhop$y.limits),unlist(volsf$y.limits),unlist(volcp$y.limits))
pred.graph2$ylab=""
print(pred.graph2)
dev.off()

trellis.device(device="pdf",file="~/MyCloud/HierViability/PredNP.pdf",width=13,height=6,col=TRUE)
pred.graph3 <- c("Hopeless Votes"=nphop,
                 "SF Ratio"=npsf,
                 "Coordination Product"=npcp,
                 layout=c(3,1,1),
                 x.same=TRUE,
                 y.same=NA)
pred.graph3$xlab <- "Observed New Party Entry (as share of M)"
pred.graph3$strip <- figure$strip
pred.graph3$y.limits <- list(unlist(nphop$y.limits),unlist(npsf$y.limits),unlist(npcp$y.limits))
pred.graph3$ylab=""
print(pred.graph3)
dev.off()

trellis.device(device="pdf",file="~/MyCloud/HierViability/PredExp.pdf",width=13,height=6,col=TRUE)
pred.graph4 <- c("Hopeless Votes"=exphop,
                 "SF Ratio"=expsf,
                 "Coordination Product"=expcp,
                 layout=c(3,1,1),
                 x.same=TRUE,
                 y.same=NA)
pred.graph4$xlab <- "Observed Number of Prior Elections"
pred.graph4$strip <- figure$strip
pred.graph4$y.limits <- list(unlist(exphop$y.limits),unlist(expsf$y.limits),unlist(expcp$y.limits))
pred.graph4$ylab=""
print(pred.graph4)
dev.off()








###Scatter Plots New parties

y.vars <- rep(wasdb.full$npshareofm,times=3)
x.vars <- c(wasdb.full$laggedvol,log(wasdb.full$magnitude),wasdb.full$experience)
cond.vars <- rep(c("Lagged Volatility","Logged District Magnitude","Prior Experience"),each=dim(wasdb.full)[1])

trellis.device(device="pdf",file="~/MyCloud/HierViability/NewPartyScatter.pdf",width=14,height=5,col=TRUE)
figure <- xyplot(y.vars~x.vars|cond.vars,
                 strip=strip.custom(bg="black",par.strip.text=list(col="white")),
                 scales=list(x="free",y="same"),
                 ylab="New Parties as Share of M",
                 xlab="",
                 layout=c(3,1,1),
                 xlim=list(c(0,1),c(0,8),c(0,100)),
                 panel=function(x,y){
                   panel.xyplot(x,y,
                                jitter.x=TRUE,
                                jitter.y=TRUE,
                                col="black",
                                pch=19,
                                cex=0.6,
                                amount=0.4)
                   panel.lmline(x,y,
                                col="white",
                                lwd=3)
                 },
                 index.cond=list(c(3,1,2))
                 )
print(figure)
dev.off()

cat("done.\n")
       

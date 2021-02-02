library(lavaan)
library(lavaanPlot)
library(semTools)
library(fields)
library(akima)
library(viridis)
library(Hmisc)
library(geosphere)
library(ape)
library(lme4)
library(nlme)
library(piecewiseSEM)
library(qgraph)
a<-read.csv('coralfishglobal.csv',header=T)
a<-a[a$cor_diversity>0,]
a<-a[a$fish_sp_div>0,]
a<-a[a$reef_fraction>0,]
reg<-unique(a$Region)
b<-data.frame(cbind(a$fish_sp_div,
                    a$cor_diversity,
                    a$fraction_30m,
                    a$isolation,
                    a$reef_fraction,
                    a$mean_Present.Surface.Temperature.Mean,
                    a$mean_Present.Surface.Temperature.Range,
                    a$min_Present.Surface.Salinity.Mean,
                    a$mean_Present.Surface.Primary.productivity.Mean,
                    a$mean_bo_ph.asc,
                    abs(a$lat)
))
ok<-which(is.finite(rowSums(b)))
b<-b[ok,]
a<-a[ok,]
for (col in 1:(dim(b)[2])){
  b[,col]<-(b[,col]-min(b[,col]))/(max(b[,col])-min(b[,col]))
}

for (r in reg){b<-cbind(b,1*(a$Reg==r))}
b<-cbind(b,a$Oblig.Fac.Linked_proportion)
b<-cbind(b,a$Region)
b<-cbind(b,a$lon)
b<-cbind(b,a$lat)


cnames<-c('fish',
               'coral',
               'fr30m',
               'isolation',
               'reef_fraction',
               'tmean',
               'trange',
               'sal',
               'pp',
               'ph',
               'abs_lat',
               'reg1',
               'reg2',
               'reg3',
               'reg4',
               'reg5',
               'lit_dep',
               'region',
               'lon',
               'lat'
)

colnames(b)<-cnames
b<-cbind(b,b$tmean**2)
b<-cbind(b,b$trange**2)
b<-cbind(b,b$sal**2)
b<-cbind(b,b$pp**2)
b<-cbind(b,b$ph**2)

colnames(b)<-c(cnames,'tmean_2','trange_2','sal_2','pp_2','ph_2')
###all models
model0_sq <- '
        fish ~ coral+ abs_lat + tmean_2 + trange_2 + isolation 
          + sal_2 + pp_2 + reef_fraction  
          + reg1 + reg2 + reg3 + reg4
        coral ~ abs_lat + tmean_2 + trange_2 
          + isolation + sal_2 + pp_2 + reef_fraction 
          + fr30m + reg1 + reg2 + reg3 + reg4
  '

model0 <- '
        fish ~ coral+ abs_lat + tmean + trange + isolation 
          + sal + pp + reef_fraction  
          + reg1 + reg2 + reg3 + reg4
        coral ~ abs_lat + tmean + trange 
          + isolation + sal + pp + reef_fraction 
          + fr30m + reg1 + reg2 + reg3 + reg4
  '

fit0_sq <- sem(model0_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

fit0 <- sem(model0, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

model1 <- '
      fish ~ coral+ abs_lat + isolation + sal + pp 
      + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4
      coral ~ abs_lat + tmean + trange + isolation + sal 
      + pp + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4
'

model1_sq <- '
      fish ~ coral+ abs_lat + isolation + sal_2 + pp_2 
      + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4
      coral ~ abs_lat + tmean_2 + trange + isolation + sal_2 
      + pp_2 + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4
'


fit1 <- sem(model1, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE)#, test = "bollen.stine")

fit1_sq <- sem(model1_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE)#, test = "bollen.stine")

# fish-coral-environmental path analysis remove tmean and fr30m from fish
model2 <- '
      fish ~ coral + abs_lat +  isolation + sal 
      + pp + reef_fraction + reg1 + reg2 + reg3 + reg4
      coral ~ abs_lat + tmean + trange + isolation 
      + sal + pp + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4
'

model2_sq <- '
      fish ~ coral + abs_lat +  isolation + sal_2 
      + pp_2 + reef_fraction + reg1 + reg2 + reg3 + reg4
      coral ~ abs_lat + tmean_2 + trange_2 + isolation 
      + sal_2 + pp_2 + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4
'

fit2 <- sem(model2, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

fit2_sq <- sem(model2_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

# SEM with environmental latent variable driven by latitude
model3 <- '
  fish ~ isolation + environment + coral + reef_fraction + reg1 + reg2 + reg3 + reg4 
  coral ~ isolation + environment + reef_fraction + reg1 + reg2 + reg3 + reg4
  environment =~ tmean + trange + sal + pp
  reef_fraction ~ environment + fr30m
  environment ~ abs_lat
  trange ~~ pp
  tmean ~~ pp
  tmean ~~ trange
  coral ~~ reef_fraction
  '

model3_sq <- '
  fish ~ isolation + environment + coral + reef_fraction + reg1 + reg2 + reg3 + reg4 
  coral ~ isolation + environment + reef_fraction + reg1 + reg2 + reg3 + reg4
  environment =~ tmean_2 + trange_2 + sal_2 + pp_2
  reef_fraction ~ environment + fr30m
  environment ~ abs_lat
  trange_2 ~~ pp_2
  tmean_2 ~~ pp_2
  tmean_2 ~~ trange_2
  coral ~~ reef_fraction
  '


fit3 <- sem(model3, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")


fit3_sq <- sem(model3_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

# SEM  with environmental and reef latent variables
model4 <- '
  fish ~ isolation + environment + coral_reef + reg1 + reg2 + reg3 + reg4
  coral ~ isolation + environment + reef_fraction + reg1 + reg2 + reg3 + reg4
  environment =~ tmean + trange + sal + pp
  coral_reef =~ coral + reef_fraction
  reef_fraction ~ environment + fr30m
  environment ~ abs_lat
  trange ~~ pp
  tmean ~~ pp
  tmean ~~ trange
  trange ~~ sal
  sal ~~ pp
  '

model4_sq <- '
  fish ~ isolation + environment + coral_reef + reg1 + reg2 + reg3 + reg4
  coral ~ isolation + environment + reef_fraction + reg1 + reg2 + reg3 + reg4
  environment =~ tmean_2 + trange_2 + sal_2 + pp
  coral_reef =~ coral + reef_fraction
  reef_fraction ~ environment + fr30m
  environment ~ abs_lat
  trange_2 ~~ pp_2
  tmean_2 ~~ pp_2
  tmean_2 ~~ trange_2
  trange_2 ~~ sal_2
  sal_2 ~~ pp_2
  '

fit4 <- sem(model4, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

fit4_sq <- sem(model4_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")


#model with no latent variables
model5 <- '
    coral ~ tmean + sal + pp + fr30m + reg1 + reg2 + reg3 + reg4
    reef_fraction ~ fr30m + coral + tmean + sal + reg1 + reg2 + reg3 + reg4
    fish ~  reef_fraction + coral + tmean + sal + pp + reg1 + reg2 + reg3 + reg4
'

model5_sq <- '
    coral ~ tmean_2 + sal_2 + pp_2 + fr30m + reg1 + reg2 + reg3 + reg4
    reef_fraction ~ fr30m + coral + tmean_2 + sal_2 + reg1 + reg2 + reg3 + reg4
    fish ~  reef_fraction + coral + tmean_2 + sal_2 + pp_2 + reg1 + reg2 + reg3 + reg4
'

fit5 <- sem(model5, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

fit5_sq <- sem(model5_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

###no latent variables no regions
model6 <- '
    coral ~ tmean + sal + pp + fr30m
    reef_fraction ~ fr30m + coral + tmean + sal
    fish ~  reef_fraction + coral + tmean + sal + pp
'
model6_sq <- '
    coral ~ tmean_2 + sal_2 + pp_2 + fr30m
    reef_fraction ~ fr30m + coral + tmean_2 + sal_2
    fish ~  reef_fraction + coral + tmean_2 + sal_2 + pp_2
'

fit6 <- sem(model6, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

fit6_sq <- sem(model6_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

###no latent variables no regions no reef_fraction
model7 <- '
    coral ~ tmean + sal + pp + fr30m
    fish ~  coral + tmean + sal + pp
'
model7_sq <- '
    coral ~ tmean_2 + sal_2 + pp_2 + fr30m
    fish ~  coral + tmean_2 + sal_2 + pp_2
'

fit7 <- sem(model7, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

fit7_sq <- sem(model7_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

model8 <- '
      coral ~ isolation + tmean + trange + sal + pp + fr30m + reg1 + reg2 + reg3 + reg4
      reef_fraction ~ isolation + fr30m + coral + tmean + trange + sal + reg1 + reg2 + reg3 + reg4
      fish ~  isolation + reef_fraction + coral + tmean + trange + sal + pp + reg1 + reg2 + reg3 + reg4
     
      coral ~~ reef_fraction
  '
model8_sq <- '
      coral ~ isolation + tmean_2 + trange_2 + sal_2 + pp_2 + fr30m + reg1 + reg2 + reg3 + reg4
      reef_fraction ~ isolation + fr30m + coral + tmean_2 + trange_2 + sal_2 + reg1 + reg2 + reg3 + reg4
      fish ~  isolation + reef_fraction + coral + tmean_2 + trange_2 + sal_2 + pp_2 + reg1 + reg2 + reg3 + reg4
     
      coral ~~ reef_fraction
  '

fit8 <- sem(model8, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE)

fit8_sq <- sem(model8_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE)

model9 <- '
      environment =~ tmean + trange + sal + pp + fr30m + reef_fraction
      history =~ isolation + tmean
    
      coral ~ environment + history + reg1 + reg2 + reg3 + reg4
      fish ~ coral + environment + history + reg1 + reg2 + reg3 + reg4
  '
model9_sq <- '
      environment =~ tmean_2 + trange_2 + sal_2 + pp_2 + fr30m + reef_fraction
      history =~ isolation + tmean_2
    
      coral ~ environment + history + reg1 + reg2 + reg3 + reg4
      fish ~ coral + environment + history + reg1 + reg2 + reg3 + reg4
  '

fit9 <- sem(model9, fixed.x=F, 
      data=b,check.gradient=F,meanstructure = TRUE,
      estimator = "ML")

fit9_sq <- sem(model9_sq, fixed.x=F, 
            data=b,check.gradient=F,meanstructure = TRUE,
            estimator = "ML")

  
compareFit(fit0,fit0_sq,fit1,fit1_sq,fit2,fit2_sq,fit3,fit3_sq,fit4,fit4_sq,
           fit5,fit5_sq,fit6,fit6_sq,fit7,fit7_sq,fit8,fit8_sq,fit9,fit9_sq,nested=F)
  
  
#summary(fit1,standardized = T, fit.measures = T, rsq = T)
#lavaanPlot(model=fit1,coefs = TRUE,covs=TRUE,stand = TRUE, sig = 1)
  
#dashed line indicates fixed parameter estimates
  
  
 library(semPlot)
# 
# lay<-matrix(c(12,1,
#               1,1,
#               0,-3,#lat  8
#               1,5,#t     1
#               4,5,#trange  2
#               -3,-3,#iso  7
#               7,5,#sal    3
#               10,5,#pp     4
#               -2,4,#reef  5
#               3,-3,#reg1 9
#               6,-3,#reg2 10
#               9,-3,#reg3 11
#               12,-3,#reg4 12
#               -4,3),#fr30 6
#             ncol=2,byrow=TRUE)
# 
# lbs<-c('FISH','CORAL','|LAT|','Tm','Tr','ISO','SAL','PP','REEF',
#   'WA','WIO','CIP','CP','30m')
# semPaths(fit0, "std",
#              residuals = FALSE,
#          intercepts=FALSE,
#          layout=lay,#'tree',
#              style = "lisrel",
#              exoCov = FALSE,
#              thresholds=FALSE,nCharNodes = 5,
#              optimizeLatRes = TRUE,
#               whatLabels = "std",curve = 100,
#              sizeMan=4,#shapeMan='circle',
#              nDigits=2,layoutSplit=TRUE,
#              curveAdjacent='<->',edge.label.cex=0.6,
#              edge.label.color='black',
#              edge.width= 1.5,
#              #edge.color='blue',
#          posCol=c("blue","red"),
#          nodeLabels=lbs,
#              normalize=FALSE,
#              label.color='black',
#              label.cex=0.8,
#              minimum=0.0,
#              width=8,
#              height=5,
#              normalize=FALSE,
#              trans=FALSE,
#              colFactor=0.5,
#              fade=T,
#              rescale=TRUE,
#              mar=c(0,6,0,6),
#          filetype='pdf',
#          filename = 'prova_sem'
#          )
# 


sc<-1
for (f in c(fit0_sq,fit1_sq,fit2_sq,fit3_sq,fit4_sq,fit5_sq,fit6_sq,fit7_sq,fit8_sq,fit9_sq)){
  
  
semPaths(f, "std", 
           residuals = FALSE, intercepts=FALSE,layout='tree',
           style = "lisrel",
           exoCov = FALSE,
           thresholds=FALSE,nCharNodes = 5,
           optimizeLatRes = TRUE,
           whatLabels = "std",curve = 1,
           sizeMan=4,
           nDigits=2,layoutSplit=FALSE,
           curveAdjacent='<->',edge.label.cex=0.6,
           edge.label.color='black',
           edge.width= 1.5,
         posCol=c("blue","red"),
           normalize=FALSE,
           label.color='black',
           label.cex=0.8,
           minimum=0.0,
           normalize=TRUE,
           trans=FALSE,
           colFactor=0.5,
           fade=TRUE,
           rescale=TRUE,filetype='pdf',
           filename = paste0(sc,'_sq'),
           mar=c(2,2,2,2)
  )
  sc<-sc+1
  
}


####
####check spatial autocorrelation
# fish_gls<-gls(fish ~ coral+ abs_lat + isolation + sal + pp 
#               + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4,
#               data=b)
# 
# 
# cor_gls<-gls(coral ~ abs_lat + tmean + trange + isolation + sal 
#              + pp + reef_fraction + fr30m + reg1 + reg2 + reg3 + reg4,
#              data=b)
# 
# 
# p.sem<-psem(cor_gls,fish_gls,data = b)
# #mod<-summary(p.sem)
# coef<-coefs(p.sem,standardize = "scale")#[,c(1,2,8)]
# 
# vario_fish <- Variogram(fish_gls, form = ~lon + lat, resType = "pearson",maxDist=500)
# vario_cor <- Variogram(cor_gls, form = ~lon + lat, resType = "pearson",maxDist=500)
# pdf('variograms.pdf')
# plot(vario_fish, smooth = F,main='fish')
# plot(vario_cor, smooth = F,main='coral')
# dev.off()


spaceCor<-corLin(form =~ lon+lat, nugget=F)
cor_gls_sp<-gls(fish ~ coral+ abs_lat + tmean + trange + isolation 
             + sal + pp + reef_fraction  
             + reg1 + reg2 + reg3 + reg4,
             correlation = spaceCor,
             data=b)

spaceCor<-corLin(form =~ lon+lat, nugget=T)
fish_gls_sp<-gls(coral ~ abs_lat + tmean + trange 
              + isolation + sal + pp + reef_fraction 
              + fr30m + reg1 + reg2 + reg3 + reg4,
              correlation = spaceCor,
              data=b)

sp.sem<-psem(cor_gls_sp,fish_gls_sp, data = b)
#mod<-summary(p.sem)
coef<-coefs(sp.sem,standardize = "scale")#[,c(1,2,8)]

fit<-fit0
pdf('fish_coral_sem_summary.pdf',useDingbats = F)
par(mfrow=c(2,2)) 
summary(fit, standardized=TRUE,fit.measures=T)
cv <- fitted(fit)$cov
cv_names<-row.names(cv)
i <- which(cv_names=='fish')
j <- which (cv_names!='fish')
coef <- solve(cv[j,j],cv[j,i])
cnames<-rownames(as.matrix(coef))
predictions <- as.matrix(b[,cnames])%*% coef
cor(b[,1],predictions)
summary(lm(b[,1]~predictions))
intercept<-lm(b[,1]~predictions)[[1]][1]
predictions<-predictions+intercept

# mod<-lm(b[,1]~predictions)
# ####check spatial autocorrelation in model residuals
# residuals<-resid(mod) 
# dists <- as.matrix(dist(cbind(a$lon, a$lat)))
# dists.inv <- 1/dists
# diag(dists.inv) <- 0
# Moran.I(residuals, dists.inv,scaled=T)

r2<-round(cor(b[,1],predictions)**2,2)
plot(b[,1]*2794+1,predictions*2794+1,
     xlab='observed fish diversity',ylab='modelled fish diversity',
     main = paste('R2 =',r2),
     col=plasma(1,alpha=0.5)[1],
     las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.8)


abline(0,1)

new_data<-b
new_data[,'coral']<-0
new_data[,'reef_fraction']<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
afl = 1-predictions_no_coral/predictions
afl[which(predictions==0)] = 0
mean(afl)
sd(afl)

fd0<-2794*new_data$fish+1
ppp<-sum(fd0*afl)/sum(fd0)


res<-data.frame(cbind(b,afl))
colnames(res)[dim(res)[2]]<-'pred_dep'

rho<-round(as.numeric(cor.test(res$lit_dep,res$pred_dep,method='spearman')[[4]]),2)
plot(res$lit_dep,res$pred_dep,
     xlab='literature dependency',ylab='modelled dependency',
     main=paste("Spearman's rho =",rho),col=plasma(1,alpha=0.5)[1],
     las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.8)

abline(0,1)

high<-which(afl>b[,1])
tot_dep<-b[,1]
tot_dep[high]<-afl[high]
mean(tot_dep) ###taking the max dependence (either from literature or from the model)

plot(b[,1]*2794+1,predictions*2794+1,
     xlab='observed fish diversity',ylab='modelled fish diversity',
     las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.8,col=plasma(3,alpha=0.5)[1])


points(b[,1]*2794+1,predictions_no_coral*2794+1,col=plasma(3,alpha=0.5)[2],pch=16,cex=0.8)
abline(0,1)
legend(0,1,c('with corals','without corals'),col=plasma(3)[1:2],pch=15)

#explore different levels of coral loss and reef area loss
x<-c()
y<-c()
z<-c()
for (cor_loss in seq(0,1,0.1)){
  for (cor_fract_loss in seq(0,1,0.1)){
    new_data<-b
    new_data[,2]<-new_data[,2]*(1-cor_loss)
    new_data[,5]<-new_data[,5]*(1-cor_fract_loss) 
    predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
    predictions_no_coral<-predictions_no_coral+intercept
    predictions[predictions<0]<-0
    predictions_no_coral[predictions_no_coral<0]<-0
    afl = 1-predictions_no_coral/predictions
    afl[which(predictions==0)] = 0
    afl <- mean(afl)
    x<-c(x,cor_loss)
    y<-c(y,cor_fract_loss)
    z<-c(z,afl)
  }}


image.plot(interp(x,y,z),col=plasma(20),
           las=1,cex.axis=1.2,cex.lab=1.2,
           main='mean fish loss',xlab='loss of coral diversity',
           ylab='loss of reef fraction')


dev.off()

###coral diversity theoretical scenarios
fd0<-2794*b$fish+1
x1<-c()
y1<-c()
y1_raw<-c()
n<-dim(b)[1]

for (rep in 1:100){
new_data<-b
to_del<-sample(1:n)
for (cor_loss in 1:n){
  new_data[to_del[cor_loss],2]<-0
  new_data[to_del[cor_loss],5]<-0
  predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
  predictions_no_coral<-predictions_no_coral+intercept
  predictions[predictions<0]<-0
  predictions_no_coral[predictions_no_coral<0]<-0
  afl = 1-predictions_no_coral/predictions
  afl[which(predictions==0)] = 0
  x1<-c(x1,cor_loss)
  glob_loss<-fd0*afl
  y1<-c(y1,100*(1-mean((fd0-glob_loss)/fd0)))
  y1_raw<-c(y1_raw,mean(glob_loss))
}
print(rep)
}

####best-worst case
###get the effect per locality
new_data<-b
new_data[,2]<-0
new_data[,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
afl = 1-predictions_no_coral/predictions
afl[which(predictions==0)] = 0
glob_loss<-fd0*afl
glob_loss_p<-afl

best<-sort(glob_loss,index.return=T)$ix
worst<-sort(glob_loss,index.return=T,decreasing = TRUE)$ix
best_p<-sort(glob_loss_p,index.return=T)$ix
worst_p<-sort(glob_loss_p,index.return=T,decreasing = TRUE)$ix


#worst - reefs obliterated from most coral rich to poorest
y2_raw<-c()
n<-dim(b)[1]
new_data<-b
for (to_del in worst){
new_data[to_del,2]<-0
new_data[to_del,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
afl = 1-predictions_no_coral/predictions
afl[which(predictions==0)] = 0
glob_loss<-fd0*afl
y2_raw<-c(y2_raw,mean(glob_loss))
}


y2<-c()
new_data<-b
for (to_del in worst_p){
new_data[to_del,2]<-0
new_data[to_del,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
afl = 1-predictions_no_coral/predictions
afl[which(predictions==0)] = 0
glob_loss<-fd0*afl
y2<-c(y2,100*(1-mean((fd0-glob_loss)/fd0)))
}

#best - reefs obliterated from coral poorest to richest 
y3_raw<-c()
n<-dim(b)[1]
new_data<-b
for (to_del in best){
new_data[to_del,2]<-0
new_data[to_del,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
afl = 1-predictions_no_coral/predictions
afl[which(predictions==0)] = 0
glob_loss<-fd0*afl
y3_raw<-c(y3_raw,mean(glob_loss))
}


y3<-c()
new_data<-b
for (to_del in best_p){
new_data[to_del,2]<-0
new_data[to_del,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
afl = 1-predictions_no_coral/predictions
afl[which(predictions==0)] = 0
glob_loss<-fd0*afl
y3<-c(y3,100*(1-mean((fd0-glob_loss)/fd0)))
}



pdf('theoretical_loss.pdf',height=5,width=10)
par(mfrow=c(1,2))
x1<-100*x1/max(x1)


y1_<-aggregate(y1~x1,FUN='mean')[,2]
y1_min<-aggregate(y1~x1,FUN='min')[,2]
y1_max<-aggregate(y1~x1,FUN='max')[,2]

x<-x1[1:n]
plot(x,y2,type='l',las=1,cex.axis=1.2,cex.lab=1.2,
  xlab='reef loss (%)',ylab = 'mean fish loss (%)',
  main='A',lwd=1.5,col='red')

lines(x,y3,lwd=1.5,col='blue')
polygon(c(rev(x), x), c(rev(y1_min),y1_max), 
     col = 'lightgrey', border = NA)
lines(x,y1_,lwd=1.5)


y1_raw_<-aggregate(y1_raw~x1,FUN='mean')[,2]
#y1_raw_sd<-aggregate(y1_raw~x1,FUN='sd')[,2]
y1_raw_max<-aggregate(y1_raw~x1,FUN='max')[,2]
y1_raw_min<-aggregate(y1_raw~x1,FUN='min')[,2]


plot(x,y2_raw, type='l',las=1,cex.axis=1.2,cex.lab=1.2,
        xlab='reef loss (%)',
   ylab = 'mean fish loss (no. spp.)',
   main='B',lwd=1.5,col='red')
lines(x,y3_raw,lwd=1.5,col='blue')

polygon(c(rev(x), x), c(rev(y1_raw_min),y1_raw_max), 
      col = 'lightgrey', border = NA)
lines(x,y1_raw_,lwd=1.5)


legend(60,200,bty='n',c('worst','random','best'),col=c('red','black','blue'),pch=15,pt.cex=1.5,title='reef loss scenario')

dev.off()

#####climatic scenarios
#non cumulative
  
fd0<-2794*b$fish+1
scen_fff <- c('ssp2_bl_risk.csv','ssp3_bl_risk.csv','ssp5_bl_risk.csv')
x<-seq(2015,2100,length.out=n)
res<-c()
for (scen in 1:3){
ssp<-read.csv(scen_fff[scen],header=F)
ssp<-ssp[,ok]
y1<-c()
n<-1032
for (mo in 1:n){
  new_data<-b
  anom<-which(ssp[mo,]>1)
  new_data[anom,2]<-0
  new_data[anom,5]<-0
  predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
  predictions_no_coral<-predictions_no_coral+intercept
  predictions[predictions<0]<-0
  predictions_no_coral[predictions_no_coral<0]<-0
  afl = 1-predictions_no_coral/predictions
  afl[which(predictions==0)] = 0
  glob_loss<-fd0*afl
  y1<-c(y1,100*(1-mean((fd0-glob_loss)/fd0)))
}
res<-cbind(res,y1)
}

pdf('climatic_non_cum.pdf',width=5,height=5)
plot(0,0,type='l',xlim=c(2015,2100),ylim=c(0,30),las=1,cex.axis=1.2,cex.lab=1.2,
   xlab='year',ylab='fish loss (%, monthly estimate)',col='white')

for (s in 1:3){
y1<-res[,s]
lines(x,y1,col=viridis(4,alpha=0.4)[s])}

for (s in 1:3){
y1<-res[,s]
y_<-aggregate(y1~round(x/10),FUN='mean')  
lines(y_[,1]*10,y_[,2],col=viridis(4)[s],lwd=2)  
}

legend(2015,30,c('SSP5-8.5','SSP3-7.0','SSP2-4.5'),rev(viridis(4)[1:3]),pt.cex=2)
dev.off()

  
  #####climatic_cumulative
pdf('climatic_cum.pdf',width=5,height=5)
  plot(0,0,type='l',xlim=c(2015,2100),ylim=c(0,45),las=1,cex.axis=1.2,cex.lab=1.2,
       xlab='year',ylab='fish loss (%, cumulative)')
  scen_fff <- c('ssp2_bl_risk.csv','ssp3_bl_risk.csv','ssp5_bl_risk.csv')
  for (scen in 1:3){
    ssp<-read.csv(scen_fff[scen],header=F)
    ssp<-ssp[,ok]
    y1<-c()
    n<-1032
    new_data<-b
    for (mo in 1:n){
      anom<-which(ssp[mo,]>1)
      new_data[anom,2]<-0
      new_data[anom,5]<-0
      predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
      predictions_no_coral<-predictions_no_coral+intercept
      predictions[predictions<0]<-0
      predictions_no_coral[predictions_no_coral<0]<-0
      afl = 1-predictions_no_coral/predictions
      afl[which(predictions==0)] = 0
      glob_loss<-2794*new_data$fish*afl+1
      y1<-c(y1,100*(1-mean((fd0-glob_loss)/fd0)))}
    
    lines(x,y1,col=viridis(4,alpha=1)[scen],lwd=2)
  }
  legend(2070,15,c('SSP5-8.5','SSP3-7.0','SSP2-4.5'),rev(viridis(4)[1:3]),pt.cex=2)
  dev.off()
  
  
  ###check coral loss through time 
  scen_fff <- c('ssp2_bl_risk.csv','ssp3_bl_risk.csv','ssp5_bl_risk.csv')
  res<-c()
  for (scen in 1:3){
    ssp<-read.csv(scen_fff[scen],header=F)
    ssp<-ssp[,ok]
    y1<-c()
    n<-1032
    for (mo in 1:n){
      new_data<-b
      anom<-which(ssp[mo,]>1)
      new_data[anom,2]<-0
      new_data[anom,5]<-0
      y1<-c(y1,100*sum(new_data[,2]==0)/dim(b)[1])
    }
    res<-cbind(res,y1)
  }
  
  
  pdf('coral_loss_anom.pdf',width=5,height=5)
  plot(0,0,type='l',xlim=c(2015,2100),ylim=c(0,65),las=1,cex.axis=1.2,cex.lab=1.2,
       xlab='year',ylab='reefs at risk (%, monthly)')
  
  for (s in 1:3){
    y1<-res[,s]
    lines(x,y1,col=viridis(4,alpha=0.4)[s])}
  
  for (s in 1:3){
    y1<-res[,s]
    y_<-aggregate(y1~round(x/10),FUN='mean')  
    lines(y_[,1]*10,y_[,2],col=viridis(4)[s],lwd=2)  
  }
  
  legend(2015,65,c('SSP5-8.5','SSP3-7.0','SSP2-4.5'),rev(viridis(4)[1:3]),pt.cex=2)
  dev.off()
  
  #################cumulative reefs at risk
  pdf('coral_loss_anom_cum.pdf',width=5,height=5)
  plot(0,0,type='l',xlim=c(2015,2100),ylim=c(0,100),las=1,cex.axis=1.2,cex.lab=1.2,
       xlab='year',ylab='reefs at risk (%, cumulative)')
  
  scen_fff <- c('ssp2_bl_risk.csv','ssp3_bl_risk.csv','ssp5_bl_risk.csv')
  for (scen in 1:3){
    ssp<-read.csv(scen_fff[scen],header=F)
    ssp<-ssp[,ok]
    y1<-c()
    n<-1032
    new_data<-b
    for (mo in 1:n){
      anom<-which(ssp[mo,]>1)
      new_data[anom,2]<-0
      new_data[anom,5]<-0
      y1<-c(y1,100*sum(new_data[,2]==0)/dim(b)[1])
    }
    
    lines(seq(2015,2100,length.out=n),y1,col=viridis(4)[scen],lwd=2)
  }
  
  legend(2015,100,c('SSP5-8.5','SSP3-7.0','SSP2-4.5'),rev(viridis(4)[1:3]),pt.cex=2)
  
  dev.off()
  


####raw div with/without corals
par(mfrow=c(1,1))
new_data<-b
new_data[,2]<-0
new_data[,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
afl = 1-predictions_no_coral/predictions
afl[which(predictions==0)] = 0
glob_loss<-new_data$fish-new_data$fish*afl

pdf('figS8.pdf',height=5,width=5)
plot(density(log(new_data$fish*2794+1)),col=plasma(3)[1],lwd=1.5,main='',
     las=1,cex.axis=1.2,cex.lab=1.2,xlab = 'log(fish spp)')
polygon(density(log(new_data$fish*2794+1)),col=plasma(3,alpha=0.2)[1],border=NA)
lines(density(log(glob_loss*2794+1)),col=plasma(3)[2],lwd=1.5)
polygon(density(log(glob_loss*2794+1)),col=plasma(3,alpha=0.2)[2],border=NA)
legend(0,0.8,c('with corals','without corals'),col=plasma(3)[1:2],pch=15)
dev.off()


# ##explore loss vs div
# 
# bb<-aggregate(res$pred_dep~round(log(res$fish*2794+1),1),FUN='mean')
# bb_sd<-aggregate(res$pred_dep~round(log(res$fish*2794+1),1),FUN='sd')[,2]
# rho<-round(as.numeric(cor.test(res$pred_dep,log(res$fish*2794+1),method='spearman')[[4]]),2)
# 
# pdf('dependency_vs_fish_diversity.pdf',useDingbats = F)
# plot(bb[,1],bb[,2],type='n',
#      ylim=c(0,1),main=paste("Spearman's rho=",rho),
#      xlab='fish diversity',ylab='coral dependency',
#      las=1, cex.axis=1.2,cex.lab=1.2)
# 
# lines(bb[,1],bb[,2],col=plasma(1)[1],lwd=1.5)
# points(bb[,1],bb[,2],col=plasma(1)[1],pch=16,cex=0.5)
# polygon(c(rev(bb[,1]),bb[,1]),c(rev(bb[,2]-bb_sd),bb[,2]+bb_sd),
#         col=plasma(1,alpha=0.3)[1],border=NA)
# 
# abline(lm(res$pred_dep~res$fish),col='darkred',lty=2,lwd=1.5)
# abline(h=median(res$pred_dep),lwd=1,lty=2)
# dev.off()


# levels(res$region)
# [1] "central Indo-Pacific"    "central Pacific"        
# [3] "Tropical Eastern Region" "western Atlantic"       
# [5] "western Indian" 
levels(res$region)<-c('CIP','CP','TEP','WA','WI')
#boxplot(res$pred_dep~res$region)

par(mfrow=c(1,1))
pdf('barplot.pdf',useDingbats=F,width=6,height=4)
t<-c()
count<-0
for (col in c(5,6,7)){
  row = c()
  c_id = 1
  for (reg in c('central Pacific','western Indian','central Indo-Pacific','western Atlantic','Tropical Eastern Region')){
    val<-as.numeric(mean(a[a$Region==reg,col]))
    if (count>0){val<-val-as.numeric(t[count,c_id])}
    row<-c(row,val)
    c_id<-c_id+1	
  }
  t<-rbind(t,row)
  count<-count+1
}

colnames(t)<-c('CP','WIP','CIP','WA','TEP')
rownames(t)<-c('obligate corallivorous','facultative corallivorous','coral associated')
bp<-barplot(t, col=colors()[c(50,89,12,6)] , border="white", space=0.1,
            las=1,cex.axis=1.2,cex.lab=1.2, xlab="",ylim=c(0,1))
pred = c()
pred_sd = c()
pred_div<-which(colnames(res)=='pred_dep')
for (reg in c('CP','WI','CIP','WA','TEP')){
  val<-as.numeric(mean(res[res$region==reg,pred_div],na.rm=T))
  val_sd<-as.numeric(sd(res[res$region==reg,pred_div],na.rm=T))
  pred<-c(pred,val)
  pred_sd<-c(pred_sd,val_sd)		
}

points(bp,pred,pch=16,cex=2)
arrows(bp, pred-pred_sd, bp, pred+pred_sd, length=0.05, angle=90, code=3)
legend(bp[3]+0.3,0.8,c(rownames(t),'predicted'),
       col=c(colors()[c(50,89,12)],'black'),
       pch=c(15,15,15,16),pt.cex=c(2,2,2,1),
       cex = 0.8)
dev.off()
  
###fig s1
r2 = round(cor(a$cor_diversity,a$fish_sp_div),2)
pdf('figS1.pdf',height=5,width = 5)
plot(a$cor_diversity,a$fish_sp_div,
     xlab='coral genera',ylab='fish species',
     main = paste('R2 =',r2),
     col=plasma(1,alpha=0.5)[1],
     las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.8)
    abline(lm(a$fish_sp_div~a$cor_diversity))  

dev.off()     
     
     
     
###data for maps
data4maps<-data.frame(cbind(a$lat,a$lon,a$cor_diversity,a$fish_sp_div,
  res$pred_dep*a$fish_sp_div,a$fish_sp_div-(res$pred_dep*a$fish_sp_div),
  res$pred_dep,res$lit_dep))

colnames(data4maps)<-c('lat','lon',
                       'coral','fish',
                       'dependent_fish_model','fish_loss_model',
                       'model_dependency','lit_dependency')
      
write.table(data4maps,file='data4maps.csv',quote=F,sep=',',row.names=F)
      

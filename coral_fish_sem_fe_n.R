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
library(semPlot)

#adj_coeff
ac<-3

prefix<-'fe_n_'
norm_to_raw<-function(norm_l_,raw_l,ac=3){
  norm_l<-norm_l_**(1/ac)
  return(norm_l*max(raw_l)-norm_l*min(raw_l)+min(raw_l))}


a<-read.csv('coralfishglobal.csv',header=T)
a<-a[a$cor_diversity>0,]
a<-a[a$fish_sp_div>0,]
a<-a[a$reef_fraction>0,]
a<-a[a$phylo_div>0,]
a<-a[a$fe_n>0,]
pdf('fish_richness_vs_functional_entities.pdf',width=5,height=5)
plot(a$fish_sp_div,a$fe_n,xlab='fish species richness',ylab='unique functional entities',
     las=1,cex.lab=1.2,cex.axis=1.2)
dev.off()





reg<-unique(a$Region)
b<-data.frame(cbind(a$fe_n**ac,
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




fd0<-norm_to_raw(b$fish,a$fe_n)

pdf('fe_n_observed_vs_model_all.pdf')
names<-paste0('fit',sort(c(0:9,0:9)),c('','_sq'))
sc<-1
for (fit in c(fit0,fit0_sq,fit1,fit1_sq,fit2,fit2_sq,
              fit3,fit3_sq,fit4,fit4_sq,
           fit5,fit5_sq,fit6,fit6_sq,fit7,
           fit7_sq,fit8,fit8_sq,fit9,fit9_sq)){

cv <- fitted(fit)$cov
cv_names<-row.names(cv)
i <- which(cv_names=='fish')
j <- which (cv_names!='fish')
coef <- solve(cv[j,j],cv[j,i])
cnames<-rownames(as.matrix(coef))
predictions <- as.matrix(b[,cnames])%*% coef
intercept<-lm(b[,1]~predictions)[[1]][1]
predictions<-predictions+intercept
raw_predictions<-norm_to_raw(predictions,a$fe_n)
new_data<-b
new_data[,'coral']<-0
new_data[,'reef_fraction']<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
raw_predictions<-norm_to_raw(predictions,a$fe_n)
raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
afl = 1-raw_predictions_no_coral/raw_predictions
r2<-round(cor(a$fe_n,raw_predictions)**2,2)

max_xy<-max(c(a$fe_n,raw_predictions))
plot(a$fe_n,raw_predictions,
     xlab='observed fish diversity',ylab='modelled fish diversity',
     main = paste(names[sc],'R2 =',r2),
     xlim = c(0,max_xy),
     ylim = c(0,max_xy),
     col=plasma(1,alpha=0.5)[1],
     las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.8)
abline(0,1)
sc<-sc+1}

dev.off()




fits<-compareFit(fit0,fit0_sq,fit1,fit1_sq,fit2,fit2_sq,
                 fit3,fit3_sq,fit4,fit4_sq,
                 fit5,fit5_sq,fit6,fit6_sq,fit7,
                 fit7_sq,fit8,fit8_sq,fit9,fit9_sq,nested=F)


fit<-fit0

semPaths(fit, "std", 
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
         filename = paste0(prefix,'best_model'),
         mar=c(2,2,2,2)
)



cv <- fitted(fit)$cov
cv_names<-row.names(cv)
i <- which(cv_names=='fish')
j <- which (cv_names!='fish')
coef <- solve(cv[j,j],cv[j,i])
cnames<-rownames(as.matrix(coef))
predictions <- as.matrix(b[,cnames])%*% coef
intercept<-lm(b[,1]~predictions)[[1]][1]
predictions<-predictions+intercept
new_data<-b
new_data[,'coral']<-0
new_data[,'reef_fraction']<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
raw_predictions<-norm_to_raw(predictions,a$fe_n)
raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
afl = 1-raw_predictions_no_coral/raw_predictions
r2<-round(cor(a$fe_n,raw_predictions)**2,2)

pdf('fish_spp_vs_predicted_fe_n.pdf',width=5,height=5)
plot(a$fish_sp_div,raw_predictions,xlab='fish diversity',
     ylab='SEM predicted functional entities',las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16,col=viridis(3,alpha=1)[1]
     )
dev.off()


sink(paste0(prefix,'model_fits.txt'))
fits
paste('R2',r2)
paste('mean dep',mean(afl))
paste('sd dep',sd(afl))
sink()



###data for maps
data4maps<-data.frame(cbind(a$lat,a$lon,a$cor_diversity,a$fe_n,
                            afl*a$fe_n,a$fe_n-afl*a$fe_n,afl))

colnames(data4maps)<-c('lat','lon',
                       'coral','fish',
                       'dependent_fish_model','fish_loss_model',
                       'model_dependency')

write.table(data4maps,file=paste0(prefix,'data4maps.csv'),quote=F,sep=',',row.names=F)



###coral diversity theoretical scenarios
nm<-read.csv('fen_null_model.csv',header=T)

nm_<-aggregate(nm$fen~nm$reef_n,FUN='mean')
nm<-cbind(nm_,nm[1:nrow(nm_),3:4])
lat_lon<-apply(cbind(b$lat,b$lon),1,paste,collapse="-")
loc_ok<-which(apply(nm[,c(3,4)],1,paste,collapse= "-") %in% lat_lon)
nm<-nm[loc_ok,] #this includes the expected functional diversity after total coral loss


x1<-c()
y1<-c()
y1_raw<-c()
y1_nm<-c()
y1_nm_raw<-c()

n<-dim(b)[1]
for (rep in 1:100){
new_data<-b
fd0_nm<-fd0
to_del<-sample(1:n)
for (cor_loss in 1:n){
  new_data[to_del[cor_loss],2]<-0
  new_data[to_del[cor_loss],5]<-0
  predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
  predictions_no_coral<-predictions_no_coral+intercept
  predictions[predictions<0]<-0
  predictions_no_coral[predictions_no_coral<0]<-0
  raw_predictions<-norm_to_raw(predictions,a$fe_n)
  raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
  afl = 1-raw_predictions_no_coral/raw_predictions
  x1<-c(x1,cor_loss)
  glob_loss<-fd0*afl #lost functional div
  y1<-c(y1,100*(1-mean((fd0-glob_loss)/fd0)))
  y1_raw<-c(y1_raw,mean(glob_loss))
  #do the same for the null model
  fd0_nm[to_del[cor_loss]]<-nm[to_del[cor_loss],2]
  y1_nm<-c(y1_nm,100*(mean((fd0-fd0_nm)/fd0)))
  y1_nm_raw<-c(y1_nm_raw,mean(fd0_nm))
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
raw_predictions<-norm_to_raw(predictions,a$fe_n)
raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
afl = 1-raw_predictions_no_coral/raw_predictions
glob_loss<-fd0*afl
glob_loss_p<-afl
best<-sort(glob_loss,index.return=T)$ix
worst<-sort(glob_loss,index.return=T,decreasing = TRUE)$ix
best_p<-sort(glob_loss_p,index.return=T)$ix
worst_p<-sort(glob_loss_p,index.return=T,decreasing = TRUE)$ix


#worst - reefs obliterated from most coral rich to poorest
n<-dim(b)[1]

new_data<-b

y2<-c()
y2_raw<-c()
y2_nm<-c()
y2_nm_raw<-c()

fd0_nm<-fd0
for (to_del in worst){
new_data[to_del,2]<-0
new_data[to_del,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
raw_predictions<-norm_to_raw(predictions,a$fe_n)
raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
afl = 1-raw_predictions_no_coral/raw_predictions
glob_loss<-fd0*afl
y2_raw<-c(y2_raw,mean(glob_loss))
y2<-c(y2,100*(1-mean((fd0-glob_loss)/fd0)))

fd0_nm[to_del]<-nm[to_del,2]
y2_nm<-c(y2_nm,100*(mean((fd0-fd0_nm)/fd0)))
y2_nm_raw<-c(y2_nm_raw,mean(fd0_nm))

}


#best - reefs obliterated from coral poorest to richest 
y3_raw<-c()
y3<-c()
y3_nm<-c()
y3_nm_raw<-c()
fd0_nm<-fd0
n<-dim(b)[1]
new_data<-b
for (to_del in best){
new_data[to_del,2]<-0
new_data[to_del,5]<-0
predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
predictions_no_coral<-predictions_no_coral+intercept
predictions[predictions<0]<-0
predictions_no_coral[predictions_no_coral<0]<-0
raw_predictions<-norm_to_raw(predictions,a$fe_n)
raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
afl = 1-raw_predictions_no_coral/raw_predictions
glob_loss<-fd0*afl
y3_raw<-c(y3_raw,mean(glob_loss))
y3<-c(y3,100*(1-mean((fd0-glob_loss)/fd0)))
fd0_nm[to_del]<-nm[to_del,2]
y3_nm<-c(y3_nm,100*(mean((fd0-fd0_nm)/fd0)))
y3_nm_raw<-c(y3_nm_raw,mean(fd0_nm))
}

y1_<-aggregate(y1~x1,FUN='mean')[,2]
y1_raw_<-aggregate(y1_raw~x1,FUN='mean')[,2]
y1_nm_<-aggregate(y1_nm~x1,FUN='mean')[,2]
y1_nm_raw_<-aggregate(y1_nm_raw~x1,FUN='mean')[,2]

n<-dim(b)[1]
x<-seq(0,100,length.out = n)

th_tab<-cbind(x,y1_,y1_raw_,y2,y2_raw,y3,y3_raw,y1_nm_,y1_nm_raw_,y2_nm,y2_nm_raw,y3_nm,y3_nm_raw)
colnames(th_tab)<-c('coral_loss','random','random_raw','worst','worst_raw','best','best_raw',
                    'random_nm','random_nm_raw','worst_nm','worst_nm_raw','best_nm','best_nm_raw')

write.table(th_tab,paste0(prefix,'th_tab.csv'),col.names=T,row.names=F,quote=F,sep=',')

#####climatic scenarios
#non cumulative
n<-1032  


scen_fff <- c('ssp2_bl_risk.csv','ssp3_bl_risk.csv','ssp5_bl_risk.csv')
x<-seq(2015,2100,length.out=n)
res<-c()
for (scen in 1:3){
  ssp<-read.csv(scen_fff[scen],header=F)
  ssp<-ssp[,ok]
  y1<-c()
  y1_nm<-c()
  for (mo in 1:n){
    fd0_nm<-fd0
    new_data<-b
    anom<-which(ssp[mo,]>1)
    new_data[anom,2]<-0
    new_data[anom,5]<-0
    predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
    predictions_no_coral<-predictions_no_coral+intercept
    predictions[predictions<0]<-0
    predictions_no_coral[predictions_no_coral<0]<-0
    raw_predictions<-norm_to_raw(predictions,a$fe_n)
    raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
    afl = 1-raw_predictions_no_coral/raw_predictions
    glob_loss<-fd0*afl
    y1<-c(y1,100*(1-mean((fd0-glob_loss)/fd0)))
    fd0_nm[anom]<-nm[anom,2]
    y1_nm<-c(y1_nm,100*(mean((fd0-fd0_nm)/fd0)))
  }
  res<-cbind(res,y1,y1_nm)
}

res_clim<-cbind(x,res)


#####climatic_cumulative
scen_fff <- c('ssp2_bl_risk.csv','ssp3_bl_risk.csv','ssp5_bl_risk.csv')
for (scen in 1:3){
    ssp<-read.csv(scen_fff[scen],header=F)
    ssp<-ssp[,ok]
    y1<-c()
    y1_nm<-c()
    new_data<-b
    fd0_nm<-fd0
    for (mo in 1:n){
      anom<-which(ssp[mo,]>1)
      new_data[anom,2]<-0
      new_data[anom,5]<-0
      predictions_no_coral <- as.matrix(new_data[,cnames])%*% coef
      predictions_no_coral<-predictions_no_coral+intercept
      predictions[predictions<0]<-0
      predictions_no_coral[predictions_no_coral<0]<-0
      raw_predictions<-norm_to_raw(predictions,a$fe_n)
      raw_predictions_no_coral<-norm_to_raw(predictions_no_coral,a$fe_n)
      afl = 1-raw_predictions_no_coral/raw_predictions
      glob_loss<-norm_to_raw(new_data$fish,a$fe_n)*afl
      y1<-c(y1,100*(1-mean((fd0-glob_loss)/fd0)))
      fd0_nm[anom]<-nm[anom,2]
      y1_nm<-c(y1_nm,100*(mean((fd0-fd0_nm)/fd0)))
      }
    res_clim<-cbind(res_clim,y1,y1_nm)
    }


colnames(res_clim)<-c('year','loss_SSP2','loss_SSP2_nm','loss_SSP3','loss_SSP3_nm','loss_SSP5','loss_SSP5_nm',
                      'cum_loss_SSP2','cum_loss_SSP2_nm','cum_loss_SSP3','cum_loss_SSP3_nm','cum_loss_SSP5','cum_loss_SSP5_nm')


write.table(res_clim,paste0(prefix,'clim_tab.csv'),col.names=T,row.names=F,quote=F,sep=',')

     
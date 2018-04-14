#Mekong basin strontium isoscape project
#last edit 3/10/18
#author: Mike Vlah
#contact: vlahm13@gmail.com

#clear console and environment
rm(list=ls()); cat('\014')

#setup ####
library(SSN)
library(RColorBrewer)

#read in SSN objects
jan = importSSN('~/git/mekong_sr/round1/sr_jan.ssn', predpts='preds',
    o.write=TRUE)

jul = importSSN('~/git/mekong_sr/round1/sr_jul.ssn', predpts='preds',
    o.write=TRUE)

#check edge attributes
names(jan@data)
names(jul@data)
#'sr_erosion' and 'sr_noErosi' are sr ratios; '..._1' are these ratios
#multiplied by watershed area; '..._2' are the same as '..._1', but
#accumuated downstream. Dividing '..._2' by RCAs should yield weighted sr
#ratios, but some of them are non-ratio at the moment.)

#observation point attributes
names(jan@obspoints@SSNPoints[[1]]@point.data)
names(jul@obspoints@SSNPoints[[1]]@point.data)

#prediction point attributes
names(jan@predpoints@SSNPoints[[1]]@point.data)
names(jul@predpoints@SSNPoints[[1]]@point.data)

#how to get S4 object slots and their contents
getSlots('SpatialStreamNetwork')

#create distance matrices; specify amongpreds=TRUE if doing block kriging?
createDistMat(jan, predpts='preds', amongpreds=TRUE)
createDistMat(jul, predpts='preds', amongpreds=TRUE)

#model and diagnostics ####

#torgegram
# torg = Torgegram(jan, 'F87Sr_86Sr', nlag=6, maxlag=1000000)
pdf(width=7, height=7, file='~/git/mekong_sr/round1/figs/torgegram_jan.pdf',
    compress=FALSE)
plot(Torgegram(jan, 'F87Sr_86Sr', nlag=8, maxlag=600000))
dev.off()

pdf(width=7, height=7, file='~/git/mekong_sr/round1/figs/torgegram_jul.pdf',
    compress=FALSE)
plot(Torgegram(jul, 'F87Sr_86Sr', nlag=8, maxlag=600000))
dev.off()

#next, model sr87/86 throughout the basin; need to figure out how to associate
#Clement data with prediction sites (right now only associated with edges,
#so can't use as predictors)

#right now using additive function as the predictor below, just to have
#something (not informative).

#list of available predictors
names(jul@obspoints@SSNPoints[[1]]@point.data)

#we'd also want to pick only tailup or taildown I reckon (can't remember which
#is which at the moment). And we can talk about whether to include
#eucidean spatial correlation.
mod_jan = glmssn(F87Sr_86Sr ~ afvArea, jan,
    CorModels=c('Spherical.tailup', 'Spherical.taildown', 'Exponential.Euclid'),
    addfunccol='afvArea')
preds_jan = predict(mod_jan, 'preds')

mod_jul = glmssn(F87Sr_86Sr ~ afvArea, jul,
    CorModels=c('Spherical.tailup', 'Spherical.taildown', 'Exponential.Euclid'),
    addfunccol='afvArea')
preds_jul = predict(mod_jul, 'preds')

# plots ####

#docs for this particular plotting method
?plot.SpatialStreamNetwork

#plot network and sample points, save color break points for later (brks).

# defpar = par(bg='gray80')

#january
pdf(width=7, height=7, file='~/git/mekong_sr/round1/figs/sites_jan.pdf',
    compress=FALSE)
brks = plot(jan, VariableName='F87Sr_86Sr', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8)
    # PredPointsID='preds', add=FALSE, addWithLegend=TRUE)
    # color.palette=brewer.pal(8, 'Reds'), breaktype='even')
dev.off()

#july
pdf(width=7, height=7, file='~/git/mekong_sr/round1/figs/sites_jul.pdf',
    compress=FALSE)
brks = plot(jul, VariableName='F87Sr_86Sr', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8)
dev.off()


#plot predictions sized by certainty
pdf(width=7, height=7, file='~/git/mekong_sr/round1/figs/predictions_jan.pdf',
    compress=FALSE)
plot(preds_jan, SEcex.max=1.3, SEcex.min=0.4, breaktype='user', brks=brks,
    asp=1, dec.dig=7)
dev.off()

pdf(width=7, height=7, file='~/git/mekong_sr/round1/figs/predictions_jul.pdf',
    compress=FALSE)
plot(preds_jul, SEcex.max=1.3, SEcex.min=0.4, breaktype='user', brks=brks,
    asp=1, dec.dig=7)
dev.off()

# par(defpar)
# head(preds$ssn.object@predpoints@SSNPoints[[1]]@point.data$F87Sr_86Sr, 30)

# woud be nice to get google maps working eventually. here's a start ####
library(RgoogleMaps)
GetMap.bbox(lonR=105.21959, latR=14.89149,
map = GetMap.bbox(center=c(14.89149,105.21959), zoom=7, maptype='terrain')
PlotOnStaticMap(map, GRAYSCALE=TRUE)

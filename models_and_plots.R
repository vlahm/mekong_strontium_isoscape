#Mekong basin strontium isoscape project
#last edit 3/26/18
#author: Mike Vlah
#contact: vlahm13@gmail.com

#clear console and environment
rm(list=ls()); cat('\014')

#setup ####
library(SSN)
library(RColorBrewer)

#read in SSN objects
jan = importSSN('~/git/mekong_sr/round2/jan.ssn', predpts='preds',
    o.write=TRUE)

jul = importSSN('~/git/mekong_sr/round2/jul.ssn', predpts='preds',
    o.write=TRUE)

#check edge attributes
names(jan@data)
names(jul@data)
#'srEro' and 'srNoE' are sr ratios with and without erosion.
#'srEro_wtd' and 'srNoE_wtd' are the same, but with upstream influence accounted for

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

#torgegram (max lag should be about half the length of the study area;
#try different bin sizes)
# torg = Torgegram(jan, 'F87Sr_86Sr', nlag=6, maxlag=1000000)
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/torgegram_jan.pdf',
    compress=FALSE)
plot(Torgegram(jan, 'F87Sr_86Sr', nlag=8, maxlag=600000))
dev.off()

pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/torgegram_jul.pdf',
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

#tailup restricts autocorrelation to flow-connected sites.
#taildown allows flow-conn and flow-unconn autocorrelation.
#euclidean ignores network topology.
#all may be important in our case, so we should try all combinations.
#the "addfunccol" argument is needed if you specify a tailup model
mod_jan_ero = glmssn(F87Sr_86Sr ~ srEro_wtd, jan,
    CorModels=c('Spherical.tailup', 'Spherical.taildown', 'Exponential.Euclid'),
    addfunccol='afvArea')
summary(mod_jan_ero)
preds_jan_ero = predict(mod_jan_ero, 'preds')

mod_jul_ero = glmssn(F87Sr_86Sr ~ srEro_wtd, jul,
    CorModels=c('Spherical.tailup', 'Spherical.taildown', 'Exponential.Euclid'),
    addfunccol='afvArea')
preds_jul_ero = predict(mod_jul_ero, 'preds')

mod_jan_noEro = glmssn(F87Sr_86Sr ~ srNoE_wtd, jan,
    CorModels=c('Spherical.tailup', 'Spherical.taildown', 'Exponential.Euclid'),
    addfunccol='afvArea')
preds_jan_noEro = predict(mod_jan_noEro, 'preds')

mod_jul_noEro = glmssn(F87Sr_86Sr ~ srNoE_wtd, jul,
    CorModels=c('Spherical.tailup', 'Spherical.taildown', 'Exponential.Euclid'),
    addfunccol='afvArea')
preds_jul_noEro = predict(mod_jul_noEro, 'preds')

# plots ####
library(plotrix)

#docs for this particular plotting method
?plot.SpatialStreamNetwork

#experimenting...
# defpar = par(bg='gray80')

# par(mfrow=c(1,3))
# plot(as.SpatialLines(jan), col="gray50")
#
# #trying to do my own color ramp; not lining up so far
# pal = colorRampPalette(c('blue', 'red'))
# xx = as.SpatialPointsDataFrame(jan)
# b = as.numeric(cut(xx$srEro_wtd, breaks=8))
# table(b)
# cl = pal(length(unique(b)))
# cl = as.character(factor(b, labels=cl))
#
# # plot(as.SpatialPoints(jan), pch=20, add=TRUE, col=cl, cex=3)
# plot(as.SpatialPoints(jan), pch=21, add=TRUE, col=cl,#col='purple', bg='cyan1',
#     cex=rescale(as.SpatialPointsDataFrame(jan)$F87Sr_86Sr, c(1, 3)))
# # plot(as.SpatialPoints(jan, data="preds"), cex=1.5, add=TRUE, pch=20, col=cl)

#plot jan sr measurements
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/jan_samples.pdf',
    compress=FALSE)
brks = plot(jan, VariableName='F87Sr_86Sr', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8,
    breaktype='user', brks=c(.7,.7025,.705,.7075,.71,.7125,.715,.7175,.72))
    # PredPointsID='preds', add=FALSE, addWithLegend=TRUE)
    # color.palette=brewer.pal(8, 'Reds'), breaktype='even')
dev.off()

#plot modeled sr with erosion (at jan sample sites)
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/jan_erosion.pdf',
    compress=FALSE)
brks = plot(jan, VariableName='srEro_wtd', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8,
    breaktype='user', brks=c(.7,.7025,.705,.7075,.71,.7125,.715,.7175,.72))
dev.off()

#plot modeled sr without erosion (at jan sample sites)
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/jan_noErosion.pdf',
    compress=FALSE)
brks = plot(jan, VariableName='srNoE_wtd', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8,
    breaktype='user', brks=c(.7,.7025,.705,.7075,.71,.7125,.715,.7175,.72))
dev.off()

#plot july sr measurements
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/jul_samples.pdf',
    compress=FALSE)
brks = plot(jul, VariableName='F87Sr_86Sr', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8,
    breaktype='user', brks=c(.7,.7025,.705,.7075,.71,.7125,.715,.7175,.72))
# PredPointsID='preds', add=FALSE, addWithLegend=TRUE)
# color.palette=brewer.pal(8, 'Reds'), breaktype='even')
dev.off()

#plot modeled sr with erosion (at jul sample sites)
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/jul_erosion.pdf',
    compress=FALSE)
brks = plot(jul, VariableName='srEro_wtd', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8,
    breaktype='user', brks=c(.7,.7025,.705,.7075,.71,.7125,.715,.7175,.72))
dev.off()

#plot modeled sr without erosion (at jul sample sites)
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/jul_noErosion.pdf',
    compress=FALSE)
brks = plot(jul, VariableName='srNoE_wtd', lwdLineCol='afvArea',
    lwdLineEx=15, lineCol='black', pch=19, cex=1.3,
    xlab='x-coordinate (m)', ylab='y-coordinate (m)', asp = 1, nclasses=8,
    breaktype='user', brks=c(.7,.7025,.705,.7075,.71,.7125,.715,.7175,.72))
dev.off()

#plot predictions sized by certainty (jan, with erosion)
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/predictions_jan_ero.pdf',
    compress=FALSE)
plot(preds_jan_ero, SEcex.max=1.3, SEcex.min=0.4, breaktype='user', brks=brks,
    asp=1, dec.dig=7)
dev.off()

#predictions for jan without erosion
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/predictions_jan_noEro.pdf',
    compress=FALSE)
plot(preds_jan_noEro, SEcex.max=1.3, SEcex.min=0.4, breaktype='user', brks=brks,
    asp=1, dec.dig=7)
dev.off()

#preds july ero
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/predictions_jul_ero.pdf',
    compress=FALSE)
plot(preds_jul_ero, SEcex.max=1.3, SEcex.min=0.4, breaktype='user', brks=brks,
    asp=1, dec.dig=7)
dev.off()

#preds july no ero
pdf(width=7, height=7, file='~/git/mekong_sr/round2/figs/predictions_jul_noEro.pdf',
    compress=FALSE)
plot(preds_jul_noEro, SEcex.max=1.3, SEcex.min=0.4, breaktype='user', brks=brks,
    asp=1, dec.dig=7)
dev.off()

# par(defpar)
# head(preds$ssn.object@predpoints@SSNPoints[[1]]@point.data$F87Sr_86Sr, 30)

# woud be nice to get google maps working eventually. here's a start ####
library(RgoogleMaps)
GetMap.bbox(lonR=105.21959, latR=14.89149,
map = GetMap.bbox(center=c(14.89149,105.21959), zoom=7, maptype='terrain')
PlotOnStaticMap(map, GRAYSCALE=TRUE)

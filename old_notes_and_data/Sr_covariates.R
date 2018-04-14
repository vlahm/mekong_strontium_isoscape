#checking for correlations between strontium isotope ratios and concentrations of various elements

x <- read.csv("C:\\Users\\Mike\\Desktop\\Grad\\Projects\\strontium_isoscape\\Sr_covariates.csv",
              stringsAsFactors=FALSE)
x <- read.csv("/home/mike/Desktop/grad/Projects/strontium_isoscape/Sr_covariates.csv",
              stringsAsFactors=FALSE)

#replacing all "<0.0000x" with "0.0000x"
ischar <- rep(NA, length(17:ncol(x)))
for(i in 17:(ncol(x))){
    ischar[i-16] <- is.character(x[1,i])
}

x2 <- x
for(i in 1:nrow(x)){
    for(j in 17:ncol(x)){
        if(ischar[j-16] & substr(x[i,j], 1, 1) == '<'){
            x2[i,j] <- substr(x[i,j], 2, nchar(x[i,j]))
        }
    }
}

#more fixer-uppering
colnames(x2)[15] <- 'sr87_86'
x2[,17:ncol(x2)] <- apply(x2[,17:ncol(x2)], 2, as.numeric)
x3 <- x2[-which(is.na(x2$sr87_86)),]

#take mean of duplicate samples
x4 <- aggregate(x3[,15:ncol(x3)], by=list(x3$Sample.Number), FUN=mean)

#normalizing by calcium
non_Ca_cols <- c(4:9,11:ncol(x4))
x4[,non_Ca_cols] <- x4[,non_Ca_cols]/x4$Ca

#correlations
as.matrix(round(cor(x4[,c(2,4:ncol(x4))])[,1], 3))
par(mfrow=c(4,3))
for(i in c(8,9,10,11,14,16,19,23,24,25,27,29)){
    mod <- lm(x4$sr87_86 ~ x4[,i])
    if(i==10){
        plot(x4[,i], x4$sr87_86, xlab='Ca', pch=20, ylab='Sr 86/86',
             main=paste('Adj. R^2 =', round(summary(mod)$adj.r.squared, 3)))
    } else {
    plot(x4[,i], x4$sr87_86, xlab=paste(colnames(x4[i]), '/Ca', sep=''), pch=20, ylab='Sr 86/86',
         main=paste('Adj. R^2 =', round(summary(mod)$adj.r.squared, 3)))
    }
    abline(mod)
}

#export
write.csv(x4, "/home/mike/Desktop/grad/Projects/strontium_isoscape/Sr_covariates_clean.csv")

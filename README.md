# GM-cofactors_and_LPD

## Figure 2A

    library(gplots)
    library(clusterSim)
    require(RColorBrewer)

    heatmap <- heatmap.2(as.matrix(X), distfun = function(x) dist(x, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),Rowv = T, Colv = T,     dendrogram ="both", col = colorRampPalette(c('darkblue', 'cyan','chartreuse1', 'chartreuse2','chartreuse3','chartreuse3','orange1', 'orange2','orange3','orange3','coral4'))(2000), key = T, keysize = 1.0, density.info = "none", trace ="none", na.color = NA, margins = c(20,20))

## Figure 2D
    library(ggplot2)
    library(ggridges)
    ggplot(data, aes(type,BMI, color= type)) + geom_boxplot(aes(colour = type)) + geom_jitter(aes(colour = type),width = 0.2) + scale_color_manual(values=c("dodgerblue4","lightskyblue3","tomato1")) + theme_classic()


## Figure 2E
    library(ggplot2)
    library(ggridges)
    ggplot(data, aes(type,Age, color= type)) + geom_boxplot(aes(colour = type)) + geom_jitter(aes(colour = type),width = 0.2) + scale_color_manual(values=c("dodgerblue4","lightskyblue3","tomato1")) + theme_classic()


## Figure 3A

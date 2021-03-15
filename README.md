# GM-cofactors_and_LPD

## Figure 2A

    library(gplots)
    library(clusterSim)
    require(RColorBrewer)

    heatmap <- heatmap.2(as.matrix(X), distfun = function(x) dist(x, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),Rowv = T, Colv = T,     dendrogram ="both", col = colorRampPalette(c('darkblue', 'cyan','chartreuse1', 'chartreuse2','chartreuse3','chartreuse3','orange1', 'orange2','orange3','orange3','coral4'))(2000), key = T, keysize = 1.0, density.info = "none", trace ="none", na.color = NA, margins = c(20,20))


## Figure 2B
    library(vegan)
    mod0 <- rda(UCX ~ 1, UCY) 
    mod1 <- rda(UCX ~ ., UCY)
    step.res <- ordiR2step(mod0, mod1, perm.max = 200)
    step.res$anova
    

## Figure 2D
    library(FSA)
    dunnTest(data$BMI,data$type)
    
    library(ggplot2)
    library(ggridges)
    ggplot(data, aes(type,BMI, color= type)) + geom_boxplot(aes(colour = type)) + geom_jitter(aes(colour = type),width = 0.2) + scale_color_manual(values=c("dodgerblue4","lightskyblue3","tomato1")) + theme_classic()


## Figure 2E
    library(FSA)
    dunnTest(data$BMI,data$type)

    library(ggplot2)
    library(ggridges)
    ggplot(data, aes(type,Age, color= type)) + geom_boxplot(aes(colour = type)) + geom_jitter(aes(colour = type),width = 0.2) + scale_color_manual(values=c("dodgerblue4","lightskyblue3","tomato1")) + theme_classic()

__________________________

## Figure 3A
    cor.test(X,Y, method = "spearman")

## Figure 3B
    library(ggplot2)
    ggplot(matrix1, aes(type,species, color= type)) + geom_boxplot(aes(colour = type)) + geom_jitter(aes(colour = type),width = 0.2) + scale_color_manual(values=c("dodgerblue4","lightskyblue3","tomato1")) + theme_classic()

## Figure 3C
    library(ggplot2)
    library(ggridges)
    
    ggplot(matrix1, aes(type,species, color= type)) + geom_boxplot(aes(colour = type)) + geom_jitter(aes(colour = type),width = 0.2) +  scale_color_manual(values=c("dodgerblue4","lightskyblue3","tomato1", "seagreen3")) + theme_classic()
    
    ggplot(matrix1, aes(x = species, y = type, fill = type, color = time)) +
    geom_density_ridges(alpha = 0.7, rel_min_height = 0.01,
                      color = "white", from =2200, to = 5000) + scale_fill_cyclical(values = c("dodgerblue4","lightskyblue3","tomato1", "seagreen3"),  guide = "legend") + geom_vline(aes(xintercept = 3595), linetype = "dashed") + theme_classic()

## Figure 3D
    library(ggplot2)
    library(ggridges)
    ggplot(C, aes(x=species, fill=cat_BMI)) + geom_density(aes(y = stat(count)),alpha=0.4) + geom_vline(aes(xintercept = 3595), linetype = "dashed") + scale_fill_manual(values=c("gray80", "lightblue", "#E69F00")) + theme_classic()


## Figure 3E
    library(vegan)
    adonis(GMsX ~ GMsY$type, method = "canberra")
    adonis(KOsX ~ GMsY$type, method = "canberra")



__________________________
## Figure 4A|B|D|E

For zOTUs & KO counts matrices
    
    #Varible selection based on RandomForest




    #PCoA based on selected variables
    
    library(vegan)
    library(ggplot2)
    mod1 <- capscale(X ~ Y$type, dist = "canberra")

    pcs1 <- as.data.frame(mod1$CA$u)
    pcsA1 <- pcs1[,(1:4)]
    rda1 = as.data.frame(mod1$CCA$w)
    rdaA1 = rda1[,(1:2)]
    matrix1 <- cbind(pcsA1,rdaA1,Y)

    sp1<-ggplot(matrix1, aes(x = CAP1, y = CAP2,colour=type)) + geom_point(size =3)
    sp1 + scale_color_manual(values = c("dodgerblue4","lightskyblue3","tomato1", "seagreen3")) + theme_classic()

## Figure 4B

## Figure 4C
## Figure 4D
## Figure 4E
## Figure 4F







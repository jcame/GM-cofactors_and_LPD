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

    echo '#!/usr/bin/env Rscript

    calm <- read.table("dataset.txt",header=TRUE)
    seed_numbers <- sample(1:100000, 1)
    set.seed(seed_numbers)


    A <- subset(calm, type == "Cluster1")
    C <- subset(calm, type == "Cluster2")


    s_A <- sample(82,57)
    s_C <- sample(67,47)


    A_train <- A[s_A,]
    A_test <- A[-s_A,]
    C_train <- C[s_C,]
    C_test <- C[-s_C,]


    train0 <- rbind(A_train,C_train)
    test0 <- rbind(A_test,C_test)
    write.table(test0, "test_set.txt", sep="\t", row.names=FALSE, quote =F)
    write.table(train0, "training_set.txt", sep="\t", row.names=FALSE, quote =F)



    calm <- read.table("training_set.txt",header=TRUE)
    train <-  calm[,-c(1)]

    seed_numbers <- sample(1:100000, 1)
    set.seed(seed_numbers)
    library(party)
    iris.cf <- cforest(type ~ ., data = train, control = cforest_unbiased(ntree = 6000))
    vi <- varimp(iris.cf, nperm = 10)
    write.table(vi, "decision.txt", sep="\t",quote =F) ' > cforest.R

    chmod 755 cforest.R
    ./cforest.R
    rm cforest.R

    TAB=$'\t'
    sed 1d decision.txt > decision1.txt
    echo "variable${TAB}MDA" > header.txt
    cat header.txt decision1.txt > Party-MDA.txt
    rm decision.txt
    rm decision1.txt
    rm header.txt

    #########################################################################

    echo '#!/usr/bin/env Rscript

    seed_numbers <- sample(1:100000, 1)
    set.seed(seed_numbers)
    library(randomForest)

    train0<- read.table("training_set.txt",header=TRUE)
    test0 <- read.table("test_set.txt",header=TRUE)
    train <-  train0[,-c(1)]
    test <-  test0[,-c(1)]


    X <-read.table("Party-MDA.txt",header=TRUE)
    Xs = X[X$MDA>0, ]
    names = as.vector(Xs$variable)

    train.1 <- train[,names] 
    type <- train$type
    train_1 <- cbind(type,train.1)

    test.1 <- test[,names] 
    type <- test$type
    test_1 <- cbind(type,test.1)

    GM.rf_1 <- randomForest(x=train_1[-1], y=train_1$type, data = train_1, importance = TRUE, proximity = TRUE, ntree= 6000, sampsize=c(47,47), nPerm=1000, iter=1000)
    p <- predict (GM.rf_1, test_1)
    print <- mean(test_1[,1]==p)
    write.table(print, "prediction.txt", sep="\t",quote =F)' > randomForest.R

    chmod 755 randomForest.R
    ./randomForest.R
    rm randomForest.R

    sed 1d prediction.txt > 1.txt
    cut -f2 1.txt > 2.txt
    mv 2.txt prediction_${dir}.txt
    rm 1.txt
    rm prediction.txt


    #########################################################################


    echo '#!/usr/bin/env Rscript

    calm <- read.table("dataset.txt",header=TRUE)
    #train <-  calm[,-c(1)]

    X <-read.table("Party-MDA.txt",header=TRUE)
    Xs = X[X$MDA>0, ]
    names = as.vector(Xs$variable)

    #train.1 <- train[,names] 
    #type <- train$type
    #train <- cbind(type,train.1)


    selected <- calm[,names] 
    info = calm[,c(1,2)]

    new_table <- cbind(info,selected)
    write.table(new_table, row.names = FALSE, "selected_variables.txt", sep="\t", quote = F)' > selected.R

    chmod 755  selected.R
    ./selected.R
    rm  selected.R


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



## Figure 4G|H

zOTUs (randomForest selected) vs LPDs

    cor <- rcorr(as.matrix(X), as.matrix(Y), type="spearman")

    correction <- (p.adjust(cor$P, method = "fdr")) < 0.05
    corrected <- correction

    corrected[corrected==TRUE] = 1

    new <-  corrected*cor$r

    new1 <- new[,-c(1:206)]
    new2 <- new1[-c(207:260),]
    cclean <- (colSums(new2, na.rm=T) != 0)
    new3 <- new2[, cclean] 
    rclean <- (rowSums(new3, na.rm=T) != 0)
    new4 <- new3[rclean,] 

    heatmap <- heatmap.2(t(new4), Rowv = T, Colv = T, dendrogram ="both", col = colorRampPalette(c('red4','red2','coral1','white','gray90','royalblue2','royalblue4'))(1000), key = T, keysize = 1.0, density.info = "none", trace ="none", na.color = NA, margins = c(8,30))


KOs (randomForest selected) vs LPDs

    cor <- rcorr(as.matrix(Z), as.matrix(Y), type="spearman")


    correction <- (p.adjust(cor$P, method = "fdr")) < 0.05
    corrected <- correction

    corrected[corrected==TRUE] = 1

    new <-  corrected*cor$r

    new1 <- new[,-c(1:55)]
    new2 <- new1[-c(56:109),]
    cclean <- (colSums(new2, na.rm=T) != 0)
    new3 <- new2[, cclean] 
    rclean <- (rowSums(new3, na.rm=T) != 0)
    new4 <- new3[rclean,] 



    heatmap <- heatmap.2(new4, Rowv = T, Colv = T, dendrogram ="both", col = colorRampPalette(c('red4','red2','white','deepskyblue1','blue4'))(1000), key = T, keysize = 1.0, density.info = "none", trace ="none", na.color = NA, margins = c(8,30))






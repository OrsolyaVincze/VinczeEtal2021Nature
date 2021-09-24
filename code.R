setwd("./VinczeEtal2021Nature")

library(ape)
library(geiger)
library(car)
library(phytools)
library(nlme)
library(phylolm) # phyloglm
library(emmeans)

# Zero-inflated phylogenetic model
ZIlogis <- function(formbin,   # model formula for binomial model
                    formlogis, # model formula for logistic model
                    databin,   # dataset for binomial GLM repsonse variable
                    datalogis, # dataset with logistic response variable
                    phyloTree, # phylogeneitc tree for logistic regression
                    phylobin,  # phylogenetic tree for binomial regression
                    l0) {      # starting value of Pagel's lambda in PGLS regression
                            out <- list()
                            out$bin <- binaryPGLMM(formbin, data = databin, phy = phylobin)
                            out$wlogis <- gls(formlogis,
                                             correlation = corPagel(l0, phy = phyloTree, fixed = FALSE, form=~Species),
                                             data = datalogis, weights = ~1/log(knownDeaths) )
                                if(out$wlogis$modelStruct[1] < 0){ # refit model with lanbda=0 if lambda converged to negative
                                        out$wlogis <- gls(formlogis,
                                             correlation = corPagel(0, phy = phyloTree, fixed = TRUE, form=~Species),
                                             data = datalogis, weights = ~1/log(knownDeaths))}
                            class(out) <- "ZILogis"
                            return(out)}

# Summary function for the zero-inflated phylogenetic model
sumZILogis <- function(x,...){
        x1 <- cbind(as.data.frame(x$bin[[2]]), as.data.frame(x$bin[3]), as.data.frame(x$bin[5]), as.data.frame(x$bin[6]))
        x1[,1:3] <- round(x1[,1:3],2)
        x1[,4] <- round(x1[,4],4)
        x1[,1] <- paste(x1[,1]," (", x1[,2],")", sep='')
        x1 <- x1[,-2]
        x1b <- cbind(row.names(x1),x1[,1:3]);
          row.names(x1b) <- NULL
          x1b[nrow(x1)+1, 4] <- paste("n=", nrow(as.data.frame(x$bin[9])))
          x1b[nrow(x1)+1, 3] <- paste("s2=", round(as.numeric(x$bin[7]),2))
          x1b[nrow(x1)+1, 2] <- paste("P_s2=", round(as.numeric(x$bin[8]),4))
          x1b[,1] <- as.character(x1b[,1])
          x1b[nrow(x1)+1, 1] <- "ModelStats"
        x2w <- as.data.frame(summary(x$wlogis)$tTable); names(x2w) <- names(x1)
          x2w[,1:3] <- round(x2w[,1:3],2)
          x2w[,4] <- round(x2w[,4],4)
          x2w[,1] <- paste(x2w[,1]," (", x2w[,2],")", sep='')
          x2w <- x2w[,-2]
          x2wb <- cbind(row.names(x2w),x2w[,1:3]);
          row.names(x2wb) <- NULL
          x2wb[nrow(x2w)+1, 2] <- paste("AIC=", round(summary(x$wlogis)$AIC,2))
          x2wb[nrow(x2w)+1, 4] <- paste("n=", length(x$wlogis$residuals))
          x2wb[nrow(x2w)+1, 3] <- paste("L=", round(unlist(summary(x$wlogis))$modelStruct.corStruct,2))
          x2wb[,1] <- as.character(x2wb[,1])
          x2wb[nrow(x2w)+1, 1] <- "ModelStats"
          names(x2wb) <- names(x1b)
          names(x1b)[1] <- ""
          x3 <- x1b[1,]; x3[,1:4] <- ""; x3[1,1] <- "Probability of zeros"
        names(x2wb)[1] <- ""
        x5 <- x3; x5[1,1] <- "Weighted logistic reg."
        X <- rbind(x3,x1b,x5,x2wb)
        names(X)[1:4] <- c('_',"b (SE)","z/t-value", "p-value")
        return(X)        }

sumGLS <- function(model){ # summary function GLS
      su <- as.data.frame(summary(model)$tTable)
      su[,1:3] <- round(su[,1:3],2); su[,4] <- round(su[,4],4)
      su[,1] <- paste(su[,1],'_(', su[,2],')', sep='')
      su <- su[,c(1,3,4)]
        x <- nrow(su)
      #su[x+1,1] <- paste("ShaP=",round(as.numeric(unlist(shapiro.test(resid(model)))[[2]][1]),4), sep="")
      su[x+1,2] <- paste("AIC=", round(AIC(model),2), sep="")
      su[x+1,3] <- paste("n=", length(resid(model)), sep='')
      su[x+1,1] <-paste('lambda=',round(as.numeric(model[[1]])[1], 2), sep='') # put lambda value in bottom right corner
      row.names(su)[nrow(su)] <- "ModelStats"
      return(su)
    }

  # data 
  data <- read.csv("data.csv")
  row.names(data) <- data$Species

  # Species-specific body mass
  data$BodyMass <- (data$MaleMeanMass + data$FemaleMeanMass)/2
  
  # phylogeny
  phy <- read.nexus("consensus_phylogeny.tre") # consensus vertlife tree
  phy <- bind.tip(phy, "Cervus_canadensis", where = which(phy$tip.label=="Cervus_elaphus"),
                   edge.length=0.5, position = 0.5)
  phy <- bind.tip(phy, "Gazella_marica", where = which(phy$tip.label=="Gazella_subgutturosa"),
                   edge.length=0.5, position = 0.5)

# Phylogenetic signal in cancer mortality
  ICM <- data$ICM; names(ICM) <- data$Species
  ICM <- ICM[!is.na(ICM)]
  # n for ICM
  length(ICM)
        #  [1] 172
  phylo <- treedata(phy,ICM, warnings = F)$phy
  phylosig(phylo, ICM, method = "lambda", test = T)
          # Phylogenetic signal lambda : 0.685124 
          # logL(lambda) : 215.992 
          # LR(lambda=0) : 33.1224 
          # P-value (based on LR test) : 8.65347e-09 

  # phylogenetic signal of Cancer prevalence ===========================
  CMR <- data$CMR; names(CMR) <- data$Species
  CMR <- CMR[!is.na(CMR)]
  # n for CMR
  length(CMR)
      # [1] 191
  phylo <- treedata(phy,CMR, warnings = F)$phy
  phylosig(phylo, CMR, method = "lambda", test = T)
      # Phylogenetic signal lambda : 0.86979 
      # logL(lambda) : 262.945 
      # LR(lambda=0) : 58.0754 
      # P-value (based on LR test) : 2.5226e-14 
  
  
#====================================================================================================================
#  Extended Data Tablex 1. - Results of linear regression models exploring order differences in cancer risk, measured as
  # a, CMR or b, ICM. Only orders represented by a minimum of two species in our database were included in these analyses.
  # Models were constructed using R function "lm" using cancer risk as a dependent variable and order as a sole predictor.
  # Values presented are estimated marginal means and their associated 95% confidence intervals. The number of species 
  # for which data was available in each order is also shown. Post-hoc test of order differences in c, CMR or d, ICM
  # across mammalian orders are shown.
  mx <- lm(CMR*100 ~ order, data[data$order%in%names(table(data$order))[which(table(data$order)>1)],])
  newdata = expand.grid(order=names(table(data$order))[which(table(data$order)>1)])
  # Extended Data Tablex 1a - Order average CMR
  cbind(newdata,
        as.data.frame(emmeans(mx,~order)))
            #             order          order    emmean        SE  df  lower.CL  upper.CL
            # 1   Artiodactyla   Artiodactyla  2.313159 0.6713050 178 0.9884187  3.637900
            # 2      Carnivora      Carnivora 10.643211 0.8452490 178 8.9752132 12.311210
            # 3     Chiroptera     Chiroptera  6.208813 1.8040781 178 2.6486801  9.768947
            # 4  Diprotodontia  Diprotodontia  5.663575 2.2095354 178 1.3033196 10.023829
            # 5 Perissodactyla Perissodactyla  5.038047 2.2095354 178 0.6777922  9.398302
            # 6       Primates       Primates  5.631310 0.8351259 178 3.9832885  7.279331
            # 7       Rodentia       Rodentia  6.077022 1.3530586 178 3.4069223  8.747122
  # Extended Data Tablex 1c - Pairwise order differences in CMR
  as.data.frame(pairs(emmeans(mx,~order)))
              #                         contrast    estimate       SE  df     t.ratio      p.value
              # 1        Artiodactyla - Carnivora -8.33005221 1.079396 178 -7.71732531 1.726741e-11
              # 2       Artiodactyla - Chiroptera -3.89565414 1.924928 178 -2.02379201 4.034288e-01
              # 3    Artiodactyla - Diprotodontia -3.35041534 2.309263 178 -1.45085893 7.731770e-01
              # 4   Artiodactyla - Perissodactyla -2.72488796 2.309263 178 -1.17998148 9.008724e-01
              # 5         Artiodactyla - Primates -3.31815076 1.071488 178 -3.09677006 3.617757e-02
              # 6         Artiodactyla - Rodentia -3.76386301 1.510436 178 -2.49190441 1.685245e-01
              # 7          Carnivora - Chiroptera  4.43439807 1.992271 178  2.22580066 2.870913e-01
              # 8       Carnivora - Diprotodontia  4.97963687 2.365691 178  2.10493994 3.543094e-01
              # 9      Carnivora - Perissodactyla  5.60516425 2.365691 178  2.36935632 2.179055e-01
              # 10           Carnivora - Primates  5.01190145 1.188226 178  4.21796974 7.645230e-04
              # 11           Carnivora - Rodentia  4.56618921 1.595372 178  2.86214614 6.914756e-02
              # 12     Chiroptera - Diprotodontia  0.54523881 2.852498 178  0.19114433 9.999958e-01
              # 13    Chiroptera - Perissodactyla  1.17076618 2.852498 178  0.41043542 9.996226e-01
              # 14          Chiroptera - Primates  0.57750338 1.987997 178  0.29049506 9.999499e-01
              # 15          Chiroptera - Rodentia  0.13179114 2.255098 178  0.05844143 1.000000e+00
              # 16 Diprotodontia - Perissodactyla  0.62552737 3.124755 178  0.20018446 9.999945e-01
              # 17       Diprotodontia - Primates  0.03226458 2.362093 178  0.01365932 1.000000e+00
              # 18       Diprotodontia - Rodentia -0.41344767 2.590910 178 -0.15957624 9.999986e-01
              # 19      Perissodactyla - Primates -0.59326280 2.362093 178 -0.25115983 9.999788e-01
              # 20      Perissodactyla - Rodentia -1.03897504 2.590910 178 -0.40100778 9.996701e-01
              # 21            Primates - Rodentia -0.44571225 1.590032 178 -0.28031647 9.999594e-01

  
  mx <- lm(ICM*100 ~ order, data [data$order%in%names(table(data$order))[which(table(data$order)>1)],])
  newdata = expand.grid(order=names(table(data$order))[which(table(data$order)>1)])
  # Extended Data Tablex 1b - Order average CMR
  cbind(newdata, as.data.frame(emmeans(mx,~order)))
              #           order          order    emmean        SE  df  lower.CL  upper.CL
              # 1   Artiodactyla   Artiodactyla  3.184424 0.7817053 160  1.640633  4.728215
              # 2      Carnivora      Carnivora 11.173098 0.9612736 160  9.274678 13.071519
              # 3     Chiroptera     Chiroptera  6.242697 2.3264295 160  1.648228 10.837166
              # 4  Diprotodontia  Diprotodontia  7.969023 2.5128311 160  3.006428 12.931617
              # 5 Perissodactyla Perissodactyla  2.418117 2.5128311 160 -2.544477  7.380711
              # 6       Primates       Primates  5.552358 1.0714748 160  3.436300  7.668415
              # 7       Rodentia       Rodentia  5.399414 1.7768399 160  1.890330  8.908498
  # Extended Data Tablex 1d - Pairwise order differences in ICM
  as.data.frame(pairs(emmeans(mx,~order)))
              #                         contrast   estimate       SE  df     t.ratio      p.value
              # 1        Artiodactyla - Carnivora -7.9886746 1.238996 160 -6.44770216 2.677097e-08
              # 2       Artiodactyla - Chiroptera -3.0582734 2.454249 160 -1.24611378 8.748351e-01
              # 3    Artiodactyla - Diprotodontia -4.7845987 2.631612 160 -1.81812447 5.378729e-01
              # 4   Artiodactyla - Perissodactyla  0.7663068 2.631612 160  0.29119290 9.999490e-01
              # 5         Artiodactyla - Primates -2.3679339 1.326319 160 -1.78534297 5.598859e-01
              # 6         Artiodactyla - Rodentia -2.2149902 1.941191 160 -1.14104691 9.143433e-01
              # 7          Carnivora - Chiroptera  4.9304012 2.517205 160  1.95868080 4.450979e-01
              # 8       Carnivora - Diprotodontia  3.2040758 2.690421 160  1.19091972 8.967189e-01
              # 9      Carnivora - Perissodactyla  8.7549814 2.690421 160  3.25413022 2.301359e-02
              # 10           Carnivora - Primates  5.6207407 1.439481 160  3.90469994 2.607203e-03
              # 11           Carnivora - Rodentia  5.7736843 2.020200 160  2.85797700 7.051598e-02
              # 12     Chiroptera - Diprotodontia -1.7263254 3.424412 160 -0.50412322 9.987710e-01
              # 13    Chiroptera - Perissodactyla  3.8245802 3.424412 160  1.11685763 9.221577e-01
              # 14          Chiroptera - Primates  0.6903395 2.561315 160  0.26952546 9.999677e-01
              # 15          Chiroptera - Rodentia  0.8432831 2.927360 160  0.28806954 9.999521e-01
              # 16 Diprotodontia - Perissodactyla  5.5509056 3.553680 160  1.56201623 7.065785e-01
              # 17       Diprotodontia - Primates  2.4166649 2.731735 160  0.88466287 9.744103e-01
              # 18       Diprotodontia - Rodentia  2.5696085 3.077577 160  0.83494532 9.809147e-01
              # 19      Perissodactyla - Primates -3.1342407 2.731735 160 -1.14734417 9.122305e-01
              # 20      Perissodactyla - Rodentia -2.9812970 3.077577 160 -0.96871566 9.599654e-01
              # 21            Primates - Rodentia  0.1529436 2.074902 160  0.07371126 1.000000e+00
              # 
  

  
#====================================================================================================================
# Extended Data Tablex 2. - Association between cancer risk and diet animal content. Models explain variance in 
  # a, CMR or b, ICM in function of body mass, life expectancy and diet animal content. Diet animal content is characterised 
  # by variables on three taxonomic level: 1. animal content (including any vertebrate or invertebrate prey); 2. invertebrate
  # or vertebrate prey; and 3. within vertebrates fish, reptiles, birds and mammals. Each diet item is coded as a two-level 
  # factor: rarely/never occurring in diet and representing the primary/secondary food item of the species. Each diet variable is added 
  # one by one to a base model containing the two significant predictors of cancer risk: body mass and life expectancy. 
  # Akaike's Information Criteria (AIC) are directly comparable among models with the same independent variable. All models are PGLS regressions.
  
  #### Extended Data Tablex 2a ---------------------------------------------------------------
  # Models of non-zero CMR in function of life expectancy, body mass and diet 
  dataCMR <- data[!data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus", "Dasyuroides_byrnei") &
                    !is.na(data$Vertebrate),]; row.names(dataCMR) <- dataCMR$Species
  phyCMR <- treedata(phy,dataCMR, warnings = FALSE)$phy
  dataCMR$occurrence <- ifelse(dataCMR$CMR==0,0,1)
  dataCMR <- dataCMR[as.character(phyCMR$tip.label), ] # organise data to follow same order as phylogeny
  # Drop tips on tree:
  phy_nonzeroCMR <- drop.tip(phyCMR, which(!dataCMR$occurrence==1))
  # logistic data:
  data_nonzeroCMR <- dataCMR[dataCMR$Species%in%phy_nonzeroCMR$tip.label, ] # organise data to follow same order as phylogeny
  #models
  DietMod0 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod1 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp)+as.factor(Animal), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod2 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp)+as.factor(Vertebrate), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod3 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp)+as.factor(Invertebrate), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod4 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp)+as.factor(Mammal), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod5 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp)+as.factor(Bird), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod6 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp)+as.factor(Herptile), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod7 <- gls(logit(CMR) ~  log(BodyMass)+log(lifeexp)+as.factor(Fish), data_nonzeroCMR,
             cor=corPagel(0,phy_nonzeroCMR, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  X <- as.data.frame(array(,c(1,3))); names(X) <- names(sumGLS(DietMod0))
  rbind(X, sumGLS(DietMod0), X, sumGLS(DietMod1),X, sumGLS(DietMod3), X, 
        sumGLS(DietMod2), X, sumGLS(DietMod7), X, sumGLS(DietMod6), X, 
        sumGLS(DietMod5), X, sumGLS(DietMod4))
            #                               Value    t-value p-value
            # 1                                <NA>       <NA>    <NA>
            #   (Intercept)              -7.13_(1.83)       -3.9   2e-04
            # log(BodyMass)            -0.04_(0.05)      -0.87   0.385
            # log(lifeexp)              0.42_(0.23)       1.83  0.0698
            # ModelStats                lambda=0.53 AIC=367.42   n=140
            # 11                               <NA>       <NA>    <NA>
            #   (Intercept)1             -7.15_(1.83)      -3.91   1e-04
            # log(BodyMass)1           -0.04_(0.05)       -0.7  0.4869
            # log(lifeexp)1             0.42_(0.23)       1.82  0.0714
            # as.factor(Animal)1        0.12_(0.23)        0.5  0.6146
            # ModelStats1               lambda=0.53 AIC=370.24   n=140
            # 12                               <NA>       <NA>    <NA>
            #   (Intercept)2             -7.12_(1.84)      -3.88   2e-04
            # log(BodyMass)2           -0.05_(0.05)      -0.99  0.3237
            # log(lifeexp)2             0.43_(0.23)       1.84  0.0682
            # as.factor(Invertebrate)1 -0.11_(0.21)      -0.53  0.5938
            # ModelStats2               lambda=0.54 AIC=370.45   n=140
            # 13                               <NA>       <NA>    <NA>
            #   (Intercept)3              -7.1_(1.77)      -4.01   1e-04
            # log(BodyMass)3           -0.02_(0.05)      -0.49  0.6232
            # log(lifeexp)3             0.41_(0.22)        1.8  0.0735
            # as.factor(Vertebrate)1    0.56_(0.22)       2.51  0.0134
            # ModelStats3                lambda=0.5 AIC=364.45   n=140
            # 14                               <NA>       <NA>    <NA>
            #   (Intercept)4             -6.94_(1.86)      -3.72   3e-04
            # log(BodyMass)4           -0.04_(0.05)       -0.9  0.3721
            # log(lifeexp)4             0.41_(0.24)       1.73  0.0864
            # as.factor(Fish)1          0.16_(0.34)       0.46   0.646
            # ModelStats4               lambda=0.52 AIC=369.54   n=140
            # 15                               <NA>       <NA>    <NA>
            #   (Intercept)5              -7.1_(1.83)      -3.89   2e-04
            # log(BodyMass)5           -0.04_(0.05)       -0.7  0.4822
            # log(lifeexp)5             0.42_(0.23)       1.81   0.072
            # as.factor(Herptile)1      0.16_(0.22)       0.71  0.4759
            # ModelStats5               lambda=0.53 AIC=370.09   n=140
            # 16                               <NA>       <NA>    <NA>
            #   (Intercept)6             -7.17_(1.86)      -3.86   2e-04
            # log(BodyMass)6           -0.04_(0.05)      -0.87  0.3873
            # log(lifeexp)6             0.43_(0.23)       1.82  0.0713
            # as.factor(Bird)1         -0.04_(0.28)      -0.14  0.8889
            # ModelStats6               lambda=0.54 AIC=370.13   n=140
            # 17                               <NA>       <NA>    <NA>
            #   (Intercept)7              -6.92_(1.7)      -4.08   1e-04
            # log(BodyMass)7           -0.06_(0.04)      -1.41  0.1622
            # log(lifeexp)7             0.44_(0.21)       2.06  0.0418
            # as.factor(Mammal)1        0.57_(0.25)       2.27  0.0245
            # ModelStats7               lambda=0.36 AIC=366.49   n=140
            
            
  ####  Extended Data Tablex 2b --------------------------------------------------------------------
  # Models of non-zero ICM in function of life expectancy, body mass and diet 
  dataICM <- data[!is.na(data$ICM) & 
                    !data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus", "Dasyuroides_byrnei") &
                    !is.na(data$Vertebrate),]
  row.names(dataICM) <- dataICM$Species
  phyICM <- treedata(phy,dataICM, warnings = FALSE)$phy
  dataICM$occurrence <- ifelse(dataICM$ICM==0,0,1)
  dataICM <- dataICM[as.character(phyICM$tip.label), ] # organise data to follow same order as phylogeny
  # Drop tips on tree:
  phy_nonzeroICM <- drop.tip(phyICM, which(!dataICM$occurrence==1))
  # logistic data:
  data_nonzeroICM <- dataICM[dataICM$Species%in%phy_nonzeroICM$tip.label, ] # organise data to follow same order as phylogeny
  #model
  DietMod0 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod1 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp)+as.factor(Animal), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod2 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp)+as.factor(Vertebrate), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod3 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp)+as.factor(Invertebrate), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod4 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp)+as.factor(Mammal), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod5 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp)+as.factor(Bird), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod6 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp)+as.factor(Herptile), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
  DietMod7 <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp)+as.factor(Fish), data_nonzeroICM,
             cor=corPagel(0,phy_nonzeroICM, form = ~Species, fixed=F), weights = ~1/log(knownDeaths))
    X <- as.data.frame(array(,c(1,3))); names(X) <- names(sumGLS(DietMod0))
    rbind(X, sumGLS(DietMod0), X, sumGLS(DietMod1),X, sumGLS(DietMod3), X, 
               sumGLS(DietMod2), X, sumGLS(DietMod7), X, sumGLS(DietMod6), X, 
               sumGLS(DietMod5), X, sumGLS(DietMod4))
            #                                 Value    t-value p-value
            # 1                                <NA>       <NA>    <NA>
            #   (Intercept)              -7.15_(2.55)      -2.81  0.0058
            # log(BodyMass)            -0.09_(0.06)      -1.48  0.1415
            # log(lifeexp)              0.51_(0.32)       1.59  0.1153
            # ModelStats                lambda=0.36 AIC=401.49   n=127
            # 11                               <NA>       <NA>    <NA>
            #   (Intercept)1             -6.81_(2.55)      -2.67  0.0085
            # log(BodyMass)1           -0.07_(0.07)      -1.11  0.2673
            # log(lifeexp)1             0.46_(0.33)        1.4  0.1653
            # as.factor(Animal)1        0.26_(0.31)       0.85  0.3981
            # ModelStats1               lambda=0.34  AIC=403.3   n=127
            # 12                               <NA>       <NA>    <NA>
            #   (Intercept)2             -7.29_(2.57)      -2.84  0.0054
            # log(BodyMass)2           -0.11_(0.07)      -1.55  0.1235
            # log(lifeexp)2             0.54_(0.33)       1.64  0.1038
            # as.factor(Invertebrate)1 -0.14_(0.29)      -0.51  0.6138
            # ModelStats2               lambda=0.37  AIC=403.9   n=127
            # 13                               <NA>       <NA>    <NA>
            #   (Intercept)3             -6.65_(2.42)      -2.75  0.0069
            # log(BodyMass)3           -0.07_(0.06)      -1.11  0.2691
            # log(lifeexp)3             0.44_(0.31)       1.41  0.1616
            # as.factor(Vertebrate)1     0.6_(0.28)        2.1   0.038
            # ModelStats3               lambda=0.27 AIC=400.29   n=127
            # 14                               <NA>       <NA>    <NA>
            #   (Intercept)4             -7.05_(2.62)      -2.69  0.0081
            # log(BodyMass)4           -0.09_(0.06)      -1.48  0.1421
            # log(lifeexp)4              0.5_(0.33)       1.51  0.1329
            # as.factor(Fish)1          0.07_(0.45)       0.16  0.8737
            # ModelStats4               lambda=0.36 AIC=403.22   n=127
            # 15                               <NA>       <NA>    <NA>
            #   (Intercept)5             -7.13_(2.55)      -2.79   0.006
            # log(BodyMass)5           -0.09_(0.07)      -1.38  0.1692
            # log(lifeexp)5             0.51_(0.33)       1.57  0.1185
            # as.factor(Herptile)1      0.06_(0.31)       0.19  0.8528
            # ModelStats5               lambda=0.36 AIC=403.99   n=127
            # 16                               <NA>       <NA>    <NA>
            #   (Intercept)6             -7.16_(2.57)      -2.78  0.0062
            # log(BodyMass)6           -0.09_(0.07)      -1.45    0.15
            # log(lifeexp)6             0.51_(0.33)       1.57  0.1183
            # as.factor(Bird)1         -0.01_(0.38)      -0.03  0.9765
            # ModelStats6               lambda=0.37 AIC=403.57   n=127
            # 17                               <NA>       <NA>    <NA>
            #   (Intercept)7             -6.73_(2.31)      -2.92  0.0042
            # log(BodyMass)7            -0.1_(0.05)      -1.89  0.0615
            # log(lifeexp)7             0.47_(0.29)       1.62  0.1081
            # as.factor(Mammal)1         0.74_(0.3)       2.47  0.0151
            # ModelStats7               lambda=0.17 AIC=400.28   n=127
    
    
#====================================================================================================================
# Extended Data Tablex 3. - Results of composite phylogenetic models exploring variation in risk of cancer related mortality.
    # Results are presented for both a, CMR and b, ICM. Binomial phylogenetic GLMMs are presented first and phylogenetic GLSs 
    # exploring variation in non-zero cancer mortality risks are presented second. For each model sample size (n) and 
    # phylogenetic inertia (s2/λ) are presented at the bottom of the results.
    #
# Data and phylogeny for binary phylogenetic model
    dataCMR <- data[!data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus", "Dasyuroides_byrnei"),] # very high leverages
    row.names(dataCMR) <- dataCMR$Species
    phyCMR <- treedata(phy,dataCMR, warnings = FALSE)$phy
    dataCMR$occurrence <- ifelse(dataCMR$CMR==0,0,1)
    dataCMR <- dataCMR[as.character(phyCMR$tip.label), ] # organise data to follow same order as phylogeny
    # Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
    phy_nonzeroCMR <- drop.tip(phyCMR, which(!dataCMR$occurrence==1))
    data_nonzeroCMR <- dataCMR[dataCMR$Species%in%phy_nonzeroCMR$tip.label, ] # organise data to follow same order as phylogeny
    #model
    modCMR <- ZIlogis(formbin = occurrence ~ +log(knownDeaths) +log(BodyMass) + log(lifeexp) ,
                      formlogis = car::logit(CMR) ~  log(BodyMass)+log(lifeexp),
                      databin = dataCMR, phylobin=phyCMR, datalogis = data_nonzeroCMR, phyloTree = phy_nonzeroCMR, l0=1)
    # Extended Data Tablex 3a.----
    sumZILogis(modCMR)
            # _                                b (SE) z/t-value p-value
            # 1    Probability of zeros                                
            # 2             (Intercept) -11.39 (5.07)     -2.25  0.0247
            # 3        log(knownDeaths)   1.52 (0.39)      3.93   1e-04
            # 4           log(BodyMass)  -0.05 (0.12)     -0.43   0.665
            # 5            log(lifeexp)   0.88 (0.57)      1.54  0.1237
            # 6              ModelStats   P_s2= 6e-04   s2= 0.1  n= 188
            # 7  Weighted logistic reg.                                
            # 8             (Intercept)  -7.11 (1.82)     -3.92   1e-04
            # 9           log(BodyMass)  -0.04 (0.05)     -0.89  0.3754
            # 10           log(lifeexp)   0.42 (0.23)      1.83  0.0688
            # 11             ModelStats   AIC= 368.98   L= 0.53  n= 141
    
    
    # Zero-inflated COMPOSITE models
    dataICM <- data[!is.na(data$ICM)& !data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus", "Dasyuroides_byrnei"),]
    row.names(dataICM) <- dataICM$Species
    phyICM <- treedata(phy,dataICM, warnings = FALSE)$phy
    dataICM$occurrence <- ifelse(dataICM$ICM==0,0,1)
    dataICM <- dataICM[as.character(phyICM$tip.label), ] # organise data to follow same order as phylogeny
    # Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
    phy_nonzeroICM <- drop.tip(phyICM, which(!dataICM$occurrence==1))
    data_nonzeroICM <- dataICM[dataICM$Species%in%phy_nonzeroICM$tip.label, ] # organise data to follow same order as phylogeny
    #model
    modICM <- ZIlogis(formbin = occurrence ~  +log(knownDeaths)+(log(BodyMass)+log(lifeexp)),
                      formlogis = car::logit(ICM) ~  (log(BodyMass)+log(lifeexp)),
                      databin = dataICM, phylobin=phyICM, datalogis = data_nonzeroICM, phyloTree = phy_nonzeroICM, l0=0)
    # Extended Data Tablex 3b.----
    sumZILogis(modICM)
            # _                                 b (SE) z/t-value p-value
            # 1    Probability of zeros                                
            # 2             (Intercept) -10.06 (5.44)     -1.85  0.0644
            # 3        log(knownDeaths)   1.46 (0.39)      3.72   2e-04
            # 4           log(BodyMass)  -0.07 (0.13)     -0.53   0.598
            # 5            log(lifeexp)   0.73 (0.63)      1.17  0.2414
            # 6              ModelStats  P_s2= 0.0065  s2= 0.09  n= 169
            # 7  Weighted logistic reg.                                
            # 8             (Intercept)  -7.27 (2.54)     -2.87  0.0049
            # 9           log(BodyMass)  -0.09 (0.06)     -1.44  0.1534
            # 10           log(lifeexp)   0.53 (0.32)      1.63  0.1049
            # 11             ModelStats   AIC= 403.93   L= 0.37  n= 128
    
    
    
#====================================================================================================================
# Supplementary Tablex 2. - Sensitivity analysis performed to test the effect of within-species sample size. Models 
    # were performed to explore the consistency of results of weighted logistic regressions of non-zero cancer
    # risks when using a different minimum threshold of number of individuals with available postmortem pathological 
    # records. Models were run with a minimum of 20, 40, 60, 80 and 100 individuals with postmortem records per 
    # species (N). Slopes and standard errors are presented for all models. While considerable decrease in the
    # number of species involved in the analyses (n) is apparent with increased threshold sample sizes, the 
    # estimated slopes remain relatively stable.
  modsICM  <- list()
  modsCMR  <- list()
  j=1
  for(k in seq(20,100, by=20)){
    dataCMR <- data[!is.na(data$CMR) & data$knownDeaths>=k & data$CMR>0 &
                      !data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus","Dasyuroides_byrnei"),]
    row.names(dataCMR) <- dataCMR$Species
    phyCMR <- treedata(phy,dataCMR, warnings = FALSE)$phy

    for(l in seq(0,1,by=0.1)){ # try different lambda starting value if model fails to converge
        m <- try(gls( car::logit(CMR) ~  log(BodyMass)+log(lifeexp), dataCMR,
                      cor=corPagel(l, phyCMR, form=~Species),
                      weights = ~1/log(knownDeaths)))
        if(!class(m)%in%"try-error"){break}}
    if(m$modelStruct<0) {m <- gls( car::logit(CMR) ~ log(BodyMass)+log(lifeexp), dataCMR,
                                   cor=corPagel(l, phyCMR, form=~Species, fixed=T),
                                   weights = ~1/log(knownDeaths)) }
    modsCMR[[j]] <- m

    dataICM <- data[!is.na(data$ICM) & data$knownDeaths>=k & data$ICM>0 &
                      !data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus","Dasyuroides_byrnei"),]
    row.names(dataICM) <- dataICM$Species
    phyICM <- treedata(phy,dataICM, warnings = FALSE)$phy

    for(l in seq(0,1,by=0.1)){ # try different lambda starting value if model fails to converge
        m <- try(gls(logit(ICM) ~  log(BodyMass)+log(lifeexp), dataICM,
                     cor=corPagel(l, phyICM, form=~Species),
                     weights = ~1/log(knownDeaths)))
        if(!class(m)%in%"try-error"){break}}
    if(m$modelStruct<0) {m <- gls(logit(ICM) ~  log(BodyMass)+log(lifeexp), dataICM,
                                  cor=corPagel(l, phyICM, form=~Species, fixed=T),
                                  weights = ~1/log(knownDeaths)) } 
    modsICM[[j]] <- m
    j=j+1}

  sumSensit <- function(model){
    su <- as.data.frame(summary(model)$tTable)[,c(1:2,4)]
    su[,1:2] <- round(su[,1:2],2)
    su[,3] <- round(su[,3],4)
    su[,1] <- paste(su[,1],' ± ', su[,2], " (", su[,3],")", sep='')
    su <- su[,c(1)]
    su[4] <- paste("n = ", length(resid(model)), sep="")
    return(su)}
  
  #### Supplementary Tablex 2a --------------------------------------------------------------------
  resCMR <- as.data.frame(cbind(c(row.names(summary(modsCMR[[1]])$tTable), "No. species"),
        sumSensit(modsCMR[[1]]),
        sumSensit(modsCMR[[2]]),
        sumSensit(modsCMR[[3]]),
        sumSensit(modsCMR[[4]]),
        sumSensit(modsCMR[[5]])))
  names(resCMR) <- c('Explanatory', 'n>=20','n>=40','n>=60','n>=80','n>=100')
  resCMR
          # Explanatory                 n>=20                n>=40                n>=60                 n>=80                n>=100
          # 1   (Intercept)  -7.11 ± 1.82 (1e-04) -4.53 ± 2.8 (0.1106) -5.4 ± 3.76 (0.1567) -3.36 ± 4.49 (0.4589) -2.55 ± 5.01 (0.6148)
          # 2 log(BodyMass) -0.04 ± 0.05 (0.3754) 0.04 ± 0.07 (0.5242) 0.04 ± 0.09 (0.6367)   0.06 ± 0.1 (0.5793)  0.02 ± 0.11 (0.8559)
          # 3  log(lifeexp)  0.42 ± 0.23 (0.0688) 0.12 ± 0.36 (0.7491) 0.25 ± 0.48 (0.6037) -0.01 ± 0.57 (0.9916) -0.09 ± 0.64 (0.8924)
          # 4   No. species               n = 141               n = 81               n = 56                n = 42                n = 32
  
  #### Supplementary Tablex 2b --------------------------------------------------------------------
  resICM <- as.data.frame(cbind(c(row.names(summary(modsCMR[[1]])$tTable), "No. species"),
        sumSensit(modsICM[[1]]),
        sumSensit(modsICM[[2]]),
        sumSensit(modsICM[[3]]),
        sumSensit(modsICM[[4]]),
        sumSensit(modsICM[[5]])))
  names(resICM) <- c('Explanatory', 'n>=20','n>=40','n>=60','n>=80','n>=100')
  resICM
          #   Explanatory                 n>=20                n>=40                 n>=60                 n>=80                n>=100
          # 1   (Intercept) -7.27 ± 2.54 (0.0049) -4.8 ± 4.03 (0.2378) -4.04 ± 5.17 (0.4388) -3.73 ± 6.35 (0.5611)  -3.2 ± 6.61 (0.6321)
          # 2 log(BodyMass) -0.09 ± 0.06 (0.1534) -0.05 ± 0.1 (0.6025) -0.07 ± 0.12 (0.5654) -0.07 ± 0.13 (0.6208) -0.15 ± 0.15 (0.3238)
          # 3  log(lifeexp)  0.53 ± 0.32 (0.1049)  0.22 ± 0.52 (0.677)  0.13 ± 0.66 (0.8395)    0.1 ± 0.8 (0.9002)  0.08 ± 0.84 (0.9291)
          # 4   No. species               n = 128               n = 77                n = 55                n = 41                n = 32
  


#====================================================================================================================
# Supplementary Tablex 3. - Sensitivity analysis performed to test for collinearity issues. In order to demonstrate 
  # the lack of collinearity in models presented in Supplementary Table 1, here we present the results of models
  # explaining variation in cancer risk with a, c body mass and b, d life expectancy entered in separate models 
  # for both a-b CMR and c-d ICM.
  
  # Data and phylogeny for binary phylogenetic model of CMR
  dataCMR <- data[!data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus", "Dasyuroides_byrnei"),]
  row.names(dataCMR) <- dataCMR$Species
  phyCMR <- treedata(phy,dataCMR, warnings = FALSE)$phy
  dataCMR$occurrence <- ifelse(dataCMR$CMR==0,0,1)
  dataCMR <- dataCMR[as.character(phyCMR$tip.label), ] # organise data to follow same order as phylogeny
  # Data and phylogeny for linear phylogenetic model of non-zero CMR
  phy_nonzeroCMR <- drop.tip(phyCMR, which(!dataCMR$occurrence==1))
  data_nonzeroCMR <- dataCMR[dataCMR$Species%in%phy_nonzeroCMR$tip.label, ] # organise data to follow same order as phylogeny
  #model
  modCMRbm <- ZIlogis(formbin = occurrence ~  log(knownDeaths) +log(BodyMass),
                   formlogis = car::logit(CMR) ~  log(BodyMass),
                   databin = dataCMR, datalogis = data_nonzeroCMR, phyloTree = phy_nonzeroCMR, phylobin=phyCMR, l0=0)
  #### Supplementary Tablex 3a --------------------------------------------------------------------
  sumZILogis(modCMRbm)
            # _       b (SE) z/t-value p-value
            # 1   Probability of zeros                               
            # 2            (Intercept) -4.09 (1.65)     -2.48  0.0131
            # 3       log(knownDeaths)  1.41 (0.38)      3.77   2e-04
            # 4          log(BodyMass)   0.06 (0.1)      0.64  0.5201
            # 5             ModelStats      P_s2= 0  s2= 0.12  n= 188
            # 6 Weighted logistic reg.                               
            # 7            (Intercept)    -4 (0.47)     -8.55       0
            # 8          log(BodyMass)  0.01 (0.04)      0.15  0.8784
            # 9             ModelStats  AIC= 369.14   L= 0.57  n= 141
  modCMRle <- ZIlogis(formbin = occurrence ~  log(knownDeaths) +log(lifeexp),
                   formlogis = car::logit(CMR) ~  log(lifeexp),
                   databin = dataCMR, datalogis = data_nonzeroCMR, phyloTree = phy_nonzeroCMR, phylobin=phyCMR, l0=0)
  #### Supplementary Tablex 3b --------------------------------------------------------------------
  sumZILogis(modCMRle)
            # _        b (SE) z/t-value p-value
            # 1   Probability of zeros                                
            # 2            (Intercept) -10.21 (4.24)     -2.41  0.0161
            # 3       log(knownDeaths)   1.49 (0.38)      3.93   1e-04
            # 4           log(lifeexp)   0.74 (0.47)      1.59   0.113
            # 5             ModelStats   P_s2= 2e-04   s2= 0.1  n= 188
            # 6 Weighted logistic reg.                                
            # 7            (Intercept)  -6.45 (1.62)     -3.97   1e-04
            # 8           log(lifeexp)    0.32 (0.2)       1.6  0.1115
            # 9             ModelStats    AIC= 363.5   L= 0.56  n= 141
  
  # Data and phylogeny for binary phylogenetic model of ICM
  dataICM <- data[!is.na(data$ICM) & !data$Species%in%c("Lagurus_lagurus", "Cricetus_cricetus", "Dasyuroides_byrnei"),]
  row.names(dataICM) <- dataICM$Species
  phyICM <- treedata(phy,dataICM, warnings = FALSE)$phy
  dataICM$occurrence <- ifelse(dataICM$ICM==0,0,1)
  dataICM <- dataICM[as.character(phyICM$tip.label), ] # organise data to follow same order as phylogeny
  # Data and phylogeny for linear phylogenetic model of non-zero ICM
  phy_nonzeroICM <- drop.tip(phyICM, which(!dataICM$occurrence==1))
  data_nonzeroICM <- dataICM[dataICM$Species%in%phy_nonzeroICM$tip.label, ] # organise data to follow same order as phylogeny
  #model
  modICMbm <- ZIlogis(formbin = occurrence ~  log(BodyMass)+log(knownDeaths),
                   formlogis = car::logit(ICM) ~  log(BodyMass),
                   databin = dataICM, datalogis = data_nonzeroICM, phyloTree = phy_nonzeroICM, phylobin=phyICM, l0=0.3)
  #### Supplementary Tablex 3c --------------------------------------------------------------------
  sumZILogis(modICMbm)
            # _                               b (SE) z/t-value p-value
            # 1   Probability of zeros                               
            # 2            (Intercept) -4.08 (1.67)     -2.45  0.0144
            # 3          log(BodyMass)   0.03 (0.1)      0.26  0.7923
            # 4       log(knownDeaths)  1.39 (0.38)      3.62   3e-04
            # 5             ModelStats  P_s2= 5e-04  s2= 0.11  n= 169
            # 6 Weighted logistic reg.                               
            # 7            (Intercept)  -3.22 (0.5)     -6.41       0
            # 8          log(BodyMass) -0.04 (0.06)     -0.67  0.5063
            # 9             ModelStats  AIC= 404.15   L= 0.39  n= 128
  
  modICMle <- ZIlogis(formbin = occurrence ~  log(lifeexp)+log(knownDeaths),
                   formlogis = car::logit(ICM) ~  log(lifeexp),
                   databin = dataICM, datalogis = data_nonzeroICM, phyloTree = phy_nonzeroICM, phylobin=phyICM, l0=0.3)
  #### Supplementary Tablex 3d --------------------------------------------------------------------
  sumZILogis(modICMle)
            # _                                b (SE) z/t-value p-value
            # 1   Probability of zeros                               
            # 2            (Intercept)  -8.45 (4.5)     -1.88  0.0603
            # 3           log(lifeexp)  0.54 (0.51)      1.06  0.2875
            # 4       log(knownDeaths)  1.42 (0.39)       3.7   2e-04
            # 5             ModelStats  P_s2= 0.002  s2= 0.09  n= 169
            # 6 Weighted logistic reg.                               
            # 7            (Intercept) -5.62 (2.27)     -2.47  0.0148
            # 8           log(lifeexp)  0.29 (0.28)      1.03  0.3029
            # 9             ModelStats  AIC= 400.26    L= 0.4  n= 128
  

  

  
#====================================================================================================================
# Supplementary Tablex 4. - Sensitivity analysis performed to confirm the consistency of the results when high leverage
  # species are included in the analyses. In order to demonstarte the consistency of the result while including the
  # three species excluded from the analyses presented in the main text due to their high leverages (i.e. Lagurus lagurus,
  # Cricetus cricetus, and Dasyuroides byrnei), here we present models exploring the association between a, CMR and b, ICM
  # and body mass as well as life expectancy on the entire species pool.
# Data and phylogeny for binary phylogenetic model of CMR
  dataCMR <- data
  row.names(dataCMR) <- dataCMR$Species
  phyCMR <- treedata(phy,dataCMR, warnings = FALSE)$phy
  dataCMR$occurrence <- ifelse(dataCMR$CMR==0,0,1)
  dataCMR <- dataCMR[as.character(phyCMR$tip.label), ] # organise data to follow same order as phylogeny
  # Data and phylogeny for linear phylogenetic model of non-zero CMR
  phy_nonzeroCMR <- drop.tip(phyCMR, which(!dataCMR$occurrence==1))
  data_nonzeroCMR <- dataCMR[dataCMR$Species%in%phy_nonzeroCMR$tip.label, ] # organise data to follow same order as phylogeny
  #model
  modCMR <- ZIlogis(formbin = occurrence ~  +log(BodyMass) + log(lifeexp)+log(knownDeaths) ,
                   formlogis = car::logit(CMR) ~  log(BodyMass)+log(lifeexp),
                   databin = dataCMR, datalogis = data_nonzeroCMR, phyloTree = phy_nonzeroCMR, phylobin=phyCMR, l0=0)
  #### Supplementary Tablex 4a --------------------------------------------------------------------
  sumZILogis(modCMR)
            # _       b (SE) z/t-value p-value
            # 1    Probability of zeros                               
            # 2             (Intercept) -8.07 (4.71)     -1.71   0.087
            # 3           log(BodyMass) -0.03 (0.12)     -0.27  0.7892
            # 4            log(lifeexp)  0.53 (0.54)      0.97  0.3302
            # 5        log(knownDeaths)  1.44 (0.38)      3.79   2e-04
            # 6              ModelStats      P_s2= 0  s2= 0.12  n= 191
            # 7  Weighted logistic reg.                               
            # 8             (Intercept) -6.06 (1.84)      -3.3  0.0012
            # 9           log(BodyMass) -0.04 (0.05)     -0.84  0.4023
            # 10           log(lifeexp)  0.29 (0.23)      1.23   0.222
            # 11             ModelStats  AIC= 387.92   L= 0.61  n= 144
  
  # Zero-inflated COMPOSITE models
  dataICM <- data[!is.na(data$ICM),]
  row.names(dataICM) <- dataICM$Species
  phyICM <- treedata(phy,dataICM, warnings = FALSE)$phy
  # Data andő phylogeny for linear phylogenetic model of non-zero ICM
  dataICM$occurrence <- ifelse(dataICM$ICM==0,0,1)
  dataICM <- dataICM[as.character(phyICM$tip.label), ] # organise data to follow same order as phylogeny
  # Data and phylogeny for linear phylogenetic model of non-zero ICM
  phy_nonzeroICM <- drop.tip(phyICM, which(!dataICM$occurrence==1))
  data_nonzeroICM <- dataICM[dataICM$Species%in%phy_nonzeroICM$tip.label, ] # organise data to follow same order as phylogeny
  #model
  modICM <- ZIlogis(formbin = occurrence ~  (log(BodyMass)+log(lifeexp))+log(knownDeaths),
                   formlogis = car::logit(ICM) ~  (log(BodyMass)+log(lifeexp)),
                   databin = dataICM, datalogis = data_nonzeroICM, phyloTree = phy_nonzeroICM, phylobin=phyICM, l0=0)
  #### Supplementary Tablex 4b --------------------------------------------------------------------
  sumZILogis(modICM)
            # _                                b (SE) z/t-value p-value
            # 1    Probability of zeros                               
            # 2             (Intercept) -6.33 (5.01)     -1.26   0.207
            # 3           log(BodyMass) -0.05 (0.13)     -0.37  0.7126
            # 4            log(lifeexp)  0.33 (0.59)      0.55  0.5791
            # 5        log(knownDeaths)  1.37 (0.38)      3.57   4e-04
            # 6              ModelStats  P_s2= 0.001  s2= 0.11  n= 172
            # 7  Weighted logistic reg.                               
            # 8             (Intercept) -5.94 (2.41)     -2.47   0.015
            # 9           log(BodyMass)  -0.1 (0.06)     -1.53  0.1282
            # 10           log(lifeexp)  0.38 (0.31)      1.23   0.221
            # 11             ModelStats  AIC= 417.41   L= 0.39  n= 131
            # 
  
  

# Sex differences in cancer mortality -------------------------------------
  # Sex-bias in CMR across Carnivora -------------------------------------
  sexspecCM <- read.csv("Sex-specific_ICM_CMR.csv")
  CMR_M <- sexspecCM$CMR_M; names(CMR_M) <- sexspecCM$Species
  CMR_F <- sexspecCM$CMR_F; names(CMR_F) <- sexspecCM$Species
  CMR_M1 <- CMR_M[!is.na(CMR_M) & !is.na(CMR_F)]
  CMR_F1 <- CMR_F[!is.na(CMR_F) & !is.na(CMR_M)]
  tr2 <- treedata(phy,CMR_M1, warnings = F)$phy
  phyl.pairedttest(tr2,CMR_M1,CMR_F1)
            # Phylogenetic paired t-test:
            #   
            #   t = 0.520704, df = 33, p-value = 0.60605
            # 
            # alternative hypothesis:
            #   true difference in means is not equal to 0
            # 
            # 95 percent confidence interval on the phylogenetic
            # difference in mean:
            #   [-0.0290284, 0.050032]
            # 
            # estimates:
            #   phylogenetic mean difference = 0.0105018 
            # sig^2 of BM model = 0.000313725 
            # lambda (fixed or estimated) = 0.0392606 
            # 
            # log-likelihood:
            #   30.4795 
            
 # Sex-bias in ICM across Carnivora -----------------------------------------
  ICM_M <- sexspecCM$ICM_M[!is.na(sexspecCM$ICM_M) & !is.na(sexspecCM$ICM_F)]
  names(ICM_M) <- sexspecCM$Species[!is.na(sexspecCM$ICM_M) & !is.na(sexspecCM$ICM_F)]
  CMR_F <- sexspecCM$ICM_F[!is.na(sexspecCM$ICM_M) & !is.na(sexspecCM$ICM_F)]
  names(CMR_F) <- sexspecCM$Species[!is.na(sexspecCM$ICM_M) & !is.na(sexspecCM$ICM_F)]
  
  ICM_M1 <- ICM_M[!is.na(ICM_M) & !is.na(CMR_F)]
  CMR_F1 <- CMR_F[!is.na(CMR_F) & !is.na(ICM_M)]
  tr2 <- treedata(phy,ICM_M1, warnings = F)$phy
  phyl.pairedttest(tr2,ICM_M1,CMR_F1)
            # 
            # Phylogenetic paired t-test:
            #   
            #   t = -0.681466, df = 27, p-value = 0.501381
            # 
            # alternative hypothesis:
            #   true difference in means is not equal to 0
            # 
            # 95 percent confidence interval on the phylogenetic
            # difference in mean:
            #   [-0.0487073, 0.0235755]
            # 
            # estimates:
            #   phylogenetic mean difference = -0.0125659 
            # sig^2 of BM model = 0.000295832 
            # lambda (fixed or estimated) = 0 
            # 
            # log-likelihood:
            #   26.2127 
            
sessionInfo()
  # R version 4.0.4 (2021-02-15)
  # Platform: x86_64-pc-linux-gnu (64-bit)
  # Running under: Debian GNU/Linux 10 (buster)
  # 
  # Matrix products: default
  # BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
  # LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
  # 
  # locale:
  # [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8      
  # [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
  # 
  # attached base packages:
  # [1] stats     graphics  grDevices utils     datasets  methods   base     
  # 
  # other attached packages:
  # [1] gdata_2.18.0     knitr_1.34       kableExtra_1.3.4 rr2_1.0.2        emmeans_1.6.3    phylolm_2.6.2    nlme_3.1-153     phytools_0.7-80  maps_3.3.0       car_3.0-11       carData_3.0-4   
  # [12] geiger_2.0.7     ape_5.5         
  # 
  # loaded via a namespace (and not attached):
  # [1] httr_1.4.2              viridisLite_0.4.0       splines_4.0.4           tmvnsim_1.0-2           gtools_3.9.2            expm_0.999-6            cellranger_1.1.0        yaml_2.2.1             
  # [9] globals_0.14.0          numDeriv_2016.8-1.1     pillar_1.6.2            lattice_0.20-45         glue_1.4.2              quadprog_1.5-8          phangorn_2.7.1          digest_0.6.27          
  # [17] rvest_1.0.1             minqa_1.2.4             colorspace_2.0-2        htmltools_0.5.2         Matrix_1.3-2            plyr_1.8.6              pkgconfig_2.0.3         listenv_0.8.0          
  # [25] haven_2.4.3             xtable_1.8-4            mvtnorm_1.1-2           webshot_0.5.2           scales_1.1.1            svglite_2.0.0           openxlsx_4.2.4          rio_0.5.27             
  # [33] lme4_1.1-27.1           tibble_3.1.4            combinat_0.0-8          ellipsis_0.3.2          mnormt_2.0.2            magrittr_2.0.1          crayon_1.4.1            readxl_1.3.1           
  # [41] estimability_1.3        evaluate_0.14           future_1.22.1           fansi_0.5.0             parallelly_1.28.1       MASS_7.3-54             xml2_1.3.2              forcats_0.5.1          
  # [49] foreign_0.8-81          tools_4.0.4             data.table_1.14.0       hms_1.1.0               lifecycle_1.0.0         stringr_1.4.0           munsell_0.5.0           plotrix_3.8-2          
  # [57] zip_2.2.0               compiler_4.0.4          systemfonts_1.0.2       clusterGeneration_1.3.7 rlang_0.4.11            grid_4.0.4              nloptr_1.2.2.2          rstudioapi_0.13        
  # [65] subplex_1.6             igraph_1.2.6            rmarkdown_2.11          boot_1.3-27             codetools_0.2-18        abind_1.4-5             deSolve_1.29            curl_4.3.2             
  # [73] R6_2.5.1                fastmap_1.1.0           future.apply_1.8.1      utf8_1.2.2              fastmatch_1.1-3         stringi_1.7.4           parallel_4.0.4          Rcpp_1.0.7             
  # [81] vctrs_0.3.8             scatterplot3d_0.3-41    xfun_0.26               coda_0.19-4            
  

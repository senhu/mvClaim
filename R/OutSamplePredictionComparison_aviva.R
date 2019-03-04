# load(naviva_data.Rdata)

#------------------------------------------------
# mainly comparing method 2 and 4
# since they allow to model each part (including covariance part) separately

dim(naviva)[1]*.8
train.index <- sample(c(1:dim(naviva)[1]), 1659, replace = FALSE)
naviva.train <- naviva[train.index,]; dim(naviva.train)
naviva.test  <- naviva[setdiff(1:dim(naviva)[1], train.index),]; dim(naviva.test)

#----------------------------------------------

#unimod1 <- glm(COST_AD  ~ VEH_FUEL_TYPE + VEH_TRANS + NCDLAST.factor + NCD_PROT + STEPBACK, family = Gamma(link="log"), data=naviva.train)
unimod1 <- glm(COST_AD  ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE, family = Gamma(link="log"), data=naviva.train)
pred1 <- predict(unimod1, newdata = naviva.test, type="response")
pred1 <- predict(unimod1, newdata = naviva.test, type="response", dispersion=MASS::gamma.dispersion(unimod1))

unimod2 <- glm(COST_TPD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE, family = Gamma(link="log"), data=naviva.train)
pred2 <- predict(unimod2, newdata = naviva.test, type="response")
pred2 <- predict(unimod2, newdata = naviva.test, type="response", dispersion=MASS::gamma.dispersion(unimod2))

plot(cbind(naviva.test$COST_AD, naviva.test$COST_TPD), main="univariate GLM predictions", pch=1, cex=.2, xlab="AD", ylab="PD")
points(pred1, pred2, pch=20, col=2)

ComparisonMetrics(naviva.test$COST_AD, pred1)
ComparisonMetrics(naviva.test$COST_TPD, pred2)

#----------------------------------------------

m1 <-  MBGR2(data   = naviva.train,
             G      = 1,
             l1     = COST_AD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
             l2     = COST_TPD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
             l3     = ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
             lp     = NULL,
             expo   = NULL,
             maxit  = 300,
             tol    = 1e-5,
             start  ="mclust.start",
             Aitken = FALSE,
             verbose= TRUE)
m2 <- MBGR2(data   = naviva.train,
            G      = 2,
            l1     = COST_AD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l2     = COST_TPD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l3     = ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            lp     = NULL, 
            expo   = NULL,
            maxit  = 300,
            tol    = 1e-5,
            start  ="mclust.start",
            Aitken = FALSE,
            verbose= TRUE)
m3 <- MBGR2(data   = naviva.train,
            G      = 2,
            l1     = COST_AD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l2     = COST_TPD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l3     = ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            lp     = ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            expo   = NULL,
            maxit  = 300,
            tol    = 1e-5,
            start  ="mclust.start",
            Aitken = FALSE,
            verbose= TRUE)
m4 <- MBGR2(data   = naviva.train,
            G      = 3,
            l1     = COST_AD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l2     = COST_TPD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l3     = ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            lp     = NULL, 
            expo   = NULL,
            maxit  = 300,
            tol    = 1e-5,
            start  ="mclust.start",
            Aitken = FALSE,
            verbose= TRUE)
m5 <- MBGR2(data   = naviva,
            G      = 3,
            l1     = COST_AD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l2     = COST_TPD ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            l3     = ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            lp     = ~ NCDLAST + NCD_PROT + NUMDRIVS + MILEAGE + MD_LIC_CAT_A + PH_PENPTS + VEH_TRANS + VEH_FUEL_TYPE,
            expo   = NULL,
            maxit  = 300,
            tol    = 1e-5,
            start  ="mclust.start",
            Aitken = FALSE,
            verbose= TRUE)

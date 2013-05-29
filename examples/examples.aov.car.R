
# exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.

data(obk.long, package = "afex")

# run univariate mixed ANOVA for the full design:
aov.car(value ~ treatment * gender + Error(id/phase*hour), 
        data = obk.long, observed = "gender")

ez.glm("id", "value", obk.long, between = c("treatment", "gender"), 
        within = c("phase", "hour"), observed = "gender")

# both calls return the same:
##                         Effect          df   MSE         F  ges     p
## 1                    treatment       2, 10 22.81    3.94 +  .20   .05
## 2                       gender       1, 10 22.81    3.66 +  .11   .08
## 3             treatment:gender       2, 10 22.81      2.86  .18   .10
## 4                        phase 1.60, 15.99  5.02 16.13 ***  .15 <.001
## 5              treatment:phase 3.20, 15.99  5.02    4.85 *  .10   .01
## 6                 gender:phase 1.60, 15.99  5.02      0.28 .003   .71
## 7       treatment:gender:phase 3.20, 15.99  5.02      0.64  .01   .61
## 8                         hour 1.84, 18.41  3.39 16.69 ***  .13 <.001
## 9               treatment:hour 3.68, 18.41  3.39      0.09 .002   .98
## 10                 gender:hour 1.84, 18.41  3.39      0.45 .004   .63
## 11       treatment:gender:hour 3.68, 18.41  3.39      0.62  .01   .64
## 12                  phase:hour 3.60, 35.96  2.67      1.18  .02   .33
## 13        treatment:phase:hour 7.19, 35.96  2.67      0.35 .009   .93
## 14           gender:phase:hour 3.60, 35.96  2.67      0.93  .01   .45
## 15 treatment:gender:phase:hour 7.19, 35.96  2.67      0.74  .02   .65


# replicating ?Anova using aov.car:
aov.car(value ~ treatment * gender + Error(id/phase*hour), 
        data = obk.long, type = 2, return = "Anova")
# in contrast to aov you do not need the within-subject factors outside Error()

# replicating ?Anova using ez.glm:
ez.glm("id", "value", obk.long, c("treatment", "gender"), 
        c("phase", "hour"), type = 2, return = "Anova")

#both return:
## Type II Repeated Measures MANOVA Tests: Pillai test statistic
##                             Df test stat approx F num Df den Df       Pr(>F)    
## (Intercept)                  1     0.970      318      1     10 0.0000000065 ***
## treatment                    2     0.481        5      2     10      0.03769 *  
## gender                       1     0.204        3      1     10      0.14097    
## treatment:gender             2     0.364        3      2     10      0.10447    
## phase                        1     0.851       26      2      9      0.00019 ***
## treatment:phase              2     0.685        3      4     20      0.06674 .  
## gender:phase                 1     0.043        0      2      9      0.82000    
## treatment:gender:phase       2     0.311        1      4     20      0.47215    
## hour                         1     0.935       25      4      7      0.00030 ***
## treatment:hour               2     0.301        0      8     16      0.92952    
## gender:hour                  1     0.293        1      4      7      0.60237    
## treatment:gender:hour        2     0.570        1      8     16      0.61319    
## phase:hour                   1     0.550        0      8      3      0.83245    
## treatment:phase:hour         2     0.664        0     16      8      0.99144    
## gender:phase:hour            1     0.695        1      8      3      0.62021    
## treatment:gender:phase:hour  2     0.793        0     16      8      0.97237    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

# ANCOVA: adding a covariate (necessary to set factorize = FALSE)
aov.car(value ~ treatment * gender + age + Error(id/phase*hour), 
        data = obk.long, observed = c("gender", "age"), factorize = FALSE)

ez.glm("id", "value", obk.long, between = c("treatment", "gender"), 
        within = c("phase", "hour"), covariate = "age", 
        observed = c("gender", "age"), factorize = FALSE)

# aggregating over one within-subjects factor (phase) with warning:
aov.car(value ~ treatment * gender + Error(id/hour), data = obk.long, observed = "gender")

ez.glm("id", "value", obk.long, c("treatment", "gender"), "hour", observed = "gender")

# runs with "numeric" factors
obk.long$hour2 <- as.numeric(as.character(obk.long$hour))

aov.car(value ~ treatment * gender + Error(id/hour2), 
        data = obk.long, type = 2,observed = c("gender"))

# only between
aov.car(value ~ treatment * gender + Error(id), 
        data = obk.long, type = 2,observed = c("gender"))
aov.car(value ~ treatment * gender + Error(id), 
        data = obk.long, type = 2, observed = c("gender"))

ez.glm("id", "value", obk.long, c("treatment", "gender"), 
        within = NULL, type = 2, print.formula = TRUE, observed = "gender")

# only within

aov.car(value ~ Error(id/phase*hour), data = obk.long, type = 2)

ez.glm("id", "value", obk.long,  NULL, c("phase", "hour"), 
        type = 2, print.formula = TRUE)

# using return = "full":

str(aov.car(value ~ Error(id/phase*hour), data = obk.long, return = "full"), 1)

## List of 4
##  $ Anova:List of 14
##   ..- attr(*, "class")= chr "Anova.mlm"
##  $ lm   :List of 11
##   ..- attr(*, "class")= chr [1:2] "mlm" "lm"
##  $ data :'data.frame':  16 obs. of  16 variables:
##  $ idata:'data.frame':  15 obs. of  2 variables:

# use args.return arguments:
aov.car(value ~ treatment * gender + Error(id/phase*hour), 
        data = obk.long, args.return = list(correction = "none", es = "pes"))

aov.car(value ~ treatment * gender + Error(id/phase*hour), 
        data = obk.long,observed = "gender", 
        args.return = list(correction = "none", MSE = FALSE))


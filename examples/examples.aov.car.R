
# exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.

data(obk.long, package = "afex")

# run univariate mixed ANCOVA for the full design:
aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long, observed = c("gender", "age"))

ez.glm("id", "value", obk.long, between = c("treatment", "gender"), within = c("phase", "hour"), covariate = "age", observed = c("gender", "age"))

# both calls return the same:
##                         Effect          df   MSE         F  ges     p
## 1                    treatment        2, 9 23.96    3.58 +  .20   .07
## 2                       gender        1, 9 23.96    3.95 +  .14   .08
## 3                          age        1, 9 23.96      0.52  .02   .49
## 4             treatment:gender        2, 9 23.96      1.28  .09   .32
## 5                        phase  1.7, 15.28  3.91 20.28 ***  .17 <.001
## 6              treatment:phase 3.39, 15.28  3.91   6.07 **  .11  .005
## 7                 gender:phase  1.7, 15.28  3.91      0.25 .002   .75
## 8                    age:phase  1.7, 15.28  3.91    3.10 +  .03   .08
## 9       treatment:gender:phase 3.39, 15.28  3.91      1.60  .03   .23
## 10                        hour 2.14, 19.23  2.48 20.52 ***  .14 <.001
## 11              treatment:hour 4.27, 19.23  2.48      0.71  .01   .60
## 12                 gender:hour 2.14, 19.23  2.48      0.71 .006   .51
## 13                    age:hour 2.14, 19.23  2.48    2.82 +  .02   .08
## 14       treatment:gender:hour 4.27, 19.23  2.48      0.59 .009   .68
## 15                  phase:hour 3.48, 31.36  2.83      0.99  .01   .42
## 16        treatment:phase:hour 6.97, 31.36  2.83      0.33 .010   .93
## 17           gender:phase:hour 3.48, 31.36  2.83      0.90  .01   .47
## 18              age:phase:hour 3.48, 31.36  2.83      0.77  .01   .54
## 19 treatment:gender:phase:hour 6.97, 31.36  2.83      0.65  .02   .71

# replicating ?Anova using aov.car:
aov.car(value ~ treatment * gender + Error(id/phase*hour), data = obk.long, type = 2, return = "Anova")
# in contrast to aov you do not need the within-subject factors outside Error()

# replicating ?Anova using ez.glm:
ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), type = 2, return = "Anova")

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


# aggregating over one within-subjects factor (phase) with warning:

aov.car(value ~ treatment * gender + age + Error(id/hour), data = obk.long, observed = c("gender", "age"))

ez.glm("id", "value", obk.long, c("treatment", "gender"), "hour", "age", observed = c("gender", "age"))


# runs with "numeric" factors
obk.long$hour2 <- as.numeric(as.character(obk.long$hour))

aov.car(value ~ treatment * gender + Error(id/hour2), data = obk.long, type = 2,observed = c("gender"))

# only between
aov.car(value ~ treatment * gender + age + Error(id), data = obk.long, type = 2,observed = c("gender", "age"))
aov.car(value ~ treatment * gender + Error(id), data = obk.long, type = 2, observed = c("gender"))

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, covariate = "age", type = 2, print.formula = TRUE, observed = c("gender", "age"))

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, type = 2, print.formula = TRUE, observed = c("gender"))

# only within

aov.car(value ~ Error(id/phase*hour), data = obk.long, type = 2)

ez.glm("id", "value", obk.long,  NULL, c("phase", "hour"), type = 2, print.formula = TRUE)

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
aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long, args.return = list(correction = "none", es = "pes"))

aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long,observed = c("gender", "age"), args.return = list(correction = "none", MSE = FALSE))


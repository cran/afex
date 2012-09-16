
# exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.

data(obk.long, package = "afex")

# run univariate mixed ANCOVA for the full design:
univ(aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long))
univ(ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), "age"))

# both calls return the same:

## $anova
##                                      SS num Df  Error SS den Df           F       Pr(>F)
## (Intercept)                 6454.236987      1 215.65658      9 269.3547893 5.152317e-08
## treatment                    171.399953      2 215.65658      9   3.5765187 7.193619e-02
## gender                        94.598340      1 215.65658      9   3.9478742 7.818280e-02
## age                           12.398975      1 215.65658      9   0.5174466 4.901885e-01
## treatment:gender              61.531858      2 215.65658      9   1.2839551 3.231798e-01
## phase                        134.586005      2  59.72439     18  20.2810632 2.448505e-05
## treatment:phase               80.604542      4  59.72439     18   6.0732385 2.826803e-03
## gender:phase                   1.634246      2  59.72439     18   0.2462681 7.843036e-01
## age:phase                     20.553392      2  59.72439     18   3.0972362 6.982439e-02
## treatment:gender:phase        21.254421      4  59.72439     18   1.6014379 2.170946e-01
## hour                         108.513510      4  47.59543     36  20.5192290 7.001584e-09
## treatment:hour                 7.547869      8  47.59543     36   0.7136275 6.779072e-01
## gender:hour                    3.746135      4  47.59543     36   0.7083708 5.915285e-01
## age:hour                      14.904567      4  47.59543     36   2.8183608 3.926421e-02
## treatment:gender:hour          6.235198      8  47.59543     36   0.5895186 7.798264e-01
## phase:hour                     9.762579      8  88.62706     72   0.9913814 4.501348e-01
## treatment:phase:hour           6.579092     16  88.62706     72   0.3340505 9.915014e-01
## gender:phase:hour              8.851396      8  88.62706     72   0.8988515 5.222336e-01
## age:phase:hour                 7.539611      8  88.62706     72   0.7656409 6.339004e-01
## treatment:gender:phase:hour   12.822199     16  88.62706     72   0.6510416 8.307936e-01
## 
## $mauchly
##                             Test statistic    p-value
## phase                         0.8217571566 0.45600959
## treatment:phase               0.8217571566 0.45600959
## gender:phase                  0.8217571566 0.45600959
## age:phase                     0.8217571566 0.45600959
## treatment:gender:phase        0.8217571566 0.45600959
## hour                          0.0966749877 0.04923980
## treatment:hour                0.0966749877 0.04923980
## gender:hour                   0.0966749877 0.04923980
## age:hour                      0.0966749877 0.04923980
## treatment:gender:hour         0.0966749877 0.04923980
## phase:hour                    0.0002379741 0.08651564
## treatment:phase:hour          0.0002379741 0.08651564
## gender:phase:hour             0.0002379741 0.08651564
## age:phase:hour                0.0002379741 0.08651564
## treatment:gender:phase:hour   0.0002379741 0.08651564
## 
## $sphericity.correction
##                                GG eps   Pr(>F[GG])    HF eps   Pr(>F[HF])
## phase                       0.8487215 8.383485e-05 1.0252867 2.448505e-05
## treatment:phase             0.8487215 5.159591e-03 1.0252867 2.826803e-03
## gender:phase                0.8487215 7.493990e-01 1.0252867 7.843036e-01
## age:phase                   0.8487215 8.073373e-02 1.0252867 6.982439e-02
## treatment:gender:phase      0.8487215 2.279698e-01 1.0252867 2.170946e-01
## hour                        0.5341747 1.302016e-05 0.7054545 8.046331e-07
## treatment:hour              0.5341747 6.010781e-01 0.7054545 6.342676e-01
## gender:hour                 0.5341747 5.137213e-01 0.7054545 5.478398e-01
## age:hour                    0.5341747 8.155027e-02 0.7054545 6.211130e-02
## treatment:gender:hour       0.5341747 6.843526e-01 0.7054545 7.263729e-01
## phase:hour                  0.4355822 4.186799e-01 0.7444364 4.402119e-01
## treatment:phase:hour        0.4355822 9.317848e-01 0.7444364 9.787985e-01
## gender:phase:hour           0.4355822 4.651930e-01 0.7444364 5.020890e-01
## age:phase:hour              0.4355822 5.395151e-01 0.7444364 5.992844e-01
## treatment:gender:phase:hour 0.4355822 7.100921e-01 0.7444364 7.878433e-01
## 
## Warning message:
## In univariate(aov.car(value ~ treatment * gender + age + Error(id/phase *  :
##   HF eps > 1 treated as 1

# To get a nicer ANOVA table use function nice.anova (see ?nice.anova):
nice.anova(ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), "age"))

##                         Effect          df   MSE         F     p
## 1                    treatment        2, 9 23.96    3.58 +   .07
## 2                       gender        1, 9 23.96    3.95 +   .08
## 3                          age        1, 9 23.96      0.52   .49
## 4             treatment:gender        2, 9 23.96      1.28   .32
## 5                        phase  1.7, 15.28  3.91 20.28 *** <.001
## 6              treatment:phase 3.39, 15.28  3.91   6.07 **  .005
## 7                 gender:phase  1.7, 15.28  3.91      0.25   .75
## 8                    age:phase  1.7, 15.28  3.91    3.10 +   .08
## 9       treatment:gender:phase 3.39, 15.28  3.91      1.60   .23
## 10                        hour 2.14, 19.23  2.48 20.52 *** <.001
## 11              treatment:hour 4.27, 19.23  2.48      0.71   .60
## 12                 gender:hour 2.14, 19.23  2.48      0.71   .51
## 13                    age:hour 2.14, 19.23  2.48    2.82 +   .08
## 14       treatment:gender:hour 4.27, 19.23  2.48      0.59   .68
## 15                  phase:hour 3.48, 31.36  2.83      0.99   .42
## 16        treatment:phase:hour 6.97, 31.36  2.83      0.33   .93
## 17           gender:phase:hour 3.48, 31.36  2.83      0.90   .47
## 18              age:phase:hour 3.48, 31.36  2.83      0.77   .54
## 19 treatment:gender:phase:hour 6.97, 31.36  2.83      0.65   .71


# replicating ?Anova using aov.car:
aov.car(value ~ treatment * gender + Error(id/phase*hour), data = obk.long, type = 2)
# in contrast to aov you do not need the within-subject factors outside Error()

# replicating ?Anova using ez.glm:
ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), type = 2)

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

aov.car(value ~ treatment * gender + age + Error(id/hour), data = obk.long)

ez.glm("id", "value", obk.long, c("treatment", "gender"), "hour", "age")


# runs with "numeric" factors
obk.long$hour2 <- as.numeric(as.character(obk.long$hour))

aov.car(value ~ treatment * gender + Error(id/hour2), data = obk.long, type = 2)

# only between
aov.car(value ~ treatment * gender + age + Error(id), data = obk.long, type = 2)
aov.car(value ~ treatment * gender + Error(id), data = obk.long, type = 2)

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, covariate = "age", type = 2, print.formula = TRUE)

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, type = 2, print.formula = TRUE)

# only within

univ(aov.car(value ~ Error(id/phase*hour), data = obk.long, type = 2))

univ(ez.glm("id", "value", obk.long,  NULL, c("phase", "hour"), type = 2, print.formula = TRUE))

# using the return argument:

str(aov.car(value ~ Error(id/phase*hour), data = obk.long, return = ""), 1)

## List of 4
##  $ Anova:List of 14
##   ..- attr(*, "class")= chr "Anova.mlm"
##  $ lm   :List of 11
##   ..- attr(*, "class")= chr [1:2] "mlm" "lm"
##  $ data :'data.frame':  16 obs. of  16 variables:
##  $ idata:'data.frame':  15 obs. of  2 variables:


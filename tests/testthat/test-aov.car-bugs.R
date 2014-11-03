
context("ANOVAs: known bugs")

test_that("another label bug (May 2014)", {
  load("../../../../bugs/another-label-error.rda")
  expect_is(suppressWarnings(ez.glm("id", "corr", d.in, between = "cond", within = c("type", "inference"), return = "Anova")), "Anova.mlm")  
})

test_that("orig label bug", {
  data(obk.long)
  obk2 <- obk.long
  levels(obk2$phase) <- c("fup test", "post-hans", "pre tenetious")
  expect_is(suppressWarnings(aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk2, factorize=FALSE, return = "Anova")), "Anova.mlm")
})

test_that("ANCOVA check bug (reported by Gang Chen), January 2013", {
  dat <- read.table('../../../../bugs/mydata.txt', header=T)
  dat$ID <- as.factor(dat$ID)
  fm <- suppressWarnings(aov.car(Value ~ Propdd00 + Group + Gender + GAS0 + MAD0 + CPD0 + Error(ID/ROI), data=dat, factorize=FALSE, return = "Anova"))
  fm0 <- suppressWarnings(aov.car(Value ~ MAD0 + CPD0 + Error(ID/ROI), data=dat, factorize=FALSE, return='full'))
  expect_is(fm, "Anova.mlm")
  expect_is(fm0, "list")
})


test_that("ANOVA: ids in multiple between.subjects conditions", {
  species<- c("a","b","c","c","b","c","b","b","a","b","c","c","a","a","b","b","a","a","b","c")
  habitat<-  c("x","x","x","y","y","y","x","x","y","z","y","y","z","z","x","x","y","y","z","z")
  mvt.rate<-c(6,5,7,8,9,4,3,5,6,9,3,6,6,7,8,9,5,6,7,8)
  ind<-as.factor(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))
  data1<-data.frame(species, habitat, mvt.rate, ind)
  # should give an error
  expect_error(ez.glm("ind", "mvt.rate", data1, within = "habitat", between = "species"), "Following ids are in more than one between subjects condition:")
})

test_that("empty factors are not causing aov.cat to choke", {
  data(sleepstudy) #Example data in lme4
  sleepstudy$Days<-factor(sleepstudy$Days)
  #Works with all factors
  expect_is(ez.glm("Subject","Reaction",sleepstudy, within="Days", return = "Anova"), "Anova.mlm")
  #If you remove a factor it fails...
  expect_is(ez.glm("Subject","Reaction",sleepstudy[sleepstudy$Days!=9,], within="Days", return = "Anova"), "Anova.mlm")
})

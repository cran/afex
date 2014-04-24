
context("Mixed: structural tests")

# note: all calls with type 2 are wrapped in suppressWarnings()!

test_that("mixed: Maxell & Delaney (2004), Table 16.4, p. 842: Type 2", {
  data(md_16.4)
  md_16.4b <- md_16.4
  md_16.4b$cog <- scale(md_16.4b$cog, scale=FALSE)
  suppressWarnings(mixed4_2 <- mixed(induct ~ cond*cog + (cog|room:cond), md_16.4b, type = 2, progress=FALSE))
  lmer4_full <- lmer(induct ~ cond*cog + (cog|room:cond), md_16.4b)
  lmer4_small <- lmer(induct ~ cond+cog + (cog|room:cond), md_16.4b)
  expect_that(fixef(mixed4_2$full.model[[2]]), equals(fixef(lmer4_full)))
  expect_that(fixef(mixed4_2$full.model[[1]]), is_equivalent_to(fixef(lmer4_small)))  
})

test_that("mixed, obk.long: type 2 and LRTs", {
  data(obk.long, package = "afex")
  suppressWarnings(t2 <- mixed(value ~ treatment*phase +(1|id), data = obk.long, method = "LRT", type = 2, progress=FALSE))
  a2.f <- lmer(value ~ treatment*phase +(1|id), data = obk.long, REML=FALSE)
  a2.h <- lmer(value ~ treatment+phase +(1|id), data = obk.long, REML=FALSE)
  a2.t <- lmer(value ~ treatment +(1|id), data = obk.long, REML=FALSE)
  a2.p <- lmer(value ~ phase +(1|id), data = obk.long, REML=FALSE)
  
  extract_anova <- function(anova) unlist(anova)[c("Df1", "Df2", "Chisq2", "Chi Df2", "Pr(>Chisq)2" )]
  
  expect_that(
    unlist(t2$anova.table[3,c(3, 2, 4:6)])
    , is_equivalent_to(
      extract_anova(anova(a2.h, a2.f))
    ))
  expect_that(
    unlist(t2$anova.table[2,c(3, 2, 4:6)])
    , is_equivalent_to(
      extract_anova(anova(a2.t, a2.h))
    ))
  expect_that(
    unlist(t2$anova.table[1,c(3, 2, 4:6)])
    , is_equivalent_to(
      extract_anova(anova(a2.p, a2.h))
    ))
})

test_that("mixed, mlmRev: type 3 and 2 LRTs; print.mixed and propagate warnings", {
  require("mlmRev")
  suppressWarnings(gm1 <- mixed(use ~ age*urban + (1 | district), family = binomial, data = Contraception, method = "LRT", progress=FALSE))
  suppressWarnings(gm2 <- mixed(use ~ age*urban + (1 | district), family = binomial, data = Contraception, method = "LRT", type = 2, progress=FALSE))
  expect_that(print(gm1), gives_warning("."))
  expect_that(print(gm2), gives_warning("."))  
})

test_that("mixed, obk.long: LMM with method = PB", {
  expect_that(mixed(value ~ treatment+phase*hour +(1|id), data = obk.long, method = "PB", args.test = list(nsim = 10), progress=FALSE), is_a("mixed"))
})

test_that("mixed, obk.long: multicore loads lme4", {
  data(obk.long, package = "afex")
  require(parallel)
  cl <- makeCluster(rep("localhost", 2)) # make cluster
  # 1. Obtain fits with multicore:
  m_mc1 <- mixed(value ~ treatment +(1|id), data = obk.long, method = "LRT", cl = cl, control = lmerControl(optCtrl=list(maxfun = 100000)), progress=FALSE)
  cl_search <- clusterEvalQ(cl, search())
  stopCluster(cl)  
  expect_that(all(vapply(cl_search, function(x) any(grepl("^package:base$", x)), NA)), is_true())
})
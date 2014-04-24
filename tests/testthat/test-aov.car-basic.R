
context("ANOVAs: replicating Maxwell & Delaney (2004, p.578)")

test_that("purely within ANOVA, return='univ': Maxell & Delaney (2004), Table 12.5 and 12.6, p. 578", {
  ### replicate results from Table 12.6
  data(md_12.1)
  # valus from table:
  f <- c(40.72, 33.77, 45.31)
  ss_num <- c(289920, 285660, 105120)
  ss_error <- c(64080, 76140, 20880)
  num_df <- c(2, 1, 2)
  den_df <- c(18, 9, 18)
  
  suppressWarnings(md_ez_r <- ez.glm("id", "rt", md_12.1, within = c("angle", "noise"), return = "univ"))
  suppressWarnings(md_car_r <- aov.car(rt ~ 1 + Error(id/angle*noise), md_12.1, return = "univ"))
  suppressWarnings(md_aov4_r <- aov4(rt ~ 1 + (angle*noise|id), md_12.1, return = "univ"))
  
  expect_that(md_ez_r, is_equivalent_to(md_car_r))
  expect_that(md_ez_r, is_equivalent_to(md_aov4_r))
  expect_that(round(md_ez_r$anova[,"F"][-1], 2), is_equivalent_to(f))
  expect_that(md_ez_r$anova[,"SS"][-1], is_equivalent_to(ss_num))  
  expect_that(md_ez_r$anova[,"Error SS"][-1], is_equivalent_to(ss_error))
  expect_that(md_ez_r$anova[,"num Df"][-1], is_equivalent_to(num_df))
  expect_that(md_ez_r$anova[,"den Df"][-1], is_equivalent_to(den_df))
})

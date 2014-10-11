
context("ANOVAs: structural tests")

test_that("dv is numeric", {
  data(obk.long)
  expect_that(aov.car(treatment ~ gender + Error(id/phase*hour), data = obk.long, observed = "gender"), throws_error("dv needs to be numeric."))
})

test_that("non Type 3 sums give warning", {
  data(obk.long)
  expect_that(aov4(value ~ treatment * gender + (phase*hour|id), data = obk.long, observed = "gender", check.contrasts = FALSE), gives_warning("contrasts"))
})

test_that("return='aov' works", {
  data(obk.long)
  data(md_12.1)
  
  # purely within
  expect_that(ez.glm("id", "rt", md_12.1, within = c("angle", "noise"), return = "aov"), is_a(c( "aovlist", "listof" )))
  expect_that(aov.car(value ~ Error(id/phase*hour), data = obk.long, return = "aov"), is_a(c( "aovlist", "listof" )))
  #purely between
  expect_that(suppressWarnings(aov.car(value ~ treatment * gender + Error(id), data = obk.long, return = "aov")), is_a(c( "aov")))
  expect_that(suppressWarnings(aov.car(value~treatment * gender + Error(id/phase*hour), data = obk.long, return = "aov")), is_a(c( "aovlist", "listof" )))
})


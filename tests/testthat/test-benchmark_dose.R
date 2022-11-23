test_that("benchmark dosing works", {
  expect_s3_class(mod1,"semibmd")
  expect_s3_class(mod2,"semibmd")
  expect_s3_class(mod3,c('simpleError','error','condition'))

  # Numbers found by setting seed as set.seed(3798) in the setup
  expect_lt(abs(get_bmd(mod1)[1] - 0.07581173),1e-06)
  expect_lt(abs(get_bmd(mod1)[2] - 0.06417901),1e-06)
  expect_lt(abs(get_bmd(mod2)[1] - 0.07636082),1e-06)
  expect_lt(abs(get_bmd(mod2)[2] - 0.06439774),1e-06)

})

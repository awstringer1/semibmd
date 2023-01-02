test_that("benchmark dosing works", {
  expect_s3_class(mod1,"semibmd")
  expect_s3_class(mod2,"semibmd")
  expect_s3_class(mod3,c('simpleError','error','condition'))

  expect_s3_class(mod1v,"semibmd")
  expect_s3_class(mod2v,"semibmd")
  expect_s3_class(mod3v,c('simpleError','error','condition'))

  expect_s3_class(mod1g,"semibmd")
  expect_s3_class(mod2g,"semibmd")
  expect_s3_class(mod3g,c('simpleError','error','condition'))


  # Numbers found by setting seed as set.seed(3798) in the setup
  expect_lt(abs(get_bmd(mod1)[1] - 0.07581173),1e-06)
  expect_lt(abs(get_bmd(mod1)[2] - 0.06417901),1e-06)
  expect_lt(abs(get_bmd(mod2)[1] - 0.07636082),1e-06)
  expect_lt(abs(get_bmd(mod2)[2] - 0.06439774),1e-06)

  # Verbose option should not affect results
  expect_lt(abs(get_bmd(mod1)[1] - get_bmd(mod1v)[1]),1e-12)
  expect_lt(abs(get_bmd(mod1)[2] - get_bmd(mod1v)[2]),1e-12)
  expect_lt(abs(get_bmd(mod2)[1] - get_bmd(mod2v)[1]),1e-12)
  expect_lt(abs(get_bmd(mod2)[2] - get_bmd(mod2v)[2]),1e-12)

  # Results from non-monotone should be close-ish
  expect_lt(abs(get_bmd(mod1)[1] - get_bmd(mod1g)[1]),1e-03)
  expect_lt(abs(get_bmd(mod1)[2] - get_bmd(mod1g)[2]),1e-02)
  expect_lt(abs(get_bmd(mod2)[1] - get_bmd(mod2g)[1]),1e-03)
  expect_lt(abs(get_bmd(mod2)[2] - get_bmd(mod2g)[2]),1e-02)
})

test_that("benchmark dosing works with TMB", {
  expect_s3_class(mod1_tmb,"semibmd")


  # Numbers found by setting seed as set.seed(3798) in the setup
  expect_lt(abs(get_bmd(mod1_tmb)[1] - 0.06900426),1e-06)
  expect_lt(abs(get_bmd(mod1_tmb)[2] - 0.03918665),1e-03)
  
  expect_lt(abs(get_bmd(mod2_tmb)[1] - 0.02598749),1e-06)
  expect_lt(abs(get_bmd(mod2_tmb)[2] - 0.01418087),1e-06)
  

  # Check BMDL calculation
  expect_length(mod1_tmb$info$bmdl_alternatives,2)
  expect_named(get_all_bmdl(mod1_tmb),c('score','delta','bmdl_bayes'))
  expect_type(get_all_bmdl(mod1),'double')

  expect_length(mod2_tmb$info$bmdl_alternatives,2)
  expect_named(get_all_bmdl(mod2_tmb),c('score','delta','bmdl_bayes'))
  expect_type(get_all_bmdl(mod2),'double')

  # Computation times
  expect_named(get_computation_times(mod1_tmb),c('model','posterior_samples','bmd_samples','bmd_estimate','bmdl_delta','bmdl_score','plot_information','total'))
  expect_length(get_computation_times(mod1_tmb),8)

  expect_named(get_computation_times(mod2_tmb),c('model','posterior_samples','bmd_samples','bmd_estimate','bmdl_delta','bmdl_score','plot_information','total'))
  expect_length(get_computation_times(mod2_tmb),8)

})

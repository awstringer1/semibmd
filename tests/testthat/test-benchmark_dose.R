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

  # Check BMDL calculation
  expect_length(mod1$info$bmdl_alternatives,1)
  expect_length(mod2$info$bmdl_alternatives,1)
  expect_named(get_all_bmdl(mod1),c('score','delta'))
  expect_named(get_all_bmdl(mod2),c('score','delta'))
  expect_type(get_all_bmdl(mod1),'double')
  expect_type(get_all_bmdl(mod2),'double')

  expect_length(get_bmd(mod1_delta),1)
  expect_length(get_bmd(mod1_none),1)
  expect_length(get_bmd(mod1_score),2)
  expect_equal(get_bmd(mod1)[1],get_bmd(mod1_delta))
  expect_equal(get_bmd(mod1)[1],get_bmd(mod1_none))
  expect_equal(get_bmd(mod1)[1],get_bmd(mod1_score)[1])

  expect_length(get_all_bmdl(mod1),2)
  expect_length(get_all_bmdl(mod1_none),0)
  expect_length(get_all_bmdl(mod1_score),1)
  expect_length(get_all_bmdl(mod1_delta),1)

  # With bootstrap
  expect_length(mod1_boot$info$bmdl_alternatives,2)
  expect_named(get_all_bmdl(mod1_boot),c('score','delta','bootstrap'))
  expect_type(get_all_bmdl(mod1_boot),'double')
  expect_length(get_all_bmdl(mod1_boot),3)

  expect_true(all(mod1_boot$info$bootstrapinfo$bootstrapsuccess))
  expect_equal(unname(get_all_bmdl(mod1_boot)['bootstrap']),unname(quantile(mod1_boot$info$bootstrapinfo$bootstrapvals,.025)))

  expect_length(mod1_bayes$info$bmdl_alternatives,2)
  expect_named(get_all_bmdl(mod1_bayes),c('score','delta','bmdl_bayes'))
  expect_type(get_all_bmdl(mod1_bayes),'double')
  expect_length(get_all_bmdl(mod1_bayes),3)

  expect_true(all(mod1_bayes$info$bootstrapinfo$bayes_boot_vals > 0))
  expect_equal(unname(get_all_bmdl(mod1_bayes)['bayes_boot']),unname(quantile(mod1_boot$info$bootstrapinfo$bayes_boot_vals,.025)))

  expect_length(mod1_bothboot$info$bmdl_alternatives,3)
  expect_named(get_all_bmdl(mod1_bothboot),c('score','delta','bootstrap','bmdl_bayes'))
  expect_type(get_all_bmdl(mod1_bothboot),'double')
  expect_length(get_all_bmdl(mod1_bothboot),4)

  # Computation times
  expect_named(get_computation_times(mod1),c('model','bmd','bmdl_score','bmdl_delta','total'))
  expect_named(get_computation_times(mod2),c('model','bmd','bmdl_score','bmdl_delta','total'))
  expect_length(get_computation_times(mod3),0)
  expect_named(get_computation_times(mod1_score),c('model','bmd','bmdl_score','total'))
  expect_named(get_computation_times(mod1_delta),c('model','bmd','bmdl_delta','total'))
  expect_named(get_computation_times(mod1_none),c('model','bmd','total'))

  expect_named(get_computation_times(mod1_boot),c('model','bmd','bmdl_score','bmdl_delta','bmdl_boot','total'))
  expect_named(get_computation_times(mod1_bayes),c('model','bmd','bmdl_score','bmdl_delta','bmdl_bayes','total'))
  expect_named(get_computation_times(mod1_bothboot),c('model','bmd','bmdl_score','bmdl_delta','bmdl_boot','bmdl_bayes','total'))

})

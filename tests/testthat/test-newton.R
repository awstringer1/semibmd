test_that("Newton's method works", {
  expect_lt(abs(bounded_newton(~.x,c(-1,1))),1e-08)
  expect_lt(abs(bounded_newton(function(x) x,c(-1,1))),1e-08)
  expect_lt(abs(bounded_newton("force",c(-1,1))),1e-08)

  expect_error(bounded_newton(~.x,1))
  expect_error(bounded_newton(~.x,c(1,2,3)))

})

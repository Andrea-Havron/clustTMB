
test_that("summary and ", {

library(clustTMB)
data("faithful")
m1 <- clustTMB(response = faithful, covariance.structure = "VVV")
l <- -m1$opt$objective
expect_equal(l, logLik(m1)[1])
expect_equal(-2*l + 2*length(m1$opt$par), AIC(m1))
expect_equal(summary(m1$sdr, "fixed"), summary(m1, "fixed"))
expect_equal(summary(m1$sdr, "fixed", p.value = TRUE), 
             summary(m1, "fixed", p.value = TRUE))
expect_equal(summary(m1$sdr, "all"), summary(m1, "all"))
expect_equal(summary(m1$sdr), summary(m1))
expect_equal(m1$opt$par, coef(m1))

})

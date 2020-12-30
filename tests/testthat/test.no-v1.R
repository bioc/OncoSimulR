inittime <- Sys.time()
cat(paste("\n Starting test.no-v1 at", date(), "\n"))

test_that("v1: error", {
    ## Do not do the capture output from oncoSimulPop,
    ## as that comes from mclapply and it is a mess.
    RNGkind("Mersenne-Twister")
    set.seed(1)
    p1 <- cbind(1L, 2L)
    expect_error(oncoSimulIndiv(p1,
                              sh = 0,
                              initSize = 1e5,
                              sampleEvery = 0.02,
                              detectionSize = 1e9,
                              model = "Exp",
                              finalTime = 2000,
                              extraTime = 3.17,
                              onlyCancer = FALSE,
                              seed = NULL),
                 "v.1 functionality has been removed",
                 fixed = TRUE)
})



cat(paste("\n Ending test.no-v1 at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

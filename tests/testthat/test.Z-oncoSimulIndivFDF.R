inittime <- Sys.time()
cat(paste("\n Starting test.Z-oncoSimulIndivFDF at", date(), "\n"))


test_that("testing output classes", {

  r <- data.frame(rfitness(2))

  colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"

  r[, "Fitness"] <- c("f_ - f_1 - f_2 - f_1_2",
                      "max(100*f_1, 10)",
                      "max(100*f_2, 10)",
                      "max((200*(f_1 + f_2) + 50*f_1_2), 1)")


  suppressWarnings(afe <- allFitnessEffects(genotFitness = r,
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel"))

  set.seed(1)

  osi <- oncoSimulIndiv(afe,
                        model = "McFL",
                        onlyCancer = FALSE,
                        finalTime = 20,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 5000,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE)

  expect_identical(class(r), "data.frame")

  expect_identical(class(afe), "fitnessEffects")


  expect_identical(class(osi), c("oncosimul", "oncosimul2"))

  if(as.character(version$major) < 4) {
  expect_identical(class(osi$Genotypes), "matrix")
  } else {
      expect_identical(class(osi$Genotypes), c("matrix", "array"))
  }
})

test_that("testing performance", {

  r <- data.frame(rfitness(2))

  colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"

  r[, "Fitness"] <- c("10*f_",
                      "10*f_1",
                      "50*f_2",
                      "200*(f_1 + f_2) + 50*f_1_2")
  ra <- r

  ra[, "Fitness"] <- c("10*n_/N",
                      "10*n_1/N",
                      "50*n_2/N",
                      "200*(n_1/N + n_2/N) + 50*n_1_2/N")


  suppressWarnings(afe <- allFitnessEffects(genotFitness = r,
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel"))

  suppressWarnings(afe_ra <- allFitnessEffects(genotFitness = ra,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs"))

  set.seed(1)

  null <- capture.output(osi <- oncoSimulIndiv(afe,
                        model = "Exp",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  set.seed(1)

  null <- capture.output(osi_ra <- oncoSimulIndiv(afe_ra,
                        model = "Exp",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  expect_output(print(osi),
                "Individual OncoSimul trajectory", fixed = TRUE)

  expect_identical(osi$NumClones, 3)

  expect_identical(osi$NumClones, osi_ra$NumClones)

  expect_true(osi$TotalPopSize >2e8)

  expect_identical(osi$TotalPopSize, osi_ra$TotalPopSize)

  expect_true(osi$geneNames[2] %in% LETTERS)

  expect_identical(osi$geneNames, osi_ra$geneNames)

  expect_false(osi$FinalTime < 1.4)

  expect_identical(osi$FinalTime, osi_ra$FinalTime)

  expect_identical(osi$NumIter, osi_ra$NumIter)

  ## It fails on Mac. Well, the key is the identity below
  skip_on_os(c("mac"))
  expect_equal(osi$NumIter, 458)
})

test_that("testing Bozic failure", {
  r <- data.frame(rfitness(2))

  colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"

  r[, "Fitness"] <- c("10*f_",
                      "10*f_1",
                      "50*f_2",
                      "200*(f_1 + f_2) + 50*f_1_2")


  suppressWarnings(afe <- allFitnessEffects(genotFitness = r,
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel"))

  set.seed(1)

  st <- capture.output(suppressWarnings(osi <- oncoSimulIndiv(afe,
                                             model = "Bozic",
                                             onlyCancer = FALSE,
                                             finalTime = 5000,
                                             verbosity = 0,
                                             mu = 1e-6,
                                             initSize = 500,
                                             keepPhylog = FALSE,
                                             seed = NULL,
                                             errorHitMaxTries = TRUE,
                                             errorHitWallTime = TRUE)))
  expect_true(st[22] == " Unrecoverable exception: Algo 2: retval not finite. Aborting. ")

  ## And test that plot gives an error
  expect_error(plot(osi), "An unrecoverable exception happened")

})

## This test is more a plotting test, but since we fix the seed
## and test plot, much as the last case above, we leave it here
## Skip on kjohnson3, arm64
if (Sys.getenv("R_PLATFORM") != "aarch64-apple-darwin20") {
  test_that("We produce an error when plotting timed out runs", {
    set.seed(1)
    data(examplesFitnessEffects)
    suppressMessages(tmp <-  oncoSimulIndiv(examplesFitnessEffects[["cbn2"]],
                                            model = "McFL",
                                            mu = 5e-6,
                                            detectionSize = 1e8,
                                            detectionDrivers = 3,
                                            sampleEvery = 0.03,
                                            max.num.tries = 10,
                                            keepEvery = 15,
                                            initSize = 20,
                                            finalTime = 1e5,
                                            max.wall.time = 0.01,
                                            onlyCancer = FALSE,
                                            keepPhylog = TRUE))
    expect_error(plot(tmp), "The simulation hit max wall time or max tries.")
  })
}


set.seed(NULL)

cat(paste("\n Ending test.Z-oncoSimulIndivFDF at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

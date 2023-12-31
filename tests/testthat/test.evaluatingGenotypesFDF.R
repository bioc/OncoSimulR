inittime <- Sys.time()
cat(paste("\n Starting test.evaluatingGenotypesFDF at", date(), "\n"))

test_that("testing single gene evaluation", {
  
  r1 <- data.frame(rfitness(2))
  
  colnames(r1)[which(colnames(r1) == "Birth")] <- "Fitness"
  r1[, "Fitness"] <- c("max(3, 2*f_)",
                      "max(1.5, 3*(f_ + f_1))",
                      "max(1.5, 3*(f_ + f_2))",
                      "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)")
  
  r2 <- rfitness(2)
  
  colnames(r2)[which(colnames(r2) == "Birth")] <- "Fitness"
  suppressWarnings(afe1 <- allFitnessEffects(genotFitness = r1, 
                            frequencyDependentFitness = TRUE, 
			                      frequencyType = "rel"))
                            #spPopSizes = c(5000, 2500, 2500, 7500))
  
  suppressWarnings(afe2 <- allFitnessEffects(genotFitness = r2, 
                            frequencyDependentFitness = FALSE))
                            #spPopSizes = c(5000, 2500, 2500, 7500))
  
  suppressWarnings(afe3 <- allFitnessEffects(genotFitness = r1, 
                            frequencyDependentFitness = TRUE, 
			    frequencyType = "rel"))
  suppressWarnings({
  evge1 <- evalGenotype(genotype = "Root", fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge2 <- evalGenotype(genotype = 0, fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge3 <- evalGenotype(genotype = "A", fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge4 <- evalGenotype(genotype = 1, fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge5 <- evalGenotype(genotype = c(1, 2), fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge6 <- evalGenotype(genotype = "A, B", fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  })
  
  expect_equal(evge1, evge2)
  
  expect_equal(evge3, evge4)
  
  expect_equal(evge5, evge6)
  suppressWarnings({
  expect_error(evalGenotype(genotype = c(0, 1), fitnessEffects = afe1,
                            spPopSizes = c(5000, 2500, 2500, 7500)), 
               "Genotype cannot contain any 0 if its length > 1")
  
  expect_error(evalGenotype(genotype = c(1, 3), fitnessEffects = afe1,
                            spPopSizes = c(5000, 2500, 2500, 7500)),
               "Genotype as vector of numbers contains genes not in fitnessEffects/mutatorEffects.")
  })
  expect_error(evalGenotype(genotype = 0, fitnessEffects = afe2), 
               "Genotype cannot be 0.")
  
  expect_error(evalGenotype(genotype = "Root", fitnessEffects = afe2), 
               "Genotype cannot be 0.")
  
  expect_error(evalGenotype(genotype = c(1, 2), fitnessEffects = afe3), 
               "You have a NULL spPopSizes")
  
  expect_error(evalGenotype(genotype = "1", fitnessEffects = afe1), 
               "You have a NULL spPopSizes")
  
  expect_error(evalGenotype(genotype = "1", fitnessEffects = afe2), 
               "Genotype contains NAs or genes not in fitnessEffects/mutatorEffects")
  
})

test_that("testing all genes evaluation", {
  
  r <- data.frame(rfitness(3))
  
  colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"
  
  r[, "Fitness"] <- c("max(2, log(1 + f_))",
                      "3",
                      "3",
                      "3", 
                      "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
                      "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
                      "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
                      "max(5, (1 + f_1_2 + f_2_3 + f_1_3)^2)")
  
  suppressWarnings(afe <- allFitnessEffects(genotFitness = r, 
                           frequencyDependentFitness = TRUE, 
			                     frequencyType = "rel"))
                           #spPopSizes = c(500, 
                            #              250, 
                             #             250, 
                              #            250, 
                               #           300,
                                #          300,
                                 #         300,
                                  #        450))
  
  genotypes <- c(0, OncoSimulR:::generateAllGenotypes(fitnessEffects = afe, 
                                                      order = FALSE, 
                                                      max = 256)$genotNums)

  suppressWarnings({
  evalGs_one_by_one <- sapply(genotypes, function(x) evalGenotype(x, afe,
                                                                  spPopSizes = c(500,
                                                                                 250, 
                                                                                 250, 
                                                                                 250, 
                                                                                 300,
                                                                                 300,
                                                                                 300,
                                                                                 450)))
  
  evalGs_all_together <- evalAllGenotypes(afe,
                                          spPopSizes = c(500,
                                                         250,
                                                         250,
                                                         250,
                                                         300,
                                                         300,
                                                         300,
                                                         450))$Fitness

  })                                     
  expect_identical(evalGs_one_by_one, evalGs_all_together)
  
})


test_that("eval single WT genotype with FDF" , {

    r <- data.frame(rfitness(2))
	
	colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"
    r[, "Fitness"] <- c("100 - n_1 - 2 * n_2 - 3 * n_1_2", 
                        "max(10*n_1, 4)", 
                        "max(10*n_2, 4)", 
                        "max((200*(n_1 + n_2) + 50*n_1_2), 1)")
						
    suppressWarnings(afe <- allFitnessEffects(genotFitness = r, 
                             frequencyDependentFitness = TRUE))


    suppressWarnings({
           ## Silencing the warnings, which are irrelevant, does not silence errors.
    ## uncomment this to see for yourself
    ## expect_true(2 == 3)
    ## expect_error(2 * 5)
    expect_equal(evalGenotype(0, afe, spPopSizes = rep(2, 4)),
                 100 - 2 - 4 - 6)

    expect_equal(evalGenotype(0, afe, spPopSizes = c(10, 10, 2, 8)),
                 100 - 10 - 4 - 24)

    expect_equal(evalGenotype(0, afe, spPopSizes = c(10, 10, 7, 8)),
                 100 - 10 - 14 - 24)


    expect_equal(evalGenotype("WT",  afe, spPopSizes = rep(2, 4)),
                 100 - 2 - 4 - 6)

    expect_equal(evalGenotype("WT",  afe, spPopSizes = c(10, 10, 2, 8)),
                 100 - 10 - 4 - 24)

    expect_equal(evalGenotype("WT", afe, spPopSizes = c(10, 10, 7, 8)),
                 100 - 10 - 14 - 24)

    
    expect_equal(evalGenotype("",  afe, spPopSizes = rep(2, 4)),
                 100 - 2 - 4 - 6)

    expect_equal(evalGenotype("",  afe, spPopSizes = c(10, 10, 2, 8)),
                 100 - 10 - 4 - 24)

    expect_equal(evalGenotype("", afe, spPopSizes = c(10, 10, 7, 8)),
                 100 - 10 - 14 - 24)


    expect_equal(evalGenotype("Root",  afe, spPopSizes = rep(2, 4)),
                 100 - 2 - 4 - 6)

    expect_equal(evalGenotype("Root",  afe, spPopSizes = c(10, 10, 2, 8)),
                 100 - 10 - 4 - 24)

    expect_equal(evalGenotype("Root", afe, spPopSizes = c(10, 10, 7, 8)),
                 100 - 10 - 14 - 24)

    
    
    ## These all do not interpret as WT but as mispelled gene
    expect_error(evalGenotype("aeiou", afe, spPopSizes = rep(10, 4)),
                 "Genotype contains NA or a gene not in fitnessEffects/mutatorEffects",
                 fixed = TRUE)
    expect_error(evalGenotype("root", afe, spPopSizes = rep(10, 4)),
                 "Genotype contains NA or a gene not in fitnessEffects/mutatorEffects",
                 fixed = TRUE)
    expect_error(evalGenotype("wt", afe, spPopSizes = rep(10, 4)),
                 "Genotype contains NA or a gene not in fitnessEffects/mutatorEffects",
                 fixed = TRUE)

 
    })
    
})


cat(paste("\n Ending test.evaluatingGenotypesFDF at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

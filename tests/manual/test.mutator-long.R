RNGkind("L'Ecuyer-CMRG") ## for the mclapplies

date()
test_that("McFL: Relative ordering of number of clones with mutator effects", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <- sample(999999, 1)
    set.seed(pseed)
    cat("\n x2: the seed is", pseed, "\n")
    pops <- 60
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "b" = 0.14,
                                           "c" = 0.16,
                                           "d" = 0.11))
    fm6 <- allMutatorEffects(noIntGenes = c("a" = 5,
                                            "b" = 10,
                                            "c" = 12,
                                            "d" = 14))
    pseed <- sample(999999, 1)
    set.seed(pseed)
    cat("\n x2a: the seed is", pseed, "\n")
    nc1 <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =250,
                        mutationPropGrowth = FALSE,
                        initSize  = 1e6, model = "McFL",
                        mc.cores = 2, seed = NULL,
                        onlyCancer = FALSE)
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "b" = 0.14,
                                           "c" = 0.16,
                                           "d" = 0.11))
    fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                            "b" = 1,
                                            "c" = 1,
                                            "d" = 1))
    pseed <- sample(999999, 1)
    set.seed(pseed)
    cat("\n x2b: the seed is", pseed, "\n")
    nc2 <- oncoSimulPop(pops, fe, muEF = fm8, finalTime =250,
                        mutationPropGrowth = FALSE,
                        initSize  = 1e6, model = "McFL",
                        mc.cores = 2, seed = NULL, 
                        onlyCancer = FALSE)
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "b" = 0.14,
                                           "c" = 0.16,
                                           "d" = 0.11))
    fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
                                            "b" = 1e-6,
                                            "c" = 1e-6,
                                            "d" = 1e-6))
    gc(); pseed <- sample(999999, 1)
    set.seed(pseed)
    cat("\n x2c: the seed is", pseed, "\n")
    nc3 <- oncoSimulPop(pops, fe, muEF = fm7, finalTime =250,
                        mutationPropGrowth = FALSE,
                        initSize  = 1e6, model = "McFL",
                        mc.cores = 2, seed = NULL, 
                        onlyCancer = FALSE)
    expect_true(median(summary(nc1)$NumClones) > median(summary(nc2)$NumClones))
    expect_true(median(summary(nc2)$NumClones) > median(summary(nc3)$NumClones))
    gc()
})
date()


date() 
test_that("Mutator, several modules differences, sample at end only, and larger effect", {
    ## we also change the effect (larger) and finalTime (shorter)
    ## but this can from time to time take veeeery long.
    pseed <- sample(99999999, 1)
    set.seed(pseed)
    cat("\n l-mmd2: the seed is", pseed, "\n")
    reps <- 10
    no <- 5e3
    ft <- 30 ## you need it large enough to get enough hits
    mu <- 1e-5
    ln <- 50 
    m1 <- 6 ## if this is too large, easy to get it to blow.
    ni <- rep(0, 3 * ln)
    gna <- paste0("a", 1:ln)
    gnb <- paste0("b", 1:ln)
    gnc <- paste0("c", 1:ln)
    names(ni) <- c(gna, gnb, gnc)
    gn1 <- paste(c(gna, gnb, gnc), collapse = ", ")
    gna <- paste(gna, collapse = ", ")
    gnb <- paste(gnb, collapse = ", ")
    gnc <- paste(gnc, collapse = ", ")
    mut1 <- allMutatorEffects(epistasis = c("A" = m1),
                              geneToModule = c("A" = gn1))
    mut2 <- allMutatorEffects(epistasis = c("A" = m1,
                                            "B" = m1,
                                            "C" = m1),
                              geneToModule = c("A" = gna,
                                               "B" = gnb,
                                               "C" = gnc))
    f1 <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(99999999, 1)
    set.seed(pseed)
    cat("\n l-mmd2a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       f1,
                       mu = mu,
                       muEF = mut1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mc.cores = 2,
                       seed = NULL, sampleEvery = ft
                       )
    gc()
    gc(); pseed <- sample(99999999, 1)
    set.seed(pseed)
    cat("\n l-mmd2b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       f1,
                       mu = mu,
                       muEF = mut2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mc.cores = 2,
                       seed = NULL, sampleEvery = ft
                       )
    gc()
    summary(b2)[, c(1:3, 8:9)]
    summary(b1)[, c(1:3, 8:9)]
    mean(mutsPerClone(b2))
    mean(mutsPerClone(b1))
    ## This is, of course, affected by sampling only at end: we do not see
    ## the many intermediate events.
    ## expect_true( median(summary(b2)$NumClones) >
    ##              median(summary(b1)$NumClones))
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
    gc()
})
date()

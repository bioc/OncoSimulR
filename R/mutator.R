## Copyright 2013-2021 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Any mutator gene must have been in fitness, even if it has fitness
## effect of 0. Yes, because otherwise the mapping of gene names to
## numerical ids becomes a PITA.


allMutatorEffects <- function(epistasis = NULL,
                              noIntGenes = NULL,
                              geneToModule = NULL,
                              keepInput = TRUE) {
    ## This is on purpose to prevent using a rT or orderEffects. Those are
    ## not tested to work with mutator.

    ## Neither do we accept a fitness landscape object either for now.
    allFitnessORMutatorEffects(
        rT = NULL,
        epistasis = epistasis,
        orderEffects = NULL,
        noIntGenes = noIntGenes,
        geneToModule = geneToModule,
        drvNames = NULL,
        keepInput = keepInput,
        calledBy = "allMutatorEffects")
}
                              

evalAllGenotypesMut <- function(mutatorEffects,
                                max = 256,
                                addwt = FALSE,
                                spPopSizes = NULL,
                                currentTime = 0
                                ) {
    evalAllGenotypesORMut(
        fmEffects = mutatorEffects,
        order = FALSE,
        max = max,
        model = "",
        spPopSizes = spPopSizes,
        calledBy_= "evalGenotypeMut",
        currentTime = currentTime
    )
}

evalGenotypeMut <- function(genotype, 
                            mutatorEffects,
                            spPopSizes = NULL,
                            verbose = FALSE,
                            echo = FALSE,
                            currentTime = 0
                            ) {
    
    if(inherits(mutatorEffects, "fitnessEffects"))
        stop("You are trying to get the mutator effects of a fitness specification. ",
             "You did not pass an object of class mutatorEffects.")
    evalGenotypeORMut(genotype = genotype,
                      fmEffects = mutatorEffects,
                      spPopSizes = spPopSizes,
                      verbose = verbose,
                      echo = echo,
                      model  = "" ,
                      calledBy_= "evalGenotypeMut",
                      currentTime = currentTime
                      )

}


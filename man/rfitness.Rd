\name{rfitness}
\alias{rfitness}
\encoding{UTF-8}

\title{Generate random fitness.}

\description{ Generate random fitness landscapes under a House of Cards,
  Rough Mount Fuji (RMF), additive (multiplicative) model, Kauffman's NK
  model, Ising model, Eggbox model and Full model}


\usage{

rfitness(g, c = 0.5, sd = 1, mu = 1, reference = "random", scale = NULL,
         wt_is_1 = c("subtract", "divide", "force", "no"),
         log = FALSE, min_accessible_genotypes = NULL,
         accessible_th = 0, truncate_at_0 = TRUE,
         K = 1, r = TRUE, i = 0, I = -1, circular = FALSE, e = 0, E = -1,
         H = -1, s = 0.1, S = -1, d = 0, o = 0, O = -1, p = 0, P = -1, 
         model = c("RMF", "Additive", "NK", "Ising", "Eggbox", "Full"),
         seed_magellan = -1)
}




\arguments{

  \item{g}{Number of genes.}

  \item{c}{The decrease in fitness of a genotype per each unit increase
    in Hamming distance from the reference genotype for the RMF model
    (see \code{reference}).}

  \item{sd}{The standard deviation of the random component (a normal
  distribution of mean \code{mu} and standard deviation \code{sd}) for
  the RMF and additive models .}

\item{mu}{The mean of the random component (a normal distribution of
mean \code{mu} and standard deviation \code{sd}) for the RMF and
additive models.}


\item{reference}{The reference genotype: in the RMF model, for the
  deterministic, additive part, this is the genotype with maximal
  fitness, and all other genotypes decrease their fitness by \code{c}
  for every unit of Hamming distance from this reference. If "random" a
  genotype will be randomly chosen as the reference. If "max" the
  genotype with all positions mutated will be chosen as the
  reference. If you pass a vector (e.g., \code{reference = c(1, 0, 1,
  0)}) that will be the reference genotype.  If "random2" a genotype
  will be randomly chosen as the reference. In contrast to "random",
  however, not all genotypes have the same probability of being chosen;
  here, what is equal is the probability that the reference genotype has
  1, 2, ..., g, mutations (and, once a number mutations is chosen, all
  genotypes with that number of mutations have equal probability of
  being the reference). }

\item{scale}{Either NULL (nothing is done) or a two- or three-element
  vector.

  If a two-element vector, fitness is re-scaled between
  \code{scale[1]} (the minimum) and \code{scale[2]} (the maximum) and,
  later, if you have selected it, \code{wt_is_1} will be enforced.

  If you pass a three element vector, fitness is re-scaled so that the
  new maximum fitness is \code{scale[1]}, the new minimum is
  \code{scale[2]} and the new wildtype is \code{scale[3]}. If you pass a
  three element vector, none of the \code{wt_is_1} options apply in this
  case, to ensure you obtain the range you want. If you want the
  wildtype to be one, pass it as the third element of the vector.

  As a consequence of using a three element vector, the amount of
  stretching/compressing (i.e., scaling) of fitness values larger than
  that of the wildtype will likely be different from the scaling of
  fitness values smaller than that of the wildtype.  In other words,
  this argument allows you to change the spread of the positive and
  negative fitness values (and you can make this difference extreme and
  make most fitness values less than wildtype be 0 by using a huge
  negative number --huge in absolute value-- for \code{scale[2]} if you
  then truncate at 0 --see \code{truncate_at_9}).

  Using a three element vector is probably the most natural way of
  changing the scale and range of fitness.

  See also \code{log} if you want the log-transformed values to respect
  the scale.
}

\item{wt_is_1}{If "divide" the fitness of all genotypes is
  divided by the fitness of the wildtype (after possibly adding a value
  to ensure no negative fitness) so that the wildtype (the genotype with
  no mutations) has fitness 1. This is a case of scaling, and it is
  applied after \code{scale}, so if you specify both
  "wt_is_1 = 'divide'" and use an argument for \code{scale} it is most
  likely that the final fitness will not respect the limits in
  \code{scale}.

  If "subtract" (the default) we shift all the fitness values (subtracting fitness of
  the wildtype and adding 1) so that the wildtype ends up with a fitness
  of 1. This is also applied after \code{scale}, so if you specify both
  "wt_is_1 = 'subtract'" and use an argument for \code{scale} it is most
  likely that the final fitness will not respect the limits in
  \code{scale} (though the distorsion might be simpler to see as just a
  shift up or down).
  
  If "force" we simply set the fitness of the wildtype to 1, without any
  divisions. This means that the \code{scale} argument would work (but
  it is up to you to make sure that the range of the scale argument
  includes 1 or you might not get what you want). Note that using this
  option can easily lead to landscapes with no accessible genotypes
  (even if you also use \code{scale}).

  If "no", the fitness of the wildtype is not modified.

  This option has no effect if you pass a three-element vector for
  \code{scale}. Using a three-element vector for \code{scale} is
  probably the most natural way of changing the scale and range of
  fitness while setting the wildtype to a value of your choice.
  
}


\item{log}{If TRUE, log-transform fitness. Actually, there are two
  cases: if \code{wt_is_1 = "no"} we simply log the fitness values;
  otherwise, we log the fitness values and add a 1, thus shifting all
  fitness values, because by decree the fitness (birth rate) of the
  wildtype must be 1.

  If you pass a three-element vector for scale, you will want to pass
  \code{exp(desired_max)}, \code{exp(desired_min)}, and
  \code{exp(desired_wildtype)} to the \code{scale} argument. (We first
  scale values in the original scale and then log them). In this case,
  we ignore whatever you passed as \code{wt_is_1}, setting \code{wt_is_1
  = "no"} to avoid modifying your requested value for the wildtype.}

\item{min_accessible_genotypes}{If not NULL, the minimum number of
  accessible genotypes in the fitness landscape. A genotype is
  considered accessible if you can reach if from the wildtype by going
  through at least one path where all changes in fitness are larger or
  equal to \code{accessible_th}. The changes in fitness are considered
  at each mutational step, i.e., at each addition of one mutation we
  compute the difference between the genotype with \code{k + 1}
  mutations minus the ancestor genotype with \code{k} mutations. Thus, a
  genotype is considered accessible if there is at least one path where
  fitness increases at each mutational step by at least
  \code{accessible_th}.

  If the condition is not satisfied, we continue generating random
  fitness landscapes with the specified parameters until the condition
  is satisfied.

  (Why check against NULL and not against zero? Because this allows you
  to count accessible genotypes even if you do not want to ensure a
  minimum number of accessible genotypes.)
}

\item{accessible_th}{The threshold for the minimal change in fitness at
  each mutation step (i.e., between successive genotypes) that allows a
  genotype to be regarded as accessible. This only applies if
  \code{min_accessible_genotypes} is larger than 0.  So if you want to
  allow small decreases in fitness in successive steps, use a small
  negative value for \code{accessible_th}.  }

\item{truncate_at_0}{If TRUE (the default) any fitness <= 0 is
  substituted by a small positive constant (a random uniform number
  between 1e-10 and 1e-9). Why? Because MAGELLAN and some plotting
  routines can have trouble (specially if you log) with values <=0. Or
  we might have trouble if we want to log the fitness. This is done
  after possibly taking logs. Noise is added to prevent creating several
  identical minimal fitness values.  Note that \code{\link{allFitnessEffects}} will remove from the table
  of genotypes any genotype with a fitness <= 1e-9, thus 
    making it a non-viable genotype during simulations. }

\item{K}{K for NK model; K is the number of loci with which each locus
  interacts, and the larger the K the larger the ruggedness of the
  landscape.}

\item{r}{For the NK model, whether interacting loci are chosen at random
  (\code{r = TRUE}) or are neighbors (\code{r = FALSE}).}
\item{i}{For de Ising model, i is the mean cost for incompatibility with which
  the genotype's fitness is penalized when in two adjacent genes, only one of 
  them is mutated.}

\item{I}{For the Ising model, I is the standard deviation for the cost 
  incompatibility (i).}
  
\item{circular}{For the Ising model, whether there is a circular arrangement, 
  where the last and the first genes are adjacent to each other.}

\item{e}{For the Eggbox model, mean effect in fitness for the neighbor
  locus +/- e.}
  
\item{E}{For the Eggbox model, noise added to the mean effect in fitness (e).}

\item{H}{For Full models, standard deviation for the House of Cards model.}

\item{s}{For Full models, mean of the fitness for the Multiplicative model.}

\item{S}{For Full models, standard deviation for the Multiplicative model.}

\item{d}{For Full models, a disminishing (negative) or increasing 
  (positive) return as the peak is approached for multiplicative model.}
  
\item{o}{For Full models, mean value for the optimum model.}

\item{O}{For Full models, standard deviation for the optimum model.}

\item{p}{For Full models, the mean production value for each non 0
  allele in the Optimum model component.}

\item{P}{For Full models, the associated stdev (of non 0 alleles) in the
Optimum model component.}



\item{model}{One of "RMF" (default) for Rough Mount Fuji, "Additive" for
 Additive model, "NK", for Kauffman's NK model, "Ising" for Ising model,
 "Eggbox" for Eggbox model or "Full" for Full models.}

\item{seed_magellan}{The seed for the random number generator in models
generated from MAGELLAN. If -1, the clock is used by MAGELLAN to
generate a seed, but you probably want to pass a seed: see Details.}
} 


\details{

  When using \code{model = "RMF"}, the model used here follows
  the Rough Mount Fuji model in Szendro et al., 2013 or Franke et al.,
  2011. Fitness is given as

  \deqn{f(i) = -c d(i, reference) + x_i}

  where \eqn{d(i, j)} is the Hamming distance between genotypes \eqn{i}
  and \eqn{j} (the number of positions that differ) and \eqn{x_i} is a
  random variable (in this case, a normal deviate of mean \code{mu}
  and standard deviation \code{sd}).

  When using \code{model = "RMF"}, setting \eqn{c = 0} we obtain a House
    of Cards model. Setting \eqn{sd = 0} fitness is given by the
    distance from the reference and if the reference is the genotype
    with all positions mutated, then we have a fully additive model
    (fitness increases linearly with the number of positions mutated),
    where all mutations have the same effect.

  More flexible additive models can be used using \code{model =
  "Additive"}. This model is like the Rough Mount Fuji model in Szendro
  et al., 2013 or Franke et al., 2011, but in this case, each locus can
  have different contributions to the fitness evaluation. This model is
  also referred to as the "multiplicative" model in the literature as it
  is additive in the log-scale (e.g., see Brouillet et al., 2015 or
  Ferretti et al., 2016). The contribution of each mutated allele to the
  log-fitness is a random deviate from a Normal distribution with
  specified mean \code{mu} and standard deviation \code{sd}, and the
  log-fitness of a genotype is the sum of the contributions of each
  mutated allele. There is no "reference" genotype in the Additive
  model.  There is no epistasis in the additve model because the effect
  of a mutation in a locus does not depend on the genetic background, or
  whether the rest of the loci are mutated or not.
  

  When using \code{model = "NK"} fitness is drawn from a uniform (0, 1)
  distribution.
  
  
  When using \code{model = "Ising"} for each pair of interacting loci, 
  there is an associated cost if both alleles are not identical 
  (and therefore 'compatible').
  
  
  When using \code{model = "Eggbox"} each locus is either high or low fitness,
  with a systematic change between each neighbor.
  
  
  When using \code{model = "Full"}, the fitness is computed with different
  parts of the previous models depending on the choosen parameters described 
  above. 
  
  
  For \code{model = "NK" | "Ising" | "Eggbox" | "Full"} the fitness
  landscape is generated by directly calling the \code{fl_generate}
  function of MAGELLAN
  (\url{http://wwwabi.snv.jussieu.fr/public/Magellan/}). See details in
  Ferretti et al. 2016, or Brouillet et al., 2015.
  

  For OncoSimulR, we often want the wildtype to have a mean of
  1. Reasonable settings when using RMF are \code{mu = 1} and \code{wt_is_1 =
  'subtract'} so that we simulate from a distribution centered in 1, and
  we make sure afterwards (via a simple shift) that the wildtype is
  actuall 1. The \code{sd} controls the standard deviation, with the
  usual working and meaning as in a normal distribution, unless \code{c}
  is different from zero. In this case, with \code{c} large, the range
  of the data can be large, specially if \code{g} (the number of genes)
  is large.

  Note that \code{\link{allFitnessEffects}} will remove from the table
  of genotypes any genotype with a fitness <= 1e-9, thus 
    making it a non-viable genotype during simulations.   


 \code{seed_magellan}: if you run code in parallel or you use sequential
 code where you generate random fitness landscapes generated by MAGELLAN
 (\code{model = "NK" | "Ising" | "Eggbox" | "Full"}) in a short time,
 MAGELLAN would likely end up using the same seed as the different calls
 would be done within the same time (within second resolution). Thus,
 especially if you are generating the same kind of fitness landscape, 
 you probably want to pass different seeds. The seed is read as a C
 long, so you should be able to use integers going from at least
 -2,147,483,647 to +2,147,483,647, in 32-bits, but probably a much
 larger range (-9,223,372,036,854,775,808 to 9,223,372,036,854,775,807)
 in 64-bits. Note, though, that some values can crash MAGELLAN (for
 example -9223372036854775806 or -2147483647). You do not want to pass
 numbers in scientific notation; for example, you could instead do
 \code{seed_magellan = format(2^40, scientific = FALSE)}. 
 
    
} 

\value{
  
  An matrix with \code{g + 1} columns. Each column corresponds to a
  gene, except the last one that corresponds to fitness. 1/0 in a gene
  column denotes gene mutated/not-mutated. (For ease of use in other
  functions, this matrix has class  "genotype_fitness_matrix".) 

  If you have specified \code{min_accessible_genotypes > 0}, the return
  object has added attributes \code{accessible_genotypes} and
  \code{accessible_th} that show the number of accessible
  genotypes under the specified  threshold.
}


\note{MAGELLAN uses its own random number generating functions; using
  \code{set.seed} does not allow to obtain the same fitness landscape
  repeatedly.}

\references{

  Szendro I.~G. et al. (2013). Quantitative analyses of empirical
fitness landscapes. \emph{Journal of Statistical Mehcanics: Theory and
  Experiment\/}, \bold{01}, P01005.

Franke, J. et al. (2011). Evolutionary accessibility of mutational
pathways. \emph{PLoS Computational Biology\/}, \bold{7}(8), 1--9.

Brouillet, S. et al. (2015). MAGELLAN: a tool to explore small fitness
landscapes. \emph{bioRxiv},
\bold{31583}. \url{http://doi.org/10.1101/031583}

Ferretti, L., Schmiegelt, B., Weinreich, D., Yamauchi, A., Kobayashi,
Y., Tajima, F., & Achaz, G. (2016). Measuring epistasis in fitness
landscapes: The correlation of fitness effects of mutations. \emph{Journal of
Theoretical Biology\/}, \bold{396}, 132--143. \url{https://doi.org/10.1016/j.jtbi.2016.01.037}

MAGELLAN web site: \url{http://wwwabi.snv.jussieu.fr/public/Magellan/}

}

\author{ Ramon Diaz-Uriarte for the RMF and general wrapping
  code. S. Brouillet, G. Achaz, S. Matuszewski, H. Annoni, and
  L. Ferreti for the MAGELLAN code. Further contributions to the
  additive model and to wrapping MAGELLAN code and documentation from
  Guillermo Gorines Cordero, Ivan Lorca Alonso, Francisco Muñoz Lopez,
  David Roncero Moroño, Alvaro Quevedo, Pablo Perez, Cristina Devesa,
  Alejandro Herrador.}

\seealso{
  
  \code{\link{oncoSimulIndiv}},
  \code{\link{plot.genotype_fitness_matrix}},
  \code{\link{evalAllGenotypes}}
  \code{\link{allFitnessEffects}}
  \code{\link{plotFitnessLandscape}}
  \code{\link{Magellan_stats}}  

}
\examples{

## Random fitness for four genes-genotypes,
## plotting and simulating an oncogenetic trajectory


## NK model
rnk <- rfitness(5, K = 3, model = "NK")
plot(rnk)
oncoSimulIndiv(allFitnessEffects(genotFitness = rnk))

## Additive model
radd <- rfitness(4, model = "Additive", mu = 0.2, sd = 0.5)
plot(radd)


\dontrun{
## Eggbox model
regg = rfitness(g=4,model="Eggbox", e = 2, E=2.4)
plot(regg)


## Ising model
ris = rfitness(g=4,model="Ising", i = 0.002, I=2)
plot(ris)


## Full model
rfull = rfitness(g=4, model="Full", i = 0.002, I=2, 
                 K = 2, r = TRUE,
                 p = 0.2, P = 0.3, o = 0.3, O = 1)
    plot(rfull)
    }
}
\keyword{ datagen }


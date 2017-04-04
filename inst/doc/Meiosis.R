## ---- include = FALSE----------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts = list(width.cutoff = 80), tidy = TRUE)

## ------------------------------------------------------------------------
library('Meiosis')
set.seed(123L) ## Seed R's rng
Meiosis::seed_rng(seed = 123L) ## Seed rng used by Meiosis

## ------------------------------------------------------------------------
n_chr <- 3L  ## number of chromosomes
L <- runif(n = n_chr, min = 100, max = 300)  ## sample length of chromosomes in cM
xoparam <- create_xoparam(L)  ## no interference, no obligate chiasma
str(xoparam)

## ------------------------------------------------------------------------
n_loci <- round(runif(n = n_chr, min = 5L, max = 10L))  ## sample number of loci
## sample positions of loci on the chromosome
positions <- lapply(seq_len(n_chr), function(i) sort(runif(n_loci[i], min = 0, max = L[i])))

## ------------------------------------------------------------------------
ind <- replicate(2L, lapply(n_loci, function(n) sample(c(0L, 1L), n, replace = TRUE)),
                 simplify = FALSE) ## simulate some genotypic data
str(ind)

p_geno <- Meiosis::cross_geno(father = ind, mother = ind, positions = positions,
                              xoparam = xoparam)
str(p_geno)

## ------------------------------------------------------------------------
f_alleles <- c(21L, 65L) ## 21 and 65 are arbitrary integers
f <- Meiosis::create_xo_founder(alleles = f_alleles, L = L)

p_xo <- Meiosis::cross_xo(father = f, mother = f, xoparam = xoparam)
str(p_xo)

## ------------------------------------------------------------------------
conv <- new(Meiosis::Converter, positions) ## create a new converter object
conv$insert_founder(f_alleles, ind) ## insert the (one and only) founder
str(conv$convert(p_xo)) ## convert the progeny

## ------------------------------------------------------------------------
n_self <- 10L  ## number of generations of selfing
n <- 30L ## number of progeny

## Genotypic representation
ind2 <- replicate(2L, lapply(n_loci, function(n) sample(c(0L, 1L), n, replace = TRUE)),
                  simplify = FALSE) # Second individual as parent.

pop <- replicate(n, Meiosis::cross_geno(ind, ind2, positions, xoparam), simplify = FALSE)
for (i in seq_len(n_self)) {
  for (j in seq_len(n)) {
    pop[[i]] <- Meiosis::cross_geno(pop[[i]], pop[[i]], positions, xoparam)
  }
}

## Segmental representation
f2 <- create_xo_founder(alleles = c(55L, 77L), L = L)
pop_xo <- replicate(n, Meiosis::cross_xo(f, f2, xoparam), simplify = FALSE)
for (i in seq_len(n_self)) {
  for (j in seq_len(n)) {
    pop_xo[[i]] <- Meiosis::cross_xo(pop_xo[[i]], pop_xo[[i]], xoparam)
  }
}


# conv$convert(pop[[1]]) ## error, because genotypic data of second founder not present
conv$insert_founder(c(55L, 77L), ind2)  ## insert second founder first
pop_geno <- lapply(pop_xo, conv$convert) ## convert whole population

## ------------------------------------------------------------------------
make_synthetic <- function(founder, n_ind, n_gen) {

  ## Cross parents
  n_founder <- length(founder)
  tmp <- combn(x = seq_len(n_founder), m = 2L)
  combinations <- split(tmp, col(tmp))
  pop_xo <- replicate(n = n_ind, simplify = FALSE, {
    pair <- unlist(sample(combinations, size = 1L))
    cross_xo(founder[[pair[1L]]], founder[[pair[2L]]], xoparam)
  })

  ## Random mating
  for (i in seq_len(n_gen)) {
    pop_xo_new <- pop_xo ## copy
    for (j in seq_len(n_ind)) {
      pair <- sample(n_ind, size = 2L, replace = TRUE)  ## selfing possible
      pop_xo_new[[j]] <- cross_xo(pop_xo[[pair[1L]]], pop_xo[[pair[2L]]], xoparam)
    }
    pop_xo <- pop_xo_new ## swap
  }
  pop_xo
}

n_founder <- 5L  ## number of founders
n_ind <- 100L  ## size of synthetic
n_gen <- 10L  # generations of random mating
alleles <- lapply(seq_len(n_founder), function(i) c(2L * i - 1L, 2L * i))
founder <- lapply(alleles, create_xo_founder, L = L)

## Create synthetic
system.time(syn <- make_synthetic(founder, n_ind, n_gen))

## ------------------------------------------------------------------------
## Calculate realized coefficients of co-ancestry.
Meiosis::realized_coancestry(f)
Meiosis::realized_coancestry(p_xo)  ## selfing progeny, expected coefficient is 0.75.

## ------------------------------------------------------------------------
library('pedigree')

## Create a simple pedigree
id <- 1:6
dam <- c(0, 0, 1, 1, 4, 4)
sire <- c(0, 0, 2, 2, 3, 5)
ped <- data.frame(id, dam, sire)

## Compute the additive genetic relationship matrix and coefficients of co-ancestry
cwd <- getwd()
tpdir <- tempdir()
setwd(tpdir)
invisible(makeA(ped, which = rep(TRUE, length(id))))
coanc <- read.table("A.txt")
setwd(cwd)

A <- matrix(NA_real_, nrow = length(id), ncol = length(id))
A[as.matrix(coanc[1:2])] <- A[as.matrix(coanc[2:1])] <- coanc[[3]]
eCoc <- A / 2  ## expected coefficients of co-ancestry

## Helper function for simulating pedigree and computing realized coefficients of co-ancestry.
sim_ped <- function() {
  f1 <- create_xo_founder(c(1L, 2L), L)
  f2 <- create_xo_founder(c(3L, 4L), L)
  i1 <- cross_xo(f1, f2, xoparam)
  i2 <- cross_xo(f1, f2, xoparam)
  i3 <- cross_xo(i1, i2, xoparam)
  i4 <- cross_xo(i2, i3, xoparam)

  tmp <- list(f1, f2, i1, i2, i3, i4)
  C <- matrix(data = NA_real_, nrow = length(id), ncol = length(id))
  for (i in seq_along(id))
    for (j in i:length(id))
      C[i, j] <- C[j, i] <- realized_coancestry(tmp[[i]], tmp[[j]])
  C
}

## Verify that, on average, the realized coefficients are equal to the expected coefficients.
n <- 1000L  ## number of replicates
rCoc_avg <- Reduce(f = `+`, x = replicate(n, sim_ped(), simplify = FALSE)) / n ##  take average
plot(as.vector(rCoc_avg), as.vector(eCoc)); abline(0, 1)

## ------------------------------------------------------------------------
## Simulate a doubled haploid individual.
str(Meiosis::dh_geno(ind, positions, xoparam))
str(conv$convert(Meiosis::dh_xo(f, xoparam)))

## ------------------------------------------------------------------------
## Calculate realized heterozygosity.
Meiosis::realized_heter(p_xo) 


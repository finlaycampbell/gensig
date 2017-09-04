##====== Sending off to the cluster =====##

## Set up parameters for cluster run
cluster.param <- function() {
  diseases <- names(create.param())
  param <- data.frame(disease = rep(diseases, each = 5))
  param$min.n <- 30
  param$n.hosts <- 100
  param$runs <- 20
  param$disease <- as.character(param$disease)
  return(param)
}

## Testing if empirical and analytical solution match
test.anlyt <- function(dis = 'tb') {

  param <- create.param()[[dis]]
  w.mean <- param$w.mean
  w.sd <- param$w.sd
  seql <- 10 #param$seql
  mu <- param$mut
  mu.transi <- mu*(2/3)
  mu.transv <- mu*(1/3)
  anylt <- w.mean
  w <- sapply(0:300, EpiEstim::DiscrSI, w.mean, w.sd)
  if(w[1] == 0) w[1] <- 1e-15

  sim <- function() {
    repeat({out <- outbreaker::simOutbreak(param$R0,
                                           n.hosts = 100,
                                           duration = 500,
                                        #dgamma(0:100, (w.mean^2)/(w.sd^2) , w.mean/(w.sd^2)),
                                           w,
                                           mu.transi = mu.transi,
                                           mu.transv = mu.transv,
                                           seq.length = seql,
                                           rate.import.case = 0)
                                           if(out$n > 10) break
    })

    ret <- out$onset - out$onset[out$ances]
    ret <- ret[!is.na(ret)]

    return(ret)
                                        #return(mean(out$nmut, na.rm = TRUE))

  }

  print(paste("Analytical solution is", anylt))
  out <- unlist(replicate(20, sim()))
  print(paste("Empirical mean is", mean(out)))
  print(paste(length(out), "cases"))

  diff <- anylt - mean(out)

  sq <- 1:20
  df <- ldply(table(out)/sum(table(out)))
  plot(df$.id[sq], df$V1[sq], type = 'l', lwd = 2, col = 'red')
  lines(sq, w[sq + 1])

  df2 <- data.frame(dis = dis, diff = diff, w.mean = w.mean, w.sd = w.sd)

  return(df2)

}

## Run phyb analysis on the cluster
run.phyb.cluster <- function(disease, min.n, n.hosts, runs) {
  
  config <- create.config(min.n = min.n,
                          n.hosts = n.hosts)
  
  param <- create.param()
  param <- list(tmp = param[[disease]])
  names(param) <- disease
  param[[disease]]$w[param[[1]]$w == 0] <- 1e-50

  store <- list()

  for(i in seq_len(runs)) {
    
    mean.gen <- param[[disease]]$w.mean
    shape.gen <- param[[disease]]$w.mean^2/param[[disease]]$w.sd^2
    mean.sample <- param[[disease]]$w.mean
    shape.sample <- param[[disease]]$w.mean^2/param[[disease]]$w.sd^2

    seql <- param[[disease]]$seql
    mut = param[[disease]]$mut

    if(seql > 1e6) {
      seql <- round(seql/100, 0)
      mut <- mut*100
    }
    
    sim.n <- 0

    while(sim.n < min.n) {
      
      sim <- sim.phybreak(obsize = NA,
                          popsize = n.hosts,
                          R0 = param[[disease]]$R0,
                          mean.gen = mean.gen,
                          shape.gen = shape.gen,
                          mean.sample = mean.sample,
                          shape.sample = shape.sample,
                          sequence.length = seql,
                          mu = mut,
                          wh.model = 2,
                          wh.slope = 1)

      ## Check for no-transmission
      if(inherits(sim, 'character')) next
      
      sim.n <- length(sim$sample.hosts)

    }

    dna.res <- phybreak(sim,
                        gen.mean = mean.gen,
                        gen.shape = shape.gen,
                        sample.mean = mean.sample,
                        sample.shape = shape.sample,
                        est.gen.mean = FALSE,
                        est.sample.mean = FALSE,
                        est.wh.slope = FALSE) %>%
      burnin.phybreak(ncycles = 1e3) %>% 
      sample.phybreak(nsample = 500, thin = 20) 

    ## for the no.dna run make all sequences the same
    nodna.sim <- sim
    for(j in seq_along(nodna.sim$sequences)) {
      nodna.sim$sequences[[j]] <- nodna.sim$sequences[[1]]
    }
    
    nodna.res <- phybreak(nodna.sim,
                          gen.mean = mean.gen,
                          gen.shape = shape.gen,
                          sample.mean = mean.sample,
                          sample.shape = shape.sample,
                          est.gen.mean = FALSE,
                          est.sample.mean = FALSE,
                          est.wh.slope = FALSE) %>%
      burnin.phybreak(ncycles = 1e3) %>% 
      sample.phybreak(nsample = 500, thin = 20) 
    
    store$param <- param
    store$config <- config
    store$phyb.sim[[i]] <- sim
    store$phyb.gensig[[i]] <- get.phyb.gensig(sim)
    store$phyb.dna.res[[i]] <- dna.res
    store$phyb.nodna.res[[i]] <- nodna.res
    store$phyb.dna.acc[[i]] <- get.phyb.acc(dna.res, sim)
    store$phyb.nodna.acc[[i]] <- get.phyb.acc(nodna.res, sim)
    store$phyb.uniq[[i]] <- get.phyb.uniq(sim)

  }

  return(store)
}

## Run analysis on the cluster
run.cluster <- function(disease, min.n, n.hosts, runs, imp, dur) {

  config <- create.config(min.n = min.n,
                          n.hosts = n.hosts,
                          imp = 0,
                          dur = 500)

  disease <- disease

  param <- create.param()
  param <- list(tmp = param[[disease]])
  names(param) <- disease
  param[[disease]]$w[param[[1]]$w == 0] <- 1e-50

  store <- list()

  for(i in seq_len(runs)) {

    gensig.result <- gensig(param = param, config = config)

    store$param <- gensig.result$param
    store$config <- gensig.result$config
    store$gensig[[i]] <- gensig.result$store[[disease]]$gensig

    sim <- gensig.result$store[[disease]]$sim
    store$sim[[i]] <- sim

    dna.outb.data <- list(w_dens = param[[disease]]$w,
                          dna = sim$dna,
                          dates = sim$onset)

    nodna.outb.data <- dna.outb.data
    nodna.outb.data$dna <- NULL

    outb.config <- list(n_iter = 1e5, sample_every = 200, find_import=FALSE,
                        move_kappa = FALSE, move_pi = FALSE, init_kappa = 1,
                        init_pi = 1, max_kappa = 1)

    dna.outb.result <- outbreaker2::outbreaker(dna.outb.data, outb.config)
    nodna.outb.result <- outbreaker2::outbreaker(nodna.outb.data, outb.config)

    store$dna.result[[i]] <- dna.outb.result
    store$nodna.result[[i]] <- nodna.outb.result

    store$dna.acc[[i]] <- get.acc(dna.outb.result, sim)
    store$nodna.acc[[i]] <- get.acc(nodna.outb.result, sim)
    store$uniq[[i]] <- get.uniq(sim)

    store$sim[[i]]$dna <- NULL

  }

  return(store)

}

## Returns the accuracy of transmission tree inference
get.acc <- function(res, sim, burnin = 1000) {

  inferred <- summary(res, burnin = burnin)$tree$from
  true <- sim$ances

  inferred[is.na(inferred)] <- 0
  true[is.na(true)] <- 0

  acc <- mean(inferred == true)

  return(acc)

}

## Returns the accuracy of transmission tree inference
get.phyb.acc <- function(res, sim) {

  mean(transtree(res, "edmonds")$infector == sim$sim.infectors)

}

## Calculate the entropy of a vector
calc.ent <- function(x) {
  x <- as.character(x)
  fk <- table(x)/sum(table(x))
  return(-sum(log(fk)*fk))
}

## Calculate the mean entropy of an outbreak run
get.ent <- function(res) {

  i <- grep('alpha', names(res))
  
  ent <- sapply(res[i], calc.ent)

  return(round(mean(ent), 2))
  
}

## Calculate the mean entropy of a phybreak object
get.phyb.ent <- function(res) {

  nsamples <- res$d$nsamples
  obsize <- res$p$obs
  samplerange <- seq_along(res$s$mu)

  ## Extract posterior distribution of ancestries (col = samples, row = id)
  mat <- res$s$nodehosts[nsamples:(nsamples + obsize - 1), samplerange]

  out <- round(mean(apply(mat, 1, calc.ent)), 2)

  return(out)
  
}

## Run gensig
gensig <- function(param = NULL, config = NULL) {

  param <- create.param(param)

  config <- create.config(config)

  store <- list()

  for(disease in names(param)) {

    sim <- run.sim(disease, param, config)

    store[[disease]][["sim"]] <- sim

    store[[disease]][["gensig"]] <- get.gensig(sim)

  }

  return(list(param = param, config = config, store = store))

}

## Describe epidemiological parameters of diseases
create.param <- function(param = NULL) {
  
  defaults <- list(
    ebola = list(R0 = 1.8 , mut = 3.10e-6 , seql = 18958   , w.mean = 14.4 , w.sd = 8.9,  dist = "gamma" ),
    sars  = list(R0 = 2.7 , mut = 1.14e-5 , seql = 29714   , w.mean = 8.7  , w.sd = 3.6,  dist = "gamma" ),
    mers  = list(R0 = 1.2 , mut = 0.25e-5 , seql = 30115   , w.mean = 10.7 , w.sd = 6.0,  dist = "gamma" ),
    ifz   = list(R0 = 1.5 , mut = 1.19e-5 , seql = 13155   , w.mean = 3.0  , w.sd = 1.5,  dist = "gamma" ),
    ##ifz.h = list(R0 = 1.5 , mut = 1.09e-5 , seql = 1701    , w.mean = 3.0  , w.sd = 1.5,  dist = "gamma" ),
    mrsa  = list(R0 = 1.3 , mut = 5.21e-9 , seql = 2842618 , w.mean = 15.6 , w.sd = 10.0, dist = "gamma" ),
    klebs = list(R0 = 2.0 , mut = 6.30e-9 , seql = 5305677 , w.mean = 62.7 , w.sd = 24.0, dist = "gamma" ),
    strep = list(R0 = 1.4 , mut = 5.44e-9 , seql = 2126652 , w.mean = 6.6  , w.sd = 1.8 , dist = "gamma" ),
    shig  = list(R0 = 1.1 , mut = 1.64e-9 , seql = 4825265 , w.mean = 7.0  , w.sd = 2   , dist = "gamma" ),
    tb    = list(R0 = 1.8 , mut = 2.36e-10, seql = 4411621 , w.mean = 324  , w.sd = 385 , dist = "gamma" ),
    cdif  = list(R0 = 1.5 , mut = 8.76e-10, seql = 4290252 , w.mean = 27.7 , w.sd = 14.9, dist = 'gamma' ))

  defaults$tb$w.sd <- defaults$tb$w.sd

  defaults$klebs$mut <- defaults$klebs$mut*7
  for(i in c("w.mean", "w.sd")) defaults$klebs[[i]] <- defaults$klebs[[i]]/7

  defaults$tb$mut <- defaults$tb$mut*7
  for(i in c("w.mean", "w.sd")) defaults$tb[[i]] <- defaults$tb[[i]]/7

  for(disease in names(defaults)) defaults[[disease]]$w <- 0
  
  if(!is.null(param)) {
    for(i in 1:length(param)) {
      if(length(param[[i]]) != 7) {
        stop("Incorrect number of parameter values provided")
      }
      if(any(!names(param[[i]]) %in% names(defaults[[1]]))) {
        stop("Incorrect parameter name provided")
      }
    }
  } else {
    param <- defaults
  }

  ## Calculate discretised distributions from the mean and standard deviations
  for(pathogen in names(param)) {
    tmp <- param[[pathogen]]
    if (tmp$dist == "weib") param[[pathogen]]$w <- discr.weib(tmp$w.mean, tmp$w.sd)
    if (tmp$dist == "gamma") param[[pathogen]]$w <- discr.gamma(tmp$w.mean, tmp$w.sd)
    
    ## Find first w that hits 0 (from decreasing slope), set the rest to 0
    ## Avoids confusion when w has lots of low but non-zero numbers
    ind <- which(c(FALSE, diff(param[[pathogen]]$w)) < 0 & param[[pathogen]]$w < 1e-10)
    if(length(ind) == 0) {
      tend <- length(param[[pathogen]]$w)
    } else {
      tend <- min(ind)
    }
    param[[pathogen]]$w[tend:length(param[[pathogen]]$w)] <- 1e-50
    
  }

  return(param)

}

## Describe parameter values for other runs
create.config <- function(...) {

  config <- list(...)
  if(length(config) == 1L && is.list(config[[1]])) {
    config <- config[[1]]
  }

  ## If no user labeller is provided, simply use the given names
  label <- function(char) return(char)

  defaults <- list(n.hosts = 200,
                   dur = 500,
                   imp = 0.05,
                   min.n = 50,
                   label = label)

  config <- modify.defaults(defaults, config)

  return(config)

}

## Returns a discretized gamma distribution
discr.gamma <- function(mean, sd) {
  w <- sapply(0:300, EpiEstim::DiscrSI, mean, sd)
  return(w)
}

## Returns a discretised weibull distribution
discr.weib <- function(mean, sd) {
  shape <- (sd/mean)^-1.086
  scale <- mean/gamma(1+1/shape)
  w <- stats::dweibull(0:300, shape = shape, scale = scale)
  return(w)
}

## Modifies default values of a function
modify.defaults <- function(default, modified) {

  not.found <- ! names(modified) %in% names(default)
  if(any(not.found)) stop(paste(paste(names(modified)[not.found], collapse = ", "),
                                "is not a valid descriptor"))

  for (i in names(modified)) default[[i]] <- modified[[i]]
  return(default)

}

## Run the simulations
run.sim <- function(disease, param, config) {

  tmp <- param[[disease]]
  
  n <- 0
  while(n < config$min.n) {
    sim <- outbreaker::simOutbreak(R0 = tmp$R0,
                                   infec.curve = tmp$w,
                                   seq.length = tmp$seql,
                                   mu.transi = tmp$mut*(2/3),
                                   mu.transv = tmp$mut*(1/3),
                                   n.hosts = config$n.hosts,
                                   duration = config$dur,
                                   rate.import.case = config$imp)

    sim$inf <- sim$onset
    
    sim$onset <- sim$onset + sample(0:300,
                                    length(sim$onset),
                                    replace = TRUE,
                                    prob = tmp$w)
    
    n <- sim$n
  }

  return(sim)

}

## Returns a function calculating the genetic signature between an id and its
## ancestor. The simulation 'sim' is enlosed as it remains unchanged, avoiding
## unnecessary passing of the sim argument
get.gensig <- function(sim) {

  ids <- sim$id[!is.na(sim$ances)]
  gensig <- sapply(ids, function(id) {
    dna <- sim$dna[c(id, sim$ances[id]),]
    gensig <- as.numeric(ape::dist.dna(dna, model="N"))
  })

  return(gensig)

}

## Returns a function calculating the genetic signature between an id and its
## ancestor, using a phybreakdata object
get.phyb.gensig <- function(sim) {

  ids <- which(sim$sim.infectors != "index")
  all.dna <- as.DNAbin(sim$sequences)
  gensig <- sapply(ids, function(id) {
    dna <- all.dna[c(id, which(sim$sample.hosts == sim$sim.infectors[id])),]
    gensig <- as.numeric(ape::dist.dna(dna, model="N"))
  })

  return(as.vector(gensig))

}

## Returns the proportion of unique sequences
get.uniq <- function(sim) {

  sim$dna <- phangorn::as.phyDat(sim$dna)
  length(unique(sim$dna))/sim$n
  
}

## Returns the proportion of unique sequences
get.phyb.uniq <- function(sim) {

  length(unique(sim$sequences))/length(sim$sample.hosts)

}

## Returns a vector of pathogen names ordered by mean genetic signature (high to low)
sort.gensig <- function(store) {

  means <- by(store$gensig$gensig, store$gensig$disease, mean)
  return(names(sort(-means)))

}

## Calculate the proportion of cases with gensig > 0
get.prop <- function(i) {
  mean(i > 0)
}

## Prop is the proportion of transmission pairs separated by at least one mutation
## Uniq is the number of unique sequences divided by the number of cases
## These numbers are directly correlated in a fully sampled tree with one introduction
## This IS assuming no reverse mutations though
prop2uniq <- function(prop, n) {

  (1 + prop*(n - 1))/n
  
}


##====== Collect cluster results =====##

## Connects to the cluster and returns the cluster object
set.up <- function(clust = 'mrc') {

  setwd("/home/fc1915/mnt/fc1915/gensig")

                                        #didewin::didewin_config_global(credentials="fc1915",cluster="fi--dideclusthn")
                                        #didewin::web_login()

  options(didehpc.cluster = "fi--dideclusthn",
          didehpc.credentials = "~/.smbcredentials")

  if(clust == "mrc") options(didehpc.cluster = "fi--didemrchnb")
  
  didehpc::didehpc_config_global(temp = didehpc::path_mapping('tmp',
                                                              '/home/fc1915/mnt/tmp',
                                                              '//fi--didef3.dide.ic.ac.uk/tmp',
                                                              'T:'))
  
  didehpc::web_login()

  ## Remove outbreaker2 if only using phybreak
  our.pkgs <- c('ggplot2', 'reshape2', 'outbreaker', 'phybreak', 'outbreaker2',
                'plyr', 'ape', 'scales', 'EpiEstim', 'magrittr', 'seedy', 'phangorn')

  pkg <- provisionr::package_sources(local = "~/mnt/fc1915/gensig/outbreaker2_1.0-0.tar.gz")
  
                                        #pkg <- provisionr::package_sources(local="~/mnt/fc1915/gensig/outbreaker2_1.0-0.tar.gz",
                                        #                                   github = 'donkeyshot/phybreak')

  our.sources <- "~/mnt/fc1915/gensig/functions.R"

  ctx <- context::context_save("contexts", packages = our.pkgs,
                               sources = our.sources,
                               package_sources = pkg)

  obj <- didehpc::queue_didehpc(ctx)

  return(obj)
}

## Create a storage vector for gensig values
create.store <- function(obj, bundle.name, dir, load = FALSE, dl = TRUE, phyb = FALSE) {

  cur.wd <- getwd()
  on.exit(setwd(cur.wd))

  store <- list(gensig = data.frame(disease = character(),
                                    gensig = numeric(),
                                    dat = character()),
                acc = data.frame(disease = character(),
                                 acc = numeric(),
                                 dat = character(),
                                 gensig = character()))

  files <- list.files(dir)
  rem <- grep("store", files)
  if(length(rem) > 0) files <- files[-grep("store", files)]
  n.files <- length(files)

  if(!phyb) adder <- add_r else adder <- add_phyb_r
  
  if(load) {
    pb <- txtProgressBar(min = 1, max = length(files), style = 3)
    for(i in seq_along(files)) {
      setTxtProgressBar(pb, i)
      load(paste0(dir, files[i]))
      store <- adder(store, r)
    }
  }

  if(dl) {
    task_bundle <- obj$task_bundle_get(bundle.name)
    ids <- task_bundle$ids
    ids <- ids[task_bundle$status() == "COMPLETE"]
    pb <- txtProgressBar(min = 1, max = length(ids), style = 3)
    for(i in seq_along(ids)) {
      setTxtProgressBar(pb, i)
      task <- obj$task_get(ids[i])
      r <- task$result()
      save(r, file = paste0(dir, "r.", n.files + i, ".RData"))
      store <- adder(store, r)
    }
  }

  for(i in names(store)) {

    store[[i]]$disease <- factor(store[[i]]$disease,
                                 levels = sort.gensig(store))

    store[[i]] <- filter(store[[i]], disease != 'ifz.h')

  }

  return(store)

}

## Adds r loaded from cluster or file to store
add_r <- function(store, r) {

  gensig <- data.frame(disease = names(r$param),
                       gensig = unlist(r$gensig))

  acc <- data.frame(disease = names(r$param),
                    dna = r$dna.acc,
                    nodna = r$nodna.acc)

  acc$improv <- with(acc, dna - nodna)
  acc$gensig <- ldply(r$gensig, mean)$V1
  acc$prop <- unlist(ldply(r$gensig, get.prop))
  ## Determine uniq empirically if it exists - otherwise do analytically from prop
  ## They should be the same in a single tree with one introduction for simOutbreak
  if(is.null(r$uniq)) {
    acc$uniq <- prop2uniq(acc$prop, sim$n)
  } else {
    acc$uniq <- r$uniq
  }
  acc$dna.ent <- sapply(seq_along(r$dna.res),
                        function(i) get.ent(r$dna.res[[i]]))
  acc$nodna.ent <- sapply(seq_along(r$nodna.res),
                          function(i) get.ent(r$nodna.res[[i]]))
  acc$w.neg <- sapply(seq_along(r$sim),
                      function(i) mean(get.w(r$sim[[i]]) < 0))
  acc$n <- sapply(seq_along(r$sim),
                  function(i) length(r$sim[[i]]$n))
  
  store$gensig <- rbind(store$gensig, gensig)
  store$acc <- rbind(store$acc, acc)

  return(store)

}

## Adds phyb_r loaded from cluster or file to store
add_phyb_r <- function(store, r) {

  gensig <- data.frame(disease = names(r$param),
                       gensig = unlist(r$phyb.gensig))

  acc <- data.frame(disease = names(r$param),
                    dna = r$phyb.dna.acc,
                    nodna = r$phyb.nodna.acc)

  acc$improv <- with(acc, dna - nodna)
  acc$gensig <- ldply(r$phyb.gensig, mean)$V1
  acc$prop <- unlist(ldply(r$phyb.gensig, get.prop))
  ## Determine uniq empirically if it exist - otherwise do analytically from prop
  ## They should be the same in a single tree with one introduction
  if(!is.null(r$phyb.uniq)) {
    acc$uniq <- r$phyb.uniq
  } else if(!is.null(r$phyb.sim)) {
    acc$uniq <- sapply(seq_along(r$phyb.sim),
                       function(i) get.phyb.uniq(r$phyb.sim[[i]]))
  }
  acc$dna.ent <- sapply(seq_along(r$phyb.dna.res),
                        function(i) get.phyb.ent(r$phyb.dna.res[[i]]))
  acc$nodna.ent <- sapply(seq_along(r$phyb.nodna.res),
                        function(i) get.phyb.ent(r$phyb.nodna.res[[i]]))
  acc$w.neg <- sapply(seq_along(r$phyb.sim),
                      function(i) mean(get.phyb.w(r$phyb.sim[[i]]) < 0))
  acc$n <- sapply(seq_along(r$phyb.sim), function(i) length(r$phyb.sim[[i]]$sample.hosts))
  
  store$gensig <- rbind(store$gensig, gensig)
  store$acc <- rbind(store$acc, acc)

  return(store)

}

## Get the mean generation time of a simOutbreak object
get.w <- function(sim) {

  out <- sim$onset - sim$onset[sim$ances]
  out <- out[!is.na(out)]
  return(out)

}

## Get the serial interval  of a phybreak simulation
get.phyb.w <- function(sim) {

  ances.inftime <- sim$sim.infection.times[match(sim$sim.infectors, sim$sample.hosts)]
  out <- sim$sim.infection.times - ances.inftime
  out <- as.vector(out[!is.na(out)])
  return(out)
  
}

## Input the index, and return the index of the ancestor [for seedy simulations]
get.seedy.ances <- function(id, sim) sim$epidata[which(sim$epidata[,1] == id), 4]

## Returns the empirical generation times of a seedy simulation
get.seedy.w <- function(sim) {

  unlist(sapply(seq_len(nrow(sim$epidata)),
                function(i) sim$epidata[i,2] - sim$epidata[get.seedy.ances(i, sim),2]))
  
}

## Returns the empirical genetic signature of a seedy simulation
get.seedy.gensig <- function(sim) {

  ## Everything is rescaled 1:N according to sim$sampledata
  ances <- unlist(sapply(sim$sampledata[,1],
                         function(i) {
                           tmp <- which(sim$sampledata[,1] == get.seedy.ances(i, sim))
                           if(length(tmp) == 0) return(NA)
                           else return(tmp)
                         }))
  
  tTree <- data.frame(i = seq_along(ances), ances = ances)

  ## The first argument gives which WGS id a given individual has
  ## The fourth argument links the WGS id to the actual sequences (args 2, 3)
  ## The resulting matrix is therefore in the order of sim$sampledata
  dist <- gd(sim$sampledata[,3],
             sim$libr,
             sim$nuc,
             sim$librstrains)

  ## Calculate genetic distance between transmission pairs
  apply(tTree, 1, function(i) dist[i[1], i[2]])

}

## print errors of a cluster bundle
get.error <- function(bundle) {
  
  ids <- bundle$ids[bundle$status() == 'ERROR']
  for(id in ids) print(obj$task_get(id)$result()$message)
  
}


##===== Auxilliary plotting functions =====##

## Create a labeller for linking variable names with disease labels
create.lab <- function(input) {

  c(ebola = "EBOV",
    sars = "SARS-CoV",
    mers = "MERS-CoV",
    ifz = "Influenza A",
    mrsa = "MRSA",
    klebs = "K. pneumoniae",
    strep = "S. pneumoniae",
    shig = "S. sonnei",
    tb = "M. tuberculosis",
    cdif = "C. difficile")

}

## Create a labeller for x and y axes from variable names
create.axlab <- function(input) {

  reference <- c(gensig = 'Average transmission divergence of outbreak',
                 improv = 'Change in accuracy of outbreak reconstruction',
                 dna = 'Accuracy of outbreak reconstruction',
                 nodna = 'Accuracy of outbreak reconstruction',
                 prop = 'Proportion of genetically distinct transmission pairs',
                 uniq = 'Number of unique sequences / outbreak size',
                 dna.ent = 'Entropy',
                 nodna.ent = 'Entropy',
                 n = 'Outbreak size')

  if(!input %in% names(reference)) {
    return(input)
  } else {
    return(reference[[input]])
  }
  
}

## Random package not downloaded properly
as.proto.list <- function(x, envir, parent, all.names = FALSE, ...,
                          funEnvir = envir, SELECT = function(x) TRUE) {
  if (missing(envir)) {
    if (missing(parent))
      parent <- parent.frame()
    envir <- if (is.proto(parent))
               parent$proto(...)
             else
               proto(parent, ...)
  }
  for (s in names(x))
    if (SELECT(x[[s]])) {
      assign(s, x[[s]], envir = envir)
      if (is.function(x[[s]]) && !identical(funEnvir, FALSE))
        environment(envir[[s]]) <- funEnvir
    }
  if (!missing(parent))
    parent.env(envir) <- parent
  as.proto.environment(envir)  # force refresh of .that and .super
}

as.proto.environment <- function(x, ...) {
  assign(".that", x, envir = x)
  assign(".super", parent.env(x), envir = x)
  structure(x, class = c("proto", "environment"))
}

as.lm <- function(object, ...) UseMethod("as.lm")

as.lm.nls <- function(object, ...) {
  if (!inherits(object, "nls")) {
    w <- paste("expected object of class nls but got object of class:",
               paste(class(object), collapse = " "))
    warning(w)
  }

  gradient <- object$m$gradient()
  if (is.null(colnames(gradient))) {
    colnames(gradient) <- names(object$m$getPars())
  }

  response.name <- if (length(formula(object)) == 2) "0" else
                                                           as.character(formula(object)[[2]])

  lhs <- object$m$lhs()
  L <- data.frame(lhs, gradient)
  names(L)[1] <- response.name

  fo <- sprintf("%s ~ %s - 1", response.name,
		paste(colnames(gradient), collapse = "+"))
  fo <- as.formula(fo, env = as.proto.list(L))

  do.call("lm", list(fo, offset = substitute(fitted(object))))

}


##===== Deprecated functions =====##

## Get the accuracy of outbreak reconstruction
get.old.acc <- function(result, outbreak) {

  id <- seq_len(outbreak$n)
  adder <- which(names(result)=="alpha_1") - 1
  samples <- length(result$step)

                                        #Determine the modal transmission network
  network <- data.frame(from=do.call(rbind, lapply(id,  function(i) ({
    modal.ances <- as.integer(names(which.max(table(result[[i+adder]]))))
    if(length(modal.ances)==0) return(NA) else return(modal.ances)
  }))),  to=id)

  import <- which(is.na(network$from))

  transmission.id <- id[!sapply(result[id+adder], function(i) any(is.na(i)))]

                                        #Determine the proportion of correctly inferred ancestries
  num.correct <-  sum(outbreak$ances==network$from, na.rm=TRUE)
  num.correct <- num.correct + sum(is.na(outbreak$ances[is.na(network$from)]))
  acc <- round(num.correct/nrow(network), 2)

  return(acc)

}

## Test the accuracy of inferring uniq from prop using simOutbreak
test.prop2uniq <- function() {

  tmp <- data.frame(true.uniq = rep(NA, 100),
                    ana.uniq = rep(NA, 100))
  n <- 1

  for(i in 1:100) {

    print(i)
    n <- 1
    while(n < 5) {
      sim <- outbreaker::simOutbreak(3, c(0, 1, 2, 3, 4, 5))
      n <- sim$n
    }
    tmp$true.uniq[i] <- length(unique(as.phyDat(sim$dna)))/sim$n

    prop <- mean(as.vector(na.omit(sim$nmut)) > 0)
    tmp$ana.uniq[i] <- prop2uniq(prop, sim$n)

  }

  plot(tmp$ana.uniq, tmp$true.uniq)

}


## ===== In development =====##

## Run phyb analysis on the cluster
run.seedy.cluster <- function(disease, min.n, n.hosts, runs, imp, dur) {
  
  config <- create.config(min.n = min.n,
                          n.hosts = n.hosts,
                          imp = imp,
                          dur = dur)

  param <- create.param()
  param <- list(tmp = param[[disease]])
  names(param) <- disease

  ## Find first w that hits 0 (from decreasing slop), set the rest to 0
  ## Avoids confusion when w has lots of low but non-zero numbers
  tend <- min(which(c(FALSE, diff(param[[1]]$w)) < 0 & param[[1]]$w < 1e-7))
  if(is.infinite(tend)) tend <- length(param[[1]]$w)
  param[[disease]]$w[tend:length(param[[disease]]$w)] <- 1e-50

  store <- list()

  mean.gen <- param[[disease]]$w.mean
  shape.gen <- param[[disease]]$w.sd^2
  mean.sample <- param[[disease]]$w.mean
  shape.sample <- param[[disease]]$w.sd^2
  
  rem.rate <- 1/param[[disease]]$w.mean
  inf.rate <- param[[disease]]$R0*rem.rate
  mut.rate  <- param[[disease]]$mut*param[[disease]]$seql
  
  for(i in seq_len(runs)) {

    ## Overlay sequence evolution on simOutbreak tTree (incredibly slow)
                                        #simo <- run.sim(disease, param, config)

                                        #simo$ances[is.na(simo$ances)] <- 0

                                        #simfix <- simfixoutbreak(ID = seq_along(simo$ances),
                                        #                         inf.times = simo$inf,
                                        #                         rec.times = simo$onset + tend,
                                        #                         inf.source = simo$ances,
                                        #                         mut.rate = mut.rate,
                                        #                         inoc.size = 1,
                                        #                         sample.times = simo$onset,
                                        #                         glen = param[[disease]]$seql)

    ## Simulate using SIR model (faster?)
    ##(but generation time then exponentially distributed)

    seql <- param[[disease]]$seql
    mut <- param[[disease]]$mut

    ## Rescale to make simulations feasible at larger sizes
    if(seql > 1e6) {
      seql <- round(seql/100, 0)
      mut <- mut*100
    }

    ## seedy::simulateoutbreak has been modified to draw sampling times
    ## from the w distribution (instead of at a fixed schedule)
    sim <- simulateoutbreak(init.sus = n.hosts,
                            inf.rate = inf.rate,
                            rem.rate = rem.rate,
                            mut.rate = seql*mut,
                            w = param[[disease]]$w,
                            init.inf = 1,
                            inoc.size = 1,
                            samples.per.time = 1,
                            samp.schedule = 'gentime',
                            mincases = min.n,
                            glen = param[[disease]]$seql,
                            full = FALSE)

    ## Everything is rescaled 1:N according to sim$sampledata
    ances <- unlist(sapply(sim$sampledata[,1],
                           function(i) {
                             tmp <- which(sim$sampledata[,1] == get.seedy.ances(i, sim))
                             if(length(tmp) == 0) return(NA)
                             else return(tmp)
                           }))
    
    tTree <- data.frame(i = seq_along(ances), ances = ances)

    ## The first argument gives which WGS-id a given individual has
    ## The fourth argument links the WGS id to the actual sequences (args 2, 3)
    ## The resulting matrix is therefore in the order of sim$sampledata
    dist <- gd(sim$sampledata[,3],
               sim$libr,
               sim$nuc,
               sim$librstrains)

    ## Calculate genetic distance between transmission pairs
    gensig <- apply(tTree, 1, function(i) dist[i[1], i[2]])

    ## Create sequences for phybreak / other analysis
                                        #ref.strain <- sample(1:4, param[[disease]]$seql, TRUE)
                                        #
                                        # dna <- librtoDNA(sim$sampledata[,3],
                                        #                  sim$libr,
                                        #                  sim$nuc,
                                        #                  ref.strain,
                                        #                  sim$librstrains,
                                        #                  strings = FALSE)
                                        # rownames(dna) <- tTree$i
                                        # dna <- as.DNAbin(dna)
                                        # onset <- sim$sampledata[,2]
                                        # 
                                        # phyb.dna.res <- phybreakdata(sequences = dna,
                                        #                              sample.times = onset) %>%
                                        #     phybreak(gen.mean = mean.gen,
                                        #              gen.shape = shape.gen,
                                        #              sample.mean = mean.sample,
                                        #              sample.shape = shape.sample,
                                        #              est.gen.mean = FALSE,
                                        #              est.sample.mean = FALSE,
                                        #              est.wh.slope = FALSE) %>%
                                        #     burnin.phybreak(ncycles = 5000) %>% 
                                        #     sample.phybreak(nsample = 5000)

    store$param <- param
    store$config <- config
    store$seedy.sim[[i]] <- sim
    store$seedy.gensig[[i]] <- gensig
                                        # store$phyb.dna.res[[i]] <- phyb.dna.res
                                        # store$phyb.dna.acc[[i]] <- mean(transtree(phyb.dna.res, "edmonds")$infector == tTree$ances)

  }

  return(store)
}

## Modified seedy::simulateoutbreak to take w as an argument
## for drawing sampling times
simulateoutbreak <- function(init.sus, inf.rate, rem.rate, mut.rate,
                             nmat=NULL, equi.pop=10000, shape=flat,  
                             init.inf=1, inoc.size=1, samples.per.time=1, samp.schedule="random", 
                             samp.freq=500, full=FALSE, mincases=1, feedback=500, glen=100000, 
                             ref.strain=NULL, w, ...) {
  
                                        # WARNINGS
  
  if (init.sus%%1!=0 || init.sus<1) {
    stop("Initial number of susceptibles must be a postive integer")
  }
  if (init.inf%%1!=0 || init.inf<1) {
    stop("Initial number of susceptibles must be a postive integer")
  }
  if (inoc.size%%1!=0 || inoc.size<1) {
    stop("Inoculum size must be a postive integer")
  }
  if (samp.freq%%1!=0 || samp.freq<1) {
    stop("samp.freq must be a postive integer")
  }
  if (feedback%%1!=0 || feedback<1) {
    stop("feedback must be a postive integer")
  }
  if (glen%%1!=0 || glen<1) {
    stop("Genome length must be a postive integer")
  }
  if (mincases%%1!=0 || mincases<1 || mincases>init.sus+init.inf) {
    stop("Minimum cases must be a postive integer not greater than init.sus+init.inf")
  }
  if (samples.per.time%%1!=0 || samples.per.time<1) {
    stop("samples.per.time must be a postive integer")
  }
  if (inf.rate<=0) {
    stop("Infection rate must be greater than zero")
  }
  if (rem.rate<=0) {
    stop("Removal rate must be greater than zero")
  }
  if (mut.rate<0 || mut.rate>=1) {
    stop("Mutation rate must be between 0 and 1")
  }
  if (!is.function(shape)) {
    stop("'shape' must be a function")
  } 
  if (!is.null(nmat)) {
    if (!is.matrix(nmat)) {
      stop("nmat must be a matrix")
    } else {
      if (nrow(nmat)!=init.sus+init.inf || ncol(nmat)!=init.sus+init.inf) {
        stop("nmat must have init.sus+init.inf rows and columns")
      } else if (sum(nmat<0)>0) {
        stop("All entries of nmat must be >= 0")
      }
    }
  }
  if (equi.pop%%1!=0 || equi.pop<=0) {
    stop("Equilibrium population size must be a postive integer")
  }
  if (!samp.schedule%in%c("random", "calendar", "individual", "gentime")) {
    stop("samp.schedule must be 'random', 'gentime', 'calendar', or 'individual'")
  }
  
#########################
  
  cat("\nSimulating outbreak:\n")
  cat("N=",equi.pop, ", b=", inoc.size, ", beta=", inf.rate, 
      ", gamma=", rem.rate, ", mu=", mut.rate, "\n\n", sep="")
  trigger <- FALSE
  
  at <- 0
  while (!trigger) { # Repeat if < mincases are infected
    newinfect <- 0
    at <- at+1
    cat("Attempt ", at, "\n", sep="")
    eff.cur.inf <- NULL
    cur.inf <- 1:init.inf # vector of infected person IDs
    cur.sus <- init.inf+(1:init.sus) # vector of susceptible person IDs
    time <- 0 # in bacterial generations
    inf.times <- rep(0,init.inf) # vector of infection times
    rec.times <- rgeom(init.inf,rem.rate)+2 # vector of removal times
    tot.inf <- init.inf
    inf.ID <- 1:init.inf
    if (samp.schedule=="random") {
      sample.times <- NULL
      for (i in 1:init.inf) {
        sample.times <- c(sample.times, sample(1:(rec.times[i]-1),1)) # vector of sampling times
      }
    } else if (samp.schedule=="gentime") {
      sample.times <- NULL
      for (i in 1:init.inf) {
        mod.w <- w[1:(rec.times[i])]
        sample.from <- time:(time+length(mod.w)-1)
        sample.times <- c(sample.times, sample(sample.from,1,replace=T,mod.w))
      }
    } else {
      sample.times <- rep(samp.freq, init.inf)
    }
    inf.source <- rep(0,init.inf) # source of infection for each individual
    if (is.null(ref.strain)) {
      ref.strain <- sample(1:4, glen, replace=T) # reference strain
    } else {
      glen <- length(ref.strain)
    }
    totcurstrains <- 1 # current list of strains
    uniquestrains <- 1 # Number of unique strain types
    
    libr <- list() # list of mutation locations for each genotype
    mut.nuc <- list() # nucleotides at mutation locations
    freq.log <- list() # List of strain frequencies for each infective
    strain.log <- list() # Strain IDs for each within host population
    
                                        # Initialize logs
    for (i in 1:init.inf) {
      if (sample.times[i]>rec.times[i]) {
        sample.times[i] <- Inf
      }
      if (i == 1) {
        libr[[i]] <- NA
        mut.nuc[[i]] <- NA      
      } else {
        libr[[i]] <- sample(glen,1)
        mut.nuc[[i]] <- sample((1:4)[-ref.strain[libr[[i]]]], 1)
      }
      freq.log[[i]] <- 1
      strain.log[[i]] <- 1
    }
    
    for (i in (init.inf+1):(init.sus+init.inf)) {
      freq.log[[i]] <- 0
      strain.log[[i]] <- 0
    }
    
    current.infected <- init.inf
    types <- 1 # Cumulative number of strain types
    
                                        #Sample logs
    if (full) {
      obs.freq <- list()
      obs.strain <- list()
      pID <- NULL
    } else {
      sampleWGS <- NULL
      samplepick <- NULL
    }
    
    sampletimes <- NULL
    sampleID <- NULL
    
    while (length(cur.inf) > 0) { # Cycle through bacterial generations until epidemic ceases
      time <- time+1
      if (time%in%rec.times) { # recovery?
        recover <- inf.ID[which(rec.times==time)] # who has recovered?
        cur.inf <- cur.inf[-which(cur.inf%in%recover)] # remove infective(s)
        for (r in 1:length(recover)) {
          strain.log[[recover[r]]] <- 0
          freq.log[[recover[r]]] <- 0
        }
        if (length(cur.inf)==0) { # If no more infectives
          cat("t=", time, ", S=", length(cur.sus), ", I=", length(cur.inf), 
              ", total genotypes=0\n", sep="")
        }
      }
      if (time%%feedback==0 && length(cur.inf)>0) { # output current status every x generations
        cat("t=", time, ", S=", length(cur.sus), ", I=", length(cur.inf), 
            ", total genotypes=", length(unique(as.numeric(unlist(strain.log)))), ", next rec time=", 
            min(rec.times[which(rec.times>time)]), sep="")
        if (length(cur.sus)>0) {
          cat("\n")
        } else {
          cat(", final removal time=", max(rec.times), "\n", sep="")
        }
      }
                                        # calculate force of infection
      
      if (is.null(nmat)) {
        curinfrate <- inf.rate*length(cur.inf)/(init.inf+init.sus)
        pinf <- rbinom(length(cur.sus),1,curinfrate)
      } else {
        inf.mat <- nmat*inf.rate/(init.inf+init.sus)
        if (length(cur.inf)>1) {
          curinfrate <- apply(inf.mat[cur.sus,cur.inf],1,sum)
        } else {
          curinfrate <- inf.mat[cur.sus,]
        }
        pinf <- rbinom(length(cur.sus),1,curinfrate)
      }
      
                                        #if (runif(1,0,1) < curinfrate) { # infection?
      if (sum(pinf)>0) {
        for (kp in 1:sum(pinf)) {
          tot.inf <- tot.inf+1
                                        # who got infected? and by whom?
          newinfect <- cur.sus[which(pinf==1)[kp]]
          
          if (length(cur.inf)==1) {
            inf.source <- c(inf.source, cur.inf)
          } else if (is.null(nmat)) {
            inf.source <- c(inf.source, sample(cur.inf,1)) # sample source at random
          } else {
            inf.source <- c(inf.source, sample(cur.inf, 1, prob=inf.mat[newinfect,cur.inf]))
          }
          
          inf.ID <- c(inf.ID, newinfect)
          cur.inf <- c(cur.inf, newinfect) # add to current infectives
          cur.sus <- cur.sus[-which(cur.sus==newinfect)] # Remove susceptible
          inf.times <- c(inf.times, time) # new infection time
          rec.times <- c(rec.times, time+max(1,rgeom(1,rem.rate))) # Don't recover today!
          if (samp.schedule == "individual") {
            sample.times <- c(sample.times, time+samp.freq)
          } else if (samp.schedule == "calendar") {
            sample.times <- c(sample.times, ceiling(time/samp.freq)*samp.freq)
          } else if (samp.schedule == "random") {
            if (rec.times[tot.inf]>time+1) {
              sample.times <- c(sample.times, sample(time:(rec.times[tot.inf]-1),1))
            } else {
              sample.times <- c(sample.times, time)
            }
          } else if (samp.schedule == "gentime") {
            if (rec.times[tot.inf]>time+1) {
              ## modify w so it doesn't exceed the recovery time
              mod.w <- w[1:(rec.times[tot.inf]-time)]
              sample.from <- time:(time+length(mod.w)-1)
              sample.times <- c(sample.times, sample(sample.from,1,replace=T,mod.w))
            } else {
              sample.times <- c(sample.times, time)
            }
          }
          if (sample.times[tot.inf]>=rec.times[tot.inf]) {
            sample.times[tot.inf] <- Inf
          }
                                        # pass on strain
          src <- inf.ID[which(inf.ID==inf.source[tot.inf])] # Source of infection
          if (length(strain.log[[src]])==1) { # if source has clonal infection
            inoc.samp <- rep(strain.log[[src]], inoc.size)
            if (0%in%inoc.samp) {
              stop("Zeroes in inoculum")
            }
          } else {
            inoc.samp <- sample(strain.log[[src]], inoc.size, 
                                prob=freq.log[[src]], replace=T) # take random sample
            if (0%in%inoc.samp) {
              stop("Zeroes in inoculum")
            }
          }
          strain.log[[newinfect]] <- unique(inoc.samp) # distinct types in new infection
          f <- numeric(length(unique(inoc.samp)))
          k <- 1
          for (i in unique(inoc.samp)) {
            f[k] <- sum(inoc.samp==i)
            k <- k+1
          }
          freq.log[[newinfect]] <- f # frequency of types
        }
      }
      
                                        # mutate existing strains for each individual
      if (is.null(eff.cur.inf)) {
        cinf <- inf.ID[which(inf.ID%in%cur.inf)]
      } else {
        cinf <- eff.cur.inf
      }
      for (i in cinf) {
        pop.size <- sum(freq.log[[i]])
        death.prob <- min(0.5 + 0.5*(pop.size-shape(time,span=rec.times[i]-inf.times[i]+1,equi.pop,...))/shape(time,span=rec.times[i]-inf.times[i]+1,equi.pop,...),1)
        if (length(freq.log)==0 || sum(is.na(freq.log))>0 || death.prob>1 || death.prob<0) {
          cat("deathprob=", death.prob, "\npop.size=", pop.size, "\nequi.pop=", equi.pop, "\nFreq.log:\n")
          if (length(freq.log)>0) {
            for (k in 1:length(freq.log)) {
              cat(freq.log[[i]][k], "\n")
            }
          }
        }
        freq.log[[i]] <- 2*rbinom(length(freq.log[[i]]), freq.log[[i]], 1-death.prob)
        if (0 %in% freq.log[[i]]) {
          zeros <- which(freq.log[[i]]==0)
          if (length(zeros)==length(freq.log[[i]])) {
            freq.log[[i]] <- 1
            strain.log[[i]] <- strain.log[[i]][1]
          } else {
            freq.log[[i]] <- freq.log[[i]][-zeros]
            strain.log[[i]] <- strain.log[[i]][-zeros]
          }
        }
        if (length(strain.log[[i]])!=length(freq.log[[i]])) {
          stop("Error")
        }
        n.mutations <- rbinom(1, sum(freq.log[[i]]), mut.rate)
        if (n.mutations > 0) {
          for (mt in 1:n.mutations) {
            types <- types+1
            if (length(strain.log[[i]])==1) {
              mutate.grp <- strain.log[[i]]
            } else {
              mutate.grp <- sample(strain.log[[i]], 1, prob=freq.log[[i]])
            }
            if (mutate.grp %in% totcurstrains) {
              mut.loc <- sample(glen, 1)
              mut.nuc[[length(totcurstrains)+1]] <- 
                mut.nuc[[which(totcurstrains==mutate.grp)]]
              if (mut.loc %in% libr[[which(totcurstrains==mutate.grp)]]) { # if mutation at existing location
                kn <- which(libr[[which(totcurstrains==mutate.grp)]]==mut.loc)
                mut.nuc[[length(totcurstrains)+1]][kn] <- sample((1:4)[-mut.nuc[[which(totcurstrains==mutate.grp)]][kn]], 1)
                libr[[length(totcurstrains)+1]] <- libr[[which(totcurstrains==mutate.grp)]]
              } else {
                mut.nuc[[length(totcurstrains)+1]] <- 
                  c(mut.nuc[[length(totcurstrains)+1]], 
                    sample((1:4)[-ref.strain[mut.loc]], 1))
                libr[[length(totcurstrains)+1]] <- 
                  c(libr[[which(totcurstrains==mutate.grp)]], mut.loc)
              }
              if (sum(is.na(mut.nuc[[length(totcurstrains)+1]]))>0) {
                mut.nuc[[length(totcurstrains)+1]] <- mut.nuc[[length(totcurstrains)+1]][-is.na(mut.nuc[[length(totcurstrains)+1]])]
                libr[[length(totcurstrains)+1]] <- libr[[length(totcurstrains)+1]][-is.na(libr[[length(totcurstrains)+1]])]
              }
              strain.log[[i]] <- c(strain.log[[i]], types)
              freq.log[[i]] <- c(freq.log[[i]], 1)
              totcurstrains <- c(totcurstrains, types)
            }
          }
        }
      }

      if(length(totcurstrains) == 0) totcurstrains <- 0
                                        # take samples, make observations
      if (time%in%sample.times) {
        smpat <- inf.ID[which(sample.times==time)]
        for (i in smpat) {
          if (full) {
            n <- length(obs.freq)+1
            obs.freq[[n]] <- freq.log[[i]]
            obs.strain[[n]] <- strain.log[[i]]
            sampleID <- c(sampleID, n)
            pID <- c(pID, i)
            sampletimes <- c(sampletimes, time)
          } else {
            for (j in 1:samples.per.time) {
              if (length(strain.log[[i]])==1) {
                pickgrp <- strain.log[[i]]
              } else {
                pickgrp <- sample(strain.log[[i]], 1, prob=freq.log[[i]])
              }
              sampleWGS <- c(sampleWGS, pickgrp)
              if (0%in%sampleWGS) {
                stop("Sampled zeroes")
              }
              sampleID <- c(sampleID, i)
              samplepick <- c(samplepick, j)
              sampletimes <- c(sampletimes, time)
            }
          }
          if (samp.schedule != "random" && rec.times[which(inf.ID==i)] > time + samp.freq) {
            sample.times[which(inf.ID==i)] <- time + samp.freq
          }
        }
      }
                                        # clean up libr etc.
      if (!full) {
        deleters <- NULL
        uniquestrains <- 0
        for (j in 1:length(totcurstrains)) {
          tottype <- 0
          for (k in cur.inf) {
            if (totcurstrains[j]%in%strain.log[[k]]) { # if strain is extant
              tottype <- tottype+1
            }
            if (sum(!strain.log[[k]]%in%totcurstrains)>0) {
              stop("Deleted sequence for observed sample")
            }
          }
          if (tottype>0) { # don't delete if still around
            uniquestrains <- uniquestrains+1
          } else if (tottype==0 && !totcurstrains[j]%in%sampleWGS) { # if not around AND not logged
            deleters <- c(deleters, j) # delete
          }
        }
        if (length(deleters)>0) {
          for (i in sort(deleters, decreasing=T)) {
            libr[[i]] <- NULL
            mut.nuc[[i]] <- NULL
          }
          deletegroup <- totcurstrains[deleters]
          totcurstrains <- totcurstrains[-deleters]
        }
      }
      if (length(cur.sus)==0) {
        eff.cur.inf <- inf.ID[which(sample.times>time)]
      }
      if (length(cur.sus)==0 && sum(sample.times>time)==0) {
        break
      }
    }
    if (tot.inf>mincases) {
      trigger <- TRUE
    } else {
      cat("Insufficient number of infections! (mincases=", mincases, ")\n", sep="")
    }
  }
  if (full) {
    return(invisible(list(epidata=cbind(inf.ID, inf.times, rec.times, inf.source), 
                          sampledata=cbind(pID, sampleID, sampletimes), obs.freq=obs.freq, obs.strain=obs.strain,
                          libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time)))
  } else {
    return(invisible(list(epidata=cbind(inf.ID, inf.times, rec.times, inf.source), 
                          sampledata=cbind(sampleID, sampletimes, sampleWGS),
                          libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time)))
  }
}

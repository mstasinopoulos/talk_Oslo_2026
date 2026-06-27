
################################################################################
################################################################################
################################################################################
################################################################################
chooseDistCV <- function(object,
                  type = c("realAll", "realline", "realplus","real0to1",
                           "counts", "binom", "extra" ), 
                 extra = NULL,  # for extra distributions to include 
                 trace = FALSE,
                  data = NULL,
              parallel = c("no", "multicore", "snow"), #
                 ncpus = 1L, #typically one would chose this to the number of available CPUs
                    cl = NULL, # An optional parallel or snow cluster for use if parallel = "snow".
                K.fold = 10,
                  rand = NULL,
                           ...)
{
## get the type of distribution
#    newData <- if(is.null(newdata)) FALSE else TRUE
       type <- match.arg(type)
   sys_call <- sys.call() 
if (is.null(data))   stop("data should be set here")
          N <- dim(data)[1]
       rand <- if (is.null(rand)) sample(K.fold , N, replace=TRUE)       
#if (is.null(rand)&&is.null(newdata)) stop("rand or newdata should be set")
       DIST <- switch(type, "realAll"=.realAll,
                 "realline"=.realline,
                 "realplus"=.realplus,
                 "real0to1"=.real0to1,
                 "counts"=.counts,
                 "binom"=.binom,
                 "extra"= extra)
if (type=="extra"&&is.null(extra)) stop("extra is not set")
if  (!is.null(extra)) DIST <- unique(c(DIST, extra))
       m0 <- object
     tgd0 <- deviance(m0)
     AiC  <- rep(NA, 1)
#--------------- PARALLEL-------------------------------------------------------
#----------------SET UP PART----------------------------------------------------
  parallel <- match.arg(parallel)
   have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) 
  {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  } 
# -------------- finish parallel------------------------------------------------     
# define the function
  fun <- function(dist)
  {
# cat(dist, "\n")
    p1 <- try(gamlssCV(formula(object, "mu"), formula(object, "sigma"),
                       formula(object, "nu"), formula(object, "tau"), rand=rand,
                                data=data, trace=F,family=dist), silent=TRUE)
    if (any(class(p1)%in%"try-error"))
    {
    p1 <- try(gamlssCV(formula(object, "mu"), formula(object, "sigma"),
                       formula(object, "nu"), formula(object, "tau"), rand=rand, 
                       data=data, trace=F,family=dist), silent=TRUE)
    }
    if (any(class(p1)%in%"try-error")) AiC <- NA else{
      AiC <- measure_of_goodness_of_fit(p1) 
      if (trace)     cat(dist, "\n", AiC, "\n") 
    }
      AiC        # output of the function
  }
################################################################################ 
# --------  parallel -----------------------------------------------------------
  MM <- if (ncpus > 1L && (have_mc || have_snow)) 
  {
    if (have_mc) 
    {# sapply(scope, fn)
      unlist(parallel::mclapply(DIST, fun, mc.cores = ncpus))
    }
    else if (have_snow) 
    {
      list(...)
      if (is.null(cl)) 
      { # make the cluster
        if (.Platform$OS.type == "windows")
        {
          cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
          clusterEvalQ(cl,pacman::p_load(gamlss)) 
          exp.data =  paste0(object$call$data)
          clusterExport(cl, c(ls(envir = .GlobalEnv), exp.data))
        } else cl <- parallel::makeForkCluster(ncpus)
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <-  unlist((parallel::parLapply(cl, DIST, fun)))
        parallel::stopCluster(cl)
        res
      } 
      else parallel::parLapply(cl, DIST, fun)# use existing cluster
    }
  }# end parallel -----
  else sapply(DIST, fun, ...) 
#-------------------------------------------------------------------------------  
   names(MM) <- DIST
#-------------------------------------------------------------------------------
  ## save it in the final model
  #m0$TGD <- MM[order(MM)]
  MM[order(MM)]
}
################################################################################
################################################################################
################################################################################
################################################################################
#source("~/Dropbox/github/GlossaryForDR/Chi2.R")
#source("~/Dropbox/2026/IWSM/analysis/Chi2.R")
# FUN <- function(object)  cal_chi2(resid(object), xvar=da$ga, breaks = c(0,17,25,34.5, 50), tail.percent = 5,  side = c("upper"))

# measure_of_goodness_of_fit <- function(object)
# {
#   cal_chi2(object$residCV, xvar=ALB$ga, breaks = c(0,17,25,34.5, 50), tail.percent = 5,  side= "upper")
# }
# measure_of_goodness_of_fit <- function(object)
# {
#   object$CV
# }
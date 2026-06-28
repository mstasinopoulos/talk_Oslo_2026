# the function is only appropriate if there is an important xvar 
# TO DO: i) what if there is no xvariable?
################################################################################
################################################################################
################################################################################
################################################################################  
cal_chi2 <- function(resid, xvar, 
                     breaks = c(0,17,25,34.5, 50),
               tail.percent = 5, 
                       side = c("upper", "lower", "both"))
{
    side <- match.arg(side)  
# if there is x   
      ff <- cut(xvar, breaks=breaks)
if (side=="upper")
  {
      In <- (resid  > qnorm(1-(tail.percent/100)))
  } 
if (side=="lower")
  {
      In <- (resid  < qnorm(tail.percent/100))
  } 
if (side=="both")
  {
      In <- (resid  < qnorm(tail.percent/100)|(resid  > qnorm(1-(tail.percent/100))))
  } 
      O <-  table(ff,In)# this is a table 
     ttm <- table(ff)
    tttm <- c(ttm, ttm)
      ll <- dim(O)[1]
       O <- as.vector(O) # we wand it as a vector
if (side=="lower")
  {# calculating the expectation 
     Pr1 <- rep((1-(tail.percent/100)),  ll)
     Pr2 <- rep((tail.percent/100),  ll)   
      Pr <- c(Pr1,Pr2)
       E <- tttm*Pr
} 
if (side=="upper")
       {# calculating the expectation 
        Pr1 <- rep((1-(tail.percent/100)),  ll)
        Pr2 <- rep((tail.percent/100),  ll)   
         Pr <- c(Pr1,Pr2)
          E <- tttm*Pr
       }        
# the value of chisq
     chi <- ((O-E)^2)/E
 chi.val <- sum(chi)
 return(chi.val)
# list(chi.val, chi, O, E)
}
################################################################################ 
################################################################################ 
################################################################################ 
################################################################################  
resid_quantile_weight <- function(resid, weights, cent = c(5,10,25,50, 75,90, 95 ), 
                       xvar, 
                     breaks = c(0,17,25,34.5, 50))
               #tail.percent = 5, 
                      # side = c("upper", "lower", "both"))
{
         N <- length(resid)
if ((missing(weights))) weights <- rep(1, length(resid))
if (missing(xvar)) # no x
{
         S <- rep(0, length(cent))
         I <- list()
  cent_val <- quantile(resid, probs = cent / 100)
         for (i in 1:length(cent)) {
             I[[i]] <- (resid <= cent_val[i])
             S[i] <- sum(weights[I[[i]]])
         }
}
         S <- (S/N)*100
  names(S) <- as.character(cent)
return(S)   
} 

#        S <- matrix(0, nrow=length(breaks), ncol=length(cent))
#        I <- list()
#       ff <- cut(xvar, breaks=breaks)
#       for (i in 1:length(cent)) {
#         I[[i]] <- (resid <= cent_val[i])
#         S[i] <- sum(weights[I[[i]]])
#   return(chi.val)
# }
################################################################################ 
################################################################################ 
################################################################################ 
################################################################################  
weighted_quantile <- function(x, weights, quantile=0.50 ) {
  ord <- order(x)
  x <- x[ord]
  w <- weights[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= quantile)[1]]
}
################################################################################
################################################################################
################################################################################
################################################################################
weighted_median_freq <- function(x, freq) {
   ord <- order(x)
     x <- x[ord]
  freq <- freq[ord]
     n <- sum(freq)
    cf <- cumsum(freq)
  if (n %% 2 == 1) {
    # odd sample size
    x[which(cf >= (n + 1) / 2)[1]]
  } else {
    # even sample size
    m1 <- x[which(cf >= n / 2)[1]]
    m2 <- x[which(cf >= n / 2 + 1)[1]]
    (m1 + m2) / 2
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
#weighted_median_freq(x, freq)

# plot(resid(m1))
# q50 <- weighted_quantile(resid(m1), weights=da$wt, quantile=0.50)
# q05 <- weighted_quantile(resid(m1), weights=da$wt, quantile=0.05)
# q95 <- weighted_quantile(resid(m1), weights=da$wt, quantile=0.95)
# abline(h=q50, col="red")
# abline(h=q05, col="yellow")
# abline(h=q95, col="yellow")
# 
# 
# ECDF(resid(m1), weights=da$wt)


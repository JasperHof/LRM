MGres <- function(fitme = NULL, data = NULL){             # Compute martingale residuals

  MGres_check(fitme, data)

  # Check the names of 'strata' columns in the data
  check_strat = strsplit(as.character(fitme$formulaList$fixed)[3],'strata')[[1]]

  if(length(check_strat) > 1){        # In this case, strata have been specified..
    name = substring(strsplit(check_strat[2], ')')[[1]][1],2)
    strats = data[,which(colnames(data) == name)]
  } else strats = rep(0, dim(data)[1])

  strat_list = unique(strats)                       # Vector of strata names
  cumhaz = rep(0, length(unique(data$subject)))     # Vector of cumulative hazards per subject
  events = as.data.frame(as.matrix(fitme$y))        # Matrix with event times and status

  for(strat in strat_list){     # For every stratum, compute the cumulative baseline hazard for all individuals within this stratum
    in_strat = which(strats == strat)

    # ---- DETERMINE THE BASELINE HAZARD FOR THIS STRATUM ----

    basehaz = rep(0,length = dim(events)[1])  # Base hazard for all possible time points

    if(dim(events)[2] > 2){                   # In case tstart is specified (gap time/calendar time model)
      eventtimes = events$stop

      for(i in 1:length(eventtimes)){
        time = eventtimes[i]
        risk.set = which(events$start < time & events$stop >= time & strats == strat)

        denom = sum(exp(fitme$linear.predictor[risk.set]))
        nom = sum(events$stop == time & events$status == 1 & strats == strat)

        basehaz[i] = nom/denom
      }
    } else {                             # If tstart not specified, gap time
      eventtimes = events$time

      for(i in 1:length(eventtimes)){
        time = eventtimes[i]
        risk.set = which(events$time >= time & strats == strat)     # No constraint on gap time

        denom = sum(exp(fitme$linear.predictor[risk.set]))
        nom = sum(events$time == time & events$status == 1 & strats == strat)

        basehaz[i] = nom/denom
      }
    }

    # ---- COMPUTE INDIVIDUAL CUMULATIVE HAZARDS FOR THIS STRATUM ----

    for(i in 1:length(unique(data$subject))){
      rows = which(data$subject == unique(data$subject)[i] & strats == strat)      # The 'rows' of individual i in the stratum

      if(length(rows) > 0){
        for(r in 1:length(rows)){
          if(dim(events)[2] > 2){
            haz = sum(basehaz[which(events$start[rows[r]] < events$stop & events$stop <= events$stop[rows[r]])])
            cumhaz[i] = cumhaz[i] + haz
          } else {
            haz = sum(basehaz[which(events$time <= events$time[rows[r]])])
            cumhaz[i] = cumhaz[i] + haz
          }
        }
      }
    }
  }

  individual_prop_haz = rep(0, length(unique(data$subject)))
  count = rep(0, length(unique(data$subject)))

  for(i in 1:length(unique(data$subject))){   # HIER GEBLEVEN !!!

    set = which(data$subject == unique(data$subject)[i])
    individual_prop_haz[i] = fitme$linear.predictor[set[1]]
    count[i] = sum(events$status[set] == 1)
  }

  lambda = cumhaz * exp(individual_prop_haz)  # Cumulative hazard
  resids = count - lambda                     # Martingale residuals
  names(resids) = unique(data$subject)

  return(resids)
}

#
# FUNCTION FOR ASSESSING NULL MODEL + SPA
#

Null_model <- function(fitme,
                       data=NULL,
                       IDs=NULL,
                       range=c(-20,20),
                       length.out = 50000,
                       ...)
{
  Call = match.call()

  ### Compute martingale residuals
  mresid = MGres(fitme, data)

  ### Check input arguments
  obj.check = check_input(data, IDs, mresid, range)

  ### Calculate empirical CGF for martingale residuals
  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)

  cumul = NULL
  print("Start calculating empirical CGF for martingale residuals...")
  c = 0
  for(i in idx1){
    c = c+1
    t = i
    e_resid = exp(mresid*t)
    M0 = mean(e_resid)
    M1 = mean(mresid*e_resid)
    M2 = mean(mresid^2*e_resid)
    K0 = log(M0)
    K1 = M1/M0
    K2 = (M0*M2-M1^2)/M0^2
    cumul = rbind(cumul, c(t, K0, K1, K2))
    if(c %% 5000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
  }

  K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2)
  K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2)
  K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2)

  re=list(resid = mresid,
          K_org_emp = K_org_emp,
          K_1_emp = K_1_emp,
          K_2_emp = K_2_emp,
          Call = Call,
          IDs = IDs)

  class(re)<-"NULL_Model"
  return(re)
}

#
# Do Linear Regression on Martingale residuals (LRM)
#

LRM = function(obj.null,
               Geno.mtx,
               missing.cutoff = 0.05,
               min.maf = 0.05,
               p.cutoff = 0.001)
{
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf)

  ### Check input
  check_input_LRM(obj.null, Geno.mtx, par.list)

  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  IDs = rownames(Geno.mtx)
  SNPs = colnames(Geno.mtx)

  #
  ### ----------------- Quality Control ---------------
  #

  ### Only select genotypes that we also have a phenotype of
  mresid = obj.null$resid
  Complete = intersect(names(mresid), IDs)
  Geno.mtx = Geno.mtx[which(rownames(Geno.mtx) %in% Complete),]

  ### Filter on call rate, minor allele frequency,
  SNP.callrate = colMeans(is.na(Geno.mtx))
  if(any(SNP.callrate > missing.cutoff)){ Geno.mtx = Geno.mtx[,-which(SNP.callrate > missing.cutoff)] }

  No_Geno <- is.na(Geno.mtx)
  Geno.mtx = na_mean(Geno.mtx)                # Impute missing values to mean for estimating MAF
  G = Geno.mtx[,which(colMeans(Geno.mtx) > 2*min.maf & colMeans(Geno.mtx < 2*(1 - min.maf)))]
  No_Geno = No_Geno[,which(colMeans(Geno.mtx) > 2*min.maf & colMeans(Geno.mtx < 2*(1 - min.maf)))]  # Also update the No_geno matrix

  mresid = mresid[which(names(mresid) %in% Complete)]
  MG = mresid[match(rownames(G), names(mresid))]

  # ------------------------------------------------------------------------------
  # ------ Couple Genotype G to martingale residuals MG and do the analysis -----------
  # ------------------------------------------------------------------------------

  print("Start linear regression..")
  print(Sys.time())

  # Define outcome ytr and intercept X
  n = dim(G)[1]
  X = as.matrix(rep(1,dim(G)[1]), nrow = dim(G)[1],ncol=1)
  ytr = as.matrix(MG)

  # Scale genotypes to have mean 0
  U3 = crossprod(X, G)
  U4 = solve(crossprod(X), U3)
  Str = G - as.matrix(X %*% U4)

  ## Compute slopes
  if(sum(No_Geno) > 0)  Str[No_Geno] <- 0                             # Exclude the entries from computation by setting them to zero
  b = crossprod(ytr, Str)/colSums(Str ^ 2)

  ## Compute P values in LRM
  S_sq = colSums(Str ^ 2)
  RSS = crossprod(ytr^2, !No_Geno) - b ^ 2 * S_sq
  sigma_hat = RSS/(n - 2)
  error = sqrt(sigma_hat/ S_sq)
  pval.mg = as.numeric(2 * pnorm(-abs (b / error)))

  # Prepare the output file
  b = as.numeric(b); error = as.numeric(error)

  outcome = cbind(SNP = colnames(G), MAF = colMeans(G)/2, Missing = colSums(is.na(G)),
                  pSPA = pval.mg, pMG = pval.mg, Beta = b, SE = error, Z = abs(b / error))
  outcome = as.data.frame(outcome, stringsAsFactors = F)
  outcome$pMG = as.numeric(outcome$pMG)
  outcome$pSPA = as.numeric(outcome$pSPA)

  #
  # ------ use SPA for SNPs with P value below cut-off point ------
  #

  low.p = which(pval.mg < p.cutoff)

  for(i in low.p){
    g = G[,i]
    g = g - mean(g)

    s = sum(g * MG)/sum(g^2)             # Statistic, which still needs to be normalized..?

    k0 = function(x) sum(obj.null$K_org_emp(g/sum(g^2) * x))
    k1 = function(x) sum(g/sum(g^2) * obj.null$K_1_emp(g/sum(g^2) * x))
    k2 = function(x) sum((g/sum(g^2))^2 * obj.null$K_2_emp(g/sum(g^2) * x))

    get_p = function(s, tail){
      k1_root = function(x) k1(x) - s
      zeta = uniroot(k1_root, c(-2000,2000))$root

      w = sign(zeta) * sqrt(2 * (zeta * s - k0(zeta)))
      v = zeta * sqrt(k2(zeta))

      pval = pnorm(w + 1/w * log(v/w), lower.tail = tail)
      return(pval)
    }

    p_SPA = get_p(abs(s), tail = F) + get_p(-abs(s), tail = T)
    outcome$pSPA[i] = p_SPA
  }

  print("Analysis Complete.")
  print(Sys.time())

  return(outcome)
}

#
# LRM for .bgen files
#

LRM.bgen <- function(bgenfile, gIDs,
                     obj.null,
                     chr = NULL,
                     missing.cutoff = 0.05,
                     min.maf = 0.05,
                     p.cutoff = 0.001){



  bgifile = paste0(bgenfile,'.bgi')

  ### Create 'myid' variable for reading .bgen files
  db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgifile)
  on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
  infos <- dplyr::collect(dplyr::tbl(db_con, "Variant"))
  infos$myid <- with(infos, paste(chromosome, position, allele1,
                                  allele2, sep = "_"))
  info = list()
  info[[1]] = infos$myid

  ### Set up chunk sizes
  size = dim(infos)[1]
  matrixsize = 5e6
  chunksize = floor(matrixsize/length(obj.null$resid))
  reps = ceiling(size/chunksize)
  outcome = NULL

  ### Analyze .bgen files per chunk
  for(r in 1:reps){

    print(paste0('Reading chunk ', r, ' of ',reps))
    indices = c(((r-1)*chunksize + 1):min(r*chunksize, size))

    if(is.null(chr)) chr = as.numeric(infos$chromosome[indices[1]])
    backfile = paste0(getwd(),'/tmpfile',chr,'_',r)

    snps = list(); snps[[1]] = infos$myid[indices]         # The names of the SNPs in this chunk
    bgen = snp_readBGEN(bgenfiles = bgenfile, list_snp_id = snps, backingfile = backfile, ncores = 1)

    file = paste0(backfile,'.rds')
    genotype <- snp_attach(file)
    G <- genotype$genotypes

    Geno.mtx = G[,]
    colnames(Geno.mtx) = infos$myid[indices]

    if(length(gIDs) != dim(Geno.mtx)[1]) stop("Length of gIDs not equal to number of samples in .bgen file.")
    rownames(Geno.mtx) = gIDs

    # Executing the linear regression
    outcome = rbind(outcome, LRM(obj.null = obj.null,
                                 Geno.mtx = Geno.mtx,
                                 missing.cutoff = missing.cutoff,
                                 min.maf = min.maf,
                                 p.cutoff = p.cutoff))

    ### Clean up directory by removing previous connections
    for(k in 1:r){
      prev.file = paste0(getwd(),'/tmpfile',chr,'_',k)
      unlink(prev.file)
      prev.file = paste0(getwd(),'/tmpfile',chr,'_',k)
      unlink(prev.file)
    }
  }
  return(outcome)
}



#
# CHECKING LRM input
#

check_input_LRM = function(obj.null, Geno.mtx, par.list)
{
  if(class(obj.null)!="NULL_Model")
    stop("obj.null should be a returned outcome from SPACox_Null_Model()")

  if(is.null(rownames(Geno.mtx))) stop("Row names of 'Geno.mtx' should be given.")
  if(is.null(colnames(Geno.mtx))) stop("Column names of 'Geno.mtx' should be given.")
  if(!is.numeric(Geno.mtx)|!is.matrix(Geno.mtx)) stop("Input 'Geno.mtx' should be a numeric matrix.")

  if(length(intersect(obj.null$IDs,rownames(Geno.mtx))) == 0) stop("None of 'IDs' are included in rownames(Geno.mtx).")
  print(paste0('In total, ', length(intersect(obj.null$IDs,rownames(Geno.mtx))),' samples with phenotype and genotype information'))

  if(!is.numeric(par.list$min.maf)|par.list$min.maf<0|par.list$min.maf>0.5) stop("Argument 'min.maf' should be a numeric value >= 0 and <= 0.5.")
  if(!is.numeric(par.list$missing.cutoff)|par.list$missing.cutoff<0|par.list$missing.cutoff>1) stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
}

#
# CHECKING INPUT NULL MODEL
#

check_input = function(data, IDs, mresid, range)
{
  if(is.null(IDs))
    stop("Argument 'IDs' is required in case of potential errors. For more information, please refer to 'Details'.")
  IDs = as.character(IDs)

  if(any(!is.element(IDs, data$subject)))
    stop("All elements in IDs should be also in data as 'subject'")

  if(anyDuplicated(IDs)!=0)
    stop("Argument 'IDs' should not have a duplicated element.")

  if(range[2]!=-1*range[1])
    stop("range[2] should be -1*range[1]")

  if(length(mresid)!=length(IDs))
    stop("length(mresid)!=length(pIDs) where mresid is the martingale residuals from coxme package.")

}

#
# FUNCTION FOR CHECKING INPUT COXME OBJECT AND DATA OBJECT
#


MGres_check = function(fitme, data){
  if(is.null(fitme)) stop('no coxme object included')
  if(class(fitme) != 'coxme') stop('object not of class coxme')
  if(is.null(data)) stop('no data object included')
  if(class(data) != 'data.frame') stop('data is not a data frame object')
  if(!'subject' %in% colnames(data)) stop('please include individuals as "subject" in dataframe')

  check_strat = strsplit(as.character(fitme$formulaList$fixed)[3],'strata')[[1]]
  if(length(check_strat) > 2){
    stop('do not include stratum/covariates with name strata')
  }

  if(length(check_strat) > 1){        # In this case, strata have been specified..
    name = substring(strsplit(check_strat[2], ')')[[1]][1],2)

    strats = data[,which(colnames(data) == name)]
    try(if(length(strats) != dim(data)[1]) stop('please include the strata once in the data frame'))
  } else {
    strats = rep(0, dim(data)[1])
  }
}

#
# Function for imputing missing values to mean
#

na_mean = function (x, option = "mean", maxgap = Inf) {
  data <- x
  if (!is.null(dim(data)[2]) && dim(data)[2] > 1) {
    for (i in 1:dim(data)[2]) {
      if (!anyNA(data[, i])) {
        next
      }
      tryCatch(data[, i] <- na_mean(data[, i], option,
                                    maxgap), error = function(cond) {
                                      warning(paste("imputeTS: No imputation performed for column",
                                                    i, "because of this", cond), call. = FALSE)
                                    })
    }
    return(data)
  }
  else {
    missindx <- is.na(data)
    if (!anyNA(data)) {
      return(data)
    }
    if (any(class(data) == "tbl")) {
      data <- as.vector(as.data.frame(data)[, 1])
    }
    if (all(missindx)) {
      stop("Input data has only NAs. Input data needs at least 1 non-NA data point for applying na_mean")
    }
    if (!is.null(dim(data)[2]) && !dim(data)[2] == 1) {
      stop("Wrong input type for parameter x")
    }
    if (!is.null(dim(data)[2])) {
      data <- data[, 1]
    }
    else if (option == "mean") {
      mean <- mean(data, na.rm = TRUE)
      data[missindx] <- mean
    }
    if (is.finite(maxgap) && maxgap >= 0) {
      rlencoding <- rle(is.na(x))
      rlencoding$values[rlencoding$lengths <= maxgap] <- FALSE
      en <- inverse.rle(rlencoding)
      data[en == TRUE] <- NA
    }
    if (!is.null(dim(x)[2])) {
      x[, 1] <- data
      return(x)
    }
    return(data)
  }
}

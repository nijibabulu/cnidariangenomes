## Adapted from
## GenomeScope: Fast Genome Analysis from Unassembled Short Reads
##
## This is the automated script for computing genome characteristics
## from a kmer histogram file, k-mer size, and readlength
##
## See https://github.com/schatzlab/genomescope fo rthe original code


## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
options(warn=-1)

## Colors for plots


## Given mean +/- stderr, report min and max value within 2 SE
###############################################################################

min_max <- function(table){
	##return (c( abs(table[1]) - 2*abs(table[2]) , abs(table[1])+ 2*abs(table[2])))
	return (c(table[1] - 2*table[2], table[1]+ 2*table[2]))
}



## Use nls to fit 4 peak model
###############################################################################

nls_4peak<-function(x, y, k, estKmercov, estLength, max_iterations, VERBOSE=F){
	model4 = NULL

    if (VERBOSE) { cat("trying nls_4peak standard algorithm\n") }

	try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                          (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                          (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                          (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                      start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

	if(class(model4) == "try-error"){
        if (VERBOSE) { cat("retrying nls_4peak with port algorithm\n") }

        try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                              (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length +
                              (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length +
                              (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length),
                          start = list(d=0, r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                          algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)
	}

	return(model4)
}


## score model by number and percent of residual errors after excluding sequencing errors
#########################################################################################

score_model<-function(kmer_hist_orig, nls, round, VERBOSE = F){
  x = kmer_hist_orig[[1]]
  y = kmer_hist_orig[[2]]

  pred=predict(nls, newdata=data.frame(x))
  model_sum=summary(nls)
  kcovfloor = floor(min_max(model_sum$coefficients['kmercov',])[[1]])

  ## Compute error rate, by counting kmers unexplained by model through first peak
  ## truncate errors as soon as it goes to zero, dont allow it to go back up
  error_xcutoff = kcovfloor
  error_xcutoff_ind = which(x==error_xcutoff)

  error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]

  first_zero = -1

  for (i in 1:error_xcutoff_ind)
  {
    if (first_zero == -1)
    {
      if (error_kmers[i] < 1.0)
      {
        first_zero = i
        if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
      }
    }
    else
    {
      error_kmers[i] = 0
    }
  }

  if (first_zero == -1)
  {
    first_zero = error_xcutoff_ind
  }

  ## The fit is residual sum of square error, excluding sequencing errors
  model_fit_all    = c(sum(as.numeric(y[first_zero:length(y)]     - pred[first_zero:length(y)])     ** 2), first_zero, x[length(y)])
  model_fit_full   = c(sum(as.numeric(y[first_zero:(5*kcovfloor)] - pred[first_zero:(5*kcovfloor)]) ** 2), first_zero, (5*kcovfloor))
  model_fit_unique = c(sum(as.numeric(y[first_zero:(3*kcovfloor)] - pred[first_zero:(3*kcovfloor)]) ** 2), first_zero, (3*kcovfloor))

  ## The score is the percentage of unexplained kmers, excluding sequencing errors
  model_fit_allscore    = c(1-sum(abs(as.numeric(y[first_zero:length(y)]     - pred[first_zero:length(y)])))     / sum(as.numeric(y[first_zero:length(y)])),     first_zero, x[length(y)])
  model_fit_fullscore   = c(1-sum(abs(as.numeric(y[first_zero:(5*kcovfloor)] - pred[first_zero:(5*kcovfloor)]))) / sum(as.numeric(y[first_zero:(5*kcovfloor)])), first_zero, (5*kcovfloor))
  model_fit_uniquescore = c(1-sum(abs(as.numeric(y[first_zero:(3*kcovfloor)] - pred[first_zero:(3*kcovfloor)]))) / sum(as.numeric(y[first_zero:(3*kcovfloor)])), first_zero, (3*kcovfloor))

  fit = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
                   full = model_fit_full,     fullscore = model_fit_fullscore,
                   unique = model_fit_unique, uniquescore = model_fit_uniquescore)

  return (fit)
}


## Pick between the two model forms, resolves ambiguity between which is the homozygous and which is the heterozygous peak
###############################################################################

eval_model<-function(kmer_hist_orig, nls1, nls2, round, VERBOSE=FALSE){
    nls1score = -1
    nls2score = -1

    ## Evaluate the score the nls1
    if (!is.null(nls1))
    {
      nls1score = score_model(kmer_hist_orig, nls1, round+0.1, VERBOSE = VERBOSE)

      if(VERBOSE){ cat(paste("nls1score$all:\t", nls1score$all[[1]], "\n"))}

      #if (VERBOSE)
      #{
        #mdir = paste(foldername, "/round", round, ".1", sep="")
        #dir.create(mdir, showWarnings=FALSE)
        #report_results(kmer_prof_orig,kmer_prof_orig, k, (list(nls1, nls1score)) , mdir)
      #}
    }
    else
    {
      if (VERBOSE) { cat("nls1score failed to converge\n") }
    }


    ## Evaluate the score of nls2
    if (!is.null(nls2))
    {
      nls2score = score_model(kmer_hist_orig, nls2, round+0.2, VERBOSE = VERBOSE)

      if(VERBOSE){ cat(paste("nls2score$all:\t", nls2score$all[[1]], "\n"))}

      #if (VERBOSE)
      #{
        #mdir = paste(foldername, "/round", round, ".2", sep="")
        #dir.create(mdir, showWarnings=FALSE)
        #report_results(kmer_prof_orig, kmer_prof_orig, k, (list(nls2, nls2score)) , mdir)
      #}
    }
    else
    {
      if (VERBOSE) { cat("nls2score failed to converge\n") }
    }

    ## Return the better of the scores
    if (!is.null(nls1))
    {
      if (!is.null(nls2))
      {
        pdiff = abs(nls1score$all[[1]] - nls2score$all[[1]]) / max(nls1score$all[[1]], nls2score$all[[1]])

        if (pdiff < SCORE_CLOSE)
        {
          het1 = summary(nls1)$coefficients['r',][[1]]
          het2 = summary(nls2)$coefficients['r',][[1]]

          if (het2 * SCORE_HET_FOLD_DIFFERENCE < het1)
          {
            if (VERBOSE) { cat(paste("returning nls1, similar score, higher het\n")) }
            return (list(nls1, nls1score))
          }
          else if (het1 * SCORE_HET_FOLD_DIFFERENCE < het2)
          {
            if (VERBOSE) { cat(paste("returning nls2, similar score, higher het\n")) }
            return (list(nls2, nls2score))
          }
        }

        if (nls1score$all[[1]] < nls2score$all[[1]])
        {
          if (VERBOSE) { cat(paste("returning nls1, better score\n")) }
          return (list(nls1, nls1score))
        }
        else
        {
          if (VERBOSE) { cat(paste("returning nls2, better score\n")) }
          return (list(nls2, nls2score))
        }
      }
      else
      {
        if (VERBOSE) { cat(paste("returning nls1, nls2 fail\n")) }
        return (list(nls1, nls1score))
      }
    }

    if (VERBOSE) { cat(paste("returning nls2 by default\n")) }
    return (list(nls2, nls2score))
}


## Wrapper function to try fitting 4 peak model with 2 forms
###############################################################################

estimate_Genome_4peak2<-function(kmer_hist_orig, x, y, k, readlength, round, MAX_ITERATIONS, VERBOSE=F){
	## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
	numofReads   = sum(as.numeric(x*y))/(readlength-k+1)
	estKmercov1  = x[which(y==max(y))][1]
	estCoverage1 = estKmercov1*readlength/(readlength-k)
	estLength1   = numofReads*readlength/estCoverage1

    if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov1, "\n")) }
	nls1    = nls_4peak(x, y, k, estKmercov1, estLength1, MAX_ITERATIONS, VERBOSE = VERBOSE)
    if (VERBOSE) { print(summary(nls1)) }

	## Second we half the max kmercoverage (typically the heterozygous peak)
	estKmercov2  = estKmercov1 / 2 ##2.5
	estCoverage2 = estKmercov2*readlength/(readlength-k)
	estLength2   = numofReads*readlength/estCoverage2

    if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov2, "\n")) }
	nls2 = nls_4peak(x, y, k, estKmercov2, estLength2, MAX_ITERATIONS, VERBOSE = VERBOSE)
    if (VERBOSE) { print(summary(nls2)) }

	return(eval_model(kmer_hist_orig, nls1, nls2, round))
}


## Format numbers
###############################################################################
bp_format<-function(num) {
  paste(formatC(round(num),format="f",big.mark=",", digits=0), "bp",sep=" ")
}

percentage_format<-function(num) {
  paste(signif(num,6)*100,"%",sep="")
}
X_format<-function(num) {
  paste(signif(num,4),"X",sep="")
}


## Report results and make plots
###############################################################################

plot_genomescope <- function(res, VERBOSE=F, title = "GenomeScope Result", subtitle = "k-mer Curve Genome Length Estimate", title_face = "bold.italic")
{

  if(is.null(res$container)) {
    stop("No convergence. Cannot report results.")
  }
  kmer_hist = res$kmer_prof
  kmer_hist_orig = res$kmer_prof_orig
  k = res$k
  container = res$container
  trim = res$trim

    x=kmer_hist_orig[[1]]
    y=kmer_hist_orig[[2]]

	#automatically zoom into the relevant regions of the plot, ignore first 15 positions
    xmax=length(x)
	start=which(y == min(y[1:trim]))
	zoomx=x[start:(xmax-1)]
	zoomy=y[start:(xmax-1)]

    ## allow for a little space above max value past the noise
	y_limit = max(zoomy[start:length(zoomy)])*1.1

	x_limit = which(y == max(y[start:length(zoomx)])) * 3

  	if (min(zoomy) > zoomy[1]){
  		x_limit=max(which(zoomy<zoomy[1])[2],600)
  	}

    if (!is.null(container[[1]]))
    {
       model_sum=summary(container[[1]])
       kcov = min_max(model_sum$coefficients['kmercov',])[1]
       x_limit = max(kcov*5.1, x_limit)
    }

    ## Uncomment this to enforce a specific number
    # x_limit=150

    ## Features to report
    het=c(-1,-1)
    total_len=c(-1,-1)
    repeat_len=c(-1,-1)
    unique_len=c(-1,-1)
    dups=c(-1,-1)
    error_rate=c(-1,-1)
    model_status="fail"

    model_fit_unique      = c(0,0,0)
    model_fit_full        = c(0,0,0)
    model_fit_all         = c(0,0,0)
    model_fit_allscore    = c(0,0,0)
    model_fit_fullscore   = c(0,0,0)
    model_fit_uniquescore = c(0,0,0)

    plot_size=2000
    font_size=1.2
    resolution=300

	if(!is.null(container[[1]]))
    {
       x=kmer_hist[[1]]
       y=kmer_hist[[2]]

       ## The model converged!
       pred=predict(container[[1]], newdata=data.frame(x))

       ## Compute the genome characteristics
       model_sum=summary(container[[1]])

       ## save the model to a file
       # capture.output(model_sum, file=paste(foldername,"/model.txt", sep=""))

       ## Identify key values
       het  = min_max(model_sum$coefficients['r',])
       dups = min_max(model_sum$coefficients['bias',])
       kcov = min_max(model_sum$coefficients['kmercov',])
       mlen = min_max(model_sum$coefficients['length',])
       md   = min_max(model_sum$coefficients['d',])

       amlen = (mlen[1] + mlen[2]) / 2
       ahet  = (het[1]  + het[2])  / 2
       amd   = (md[1]   + md[2])   / 2
       akcov = (kcov[1] + kcov[2]) / 2
       adups = (dups[1] + dups[2]) / 2

       ## Compute error rate, by counting kmers unexplained by model through first peak
       ## truncate errors as soon as it goes to zero, dont allow it to go back up
       error_xcutoff = floor(kcov[1])
       error_xcutoff_ind = which(x==error_xcutoff)

       error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]

       first_zero = -1

       for (i in 1:error_xcutoff_ind)
       {
         if (first_zero == -1)
         {
           if (error_kmers[i] < 1.0)
           {
             first_zero = i
             if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
           }
         }
         else
         {
           error_kmers[i] = 0
         }
       }

       if (first_zero == -1)
       {
         first_zero = error_xcutoff_ind
       }

       ## Rather than "0", set to be some very small number so log-log plot looks okay
       error_kmers = pmax(error_kmers, 1e-10)

      total_error_kmers = sum(as.numeric(error_kmers * x[1:error_xcutoff_ind]))
       total_kmers = sum(as.numeric(x*y))

       f1 <- function(x){
             i=seq(1,k)
             h=(1-x)^(k-i)*x^i*choose(k,i)
             sum(h)*total_kmers-total_error_kmers
       }

       error_rate_root = try(uniroot(f1, c(0,1))$root)

       if (class(error_rate_root) == "try-error")
       {
         error_rate  = c(total_error_kmers/total_kmers/k, total_error_kmers/total_kmers/k)
       }
       else
       {
          error_rate  = c(error_rate_root, error_rate_root)
       }

       total_len = (total_kmers-total_error_kmers)/(2*kcov)

       ## find kmers that fit the 2 peak model (no repeats)
       unique_hist <- (2 * (1 - amd) * (1 - (1 - ahet)^k))                         * dnbinom(x, size = akcov     / adups, mu = akcov)     * amlen +
                      ((amd * (1 - (1 - ahet)^k)^2) + (1 - amd) * ((1 - ahet)^k))  * dnbinom(x, size = akcov * 2 / adups, mu = akcov * 2) * amlen

       unique_kmers = sum(as.numeric(x*unique_hist))
       repeat_kmers = total_kmers - unique_kmers - total_error_kmers

       repeat_len=repeat_kmers/(2*kcov)
       unique_len=unique_kmers/(2*kcov)

       score = container[[2]]

       model_fit_allscore    = score$allscore
       model_fit_fullscore   = score$fullscore
       model_fit_uniquescore = score$uniquescore

       model_fit_all    = score$all
       model_fit_full   = score$full
       model_fit_unique = score$unique

       residual = y - pred

       ## Finish Log plot
       #title(paste("\nlen:",  prettyNum(total_len[1], big.mark=","),
                   #"bp",
                   #" uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                   #"% ",
                   #" het:",  format(100*ahet, digits=3),
                   #"%",
                   #" kcov:", format(akcov, digits=3),
                   #" err:",   format(100*error_rate[1], digits=3),
                   #"% ",
                   #" dup:",  format(adups, digits=3),
                   #"% ",
                   #" k:",   format(k, digits=3),
                   #sep=""),
                   #cex.main=.85)

       ## Mark the modes of the peaks
       COLOR_BGCOLOR  = "white"
       COLOR_HIST     = "#56B4E9"
       COLOR_4PEAK    = "black"
       COLOR_2PEAK    = "#F0E442"
       COLOR_ERRORS   = "red"
       COLOR_KMERPEAK = "black"
       COLOR_RESIDUAL = "purple"
       COLOR_COVTHRES = "red"
       LWD = 1.2

       kmer_observed <-
         kmer_hist_orig %>%
         tibble::as_tibble() %>%
         dplyr::rename(Coverage = V1, Frequency = V2)

       models <-
         tibble::tibble(Unique = unique_hist,
                       FullModel = pred) %>%
         dplyr::mutate(Error = c(error_kmers, rep(NA, nrow(.) - length(error_kmers)))) %>%
         tibble::rowid_to_column("Coverage") %>%
         tidyr::pivot_longer(-Coverage, names_to = "Models", values_to = "Frequency") #%>%
         #dplyr::mutate(Model = factor(Model, levels = c("Error", "Unique", "FullModel")))

       peaks <-
         tibble::tibble(name = c("Homozygous", "Heterozygous", "3x", "4x"),
                        Coverage = akcov * seq(4))

       ## Finish Linear Plot
       #title(paste("\nlen:",  prettyNum(total_len[1], big.mark=","),
       #"bp",
       #" uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
       #"% ",
       #" het:",  format(100*ahet, digits=3),
       #"%",
       #" kcov:", format(akcov, digits=3),
       #" err:",   format(100*error_rate[1], digits=3),
       #"% ",
       #" dup:",  format(adups, digits=3),
       #"% ",
       #" k:",   format(k, digits=3),
       #sep=""),
       #cex.main=.85)

       caption <- stringr::str_glue("k = {k}, Homozygous Coverage = {format(round(akcov))}\nHeterozygosity = {format(100*ahet, digits=3)}%\nLength Estimate = {prettyNum(mean(total_len), big.mark = ',')} ")

       p <- ggplot2::ggplot(kmer_observed, ggplot2::aes(x = Coverage, y = Frequency)) +
         ggplot2::geom_col(ggplot2::aes(fill = "Observed"), width = 1) +
         ggplot2::coord_cartesian(xlim = c(0, x_limit), ylim = c(0, y_limit)) +
         ggplot2::geom_line(ggplot2::aes(color = Models), alpha = 0.7, data = models,  lwd = LWD, inherit.aes = T) +
         ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0,.05))) +
         ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0)) +
         ggplot2::scale_color_manual(values = c(COLOR_ERRORS,    COLOR_2PEAK, COLOR_4PEAK)) +
         ggplot2::scale_fill_manual(values = COLOR_HIST) +
         ggplot2::scale_linetype_manual(values = "dotted") +
         ggplot2::geom_vline(ggplot2::aes(xintercept = Coverage, linetype = "Peaks"), data = peaks) +
         ggplot2::guides(fill = ggplot2::guide_legend(title = NULL),
                         linetype = ggplot2::guide_legend(title = NULL)) +
         ggplot2::theme_bw() +
         ggplot2::labs(title = title,
                       subtitle = subtitle,
                       caption = caption) +
         ggplot2::theme(panel.grid = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(face = title_face, hjust = 0.5),
                        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(0.7)),
                        plot.caption  = ggplot2::element_text( hjust = 0.5))
       return(p)
	}
    ## Finalize the progress
    #progressFilename=paste(foldername,"/progress.txt",sep="")

}




## Main program starts here
###############################################################################

## Number of rounds before giving up
#NUM_ROUNDS=4

## Coverage steps to trim off between rounds
#START_SHIFT=5

## Typical cutoff for sequencing error
#TYPICAL_ERROR = 15

## Max rounds on NLS
#MAX_ITERATIONS=20

## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
#SCORE_CLOSE = 0.20

## Overrule heterozygosity if there is a large difference in het rate
#SCORE_HET_FOLD_DIFFERENCE = 10

## Print out VERBOSEging messages (0/1)
#VERBOSE = 0

genomescope <- function(histfile, k, readlength, maxCovGenomeLen = -1, NUM_ROUNDS = 4, START_SHIFT = 5, TYPICAL_ERROR = 15,
                        MAX_ITERATIONS = 20, SCORE_CLOSE = 0.20, SCORE_HET_FOLD_DIFFERENCE = 10, VERBOSE = F) {

    ## values for testing
    #histfile <- "~/build/genomescope/simulation/simulation_results/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"
    #k <- 21
    #readlength <- 100
    #foldername <- "~/build/genomescope/simulation/simulation_analysis/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"

    if (k > readlength) { stop("K cannot be greater than readlength") }

    cat(paste("GenomeScope analyzing ", histfile, " k=", k, " readlen=", readlength,  "\n", sep=""))

	#dir.create(foldername, showWarnings=FALSE)

	kmer_prof <- read.csv(file=histfile,sep=" ", header=FALSE)

    minkmerx = 1;
    if (kmer_prof[1,1] == 0) {
        if (VERBOSE) { cat("Histogram starts with zero, reseting minkmerx\n");  }
        minkmerx = 2;
    }

	kmer_prof <- kmer_prof[c(minkmerx:(length(kmer_prof[,2])-1)),] #get rid of the last position
    kmer_prof_orig <- kmer_prof

    ## Initialize the status
    #progressFilename <- paste(foldername,"/progress.txt",sep="")
	#cat("starting", file=progressFilename, sep="\n")

    ## try to find the local minimum between errors and the first (heterozygous) peak
    start <- which(kmer_prof[,2]==min(kmer_prof[1:TYPICAL_ERROR,2]))

    maxCovIndex = -1

    ## Figure out which kmers to exclude, if any
    if(maxCovGenomeLen == -1){
        maxCovIndex <- length(kmer_prof[,1])
    }
    else
    {
        ## Figure out the index we should use for this coverage length
        x <- kmer_prof[,1]
        maxCovIndex <- length(x[x<=maxCovGenomeLen])
    }

    if (VERBOSE) { cat(paste("using maxCovGenomeLen:", maxCovGenomeLen, " with index:", maxCovIndex, "trying 4peak model... \n")) }

    ## terminate after NUM_ROUND iterations, store best result so far in container
	round <- 0
	best_container <- list(NULL,0)

	while(round < NUM_ROUNDS)
    {
        cat(paste("round", round, "trimming to", start, "trying 4peak model... ")) #, file=progressFilename, sep="", append=TRUE)
        if (VERBOSE) { cat(paste("round", round, "trimming to", start, "trying 4peak model... \n")) }

        ## Reset the input trimming off low frequency error kmers
        kmer_prof=kmer_prof_orig[1:maxCovIndex,]
        x <- kmer_prof[start:maxCovIndex,1]
        y <- kmer_prof[start:maxCovIndex,2]

        model_4peaks <- estimate_Genome_4peak2(kmer_prof, x, y, k, readlength, round, MAX_ITERATIONS = MAX_ITERATIONS, VERBOSE = VERBOSE)

        if (!is.null(model_4peaks[[1]])) {
          cat(paste("converged. score: ", model_4peaks[[2]]$all[[1]]), sep="\n") #, file=progressFilename,  append=TRUE)

          #if (VERBOSE)
          #{
            #mdir = paste(foldername, "/round", round, sep="")
	        #dir.create(mdir, showWarnings=FALSE)
	        #report_results(kmer_prof,kmer_prof_orig, k, model_4peaks, mdir)
          #}
        } else {
          cat(paste("unconverged \n"), sep="\n") #, file=progressFilename, append=TRUE)
        }

		#check if this result is better than previous
        if (!is.null(model_4peaks[[1]]))
        {
          if (is.null(best_container[[1]]))
          {
            if (VERBOSE) { cat(paste("no previous best, updating best")) }
            best_container = model_4peaks
            best_start = start
          }
          else
          {
            pdiff = abs(model_4peaks[[2]]$all[[1]] - best_container[[2]]$all[[1]]) / max(model_4peaks[[2]]$all[[1]], best_container[[2]]$all[[1]])

            if (pdiff < SCORE_CLOSE)
            {
              hetm = summary(model_4peaks[[1]])$coefficients['r',][[1]]
              hetb = summary(best_container[[1]])$coefficients['r',][[1]]

              if (hetb * SCORE_HET_FOLD_DIFFERENCE < hetm)
              {
                if (VERBOSE) { cat(paste("model has significantly higher heterozygosity but similar score, overruling")) }
                best_container = model_4peaks
                best_start = start
              }
              else if (hetm * SCORE_HET_FOLD_DIFFERENCE < hetb)
              {
                if (VERBOSE) { cat(paste("previous best has significantly higher heterozygosity and similar score, keeping")) }
              }
              else if (model_4peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
              {
                if (VERBOSE) { cat(paste("score is marginally better but het rate is not extremely different, upating")) }
                best_container = model_4peaks
                best_start = start
              }
            }
            else if (model_4peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
            {
              if (VERBOSE) { cat(paste("score is significantly better, upating")) }
              best_container = model_4peaks
              best_start = start
            }
          }
        }

        ## Ignore a larger number of kmers as errors
        start <- start + START_SHIFT
		round <- round + 1
	}

	list(kmer_prof = kmer_prof, kmer_prof_orig = kmer_prof_orig, k = k, container = best_container, trim = start)
}

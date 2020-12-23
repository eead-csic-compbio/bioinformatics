#' @title Scatter plot of background versus query k-mer frequencies. 
#' @author Najla Nksouri <najlaksouri@gmail.com> 
#' Revised by Jacques van Helden <Jacques.van-Helden@univ-amu.fr>
#' @description Draw a scatter plot comparing background and query k-mer frequencies for a given size. 
#' @param kmer.stat.table a table with k-mer statistics produced by kmer.stats()
#' @param k oligonucleotide (k-mer) size
#' @param sig.threshold=2 Threshold on significance for the volcano plot. 
#' @param volcano.ymax=350
#' @return a data.frame with merged query and background k-mer frequencies + derived statistics
#' @examples
#' ## Plot query versus background 4-mer frequencies
#' plot.kmers(k=8)
#' 
plot.kmers <- function(kmer.stat.table,
                       k,
                       sig.threshold = 2,
                       log2FC.threshold = 1, ## Enrichment threshold
                       volcano.ymax = 350) {
  message(k, "-mer frequency plot")

  ## Define dot colors  
  dot.colors <- ifelse(
    (kmer.stat.table$sig > sig.threshold & kmer.stat.table$log2.ratio > 0.1), "darkred",
         ifelse((kmer.stat.table$sig > sig.threshold & kmer.stat.table$log2.ratio < -0.1),
                "darkblue", "darkgray"))

  ## Draw a scatter plot to compare k-mer frequencies
  max.freq <- max(kmer.stat.table[, c("obs_freq_query","obs_freq_background")])
  plot(x = kmer.stat.table$obs_freq_background,
       y = kmer.stat.table$obs_freq_query,
       
       ## Set axis limits
       xlim = c(0, max.freq),
       ylim = c(0, max.freq),
       
       ## Set labels
       xlab = paste("background ", k,"-mer frequencies", sep = ""),
       ylab = paste("query ", k,"-mer frequencies", sep = ""),
       main = paste("background vs query ", k,"-mer frequencies", sep = ""),
       
       col = dot.colors,
       
       ## Insert grid
       panel.first=grid(),
       
       ## Change symbol
       pch=1
  )
  ## Insert diagonal
  abline(0, 1, col="black", lwd=1, lty = 2)
  
  ## MA plot : log-ratio as a function of the log counts
  ## Note: usually, the MA plot represent the log-ratio versus log(mean) of the two measurements.
  ## However, in our case the query sequences are a subset of the background sequences. 
  ## We will thus plot the log(background freq) on the X axis
  plot(x = log2(kmer.stat.table$obs_freq_background),
       y = kmer.stat.table$log2.ratio,
       
       ## Set labels
       xlab = paste("log2(background)"),
       ylab = paste("log2(query/background)"),
       main = paste(k,"-mer frequencies", sep = ""),
       
       col = dot.colors,
       
       ## Insert grid
       panel.first=grid(),
       
       ## Change symbol
       pch=1
  )
  ## Insert diagonal
  abline(h=0, col="black", lwd=1, lty = 2)
  
  
  ## Draw a volcano plot
  ##
  ## Abcsissa indicates the effect,i.e. the log2-ratio querym versus background. 
  ## Ordinate indicates the signficance as computed above. 
  ##
  ## Mark in red the over-represented k-mers with sig > given threshold.
  ## and in blue the under-represented k-mers with sig > given threshold.
  
  ## Set axis limit 
  x.lim <- abs(max(kmer.stat.table$log2.ratio))
  ## Plotting
  with(kmer.stat.table,plot(x=log2.ratio, 
                        y=sig, 
                        pch = 1,
                        xlim = c(-x.lim, x.lim),
                        ylim = c (0, volcano.ymax),
                        
                        ## Set labels
                        xlab = paste("Log2.ratio ", k,"-mer frequencies", sep = ""),
                        ylab = paste("Significance ", k,"-mer frequencies", sep = ""),
                        main = paste("Volcano Plot ", k,"-mer frequencies", sep = ""),
                        # Insert grid
                        panel.first=grid(),
                        
                        #set the colors 
                        col = dot.colors))
  
  ## Denote significance exceeding the trhreshold by 
  above.ymax <- kmer.stat.table$sig > volcano.ymax
  with(kmer.stat.table[above.ymax,],
       points(x=log2.ratio, y=rep(x = volcano.ymax, length.out=sum(above.ymax)), pch=17, col="red",type="p"))
  
  ## insert diagonals
  abline(h=sig.threshold, col= "black", lty = 2)
  abline(v=0, col= "black", lty = 2)
  abline(v=-log2FC.threshold, col= "black", lty = 2)
  abline(v=log2FC.threshold, col= "black", lty = 2)
}


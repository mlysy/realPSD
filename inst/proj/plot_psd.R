#' Plot log-log scale ggplot2 figure of periodogram with ticks on x-y axis
#' @param plotdata Data tibble/frame with two columns, the first colum is the x-axis frequency series, the second column is the y-axis periodogram values
#' @param xlog_scale If TRUE, x-axis will be converted to log scale
#' @param ylog_scale If TRUE, y-axis will be converted to log scale
#' @param show_ticks If TRUE, show ticks on the axis
#' @import scales

plot_psd <- function(plotdata, xlog_scale = TRUE, ylog_scale = TRUE, show_ticks = TRUE, subplot = FALSE, nbreaks) {
  ans <- ggplot(plotdata, aes(x = x, y = y)) +
    geom_line(size = 0.3, color = "deepskyblue1") +
    theme(axis.text = element_text(size=10),
          axis.title = element_text(size=10,face="plain"))
  if(xlog_scale) {
    if(missing(nbreaks)) {
      ans <- scale_x_log10()
    } else {
      ans <- scal
    }
  }
  ## if(xlog_scale) {
  ##   if(subplot) {
  ##     ans <- ans +
  ##       scale_x_continuous(trans = log10_trans(),
  ##         breaks = base_breaks(),
  ##         labels = function(x) round(x/1000))
  ##       # theme(panel.grid.minor = element_blank())
  ##   } else {
  ##     ans <- ans + scale_x_log10(
  ##     breaks = trans_breaks("log10", function(x) 10^x),
  ##     labels = trans_format("log10",
  ##       function(x) fancy_scientific(10^(x-3)) ))
  ##   }
  ## }
  ## if(ylog_scale) {
  ##   ans <- ans + scale_y_log10(
  ##     breaks = trans_breaks("log10", function(x) 10^x),
  ##     labels = trans_format("log10", math_format(10^.x)))
  ##     # labels = fancy_scientific)
  ##   # ans <- ans + scale_y_continuous(trans = "log10", breaks = 1:10)
  ## }
  if(show_ticks) ans <- ans + annotation_logticks(size = 0.3) # display log ticks on bottom and left
  ans
}

# helper functions may be useful
base_breaks <- function(n = 10){
    function(x) {
        ticks <- axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
}

fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "%*%10^", l)
     # return this as an expression
     parse(text=l)
}
# we can also try base::prettyNum

log_breaks = function(maj, radix=10) {
  function(x) {
    minx         = floor(min(logb(x,radix), na.rm=T)) - 1
    maxx         = ceiling(max(logb(x,radix), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    if (maj) {
      breaks = major_breaks
    } else {
      steps = logb(1:(radix-1),radix)
      breaks = rep(steps, times=n_major) +
               rep(major_breaks, each=radix-1)
    }
    radix^breaks
  }
}

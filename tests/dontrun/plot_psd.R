# plot PSD, 1st and 2nd eigenmodes from original data
require(tidyverse)
require(realPSD)
require(scales) # for log scale plot
require(ggpubr) # for ggarrange
data <- read_csv("data/TR400PB.csv")
ytime <- data$ytime # get time-domain measurements
SF <- 5e6 # sampling frequency 5MHz
T_s <- 5 # total time (second)
system.time(
psd_data <- periodogram(yTime = ytime, SF_s = SF, T_s = T_s) # get data (x and y) in the freq domain
)
xPSD <- psd_data$xFreq
yPSD <- psd_data$yFreq
psd_data <- psd_data %>% as_tibble()
# ggplot(data = psd_data, aes(x = xFreq, y = yFreq)) + geom_line()
# binning, compress the frequency
bin_size <- 100
xPSD_b <- binning(xPSD, bin_size, "mean")
yPSD_b <- binning(yPSD, bin_size, "mean")
# convert the frequency (x axis) to kHz
xPSD_b <- xPSD_b / 1000
# plot the full range PSD
psd_data_b <- as_tibble(cbind(xPSD_b, yPSD_b))
plot_psd <- function(plotdata) {
  ggplot(plotdata, aes(x = xPSD_b, y = yPSD_b)) + 
  geom_line(size = 0.18, color = "#29B6F6") + 
  xlab(expression(Frequency (kHz))) +
  ylab(expression(PSD (m^2/Hz))) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() # display log ticks on bottom and left
}
plot_full <- plot_psd(psd_data_b)
plot_full
# plot of first eigenmode
# find the frequency range
frange1 <- c(100:159)
plot_eigen1 <- plot_psd(psd_data_b %>% slice(frange1))
plot_eigen1
# plot of second eigenmode
frange2 <- c(280, 320)
plot_eigen2 <- plot_psd(psd_data_b %>% slice(frange2))
# arrange them into one plot
ggarrange(plot_full, 
  ggarrange(plot_eigen1, plot_eigen2, ncol = 2, 
    labels = c("(b) 1st Eigenmode", "2nd Eigenmode")),
  nrow = 2,
  labels = "(a) Periodogram"
)
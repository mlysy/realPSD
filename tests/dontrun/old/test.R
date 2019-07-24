# test
install.packages("realPSD_1.0.tar.gz", type = "source", repos = NULL, INSTALL_opts = "--install-tests")
require(tidyverse)
require(parallel)
# read fitted data
data_path_fit <- "show_fit"
# tmp <- readRDS(file = file.path(data_path_fit, paste0("show_fit_", 1, ".rds")))
nfit <- nrow(fit_descr)
system.time(
for(ii in 1:nfit) {
  assign(paste0("fit_", ii), 
        readRDS(file.path(data_path_fit, paste0("show_fit_", ii, ".rds"))))
}
)

fit_list <- list()
system.time(
for(ii in 1:nfit) {
  fit_list[[ii]] <- readRDS(file.path(data_path_fit, paste0("show_fit_", ii, ".rds"))) 
}
)

fit_data <- do.call(rbind, fit_list)
fit_data <- cbind(fit_data, fit_descr[,c("Q_level", "method")])

# then manipulate the data frame using tidyverse toolbox
fit_data <- fit_data %>% as_tibble() %>%
  mutate(Q_level = factor(Q_level,  # convert the column Q_level into a factor
    levels = c(1,10,100,500), labels = c("Q = 1", "Q = 10", "Q = 100", "Q = 500"))) %>% 
  mutate(method = factor(method, ordered = FALSE)) # factor column method 
# get a new dataset with ratios instead of fitted values
ratio_data <- fit_data %>% 
  mutate(f0_hat = f0_hat/f0) %>%
  mutate(k_hat = k_hat/k) %>%
  mutate(Aw_hat = Aw_hat/Aw) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 1", Q_hat/Q_vec[1], Q_hat)) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 10", Q_hat/Q_vec[2], Q_hat)) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 100", Q_hat/Q_vec[3], Q_hat)) %>%
  mutate(Q_hat = ifelse(Q_level == "Q = 500", Q_hat/Q_vec[4], Q_hat))

# # covnert the estimated parameters to ratios
# # Q = 1
# fit_Q1_lp <- fit_Q1_lp %*% diag(c(1/f0, 1/Q1, 1/k))
# fit_Q1_nls <- fit_Q1_nls %*% diag(c(1/f0, 1/Q1, 1/k))
# fit_Q1_mle <- fit_Q1_mle %*% diag(c(1/f0, 1/Q1, 1/k))
# colnames(fit_Q1_lp) <- c("f0", "Q", "k")
# colnames(fit_Q1_nls) <- c("f0", "Q", "k")
# colnames(fit_Q1_mle) <- c("f0", "Q", "k")
# # Q = 10
# fit_Q10_lp <- fit_Q10_lp %*% diag(c(1/f0, 1/Q10, 1/k))
# fit_Q10_nls <- fit_Q10_nls %*% diag(c(1/f0, 1/Q10, 1/k))
# fit_Q10_mle <- fit_Q10_mle %*% diag(c(1/f0, 1/Q10, 1/k))
# colnames(fit_Q10_lp) <- c("f0", "Q", "k")
# colnames(fit_Q10_nls) <- c("f0", "Q", "k")
# colnames(fit_Q10_mle) <- c("f0", "Q", "k")
# # Q = 100
# fit_Q100_lp <- fit_Q100_lp %*% diag(c(1/f0, 1/Q100, 1/k))
# fit_Q100_nls <- fit_Q100_nls %*% diag(c(1/f0, 1/Q100, 1/k))
# fit_Q100_mle <- fit_Q100_mle %*% diag(c(1/f0, 1/Q100, 1/k))
# colnames(fit_Q100_lp) <- c("f0", "Q", "k")
# colnames(fit_Q100_nls) <- c("f0", "Q", "k")
# colnames(fit_Q100_mle) <- c("f0", "Q", "k")
# # Q = 500
# fit_Q500_lp <- fit_Q500_lp %*% diag(c(1/f0, 1/Q500, 1/k))
# fit_Q500_nls <- fit_Q500_nls %*% diag(c(1/f0, 1/Q500, 1/k))
# fit_Q500_mle <- fit_Q500_mle %*% diag(c(1/f0, 1/Q500, 1/k))
# colnames(fit_Q500_lp) <- c("f0", "Q", "k")
# colnames(fit_Q500_nls) <- c("f0", "Q", "k")
# colnames(fit_Q500_mle) <- c("f0", "Q", "k")
# # merge datasets
# # Q = 1
# fit_Q1_lp <- as_tibble(fit_Q1_lp) %>% add_column(method = "LP")
# fit_Q1_nls <- as_tibble(fit_Q1_nls) %>% add_column(method = "NLS")
# fit_Q1_mle <- as_tibble(fit_Q1_mle) %>% add_column(method = "MLE")
# fit_Q1 <- bind_rows(fit_Q1_lp, fit_Q1_nls, fit_Q1_mle)
# fit_Q1 <- fit_Q1 %>% add_column(level = "Q = 1")
# # Q = 10
# fit_Q10_lp <- as_tibble(fit_Q10_lp) %>% add_column(method = "LP")
# fit_Q10_nls <- as_tibble(fit_Q10_nls) %>% add_column(method = "NLS")
# fit_Q10_mle <- as_tibble(fit_Q10_mle) %>% add_column(method = "MLE")
# fit_Q10 <- bind_rows(fit_Q10_lp, fit_Q10_nls, fit_Q10_mle) %>% 
#   add_column(level = "Q = 10")
# # Q = 100
# fit_Q100_lp <- as_tibble(fit_Q100_lp) %>% add_column(method = "LP")
# fit_Q100_nls <- as_tibble(fit_Q100_nls) %>% add_column(method = "NLS")
# fit_Q100_mle <- as_tibble(fit_Q100_mle) %>% add_column(method = "MLE")
# fit_Q100 <- bind_rows(fit_Q100_lp, fit_Q100_nls, fit_Q100_mle) %>% 
#   add_column(level = "Q = 100")
# # Q = 500
# fit_Q500_lp <- as_tibble(fit_Q500_lp) %>% add_column(method = "LP")
# fit_Q500_nls <- as_tibble(fit_Q500_nls) %>% add_column(method = "NLS")
# fit_Q500_mle <- as_tibble(fit_Q500_mle) %>% add_column(method = "MLE")
# fit_Q500 <- bind_rows(fit_Q500_lp, fit_Q500_nls, fit_Q500_mle) %>% 
#   add_column(level = "Q = 500")
# # combine all the datasets together
# fit <- bind_rows(fit_Q1, fit_Q10, fit_Q100, fit_Q500)
# boxplot
# Q_hat / Q
ggplot(ratio_data, aes(x = Q_level, y = Q_hat, fill = method)) + geom_boxplot()
# k_hat / k
ggplot(ratio_data, aes(x = Q_level, y = k_hat, fill = method)) + geom_boxplot()
# f0_hat / f0
ggplot(ratio_data, aes(x = Q_level, y = f0_hat, fill = method)) + geom_boxplot()

library(grid)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Dummy data for the plot
	y <- exp(seq(1,10,.1))
	x <- 1:length(y)
	data <- data.frame(x = x, y = y)
tikz(file = "plot_test.tex", width = 5, height = 5)
	#Simple plot of the dummy data using LaTeX elements
	plot <- ggplot(data, aes(x = x, y = y)) + 
		geom_line() +
		#Space does not appear after Latex
		ggtitle( paste("Fancy \\LaTeX ", "\\hspace{0.01cm} title")) +
		labs( x = "$x$ = Time", y = "$\\Phi$ = Innovation output") +
		theme_bw()
	#This line is only necessary if you want to preview the plot right after compiling
	#Necessary to close or the tikxDevice .tex file will not be written
	dev.off()

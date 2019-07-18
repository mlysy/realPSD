# test
tmp <- readRDS("fit_Q1_lp_20.rds")

# covnert the estimated parameters to ratios
# Q = 1
fit_Q1_lp <- fit_Q1_lp %*% diag(c(1/f0, 1/Q1, 1/k))
fit_Q1_nls <- fit_Q1_nls %*% diag(c(1/f0, 1/Q1, 1/k))
fit_Q1_mle <- fit_Q1_mle %*% diag(c(1/f0, 1/Q1, 1/k))
colnames(fit_Q1_lp) <- c("f0", "Q", "k")
colnames(fit_Q1_nls) <- c("f0", "Q", "k")
colnames(fit_Q1_mle) <- c("f0", "Q", "k")
# Q = 10
fit_Q10_lp <- fit_Q10_lp %*% diag(c(1/f0, 1/Q10, 1/k))
fit_Q10_nls <- fit_Q10_nls %*% diag(c(1/f0, 1/Q10, 1/k))
fit_Q10_mle <- fit_Q10_mle %*% diag(c(1/f0, 1/Q10, 1/k))
colnames(fit_Q10_lp) <- c("f0", "Q", "k")
colnames(fit_Q10_nls) <- c("f0", "Q", "k")
colnames(fit_Q10_mle) <- c("f0", "Q", "k")
# Q = 100
fit_Q100_lp <- fit_Q100_lp %*% diag(c(1/f0, 1/Q100, 1/k))
fit_Q100_nls <- fit_Q100_nls %*% diag(c(1/f0, 1/Q100, 1/k))
fit_Q100_mle <- fit_Q100_mle %*% diag(c(1/f0, 1/Q100, 1/k))
colnames(fit_Q100_lp) <- c("f0", "Q", "k")
colnames(fit_Q100_nls) <- c("f0", "Q", "k")
colnames(fit_Q100_mle) <- c("f0", "Q", "k")
# Q = 500
fit_Q500_lp <- fit_Q500_lp %*% diag(c(1/f0, 1/Q500, 1/k))
fit_Q500_nls <- fit_Q500_nls %*% diag(c(1/f0, 1/Q500, 1/k))
fit_Q500_mle <- fit_Q500_mle %*% diag(c(1/f0, 1/Q500, 1/k))
colnames(fit_Q500_lp) <- c("f0", "Q", "k")
colnames(fit_Q500_nls) <- c("f0", "Q", "k")
colnames(fit_Q500_mle) <- c("f0", "Q", "k")
# merge datasets
# Q = 1
fit_Q1_lp <- as_tibble(fit_Q1_lp) %>% add_column(method = "LP")
fit_Q1_nls <- as_tibble(fit_Q1_nls) %>% add_column(method = "NLS")
fit_Q1_mle <- as_tibble(fit_Q1_mle) %>% add_column(method = "MLE")
fit_Q1 <- bind_rows(fit_Q1_lp, fit_Q1_nls, fit_Q1_mle)
fit_Q1 <- fit_Q1 %>% add_column(level = "Q = 1")
# Q = 10
fit_Q10_lp <- as_tibble(fit_Q10_lp) %>% add_column(method = "LP")
fit_Q10_nls <- as_tibble(fit_Q10_nls) %>% add_column(method = "NLS")
fit_Q10_mle <- as_tibble(fit_Q10_mle) %>% add_column(method = "MLE")
fit_Q10 <- bind_rows(fit_Q10_lp, fit_Q10_nls, fit_Q10_mle) %>% 
  add_column(level = "Q = 10")
# Q = 100
fit_Q100_lp <- as_tibble(fit_Q100_lp) %>% add_column(method = "LP")
fit_Q100_nls <- as_tibble(fit_Q100_nls) %>% add_column(method = "NLS")
fit_Q100_mle <- as_tibble(fit_Q100_mle) %>% add_column(method = "MLE")
fit_Q100 <- bind_rows(fit_Q100_lp, fit_Q100_nls, fit_Q100_mle) %>% 
  add_column(level = "Q = 100")
# Q = 500
fit_Q500_lp <- as_tibble(fit_Q500_lp) %>% add_column(method = "LP")
fit_Q500_nls <- as_tibble(fit_Q500_nls) %>% add_column(method = "NLS")
fit_Q500_mle <- as_tibble(fit_Q500_mle) %>% add_column(method = "MLE")
fit_Q500 <- bind_rows(fit_Q500_lp, fit_Q500_nls, fit_Q500_mle) %>% 
  add_column(level = "Q = 500")
# combine all the datasets together
fit <- bind_rows(fit_Q1, fit_Q10, fit_Q100, fit_Q500)
# boxplot
# Q_hat / Q
ggplot(fit, aes(x = level, y = Q, fill = method)) + geom_boxplot()
# k_hat / k
ggplot(fit, aes(x = level, y = k, fill = method)) + geom_boxplot()
# f0_hat / f0
ggplot(fit, aes(x = level, y = f0, fill = method)) + geom_boxplot()


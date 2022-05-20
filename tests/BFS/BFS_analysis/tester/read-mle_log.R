
fdir <- "C:/Users/livse301/Documents/GitHub/sigex/tests/BFS/BFS_analysis/tester"
x <- readLines(file.path(fdir, "mle_log.txt"))

n <- length(x) # includes start/end time on first/last lines

x1 <- x[seq(2, (n-2), 2)]
x2 <- stringr::str_sub(x1, start= 7)
x3 <- strsplit(x = x2, split = ", ")
x4 <- unlist(x3)
x5 <- as.numeric(x4)

psiMat <- matrix(x5, ncol = 6, byrow = TRUE)

df <- data.frame(psiMat)
colnames(df) <- c("", "", "", "", "", "")


# ---- Get lik values from file ----

x1 <- x[seq(3, (n-1), 2)]
x2 <- stringr::str_sub(x1, start= 7)
x4 <- unlist(x2)
x5 <- as.numeric(x4)

df$lik = x5


summary(df$lik)




plot(df$lik, type = 'l')


library(tidyverse)

d

df %>%
  mutate(idx = 1:nrow(df)) %>%
  filter(lik < 100) %>%
  pivot_longer(cols = -idx) %>%
  ggplot(aes(y = value, x = idx )) +
  geom_line() +
  facet_wrap(~name, scales = "free_y")


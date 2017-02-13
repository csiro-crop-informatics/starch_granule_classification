library(dplyr)
library(ggplot2)
library(tidyr)


range <- seq(from = -0.5, to = 5, by = 0.01)
mub <- log(5.2)
sigb <- sqrt(0.5)
distb <- (1/(sigb * sqrt(2 * pi))) * exp(1)^-(((range - mub)^2)/(2 * sigb^2))

mua <- log(20.9)
siga <- sqrt(0.2)
dista <- (1/(siga * sqrt(2 * pi))) * exp(1)^-(((range - mua)^2)/(2 * siga^2))

distdf <- data.frame(range = range, dista, distb = 0.3 * distb)
distdf$sum <- with(distdf, dista + distb)

dist.g <- distdf %>% 
  gather(distribution, value, -range)

ribbondat <- dist.g %>% 
  filter((range > log(10) & distribution == "distb") |
         (range < log(10) & distribution == "dista"))

distplot <- ggplot(dist.g, aes(range, value)) + 
  geom_ribbon(data = ribbondat, aes(ymax = value, fill = distribution), ymin = 0, alpha = 0.5) +
  geom_line(aes(group = distribution, colour = distribution, size = distribution)) +
  geom_vline(xintercept = log(10), linetype = "dotted") + 
  annotate("text", x = 2.1, y = 0.6, label = "10~mu*m~classification~threshold", parse = T, angle = 90, size = 2.5) +
  scale_fill_manual(values = c("#ABABAB", "#333333"), guide = guide_legend(order = 3), name = "Misclassifications", labels = c("A granules classified as B", "B granules classified as A")) +
  scale_colour_manual(values = c("red", "blue", "black"), guide = guide_legend(order = 1), name = "Distributions", labels = c("A granules", "B granules", "Total distribution")) +
  scale_size_manual(values = c(0.5, 0.5, 1), guide = guide_legend(order = 1), name = "Distributions", labels = c("A granules", "B granules", "Total distribution")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) + scale_x_continuous(expand = c(0,0)) +
  labs(x = expression(paste("Log of particle size diameter (", mu, "m)")), y = "Probability Density") +
  theme(panel.background = element_rect(fill = NA, colour = "black"), 
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black"))
distplot
ggsave("figures/distributionExample.png", distplot, type = "cairo-png", width = 7, height = 4)

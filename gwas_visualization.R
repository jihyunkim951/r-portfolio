library(ggplot2)
library(dplyr)

set.seed(42)
n_snps <- 5000

gwas_data <- data.frame(
  SNP = paste0("rs", sample(1000000:9999999, n_snps)),
  CHR = sample(1:22, n_snps, replace = TRUE),
  BP  = sample(1:250000000, n_snps, replace = TRUE),
  P   = runif(n_snps, min = 1e-10, max = 1)
)

gwas_data$P[sample(1:n_snps, 10)] <- runif(10, 1e-10, 1e-7)

gwas_data <- gwas_data %>%
  arrange(CHR, BP) %>%
  mutate(
    logP    = -log10(P),
    CHR_col = ifelse(CHR %% 2 == 0, "even", "odd")
  )

manhattan_plot <- ggplot(gwas_data, aes(x = seq_along(SNP), y = logP, color = CHR_col)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
  scale_color_manual(values = c("even" = "#4A90D9", "odd" = "#2C3E6B")) +
  labs(title = "Manhattan Plot (Simulated GWAS Data)",
       subtitle = "Red line: genome-wide significance (p < 5×10⁻⁸)",
       x = "Genomic Position", y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

n   <- nrow(gwas_data)
exp <- -log10(seq(1/n, 1, length.out = n))
obs <- sort(gwas_data$logP, decreasing = TRUE)

qq_plot <- ggplot(data.frame(expected = exp, observed = obs),
                  aes(x = expected, y = observed)) +
  geom_point(size = 0.8, color = "#4A90D9", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "QQ Plot (Simulated GWAS Data)",
       x = "Expected -log10(p)", y = "Observed -log10(p)") +
  theme_minimal()

ggsave("manhattan_plot.png", manhattan_plot, width = 12, height = 4, dpi = 150)
ggsave("qq_plot.png",        qq_plot,        width = 5,  height = 5, dpi = 150)

cat("완료! manhattan_plot.png, qq_plot.png 저장됨\n")
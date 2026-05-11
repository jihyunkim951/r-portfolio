# =============================================
# GWAS Result Visualization
# Data: simulated GWAS summary statistics
# Author: 본인 이름
# Date: 2025-05
# =============================================

# 1. 패키지 설치 & 로드 ---------------------
# 처음 한 번만 실행
install.packages(c("ggplot2", "dplyr"))

library(ggplot2)
library(dplyr)

# 2. 시뮬레이션 데이터 만들기 ---------------
set.seed(42)
n_snps <- 5000

gwas_data <- data.frame(
  SNP     = paste0("rs", 1:n_snps),
  CHR     = sample(1:22, n_snps, replace = TRUE),
  BP      = sample(1:250000000, n_snps, replace = TRUE),
  P       = c(runif(4950, 0.0001, 1),       # 일반 SNP
              runif(50,   1e-10, 1e-7))      # 유의미한 SNP (신호)
)

# 3. Manhattan Plot -------------------------
# 염색체별 x축 위치 계산
gwas_plot <- gwas_data %>%
  arrange(CHR, BP) %>%
  mutate(
    BP_cum = BP + ave(BP, CHR, FUN = function(x) max(x)) %>% 
      cumsum() %>% lag(default = 0),
    log_p  = -log10(P)
  )

# 염색체 중간 위치 (x축 라벨용)
chr_labels <- gwas_plot %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum))

# 색상 (홀짝 염색체 번갈아)
chr_colors <- rep(c("#4E79A7", "#A0CBE8"), 11)

ggplot(gwas_plot, aes(x = BP_cum, y = log_p, color = factor(CHR))) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(5e-8), color = "red",    linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -log10(1e-5), color = "orange", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = chr_colors) +
  scale_x_continuous(label = chr_labels$CHR, breaks = chr_labels$center) +
  labs(
    title    = "Manhattan Plot — Simulated GWAS",
    subtitle = "Red: genome-wide significance (p < 5×10⁻⁸) | Orange: suggestive (p < 1×10⁻⁵)",
    x        = "Chromosome",
    y        = "-log₁₀(p-value)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 7)
  )

ggsave("manhattan_plot.png", width = 12, height = 5, dpi = 300)

# 4. QQ Plot --------------------------------
observed  <- sort(-log10(gwas_data$P))
expected  <- sort(-log10(ppoints(nrow(gwas_data))))

qq_df <- data.frame(expected, observed)

ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(size = 0.8, color = "#4E79A7", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.8) +
  labs(
    title    = "QQ Plot — Simulated GWAS",
    subtitle = "Points above the line indicate true association signals",
    x        = "Expected -log₁₀(p)",
    y        = "Observed -log₁₀(p)"
  ) +
  theme_classic()

ggsave("qq_plot.png", width = 6, height = 6, dpi = 300)
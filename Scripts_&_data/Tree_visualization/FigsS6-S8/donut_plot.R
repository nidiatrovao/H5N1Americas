# Donut Plot: Flyway Root Probabilities
# Install ggplot2 if needed: install.packages("ggplot2")

library(ggplot2)

# Load data
df <- read.csv("your_root_prob.csv")

# Calculate label positions (midpoint of each slice)
df <- df[order(-df$probability), ]  # sort descending
df$fraction <- df$probability / sum(df$probability)
df$ymax <- cumsum(df$fraction)
df$ymin <- c(0, head(df$ymax, -1))
df$label_pos <- (df$ymax + df$ymin) / 2
df$label <- paste0(df$flyway, "\n", scales::percent(df$probability, accuracy = 0.1))

# Color palette
flyway_colors <- c(
  "Pacific"     = "#f46464",
  "Central"     = "#e8a989",
  "Mississippi" = "#87c9a9",
  "Atlantic"    = "#005895"
)

# Build donut plot
p <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.2, fill = flyway)) +
  geom_rect(color = "white", linewidth = 0) +
  geom_text(
    aes(x = 4.6, y = label_pos, label = label),
    size = 3.8,
    lineheight = 1.2
  ) +
  coord_polar(theta = "y") +
  xlim(c(0.5, 5)) +
  scale_fill_manual(values = flyway_colors) +
  labs(
    title = "Root State Probability",
    fill  = "Flyway"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 10)),
    legend.position = "none"
  )

# Save output
ggsave("D11_donut_plot.png", plot = p, width = 7, height = 6, dpi = 150)
message("Saved donut_plot.png")

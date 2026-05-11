# =============================================================================
# Dot-Range Plot — Between-Group Pairwise Similarity
# =============================================================================

# install.packages(c("ggplot2", "dplyr", "readr"))

library(ggplot2)
library(dplyr)
library(readr)

# =============================================================================
# CONFIG
# =============================================================================

INPUT_FILE  <- "NS1_grp_pairwise_similarity_results.csv"
OUTPUT_FILE <- "between_group_similarity.pdf"

Y_MIN <- 0.70   # adjust to zoom y-axis
Y_MAX <- 1.00

# =============================================================================

df <- read_csv(INPUT_FILE, show_col_types = FALSE) %>%
  filter(comparison_type == "between") %>%
  mutate(label = paste0(group1, "\nvs\n", group2)) #%>%
  mutate(label = case_when(
    group1 == "2009 H1N1 pandemic" & group2 == "H5 2.3.4.4b"      ~ "2009 H1N1 pandemic\nvs\nH5 2.3.4.4b",
    group1 == "2009 H1N1 pandemic" & group2 == "D1.1"  ~ "D1.1\nvs\n2009 H1N1 pandemic",
    group1 == "D1.1" & group2 == "H5 2.3.4.4b"  ~ "D1.1\nvs\nH5 2.3.4.4b",
    TRUE ~ label
  )) %>%
  mutate(label = factor(label, levels = c("2009 H1N1 pandemic\nvs\nH5 2.3.4.4b", "D1.1\nvs\nH5 2.3.4.4b","D1.1\nvs\n2009 H1N1 pandemic")))

ggplot(df, aes(x = label, y = mean_similarity)) +
  # min-max range (thin)
  geom_linerange(aes(ymin = min_similarity, ymax = max_similarity),
                 linewidth = 1, color = "#185FA5", alpha = 0.4) +
  # SD range (thick)
  geom_linerange(aes(ymin = mean_similarity - sd_similarity,
                     ymax = mean_similarity + sd_similarity),
                 linewidth = 4, color = "#185FA5", alpha = 0.7) +
  # mean dot
  geom_point(size = 4, shape = 21, fill = "white",
             color = "#185FA5", stroke = 1.5) +
  scale_y_continuous(
    limits = c(Y_MIN, Y_MAX),
    oob    = scales::squish,
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  coord_flip() +
  labs(
    title    = "Between-Group Pairwise Similarity",
    #subtitle = "Dot = mean  |  Thick bar = +/-SD  |  Thin bar = min-max",
    x = NULL,
    y = "Similarity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey50", size = 10),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor   = element_blank()
  )

ggsave(OUTPUT_FILE, width = 7, height = 4)
cat("Plot saved to:", OUTPUT_FILE, "\n")

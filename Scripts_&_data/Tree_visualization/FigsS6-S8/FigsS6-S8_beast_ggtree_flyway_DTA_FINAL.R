library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)
library(coda)  # for HPDinterval

# input files
TREE_FILE <- "your_mcc.tre"
LOG_FILE  <- "your_BEAST.log"

BURNIN_FRAC <- 0.10
MRSD <- "2025-12-14"  # MRSD = Most recent sample date
MRSD_DATE <- as.Date(MRSD)

# check names(log_df) if this is wrong - could be different between BEAST versions
ROOT_COL <- "rootHeight"

# flyway colors
TRAIT_COLORS <- c(
  "Pacific"     = "#f46464",
  "Central"     = "#e8a989",
  "Mississippi" = "#87c9a9",
  "Atlantic"    = "#005895"
)

# how tall to draw the root age density (relative to tree height, 0-1)
DENSITY_HEIGHT <- 0.3


# read tree and check what annotation fields are available
beast_tree <- read.beast(TREE_FILE)
get.fields(beast_tree)
# look for: flyway, height_0.95_HPD, posterior, etc.

HPD_FIELD <- "height_0.95_HPD"


# read log, strip burnin
log_raw  <- read.table(LOG_FILE, header = TRUE, sep = "\t", comment.char = "#")
n_burnin <- floor(nrow(log_raw) * BURNIN_FRAC)
log_df   <- log_raw[(n_burnin + 1):nrow(log_raw), ]

# quick sanity check
names(log_df)[1:20]

if (!ROOT_COL %in% names(log_df)) {
  stop(paste("Can't find", ROOT_COL, "in log — available cols:",
             paste(names(log_df), collapse = ", ")))
}

# convert root height (years before MRSD) to calendar dates
root_height_samples <- log_df[[ROOT_COL]]
root_date_samples   <- MRSD_DATE - round(root_height_samples * 365.25)

hpd <- HPDinterval(as.mcmc(as.numeric(root_date_samples)), prob = 0.95)
hpd_lo_date <- as.Date(hpd[1, "lower"], origin = "1970-01-01")
hpd_hi_date <- as.Date(hpd[1, "upper"], origin = "1970-01-01")
med_date    <- as.Date(median(as.numeric(root_date_samples)), origin = "1970-01-01")

cat(sprintf("Median TMRCA: %s\n95%% HPD: %s to %s\n",
            format(med_date, "%b %Y"),
            format(hpd_lo_date, "%b %Y"),
            format(hpd_hi_date, "%b %Y")))


# build density for the root date distribution
dens <- density(as.numeric(root_date_samples))
dens_df <- data.frame(
  x = as.Date(dens$x, origin = "1970-01-01"),
  y = dens$y
)
hpd_region <- subset(dens_df, x >= hpd_lo_date & x <= hpd_hi_date)

# normalize so we can scale to tree coords later
dens_df$y_norm    <- dens_df$y / max(dens_df$y)
hpd_region$y_norm <- hpd_region$y / max(dens_df$y)


# build the base tree first so we can grab the y range
# (need this before we can position the density ribbon correctly)
p_tree_base <- ggtree(
  beast_tree,
  mrsd    = MRSD,
  as.Date = TRUE,
  aes(color = flyway),
  right = TRUE
)

tree_built <- ggplot_build(p_tree_base)
y_range <- tree_built$layout$panel_params[[1]]$y.range
y_min   <- y_range[1]
y_span  <- diff(y_range)

dens_df$y_tree    <- y_min + dens_df$y_norm    * y_span * DENSITY_HEIGHT
hpd_region$y_tree <- y_min + hpd_region$y_norm * y_span * DENSITY_HEIGHT


# assemble final plot — density goes in first so it sits behind the branches
final_plot <- p_tree_base +

  # full density curve
  geom_ribbon(
    data = dens_df,
    aes(x = x, ymin = y_min, ymax = y_tree),
    fill      = "#EEEDFE",
    color     = "#7F77DD",
    linewidth = 0.6,
    alpha     = 0.7,
    inherit.aes = FALSE
  ) +

  # shaded HPD region
  geom_ribbon(
    data = hpd_region,
    aes(x = x, ymin = y_min, ymax = y_tree),
    fill  = "#7F77DD",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +

  # commented these out for now — a bit cluttered with both lines and the ribbon
  #geom_vline(xintercept = c(hpd_lo_date, hpd_hi_date),
  #           linetype = "dashed", color = "#534AB7", linewidth = 0.5) +
  #geom_vline(xintercept = med_date, color = "#3C3489", linewidth = 0.8) +

  coord_cartesian(clip = 'off', ylim = c(y_min - 4, NA)) +

  theme_tree2() +
  scale_x_date(date_labels = "%b\n%Y") +

  scale_color_manual(name = "Flyway", values = TRAIT_COLORS, na.value = "gray60") +

  geom_tippoint(aes(fill = flyway), shape = 21, colour = "black",
                show.legend = FALSE, stroke = 0.25) +

  scale_fill_manual(values = TRAIT_COLORS, na.value = "gray60") +

  # node HPD bars — leaving off for now, gets messy on a big tree
  #geom_range(range = HPD_FIELD, color = "gray50", alpha = 0.4, size = 1.5) +

  # posterior support labels (only strong nodes)
  #geom_nodelab(aes(label = ifelse(posterior >= 0.95, round(posterior, 2), "")),
  #             size = 2.5, hjust = -0.1, color = "gray30") +

  #geom_tiplab(size = 2.5, offset = 0.5) +

  annotate("text",
           x     = med_date,
           y     = y_min + y_span * DENSITY_HEIGHT * 1.3,
           label = paste0("median TMRCA: ", format(med_date, "%Y-%m-%d")),
           hjust = 0.5, vjust = 0, size = 3.5, color = "#3C3489") +

  annotate("text",
           x     = med_date,
           y     = y_min + y_span * DENSITY_HEIGHT * 1.25,
           label = paste0("[95% HPDI: ", format(hpd_lo_date, "%Y-%m-%d"),
                          " to ", format(hpd_hi_date, "%Y-%m-%d"), "]"),
           hjust = 0.5, vjust = 1.4, size = 3.5, color = "#3C3489") +

  labs(title = "D1.1") +

  theme(
    plot.title = element_text(hjust = 0, vjust = 0.9, size = 18),
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.8),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    axis.text.x  = element_text(size = 9)
    #plot.margin = margin(10, 20, 10, 10)
  ) +

  guides(color = guide_legend(override.aes = list(linewidth = 2)))


# save outputs
ggsave("beast_discrete_trait_tree.pdf", plot = final_plot,
       width = 12, height = 10, units = "in", device = cairo_pdf)

ggsave("beast_discrete_trait_tree.png", plot = final_plot,
       width = 12, height = 10, units = "in", dpi = 300)

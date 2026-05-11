#### Pairwise Similarity Between Groups ####

library(Biostrings)
library(dplyr)
library(readr)

# =============================================================================
# CONFIG — edit these before running
# =============================================================================

FASTA_FILE  <- "alignment.fasta"   # path to your aligned FASTA file
CSV_FILE    <- "N1-pairwise-sim-groups.csv"         # path to your groups CSV file
ID_COL      <- "seq_id"                  # column in CSV with sequence IDs
GROUP_COL   <- "group"                   # column in CSV with group labels
OUTPUT_FILE <- "N1_pairwise_similarity_results_D11_v2.csv"  # where to save results

# 
# CSV FORMAT REMINDER
# Your CSV needs at least two columns, e.g.:
#
#   seq_id, group
#   EPI_ISL_11628213, avian
#   EPI_ISL_14098916, mammal
#
# The seq_id should match the FASTA header up to the first "|" or space.
# 



# 1. Load FASTA

cat("Reading FASTA:", FASTA_FILE, "\n")
seqs <- readDNAStringSet(FASTA_FILE)

# Use only the part of the header before the first "|" or space as the ID
seq_ids    <- sub("[| ].*", "", names(seqs))
names(seqs) <- seq_ids
cat("  Loaded", length(seqs), "sequences, alignment width =", width(seqs)[1], "bp\n")



# 2. Load group CSV

cat("Reading groups CSV:", CSV_FILE, "\n")
groups_df <- read_csv(CSV_FILE, show_col_types = FALSE)

# Validate columns
if (!ID_COL %in% colnames(groups_df)) {
  stop("Column '", ID_COL, "' not found in CSV. Available: ",
       paste(colnames(groups_df), collapse = ", "))
}
if (!GROUP_COL %in% colnames(groups_df)) {
  stop("Column '", GROUP_COL, "' not found in CSV. Available: ",
       paste(colnames(groups_df), collapse = ", "))
}

groups_df <- groups_df %>%
  select(seq_id = all_of(ID_COL), group = all_of(GROUP_COL)) %>%
  filter(!is.na(seq_id), !is.na(group))

# Match sequences present in both files
common_ids <- intersect(seq_ids, groups_df$seq_id)
if (length(common_ids) == 0) {
  stop("No IDs matched between FASTA headers and CSV column '", ID_COL, "'.\n",
       "  Example FASTA IDs: ", paste(head(seq_ids, 3), collapse = ", "), "\n",
       "  Example CSV IDs:   ", paste(head(groups_df$seq_id, 3), collapse = ", "))
}

n_fasta_only <- length(seq_ids) - length(common_ids)
n_csv_only   <- nrow(groups_df) - length(common_ids)
if (n_fasta_only > 0) warning(n_fasta_only, " FASTA sequence(s) not in CSV - excluded.")
if (n_csv_only   > 0) warning(n_csv_only,   " CSV row(s) not in FASTA - excluded.")

seqs      <- seqs[common_ids]
groups_df <- groups_df %>% filter(seq_id %in% common_ids)

group_levels <- sort(unique(groups_df$group))
cat("  Matched", length(common_ids), "sequences across",
    length(group_levels), "groups:", paste(group_levels, collapse = ", "), "\n\n")



# 3. Convert alignment to character matrix

cat("Converting alignment to matrix ...\n")
seq_matrix <- as.matrix(seqs)

# Remove columns that are entirely gaps across all sequences
gap_only_cols <- apply(seq_matrix, 2, function(col) all(col == "-"))
seq_matrix    <- seq_matrix[, !gap_only_cols, drop = FALSE]
cat("  Retained", ncol(seq_matrix), "non-gap-only alignment positions\n\n")



# 4. Pairwise similarity function

# Proportion of identical sites, ignoring positions where either sequence
# has a gap (-) or ambiguous base (N, ?).
calc_similarity <- function(s1, s2) {
  valid <- s1 != "-" & s2 != "-" & s1 != "N" & s2 != "N" & s1 != "?" & s2 != "?"
  n_valid <- sum(valid)
  if (n_valid == 0) return(NA_real_)
  sum(s1[valid] == s2[valid]) / n_valid
}


# 5. Compute within- and between-group average pairwise similarities

n_groups <- length(group_levels)
cat("Computing pairwise similarities ...\n")

# Index sequences by group
group_indices <- lapply(group_levels, function(g) {
  ids <- groups_df$seq_id[groups_df$group == g]
  which(rownames(seq_matrix) %in% ids)
})
names(group_indices) <- group_levels

results <- vector("list", n_groups * (n_groups + 1) / 2)
k <- 1

for (i in seq_len(n_groups)) {
  for (j in i:n_groups) {
    g1   <- group_levels[i]
    g2   <- group_levels[j]
    idx1 <- group_indices[[g1]]
    idx2 <- group_indices[[g2]]

    if (i == j) {
      # Within-group: all unique pairs
      if (length(idx1) < 2) {
        sim_vals <- NA_real_
        n_pairs  <- NA_integer_
      } else {
        pairs    <- combn(idx1, 2)
        n_pairs  <- ncol(pairs)
        sim_vals <- vapply(seq_len(n_pairs), function(p)
          calc_similarity(seq_matrix[pairs[1, p], ], seq_matrix[pairs[2, p], ]),
          numeric(1))
      }
    } else {
      # Between-group: all cross-pairs
      n_pairs  <- length(idx1) * length(idx2)
      sim_vals <- as.vector(outer(idx1, idx2, Vectorize(function(a, b)
        calc_similarity(seq_matrix[a, ], seq_matrix[b, ]))))
    }

    results[[k]] <- tibble(
      group1          = g1,
      group2          = g2,
      comparison_type = ifelse(i == j, "within", "between"),
      n_seq_g1        = length(idx1),
      n_seq_g2        = if (i == j) length(idx1) else length(idx2),
      n_pairs         = n_pairs,
      mean_similarity = round(mean(sim_vals, na.rm = TRUE), 6),
      sd_similarity   = round(sd(sim_vals,   na.rm = TRUE), 6),
      min_similarity  = round(min(sim_vals,  na.rm = TRUE), 6),
      max_similarity  = round(max(sim_vals,  na.rm = TRUE), 6)
    )

    cat(sprintf("  %s vs %s - mean similarity: %.4f\n",
                g1, g2, mean(sim_vals, na.rm = TRUE)))
    k <- k + 1
  }
}

results_df <- bind_rows(results)


write_csv(results_df, OUTPUT_FILE)
print(as.data.frame(results_df), row.names = FALSE)

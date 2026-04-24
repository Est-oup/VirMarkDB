library(tidyverse)

dir.create(OUT_LOGS, recursive = TRUE, showWarnings = FALSE)

genome_map_tsv <- "output/config/genome_marker_map.tsv"

# Plot theme
plot_theme <- theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

# Add metadata to ranked hits
tblout_ranked_log <- tblout_ranked %>%
  left_join(
    marker_map %>%
      select(marker_group_id, group_id, marker) %>%
      distinct(),
    by = "marker_group_id"
  ) %>%
  left_join(
    manifest %>%
      select(virus_id, Virus_names, ICTV_ID, any_of(taxonomy_all)),
    by = "virus_id"
  ) %>%
  group_by(virus_id, marker_group_id) %>%
  mutate(
    n_orfs_before = n(),
    keep = copy_rank == 1 |
      (
        score_ratio >= default_threshold_score &
        length_ratio >= default_threshold_length
      )
  ) %>%
  ungroup()

# Keep selected ORFs as is
selected_orfs_log <- selected_orfs

# Log by genome and marker
log_genome_marker <- tblout_ranked_log %>%
  group_by(virus_id, Virus_names, marker_group_id, group_id, marker) %>%
  summarise(
    n_orfs_before = first(n_orfs_before),
    n_orfs_after = sum(keep),
    best_orf_name = first(orf_name),
    best_score = first(best_score),
    best_orf_length = first(best_orf_length),
    orf_names_before = str_c(orf_name, collapse = ";"),
    orf_names_after = str_c(orf_name[keep], collapse = ";"),
    multicopy_before = first(n_orfs_before) > 1,
    multicopy_after = sum(keep) > 1,
    .groups = "drop"
  )

write_tsv(
  log_genome_marker,
  file.path(OUT_LOGS, "log_genome_marker_before_after.tsv")
)

# Multicopy genomes only
log_duplicated_genomes <- log_genome_marker %>%
  filter(multicopy_before | multicopy_after)

write_tsv(
  log_duplicated_genomes,
  file.path(OUT_LOGS, "log_duplicated_genomes.tsv")
)

# Summary by marker
log_marker_summary <- log_genome_marker %>%
  group_by(marker_group_id, group_id, marker) %>%
  summarise(
    genomes_with_hit = n(),
    genomes_multicopy_before = sum(multicopy_before, na.rm = TRUE),
    genomes_multicopy_after = sum(multicopy_after, na.rm = TRUE),
    total_orfs_before = sum(n_orfs_before, na.rm = TRUE),
    total_orfs_after = sum(n_orfs_after, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  log_marker_summary,
  file.path(OUT_LOGS, "log_marker_summary.tsv")
)

# Shared ORFs across marker groups
orf_shared_log <- bind_rows(
  tblout_ranked_log %>%
    group_by(virus_id, orf_name) %>%
    summarise(
      stage = "before",
      n_marker_groups = n_distinct(marker_group_id),
      marker_group_ids = str_c(sort(unique(marker_group_id)), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(n_marker_groups > 1),

  selected_orfs_log %>%
    group_by(virus_id, orf_name) %>%
    summarise(
      stage = "after",
      n_marker_groups = n_distinct(marker_group_id),
      marker_group_ids = str_c(sort(unique(marker_group_id)), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(n_marker_groups > 1)
)

write_tsv(
  orf_shared_log,
  file.path(OUT_LOGS, "log_orf_shared_across_markers.tsv")
)

# Raw ORF lengths
orf_length_raw <- bind_rows(
  tblout_ranked_log %>%
    transmute(
      stage = "before",
      virus_id, Virus_names, marker_group_id, group_id, marker,
      orf_name, orf_length
    ),
  selected_orfs_log %>%
    transmute(
      stage = "after",
      virus_id, Virus_names, marker_group_id, group_id, marker,
      orf_name, orf_length
    )
)

write_tsv(
  orf_length_raw,
  file.path(OUT_LOGS, "log_orf_length_raw.tsv")
)

# ORF length summary
orf_length_summary <- orf_length_raw %>%
  group_by(stage, marker_group_id, group_id, marker) %>%
  summarise(
    n_orfs = n(),
    min_length = min(orf_length, na.rm = TRUE),
    median_length = median(orf_length, na.rm = TRUE),
    mean_length = mean(orf_length, na.rm = TRUE),
    max_length = max(orf_length, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  orf_length_summary,
  file.path(OUT_LOGS, "log_orf_length_summary.tsv")
)

# Missing expected markers
genome_marker_map <- read_tsv(genome_map_tsv, show_col_types = FALSE) %>%
  mutate(
    marker_group_id = if_else(
      is.na(marker_group_id) | marker_group_id == "",
      str_c(group_id, "__", marker),
      marker_group_id
    )
  ) %>%
  distinct(virus_id, marker_group_id, group_id, marker)

missing_marker_log <- genome_marker_map %>%
  left_join(
    tblout_ranked_log %>%
      select(virus_id, marker_group_id) %>%
      distinct() %>%
      mutate(detected_before = TRUE),
    by = c("virus_id", "marker_group_id")
  ) %>%
  left_join(
    selected_orfs_log %>%
      select(virus_id, marker_group_id) %>%
      distinct() %>%
      mutate(detected_after = TRUE),
    by = c("virus_id", "marker_group_id")
  ) %>%
  left_join(
    manifest %>%
      select(virus_id, Virus_names, ICTV_ID, any_of(taxonomy_all)),
    by = "virus_id"
  ) %>%
  mutate(
    detected_before = coalesce(detected_before, FALSE),
    detected_after = coalesce(detected_after, FALSE),
    missing_before = !detected_before,
    missing_after = !detected_after
  )

write_tsv(
  missing_marker_log,
  file.path(OUT_LOGS, "log_missing_expected_markers.tsv")
)

missing_marker_summary <- missing_marker_log %>%
  group_by(marker_group_id, group_id, marker) %>%
  summarise(
    n_expected_genomes = n(),
    n_missing_before = sum(missing_before, na.rm = TRUE),
    n_missing_after = sum(missing_after, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  missing_marker_summary,
  file.path(OUT_LOGS, "log_missing_expected_markers_summary.tsv")
)

# Global summary
database_summary <- tibble(
  n_viruses_total = n_distinct(manifest$virus_id),
  n_viruses_with_at_least_1_marker = n_distinct(selected_orfs_log$virus_id),
  n_marker_groups = n_distinct(selected_orfs_log$marker_group_id),
  n_selected_orfs_total = nrow(selected_orfs_log),
  n_viruses_multicopy_total = n_distinct(log_genome_marker$virus_id[log_genome_marker$multicopy_after]),
  n_shared_orfs_across_markers = nrow(
    orf_shared_log %>%
      filter(stage == "after")
  )
)

write_tsv(
  database_summary,
  file.path(OUT_LOGS, "database_summary.tsv")
)

# Summary for README
marker_summary <- selected_orfs_log %>%
  group_by(marker_group_id, group_id, marker, virus_id) %>%
  summarise(
    n_copies = n(),
    .groups = "drop"
  ) %>%
  group_by(marker_group_id, group_id, marker) %>%
  summarise(
    n_viruses_with_hit = n(),
    n_selected_orfs = sum(n_copies),
    n_viruses_multicopy = sum(n_copies > 1),
    mean_copies_per_positive_virus = mean(n_copies),
    max_copies_in_one_virus = max(n_copies),
    .groups = "drop"
  ) %>%
  left_join(
    orf_length_summary %>%
      filter(stage == "after") %>%
      select(marker_group_id, median_length, min_length, max_length),
    by = "marker_group_id"
  )

write_tsv(
  marker_summary,
  file.path(OUT_LOGS, "marker_summary.tsv")
)

# Marker coverage per genome
marker_coverage_distribution <- selected_orfs_log %>%
  distinct(virus_id, marker_group_id) %>%
  count(virus_id, name = "n_markers_detected") %>%
  count(n_markers_detected, name = "n_viruses")

write_tsv(
  marker_coverage_distribution,
  file.path(OUT_LOGS, "marker_coverage_distribution.tsv")
)

# Covered genomes
covered_viruses <- selected_orfs_log %>%
  distinct(virus_id) %>%
  mutate(is_covered = TRUE)

# Global taxonomy coverage
make_taxonomy_coverage <- function(taxo_rank_plot, top_n = 20) {
  taxonomy_coverage_summary <- manifest %>%
    select(virus_id, any_of(taxonomy_all)) %>%
    left_join(covered_viruses, by = "virus_id") %>%
    mutate(
      is_covered = coalesce(is_covered, FALSE),
      taxo_value = .data[[taxo_rank_plot]],
      taxo_value = replace_na(taxo_value, "NA")
    ) %>%
    group_by(taxo_value) %>%
    summarise(
      n_viruses_total = n(),
      n_viruses_covered = sum(is_covered, na.rm = TRUE),
      n_viruses_not_covered = sum(!is_covered, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(taxonomy_group = taxo_value) %>%
    arrange(desc(n_viruses_total), taxonomy_group)

  write_tsv(
    taxonomy_coverage_summary,
    file.path(
      OUT_LOGS,
      str_c("taxonomy_coverage_summary_", str_to_lower(taxo_rank_plot), ".tsv")
    )
  )

  taxonomy_coverage_plot <- taxonomy_coverage_summary %>%
    slice_head(n = top_n)

  p_taxonomy_coverage <- taxonomy_coverage_plot %>%
    pivot_longer(
      cols = c(n_viruses_covered, n_viruses_not_covered),
      names_to = "status",
      values_to = "n_viruses"
    ) %>%
    mutate(
      status = recode(
        status,
        n_viruses_covered = "covered",
        n_viruses_not_covered = "not_covered"
      ),
      taxonomy_group = factor(
        taxonomy_group,
        levels = rev(taxonomy_coverage_plot$taxonomy_group)
      )
    ) %>%
    ggplot(aes(x = taxonomy_group, y = n_viruses, fill = status)) +
    geom_col(position = "stack") +
    coord_flip() +
    labs(
      x = taxo_rank_plot,
      y = "Number of viruses",
      title = str_c("Taxonomy coverage by ", taxo_rank_plot)
    ) +
    plot_theme

  ggsave(
    file.path(
      OUT_LOGS,
      str_c("plot_taxonomy_coverage_", str_to_lower(taxo_rank_plot), ".png")
    ),
    p_taxonomy_coverage,
    width = 10,
    height = 7
  )
}

# Taxonomy coverage by marker
make_taxonomy_coverage_by_marker <- function(taxo_rank_plot, top_n = 12) {
  taxonomy_coverage_summary <- genome_marker_map %>%
    left_join(
      manifest %>%
        select(virus_id, any_of(taxonomy_all)),
      by = "virus_id"
    ) %>%
    left_join(
      selected_orfs_log %>%
        select(virus_id, marker_group_id) %>%
        distinct() %>%
        mutate(is_covered = TRUE),
      by = c("virus_id", "marker_group_id")
    ) %>%
    mutate(
      is_covered = coalesce(is_covered, FALSE),
      taxo_value = .data[[taxo_rank_plot]],
      taxo_value = replace_na(taxo_value, "NA")
    ) %>%
    group_by(marker_group_id, group_id, marker, taxo_value) %>%
    summarise(
      n_viruses_total = n(),
      n_viruses_covered = sum(is_covered, na.rm = TRUE),
      n_viruses_not_covered = sum(!is_covered, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(taxonomy_group = taxo_value) %>%
    arrange(marker, desc(n_viruses_total), taxonomy_group)

  write_tsv(
    taxonomy_coverage_summary,
    file.path(
      OUT_LOGS,
      str_c("taxonomy_coverage_summary_", str_to_lower(taxo_rank_plot), "_by_marker.tsv")
    )
  )

  taxonomy_coverage_plot <- taxonomy_coverage_summary %>%
    group_by(marker) %>%
    slice_max(order_by = n_viruses_total, n = top_n, with_ties = FALSE) %>%
    ungroup()

  p_taxonomy_coverage <- taxonomy_coverage_plot %>%
    pivot_longer(
      cols = c(n_viruses_covered, n_viruses_not_covered),
      names_to = "status",
      values_to = "n_viruses"
    ) %>%
    mutate(
      status = recode(
        status,
        n_viruses_covered = "covered",
        n_viruses_not_covered = "not_covered"
      )
    ) %>%
    ggplot(aes(x = reorder(taxonomy_group, n_viruses_total), y = n_viruses, fill = status)) +
    geom_col(position = "stack") +
    coord_flip() +
    facet_wrap(~marker, scales = "free_y") +
    labs(
      x = taxo_rank_plot,
      y = "Number of viruses",
      title = str_c("Taxonomy coverage by ", taxo_rank_plot, " and marker")
    ) +
    plot_theme

  ggsave(
    file.path(
      OUT_LOGS,
      str_c("plot_taxonomy_coverage_", str_to_lower(taxo_rank_plot), "_by_marker.png")
    ),
    p_taxonomy_coverage,
    width = 14,
    height = 10
  )
}

make_taxonomy_coverage("Family")
make_taxonomy_coverage("Genus")
make_taxonomy_coverage_by_marker("Family")
make_taxonomy_coverage_by_marker("Genus")

# Multicopy plot
plot_multicopy_data <- log_marker_summary %>%
  select(marker_group_id, marker, genomes_multicopy_before, genomes_multicopy_after) %>%
  pivot_longer(
    cols = c(genomes_multicopy_before, genomes_multicopy_after),
    names_to = "stage",
    values_to = "n_genomes"
  ) %>%
  mutate(
    stage = recode(
      stage,
      genomes_multicopy_before = "before",
      genomes_multicopy_after = "after"
    )
  )

p_multicopy_before <- plot_multicopy_data %>%
  filter(stage == "before") %>%
  ggplot(aes(x = reorder(marker_group_id, n_genomes), y = n_genomes)) +
  geom_col(fill = "#F8766D") +
  coord_flip() +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = "Marker group",
    y = "Genomes with multicopy",
    title = "Multicopy genomes before filtering"
  ) +
  plot_theme +
  theme(legend.position = "none")

ggsave(
  file.path(OUT_LOGS, "plot_multicopy_before.png"),
  p_multicopy_before,
  width = 12,
  height = 8
)

p_multicopy_after <- plot_multicopy_data %>%
  filter(stage == "after") %>%
  ggplot(aes(x = reorder(marker_group_id, n_genomes), y = n_genomes)) +
  geom_col(fill = "#00BFC4") +
  coord_flip() +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = "Marker group",
    y = "Genomes with multicopy",
    title = "Multicopy genomes after filtering"
  ) +
  plot_theme +
  theme(legend.position = "none")

ggsave(
  file.path(OUT_LOGS, "plot_multicopy_after.png"),
  p_multicopy_after,
  width = 12,
  height = 8
)

# ORF length plots
p_lengths_before <- orf_length_raw %>%
  filter(stage == "before") %>%
  ggplot(aes(x = marker_group_id, y = orf_length)) +
  geom_boxplot(outlier.size = 0.5, fill = "#F8766D") +
  coord_flip() +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = "Marker group",
    y = "ORF length",
    title = "ORF length distributions before filtering"
  ) +
  plot_theme +
  theme(legend.position = "none")

ggsave(
  file.path(OUT_LOGS, "plot_orf_length_boxplot_before.png"),
  p_lengths_before,
  width = 12,
  height = 8
)

p_lengths_after <- orf_length_raw %>%
  filter(stage == "after") %>%
  ggplot(aes(x = marker_group_id, y = orf_length)) +
  geom_boxplot(outlier.size = 0.5, fill = "#00BFC4") +
  coord_flip() +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = "Marker group",
    y = "ORF length",
    title = "ORF length distributions after filtering"
  ) +
  plot_theme +
  theme(legend.position = "none")

ggsave(
  file.path(OUT_LOGS, "plot_orf_length_boxplot_after.png"),
  p_lengths_after,
  width = 12,
  height = 8
)

# Missing markers plots
plot_missing_data <- missing_marker_summary %>%
  select(marker_group_id, marker, n_missing_before, n_missing_after) %>%
  pivot_longer(
    cols = c(n_missing_before, n_missing_after),
    names_to = "stage",
    values_to = "n_genomes"
  ) %>%
  mutate(
    stage = recode(
      stage,
      n_missing_before = "before",
      n_missing_after = "after"
    )
  )

p_missing_before <- plot_missing_data %>%
  filter(stage == "before") %>%
  ggplot(aes(x = reorder(marker_group_id, n_genomes), y = n_genomes)) +
  geom_col(fill = "#F8766D") +
  coord_flip() +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = "Marker group",
    y = "Genomes with missing marker",
    title = "Missing expected markers before filtering"
  ) +
  plot_theme +
  theme(legend.position = "none")

ggsave(
  file.path(OUT_LOGS, "plot_missing_expected_markers_before.png"),
  p_missing_before,
  width = 12,
  height = 8
)

p_missing_after <- plot_missing_data %>%
  filter(stage == "after") %>%
  ggplot(aes(x = reorder(marker_group_id, n_genomes), y = n_genomes)) +
  geom_col(fill = "#00BFC4") +
  coord_flip() +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = "Marker group",
    y = "Genomes with missing marker",
    title = "Missing expected markers after filtering"
  ) +
  plot_theme +
  theme(legend.position = "none")

ggsave(
  file.path(OUT_LOGS, "plot_missing_expected_markers_after.png"),
  p_missing_after,
  width = 12,
  height = 8
)
library(tidyverse)

genome_map_tsv <- "output/config/genome_marker_map.tsv"

dir.create(OUT_LOGS, recursive = TRUE, showWarnings = FALSE)

# Plot theme
plot_theme <- theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

# Helpers
write_log <- function(x, file_name) {
  write_tsv(x, file.path(OUT_LOGS, file_name))
}

save_log_plot <- function(plot, file_name, width = 12, height = 8) {
  ggsave(
    file.path(OUT_LOGS, file_name),
    plot,
    width = width,
    height = height
  )
}

add_metadata <- function(df) {
  df %>%
    select(
      -any_of(c(
        "group_id", "marker", "Virus_names", "ICTV_ID",
        taxonomy_all
      ))
    ) %>%
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
    )
}

plot_stage_bar <- function(df, stage_value, file_name, title, y_label) {
  p <- df %>%
    filter(stage == stage_value) %>%
    ggplot(aes(x = reorder(marker_group_id, n_genomes), y = n_genomes)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~marker, scales = "free_y") +
    labs(
      x = "Marker group",
      y = y_label,
      title = title
    ) +
    plot_theme +
    theme(legend.position = "none")

  save_log_plot(p, file_name)
}

plot_length_boxplot <- function(df, stage_value, file_name, title) {
  p <- df %>%
    filter(stage == stage_value) %>%
    ggplot(aes(x = marker_group_id, y = orf_length)) +
    geom_boxplot(outlier.size = 0.5) +
    coord_flip() +
    facet_wrap(~marker, scales = "free_y") +
    labs(
      x = "Marker group",
      y = "ORF length",
      title = title
    ) +
    plot_theme +
    theme(legend.position = "none")

  save_log_plot(p, file_name)
}

# Add metadata to ranked hits
tblout_ranked_log <- tblout_ranked %>%
  add_metadata() %>%
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

# Keep selected ORFs
selected_orfs_log <- selected_orfs %>%
  add_metadata()

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

write_log(log_genome_marker, "log_genome_marker_before_after.tsv")

# Multicopy genomes only
log_duplicated_genomes <- log_genome_marker %>%
  filter(multicopy_before | multicopy_after)

write_log(log_duplicated_genomes, "log_duplicated_genomes.tsv")

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

write_log(log_marker_summary, "log_marker_summary.tsv")

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

write_log(orf_shared_log, "log_orf_shared_across_markers.tsv")

# ORF lengths
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

write_log(orf_length_raw, "log_orf_length_raw.tsv")

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

write_log(orf_length_summary, "log_orf_length_summary.tsv")

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

write_log(missing_marker_log, "log_missing_expected_markers.tsv")

missing_marker_summary <- missing_marker_log %>%
  group_by(marker_group_id, group_id, marker) %>%
  summarise(
    n_expected_genomes = n(),
    n_missing_before = sum(missing_before, na.rm = TRUE),
    n_missing_after = sum(missing_after, na.rm = TRUE),
    .groups = "drop"
  )

write_log(missing_marker_summary, "log_missing_expected_markers_summary.tsv")

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

write_log(database_summary, "database_summary.tsv")

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

write_log(marker_summary, "marker_summary.tsv")

# Family x marker summary
family_marker_summary <- selected_orfs_log %>%
  mutate(
    Family = replace_na(Family, "Unclassified"),
    Genus = replace_na(Genus, "Unclassified"),
    Species = replace_na(Species, "Unclassified")
  ) %>%
  group_by(Family, marker_group_id, group_id, marker, virus_id) %>%
  summarise(
    n_copies = n(),
    Genus = first(Genus),
    Species = first(Species),
    .groups = "drop"
  ) %>%
  group_by(Family, marker_group_id, group_id, marker) %>%
  summarise(
    n_virus_genomes = n_distinct(virus_id),
    n_selected_orfs = sum(n_copies),
    n_multicopy_viruses = sum(n_copies > 1),
    n_genera = n_distinct(Genus),
    n_species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  arrange(Family, group_id, marker, marker_group_id)

write_log(family_marker_summary, "family_marker_summary.tsv")

# Wide version for README
family_marker_summary_wide <- family_marker_summary %>%
  mutate(marker_label = str_remove(marker_group_id, "_.*")) %>%
  select(Family, marker_label, n_virus_genomes) %>%
  pivot_wider(
    names_from = marker_label,
    values_from = n_virus_genomes,
    values_fill = 0
  ) %>%
  arrange(Family) %>%
  bind_rows(
    summarise(., across(where(is.numeric), sum), Family = "Total genomes")
  )

write_log(family_marker_summary_wide, "family_marker_summary_wide.tsv")


family_marker_summary_wide_md <- knitr::kable(
  family_marker_summary_wide,
  format = "pipe"
)

write_lines(
  family_marker_summary_wide_md,
  "output/hmm/logs/family_marker_summary_wide.md"
)


# Family x marker figure
p_family_marker_summary <- family_marker_summary %>%
  mutate(
    Family = fct_reorder(Family, n_virus_genomes, .fun = sum),
    marker_group_id = fct_reorder(marker_group_id, n_virus_genomes, .fun = sum)
  ) %>%
  ggplot(aes(x = marker_group_id, y = Family, fill = n_virus_genomes)) +
  geom_tile() +
  labs(
    x = "Marker group",
    y = "Family",
    fill = "Virus genomes",
    title = "Marker diversity by viral family"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_log_plot(
  p_family_marker_summary,
  "plot_family_marker_summary.png",
  width = 14,
  height = 9
)


# Marker coverage per genome
marker_coverage_distribution <- selected_orfs_log %>%
  distinct(virus_id, marker_group_id) %>%
  count(virus_id, name = "n_markers_detected") %>%
  count(n_markers_detected, name = "n_viruses")

write_log(marker_coverage_distribution, "marker_coverage_distribution.tsv")

# Covered genomes
covered_viruses <- selected_orfs_log %>%
  distinct(virus_id) %>%
  mutate(is_covered = TRUE)

# Global taxonomy coverage
make_taxonomy_coverage <- function(taxo_rank_plot, top_n = 20) {
  taxo_name <- str_to_lower(taxo_rank_plot)

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

  write_log(
    taxonomy_coverage_summary,
    str_c("taxonomy_coverage_summary_", taxo_name, ".tsv")
  )

  taxonomy_coverage_plot <- taxonomy_coverage_summary %>%
    slice_head(n = top_n)

  p <- taxonomy_coverage_plot %>%
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

  save_log_plot(
    p,
    str_c("plot_taxonomy_coverage_", taxo_name, ".png"),
    width = 10,
    height = 7
  )
}

# Taxonomy coverage by marker
make_taxonomy_coverage_by_marker <- function(taxo_rank_plot, top_n = 12) {
  taxo_name <- str_to_lower(taxo_rank_plot)

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

  write_log(
    taxonomy_coverage_summary,
    str_c("taxonomy_coverage_summary_", taxo_name, "_by_marker.tsv")
  )

  taxonomy_coverage_plot <- taxonomy_coverage_summary %>%
    group_by(marker) %>%
    slice_max(order_by = n_viruses_total, n = top_n, with_ties = FALSE) %>%
    ungroup()

  p <- taxonomy_coverage_plot %>%
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

  save_log_plot(
    p,
    str_c("plot_taxonomy_coverage_", taxo_name, "_by_marker.png"),
    width = 14,
    height = 10
  )
}

make_taxonomy_coverage("Family")
make_taxonomy_coverage("Genus")
make_taxonomy_coverage_by_marker("Family")
make_taxonomy_coverage_by_marker("Genus")

# Multicopy plots
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

plot_stage_bar(
  plot_multicopy_data,
  "before",
  "plot_multicopy_before.png",
  "Multicopy genomes before filtering",
  "Genomes with multicopy"
)

plot_stage_bar(
  plot_multicopy_data,
  "after",
  "plot_multicopy_after.png",
  "Multicopy genomes after filtering",
  "Genomes with multicopy"
)

# ORF length plots
plot_length_boxplot(
  orf_length_raw,
  "before",
  "plot_orf_length_boxplot_before.png",
  "ORF length distributions before filtering"
)

plot_length_boxplot(
  orf_length_raw,
  "after",
  "plot_orf_length_boxplot_after.png",
  "ORF length distributions after filtering"
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

plot_stage_bar(
  plot_missing_data,
  "before",
  "plot_missing_expected_markers_before.png",
  "Missing expected markers before filtering",
  "Genomes with missing marker"
)

plot_stage_bar(
  plot_missing_data,
  "after",
  "plot_missing_expected_markers_after.png",
  "Missing expected markers after filtering",
  "Genomes with missing marker"
)
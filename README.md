# MultiSiteMicrobiome
Multi_Bodysite Microbiome Dynamics
# ==============================================================================
# MULTI-BODYSITE MICROBIOME DIVERSITY: MOTHER VS NEONATE
# ==============================================================================
suppressMessages({
  library("tidyverse")
  library("phyloseq")
  library("rstatix")
  library("vegan")
  library("picante")
  library("ggpubr")
  library("patchwork")
})

# --- 1. GLOBAL SETTINGS & DATA COHESION ---
ps_global <- readRDS("ps_gtdb_final.rds")

# Target Body Sites mapped to their specific figure tags
body_sites  <- c("rectal", "vagina", "nasal", "oral", "skin")
site_titles  <- c("rectal" = "Rectal", "vagina" = "Vaginal", "nasal" = "Nasal", "oral" = "Oral", "skin" = "Skin")

# Sequential lettering: A-E for Alpha (Row 1), F-J for Beta (Row 2)
site_alpha_tags <- c("rectal" = "(A)", "vagina" = "(B)", "nasal" = "(C)", "oral" = "(D)", "skin" = "(E)")
site_beta_tags  <- c("rectal" = "(F)", "vagina" = "(G)", "nasal" = "(H)", "oral" = "(I)", "skin" = "(J)")

# Unified plot aesthetics
global_colors <- c("mother" = "#E41A1C", "neonate" = "#377EB8")
order_groups  <- c("mother", "neonate")

# Containers for assembled site-specific subplots
main_plots_list <- list()
supp_plots_list <- list()

# --- 2. UNIVERSAL PLOTTING FUNCTIONS ---

# Generalized Alpha Diversity Plotter
create_alpha_plot <- function(df, metric, title, tag) {
  plot_df <- df %>% drop_na(all_of(metric))
  
  n_moth <- sum(plot_df$Group == "mother")
  n_neon <- sum(plot_df$Group == "neonate")
  
  x_labs <- c(
    "mother" = paste0("Mother\n(n=", n_moth, ")"),
    "neonate" = paste0("Neonate\n(n=", n_neon, ")")
  )
  
  res <- plot_df %>% wilcox_test(as.formula(paste(metric, "~ Group")))
  p_lab <- if(res$p < 0.001) "p < 0.001" else paste0("p = ", round(res$p, 4))
  stats_sub <- paste0("Wilcoxon: ", p_lab)
  
  ggboxplot(plot_df, x = "Group", y = metric, fill = "Group", palette = global_colors,
            add = "jitter", add.params = list(size = 4, alpha = 0.4), outlier.shape = NA) +
    labs(title = paste(tag, title), subtitle = stats_sub, x = "") +
    scale_x_discrete(labels = x_labs) + 
    theme_bw(base_size = 16) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 13, face = "bold.italic"),
          legend.position = "none")
}

# Generalized Beta Diversity Plotter
create_beta_plot <- function(ps_obj, dist_method, title_label, label_tag, show_legend = FALSE) {
  dist_mat <- phyloseq::distance(ps_obj, method = dist_method)
  curr_meta <- data.frame(sample_data(ps_obj))
  
  perm_res <- adonis2(dist_mat ~ mother_neonate, data = curr_meta)
  r2_val   <- round(perm_res$R2[1], 3)
  p_val    <- perm_res$`Pr(>F)`[1]
  p_lab    <- if(p_val < 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
  stats_sub <- paste0("PERMANOVA: R² = ", r2_val, " | ", p_lab)
  
  pcoa <- cmdscale(dist_mat, eig = TRUE)
  var_exp <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
  df_plot <- data.frame(Axis.1 = pcoa$points[,1], Axis.2 = pcoa$points[,2], 
                        Group = factor(curr_meta$mother_neonate, levels = order_groups))
  
  n_labs <- c(paste0("Mother (n=", sum(df_plot$Group == "mother"), ")"), 
              paste0("Neonate (n=", sum(df_plot$Group == "neonate"), ")"))
  
  p <- ggplot(df_plot, aes(Axis.1, Axis.2, fill = Group)) +
    stat_ellipse(aes(color = Group), linewidth = 1, show.legend = FALSE) +
    geom_point(shape = 21, size = 5, color = "black", stroke = 1) +
    scale_fill_manual(values = global_colors, labels = n_labs) +
    scale_color_manual(values = global_colors) +
    labs(title = paste(label_tag, title_label), subtitle = stats_sub,
         x = paste0("PCoA 1 (", var_exp[1], "%)"), y = paste0("PCoA 2 (", var_exp[2], "%)"),
         fill = "Participant Group") +
    theme_bw(base_size = 16) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 13, face = "bold.italic"))
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  return(p)
}

# --- 3. PROCESSING LOOP ACROSS ALL BODY SITES ---
for (site in body_sites) {
  message("Processing body site: ", site)
  
  # Subset and clean data contextually
  ps.site <- ps_global %>%
    subset_samples(all_bodysite == site) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    subset_taxa(Kingdom == "Bacteria" & !is.na(Phylum))
  
  # Check if sample pool contains both groups to prevent math breaks
  meta_check <- data.frame(sample_data(ps.site))
  if(length(unique(meta_check$mother_neonate)) < 2) next
  
  # Matrix translation
  otu_mat <- as.matrix(otu_table(ps.site))
  if (taxa_are_rows(ps.site)) { otu_mat <- t(otu_mat) }
  
  # Math Engines
  shannon_vec <- diversity(otu_mat, index = "shannon")
  obs_vec     <- specnumber(otu_mat)
  pielou_vec  <- shannon_vec / log(obs_vec)
  pielou_vec[is.infinite(pielou_vec)] <- 0 
  faith_pd    <- pd(otu_mat, phy_tree(ps.site), include.root = FALSE)
  
  # Consolidate Alpha Table
  adiv_df <- data.frame(
    Shannon  = shannon_vec, 
    Observed = obs_vec, 
    FaithPD  = faith_pd$PD, 
    Pielou   = pielou_vec,
    Group    = factor(meta_check$mother_neonate, levels = order_groups)
  )
  
  # Extract correct tags for this loop step
  alpha_tag    <- site_alpha_tags[site]
  beta_tag     <- site_beta_tags[site]
  current_name <- site_titles[site]
  
  # --- MAIN FIGURE SUBPLOTS ---
  p_alpha <- create_alpha_plot(adiv_df, "Shannon", paste(current_name, "Shannon"), alpha_tag)
  p_beta  <- create_beta_plot(ps.site, "bray", paste(current_name, "Bray-Curtis PCoA"), beta_tag)
  
  # Stack Alpha directly on top of its respective Beta plot
  site_main_composite <- (p_alpha / p_beta) + plot_layout(heights = c(1, 1.1))
  main_plots_list[[site]] <- site_main_composite
  
  # --- SUPPLEMENTARY SUBPLOTS ---
  s1 <- create_alpha_plot(adiv_df, "Observed", "Observed Richness", "(A)")
  s2 <- create_alpha_plot(adiv_df, "FaithPD", "Faith's PD", "(B)")
  s3 <- create_alpha_plot(adiv_df, "Pielou", "Pielou's Evenness", "(C)")
  s4 <- create_beta_plot(ps.site, "wunifrac", "Weighted UniFrac", "(D)")
  s5 <- create_beta_plot(ps.site, "unifrac", "Unweighted UniFrac", "(E)")
  s6 <- create_beta_plot(ps.site, "jaccard", "Jaccard Index", "(F)", show_legend = TRUE)
  
  supp_plots_list[[site]] <- (s1 | s2 | s3) / (s4 | s5 | s6) + 
    plot_annotation(subtitle = paste(current_name, "Microbiome Matrix Metrics"),
                    theme = theme(plot.subtitle = element_text(size = 22, face = "bold.italic")))
}

# --- 4. COMPOSITING THE GRAND FIGURES WITH PATCHWORK ---

# Dynamic multi-line caption building for publication layout clarity
fig_subtitle <- paste(
  "Mother cohorts are marked in Red, Neonate cohorts are marked in Blue."
)

# Main Composite Layout
main_fig <- (main_plots_list[["rectal"]] | main_plots_list[["vagina"]] | 
               main_plots_list[["nasal"]] | main_plots_list[["oral"]] | 
               main_plots_list[["skin"]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Figure 4. Multi-Site Comparative Microbiome Diversity Profiles",
                  subtitle = fig_subtitle,
                  theme = theme(plot.title = element_text(size = 26, face = "bold"),
                                plot.subtitle = element_text(size = 15, face = "plain", margin = margin(t = 5, b = 20)),
                                legend.position = "bottom"))

# Supplementary Outputs Configuration
supp_fig <- (supp_plots_list[["rectal"]] /  supp_plots_list[["vagina"]] / 
               supp_plots_list[["nasal"]] /  supp_plots_list[["oral"]] / 
               supp_plots_list[["skin"]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Supplementary Figure 4. Comprehensive Multi-Site Diversity Metrics Across All Ecosystems",
                  theme = theme(plot.title = element_text(size = 30, face = "bold")))

# --- 5. EXPORT SYSTEM ---
ggsave("Figure_4_All_Bodysites_Main_Diversity.png", main_fig, width = 30, height = 15, dpi = 300)
ggsave("Figure_S4_All_Bodysites_Supp_Diversity.png", supp_fig, width = 22, height = 45, dpi = 200, limitsize = FALSE)

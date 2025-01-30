#!/usr/bin/env Rscript

scale_fill_scico_mid <- function(..., mid = 0, alpha = NULL, begin = 0, end = 1, direction = 1, reverse = TRUE ,palette = "broc") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for this functionality", call. = FALSE)
  }
  force(mid)
  ggplot2::continuous_scale(
    aesthetics = "fill", 
    scale_name = "gradient2",
    palette = scales::gradient_n_pal(
      colours = scico(256, alpha, begin, end, direction, palette), 
      values = NULL, space = "Lab"),
    guide="colourbar",
    rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
      scales::rescale_mid(x, to, from, mid)
    },
    ...
  )
}

compute_difference <- function(test, control, sumstats) {
    
    df1 <- sumstats %>% filter(condition == test) %>% dplyr::select(-condition)
    df2 <- sumstats %>% filter(condition == control) %>% dplyr::select(-condition)
    
    df <- inner_join(df1, df2,
                     by = c("chunk", "pos", "aa"))
    
    new_stats <- df %>%
        mutate(estimate = estimate.x - estimate.y,
               std.error = sqrt(std.error.x^2 + std.error.y^2)) %>%
        dplyr::select(chunk, pos, aa, estimate, std.error) %>%
        ungroup()
    
    new_stats$condition = paste0(test, " - ", control)
    
    return(new_stats)
    
}

plot_coverage <- function(bcs, wt_aa, sample_id, destdir) {
    
    bcs_subset <- bcs %>%
        filter(sample == sample_id)
    
    allpos <- expand_grid("pos" = unique(bcs$pos),
                     "mut_aa" = unique(bcs$mut_aa)) %>%
        filter(!paste0(pos, mut_aa) %in% paste0(wt_aa$pos, wt_aa$wt_aa))

    bcs_subset_expanded <- allpos %>%
        left_join(bcs_subset) %>%
        mutate(n = ifelse(is.na(n), 0, n))
    
    the_plot <- ggplot(bcs_subset_expanded) +
        geom_line(aes(x = as.numeric(pos), y = log10(n))) +
        theme_bw(base_size = 16) +
        ggtitle(paste0("Sample ", sample_id)) +
        facet_wrap(~mut_aa, ncol = 1) +
        xlab("TYK2 position") + ylab("log10(number of unique barcodes)")
    
    ggsave(str_c("notebooks/coverage-plots/", destdir, "-sample", sample_id, ".pdf"),
           the_plot, width = 10, height = 20)
    
    return(the_plot)
}

make_profile <- function(poss, df){
    
    p1 <- df %>%
        filter(pos == poss) %>%
        ggplot() +
            geom_pointrange(aes(x = aa, y = estimate, color = condition,
                                ymin = estimate - 2*std.error,
                                ymax = estimate + 2*std.error),
                            position = position_dodge(width = 0.4)) +
            theme_pubr(base_size = 16) +
            geom_hline(yintercept = 0) + ggtitle(poss) +
            scale_color_manual(values = c("red", "black"))

    return(p1)
}

generate_resamples <- function(mean_vec, se_vec, metadata, num){
    
    df <- MASS::mvrnorm(n = num,
                  mu = mean_vec,
                  Sigma = diag(se_vec^2)) %>%
        as_tibble() %>%
        mutate(n = row_number()) %>%
        rename(c("0.125" = "V1", "0.375" = "V2", "0.625" = "V3", "0.875" = "V4")) %>%
        pivot_longer(names_to = "bin", values_to = "value", `0.125`:`0.875`)

    df_score <- df %>%
        group_by(n) %>% 
        mutate(bin = as.numeric(bin),
               value = (value - min(value)) / sum(value - min(value))) %>%
        summarize(score = sum(bin*value)) 
    
    return(bind_cols(tibble(
        "score_mean" = mean(df_score$score),
        "score_se" = sd(df_score$score)),
        "chunk" = metadata[1],
        "pos" = metadata[2],
        "aa" = metadata[3]))

}

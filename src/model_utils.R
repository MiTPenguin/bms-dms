#!/usr/bin/env Rscript

library(magrittr)
library(data.table)
library(future)
library(furrr)
library(broom)
library(broom.mixed)
library(glmmTMB)
library(emmeans)
library(future.callr)
library(tidyverse)

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

select_complete_mutants <- function(df, cvar) {

    expected_conditions <- df %>%
        ungroup() %>%
        select(all_of(cvar)) %>%
        distinct() %>%
        unlist() %>%
        length()
    
    valid_aa <- df %>%
        ungroup() %>%
        rename("cvar" = {{cvar}}) %>%
        filter(mut_aa != "WT") %>%
        group_by(mut_aa) %>%
        summarize(num_cond_found = length(unique(cvar))) %>%
        filter(num_cond_found == expected_conditions) %>%
        select(mut_aa) %>% unlist()

    valid_aa <- c(as.character(valid_aa), "WT")
    
    df_filt <- df %>%
        ungroup() %>%
        filter(mut_aa %in% valid_aa)

    return(df_filt)

}

generate_formula <- function(keyword = "global") {

    if (keyword == "global") {
        formula <- as.formula(count ~ -1 + condition + condition:mut_aa + (1 | barcode) + offset(stop_counts))
    } else if (keyword == "interaction") {
        formula <- as.formula(count ~ -1 + condition + mut_aa + condition:mut_aa + (1 | barcode) + offset(stop_counts))
    } else if (keyword == "flow") {
        formula <- as.formula(count ~ -1 + condition_conc + condition_conc:mut_aa + (1 | barcode))
    } else if (keyword == "drc") {
        formula <- as.formula(count ~ -1 + condition + condition:mut_aa + (1 | condition:barcode) + offset(stop_counts))
    } else {
        stop("Invalid model type")
    }
    return(formula)
}

rand_effect <- function(data, mod_path, formula) {

    # Model fitting parameters:
    ## L-BFGS-B is a particularly memory efficient quasi-Newton method

    ## start is the starting value for the model, and given the log link
    ### this value closely matches the end of the optimization
    ### for most positions and seemed to make sense

    ## rel.tol is the relative tolerance for convergence
    ## pgtol is the tolerance for the projected gradient
    ### Both of these are reduced to ensure the model converges, and
    ### decreasing them did not obviously alter model fit, though
    ### a formal grid search was not performed

    to_return <- tryCatch({

        if (any(grepl("condition_conc", formula))){
            cond_var <- "condition_conc"
        } else {
            cond_var <- "condition"
        }

        data_filt <- select_complete_mutants(data, cond_var)

        mod <- glmmTMB(formula = formula,
            start = -1,
            REML = FALSE,
            control = glmmTMBControl(optimizer = optim,
                profile = TRUE,
                optArgs = list(method = "L-BFGS-B",
                    pgtol = 0,
                    rel.tol = 0.1)),
            data = data_filt,
            sparseX = c(cond=TRUE),
            family = nbinom2)

        saveRDS(mod, file = str_c(mod_path, ".RDS"))

        coefs <- broom.mixed::tidy(mod) %>%
            unnest_longer(term) %>%
            mutate(dispersion = sigma(mod),
                estimate = estimate / log(2),
                std.error = std.error / log(2)) %>%
            rename("log2FoldChange" = "estimate",
                "log2StdError" = "std.error")

        marginals <- broom::tidy(emmeans(mod, as.formula(str_c("~ mut_aa + ", cond_var)))) %>%
            mutate(estimate = estimate / log(2),
                std.error = std.error / log(2)) %>%
            rename("log2Marginal" = "estimate",
                "log2MarginalError" = "std.error")
        
        remove(data_filt, mod)

        list("coefs"=coefs, "marginals"=marginals)

        }, error = function(e) {

            message("Error fitting model")
            NULL

        }
    )

    return(to_return)
}

rand_effect_wrap <- function(mapped_counts, form, nestvars, model_output_path, doubles, nworkers) {

    plan(callr, workers = nworkers)
    allvars <- names(mapped_counts %>% select(-pos))

    wt_df <- mapped_counts %>%
        filter(mut_aa == "WT") %>%
        select(-pos) %>%
        nest(wt = all_of(setdiff(allvars, nestvars)))

    nested_counts <- mapped_counts %>%
        filter(mut_aa != "WT") %>%
        nest(data = all_of(setdiff(allvars, c(nestvars, "pos"))))

    joined_counts <- inner_join(nested_counts, wt_df) %>%
        mutate(df = map2(data, wt, bind_rows),
            name = str_c("chunk", chunk, "pos", pos, sep = "_"),
            full_path = str_c(model_output_path, name, sep = ""))

    fit_df <- joined_counts %>%
        mutate(sumstats = future_map2(df, full_path, rand_effect,
            formula = form,
            .options = furrr_options(seed = .Random.seed, scheduling = FALSE)))

    sumstats_wide <- fit_df %>%
        select(all_of(nestvars), pos, sumstats) %>%
        unnest_wider(sumstats)

  return(sumstats_wide)
}
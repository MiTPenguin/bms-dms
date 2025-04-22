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
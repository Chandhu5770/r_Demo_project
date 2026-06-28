# UAV Remote Sensing Tutorial: Vegetation Indices and SPAD Prediction
# -----------------------------------------------------------------------------
# Purpose
#   This script is a teaching-friendly version of the UAV/SPAD workflow.
#   It shows students how to:
#     1. read field boundaries, subplot polygons, a multispectral UAV raster, and SPAD data;
#     2. compute vegetation indices;
#     3. extract raster statistics for sampling plots;
#     4. select simple predictors using correlation;
#     5. fit and compare simple linear models with repeated cross-validation;
#     6. map predicted SPAD across the field;
#     7. export tables, figures, rasters, and vectors.
#
# Notes for students
#   - SPAD is used as the response variable in this tutorial.
#   - The model predicts SPAD, not biomass. If your response variable is biomass,
#     rename the CSV column and update the code consistently.
#   - We use base lm() and caret for cross-validation to keep the modelling simple.
#   - This is a tutorial workflow. Real research projects need additional checks,
#     especially calibration, spatial validation, prediction uncertainty, and
#     extrapolation beyond observed field data.
# -----------------------------------------------------------------------------

# 1. Packages -----------------------------------------------------------------
# Install these packages once if needed:
# install.packages(c("terra", "sf", "tidyverse", "caret", "viridis", "janitor", "broom"))

required_packages <- c(
  "terra",      # raster reading, processing, prediction, export
  "sf",         # vector data reading and coordinate reference systems
  "tidyverse",  # data wrangling and ggplot2 plotting
  "caret",      # simple train/test split and repeated cross-validation
  "viridis",    # colour scale for raster maps
  "janitor",    # clean column names
  "broom"       # tidy summaries of lm() objects
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the missing packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Reproducibility: the same random seed gives the same train/test split and CV folds.
set.seed(26)

# 2. Project paths ------------------------------------------------------------
# Change this path if your project folder is elsewhere.
# If this directory does not exist, the script uses the current working directory.
project_dir <- "/home/pc4dl/UAV_R/"
if (dir.exists(project_dir)) {
  setwd(project_dir)
}

# Input files. Keep the data in a folder called "data" inside the project folder.
field_file <- "data/am_sande_field.gpkg"
uav_file   <- "data/ms_image.tif"
plots_file <- "data/am_sande_subplots.gpkg"
spad_file  <- "data/am_sande_spad.csv"

# Output folders. The script creates these automatically.
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/rasters", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/vectors", recursive = TRUE, showWarnings = FALSE)

# A small helper function to stop early if an input file is missing.
check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
}

invisible(lapply(c(field_file, uav_file, plots_file, spad_file), check_file))

# 3. Read the field boundary --------------------------------------------------
field_poly <- st_read(field_file, quiet = TRUE) %>%
  st_make_valid()

print(field_poly)
print(st_crs(field_poly))
print(st_bbox(field_poly))

p_field <- ggplot(field_poly) +
  geom_sf(fill = "grey90", colour = "black", linewidth = 0.8) +
  labs(title = "Field boundary") +
  theme_minimal()

print(p_field)

# 4. Read the UAV multispectral image ----------------------------------------
uav <- rast(uav_file)

# Adjust these names if your raster band order is different.
# This tutorial assumes: blue, green, red, red edge, near-infrared.
stopifnot(nlyr(uav) >= 5)
names(uav)[1:5] <- c("blue", "green", "red", "red_edge", "nir")

print(uav)
print(res(uav))
print(ext(uav))
print(crs(uav))
print(terra::minmax(uav))

# Base R raster plots are useful for quick visual checks.
plot(uav[[1:5]], nc = 2, main = names(uav)[1:5])

# True-colour composite: red-green-blue.
plotRGB(
  uav,
  r = which(names(uav) == "red"),
  g = which(names(uav) == "green"),
  b = which(names(uav) == "blue"),
  stretch = "hist",
  main = "True-colour composite"
)

# False-colour composite: NIR-red-green.
# Healthy vegetation often appears bright in this display.
plotRGB(
  uav,
  r = which(names(uav) == "nir"),
  g = which(names(uav) == "red"),
  b = which(names(uav) == "green"),
  stretch = "hist",
  main = "False-colour composite: NIR-Red-Green"
)

# 5. Clip the UAV image to the field boundary --------------------------------
# The field polygon and UAV raster must use the same coordinate reference system.
field_poly_uav <- st_transform(field_poly, terra::crs(uav))
field_vect <- vect(field_poly_uav)

# Crop first to the bounding box, then mask to the exact field shape.
uav_field <- terra::crop(uav, field_vect) %>%
  terra::mask(field_vect)

print(uav_field)

plotRGB(
  uav_field,
  r = "red", g = "green", b = "blue",
  stretch = "lin",
  main = "UAV image clipped to field boundary"
)
plot(field_vect, add = TRUE, border = "yellow", lwd = 2)

# 6. Compute vegetation indices ----------------------------------------------
# The small value eps avoids division by zero.
eps <- 1e-10

ndvi <- (uav_field$nir - uav_field$red) /
  (uav_field$nir + uav_field$red + eps)

ndre <- (uav_field$nir - uav_field$red_edge) /
  (uav_field$nir + uav_field$red_edge + eps)

gndvi <- (uav_field$nir - uav_field$green) /
  (uav_field$nir + uav_field$green + eps)

# EVI is included for teaching. It is most interpretable when bands are properly
# calibrated/scaled reflectance values.
evi <- 2.5 * (uav_field$nir - uav_field$red) /
  (uav_field$nir + 6 * uav_field$red - 7.5 * uav_field$blue + 1 + eps)

names(ndvi)  <- "ndvi"
names(ndre)  <- "ndre"
names(gndvi) <- "gndvi"
names(evi)   <- "evi"

vi_stack <- c(ndvi, ndre, gndvi, evi)

# The predictor stack includes original bands and vegetation indices.
predictor_stack <- c(
  uav_field[[c("blue", "green", "red", "red_edge", "nir")]],
  vi_stack
)

print(predictor_stack)

plot(vi_stack, nc = 2, col = viridis::viridis(100), main = toupper(names(vi_stack)))


# 7. Create simple NDRE management zones -------------------------------------
# Here we split NDRE into low, medium, and high zones using tertiles.
# These are data-driven teaching zones, not agronomic recommendation thresholds.
ndre_quantiles <- terra::global(ndre, quantile, probs = c(0.33, 0.67), na.rm = TRUE)
q33 <- as.numeric(ndre_quantiles[1, 1])
q67 <- as.numeric(ndre_quantiles[1, 2])

zone_matrix <- matrix(
  c(-Inf, q33, 1,
    q33, q67, 2,
    q67, Inf, 3),
  ncol = 3,
  byrow = TRUE
)

ndre_zones <- terra::classify(ndre, zone_matrix)
names(ndre_zones) <- "ndre_zone"
levels(ndre_zones) <- data.frame(
  value = 1:3,
  zone = c("Low", "Medium", "High")
)

plot(ndre_zones, main = "NDRE management zones")

# 8. Read field-sampling polygons and SPAD data ------------------------------
plots <- st_read(plots_file, quiet = TRUE) %>%
  st_make_valid() %>%
  janitor::clean_names()

# The polygon file must contain p_id because p_id is used to join with SPAD data.
stopifnot("p_id" %in% names(plots))
stopifnot(!anyDuplicated(plots$p_id))

plots <- st_transform(plots, terra::crs(uav_field))
print(plots)

plotRGB(
  uav_field,
  r = "red", g = "green", b = "blue",
  stretch = "hist",
  main = "Ground-sampling polygons"
)
plot(vect(plots), add = TRUE, border = "yellow", lwd = 2)

spad_data <- read_csv(spad_file, show_col_types = FALSE) %>%
  janitor::clean_names()

# The CSV must have p_id and spad columns.
stopifnot(all(c("p_id", "spad") %in% names(spad_data)))
print(spad_data)

# 9. Extract raster statistics for each subplot ------------------------------
# Each field subplot covers many pixels. We summarize those pixels using simple
# statistics such as mean, median, quartiles, and standard deviation.
extract_one_stat <- function(r, polygons, fun, suffix, ...) {
  out <- terra::extract(
    r,
    terra::vect(polygons),
    fun = fun,
    na.rm = TRUE,
    ...
  ) %>%
    tibble::as_tibble() %>%
    dplyr::select(-ID)

  names(out) <- paste0(names(out), "_", suffix)
  out
}

q25_fun <- function(x, na.rm = TRUE) {
  stats::quantile(x, probs = 0.25, na.rm = na.rm, names = FALSE)
}

q75_fun <- function(x, na.rm = TRUE) {
  stats::quantile(x, probs = 0.75, na.rm = na.rm, names = FALSE)
}

plot_attributes <- plots %>%
  st_drop_geometry() %>%
  dplyr::select(p_id)

stats_table <- dplyr::bind_cols(
  plot_attributes,
  extract_one_stat(predictor_stack, plots, mean,    "mean"),
  extract_one_stat(predictor_stack, plots, median,  "median"),
  extract_one_stat(predictor_stack, plots, min,     "min"),
  extract_one_stat(predictor_stack, plots, max,     "max"),
  extract_one_stat(predictor_stack, plots, q25_fun, "q25"),
  extract_one_stat(predictor_stack, plots, q75_fun, "q75"),
  extract_one_stat(predictor_stack, plots, sd,      "sd")
)

print(stats_table)
readr::write_csv(stats_table, "outputs/tables/uav_polygon_statistics.csv")

missing_summary <- stats_table %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.x)))) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "n_missing") %>%
  dplyr::arrange(dplyr::desc(n_missing))

print(missing_summary)
readr::write_csv(missing_summary, "outputs/tables/missing_values_summary.csv")

# 10. Join UAV statistics with SPAD data -------------------------------------
analysis_data <- spad_data %>%
  dplyr::inner_join(stats_table, by = "p_id")

# These checks show whether any p_id values failed to match.
spad_without_raster_stats <- dplyr::anti_join(spad_data, stats_table, by = "p_id")
stats_without_spad <- dplyr::anti_join(stats_table, spad_data, by = "p_id")

print(spad_without_raster_stats)
print(stats_without_spad)
print(analysis_data)

readr::write_csv(analysis_data, "outputs/tables/analysis_dataset.csv")

# 11. Correlation analysis ----------------------------------------------------
# Correlation helps identify promising variables, but it is not model validation.
predictor_names <- analysis_data %>%
  dplyr::select(-dplyr::any_of(c("p_id", "spad"))) %>%
  dplyr::select(where(is.numeric)) %>%
  names()

safe_correlation_row <- function(variable_name) {
  x <- analysis_data[[variable_name]]
  y <- analysis_data$spad
  ok <- stats::complete.cases(x, y)

  if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(tibble::tibble(
      variable = variable_name,
      correlation = NA_real_,
      abs_correlation = NA_real_,
      p_value = NA_real_,
      n = sum(ok)
    ))
  }

  test <- stats::cor.test(x[ok], y[ok], method = "pearson")

  tibble::tibble(
    variable = variable_name,
    correlation = unname(test$estimate),
    abs_correlation = abs(unname(test$estimate)),
    p_value = test$p.value,
    n = sum(ok)
  )
}

correlation_table <- purrr::map_dfr(predictor_names, safe_correlation_row) %>%
  dplyr::arrange(dplyr::desc(abs_correlation))

print(correlation_table)
readr::write_csv(correlation_table, "outputs/tables/spad_correlations.csv")

# Select the three strongest variables for simple teaching models.
top3_variables <- correlation_table %>%
  dplyr::filter(!is.na(abs_correlation)) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::pull(variable)

if (length(top3_variables) == 0) {
  stop("No usable predictors were found for modelling.", call. = FALSE)
}

print(top3_variables)

p_corr <- correlation_table %>%
  dplyr::filter(!is.na(correlation)) %>%
  dplyr::slice_head(n = min(20, nrow(correlation_table))) %>%
  dplyr::mutate(variable = forcats::fct_reorder(variable, correlation)) %>%
  ggplot(aes(correlation, variable)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = "Variables most strongly correlated with SPAD",
    x = "Pearson correlation coefficient",
    y = NULL
  ) +
  theme_minimal()

print(p_corr)

# 12. K-fold cross-validation model comparison -------------------------------

# We are NOT making a separate train/test split here.
# Instead, all model comparison is based on k-fold cross-validation.

# Choose 5-fold CV by default.
# For very small datasets, you can change this to 4.
cv_folds <- 5

# Safety checks for small datasets
if (nrow(model_data) < cv_folds) {
  cv_folds <- min(4, nrow(model_data))
}

if (cv_folds < 2) {
  stop("There are too few observations for cross-validation.", call. = FALSE)
}

cat("Number of plots used for modelling:", nrow(model_data), "\n")
cat("Number of CV folds:", cv_folds, "\n")

# -------------------------------------------------------------------------
# Candidate formulas
# -------------------------------------------------------------------------
# We fit one simple linear model for each of the top three UAV variables.
# Example:
#   spad ~ ndre_mean
#   spad ~ gndvi_q25
#   spad ~ nir_median

candidate_formulas <- purrr::map(
  top3_variables,
  ~ stats::as.formula(paste("spad ~", .x))
)

names(candidate_formulas) <- paste0("LM_", top3_variables)

print(candidate_formulas)

# -------------------------------------------------------------------------
# Create the same CV folds for every model
# -------------------------------------------------------------------------
# This is important because every model should be evaluated on the same folds.
# Otherwise, the comparison is less fair.

set.seed(26)

cv_index <- caret::createFolds(
  model_data$spad,
  k = cv_folds,
  returnTrain = TRUE
)

train_control <- caret::trainControl(
  method = "cv",
  number = cv_folds,
  index = cv_index,
  savePredictions = "final"
)

# -------------------------------------------------------------------------
# Fit all candidate models using k-fold CV
# -------------------------------------------------------------------------

caret_models <- purrr::imap(
  candidate_formulas,
  function(formula_i, model_name) {
    caret::train(
      formula_i,
      data = model_data,
      method = "lm",
      trControl = train_control,
      metric = "RMSE"
    )
  }
)

# -------------------------------------------------------------------------
# Summarise CV results
# -------------------------------------------------------------------------

get_caret_metric <- function(result_row, metric_name) {
  if (metric_name %in% names(result_row)) {
    return(as.numeric(result_row[[metric_name]][1]))
  }
  NA_real_
}

cv_summary <- purrr::imap_dfr(
  caret_models,
  function(model_i, model_name) {
    result_i <- model_i$results[1, ]
    
    tibble::tibble(
      model = model_name,
      formula = paste(deparse(candidate_formulas[[model_name]]), collapse = " "),
      CV_RMSE = get_caret_metric(result_i, "RMSE"),
      CV_MAE  = get_caret_metric(result_i, "MAE"),
      CV_R2   = get_caret_metric(result_i, "Rsquared")
    )
  }
) %>%
  dplyr::mutate(
    # Low RMSE is good, so smaller rank is better.
    rank_RMSE = dplyr::min_rank(CV_RMSE),
    
    # High R² is good, so we rank the negative of R².
    rank_R2 = dplyr::min_rank(dplyr::desc(CV_R2)),
    
    # Simple combined score.
    # The best model has low RMSE and high R².
    selection_score = rank_RMSE + rank_R2
  ) %>%
  dplyr::arrange(selection_score, CV_RMSE, dplyr::desc(CV_R2))

print(cv_summary)

readr::write_csv(
  cv_summary,
  "outputs/tables/cross_validation_summary.csv"
)

# -------------------------------------------------------------------------
# Select the best model
# -------------------------------------------------------------------------
# We choose the model with the best combined CV ranking:
#   - lower CV_RMSE is better
#   - higher CV_R2 is better

best_model_name <- cv_summary$model[1]
best_formula <- candidate_formulas[[best_model_name]]
best_caret_model <- caret_models[[best_model_name]]

cat("Best model selected by CV:", best_model_name, "\n")
cat("Best formula:", deparse(best_formula), "\n")

# 13. Plot cross-validated results for the best model only -------------------

# The best model was selected in Section 12 using CV_RMSE and CV_R2.
# Here we plot only the cross-validated predictions from that best model.

cv_predictions_best <- best_caret_model$pred %>%
  dplyr::mutate(
    model = best_model_name,
    p_id = model_data$p_id[rowIndex]
  ) %>%
  dplyr::select(
    p_id,
    model,
    obs,
    pred,
    rowIndex,
    Resample
  ) %>%
  dplyr::arrange(rowIndex)

print(cv_predictions_best)

readr::write_csv(
  cv_predictions_best,
  "outputs/tables/best_model_cv_predictions.csv"
)

# Extract the CV summary for the selected best model
best_cv_summary <- cv_summary %>%
  dplyr::filter(model == best_model_name)

print(best_cv_summary)

# -------------------------------------------------------------------------
# Observed vs predicted plot for the best model
# -------------------------------------------------------------------------

p_cv_obs_pred <- ggplot(cv_predictions_best, aes(obs, pred)) +
  geom_point(size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(
    title = "Observed versus predicted SPAD: best CV model",
    subtitle = paste0(
      best_model_name,
      " | CV RMSE = ", round(best_cv_summary$CV_RMSE, 2),
      " | CV R² = ", round(best_cv_summary$CV_R2, 2)
    ),
    x = "Observed SPAD",
    y = "Predicted SPAD"
  ) +
  theme_minimal()

print(p_cv_obs_pred)

# -------------------------------------------------------------------------
# Residual plot for the best model
# -------------------------------------------------------------------------

p_cv_resid <- cv_predictions_best %>%
  dplyr::mutate(residual = obs - pred) %>%
  ggplot(aes(pred, residual)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    title = "Residual diagnostics: best CV model",
    subtitle = best_model_name,
    x = "Predicted SPAD",
    y = "Observed minus predicted"
  ) +
  theme_minimal()

print(p_cv_resid)
# 14. Select the best model and refit on all complete data -------------------
# We select the model with the lowest cross-validated RMSE.
best_model_name <- cv_summary %>%
  dplyr::slice_min(CV_RMSE, n = 1, with_ties = FALSE) %>%
  dplyr::pull(model)

best_formula <- candidate_formulas[[best_model_name]]
best_model_terms <- attr(stats::terms(best_formula), "term.labels")

cat("Best model:", best_model_name, "\n")
print(best_formula)

final_training_data <- analysis_data %>%
  dplyr::select(p_id, spad, dplyr::all_of(best_model_terms)) %>%
  tidyr::drop_na()

final_lm <- stats::lm(best_formula, data = final_training_data)

print(summary(final_lm))

coef_table <- broom::tidy(final_lm, conf.int = TRUE)
model_fit_table <- broom::glance(final_lm)

print(coef_table)
print(model_fit_table)

# 15. Apply the best model to the whole field image --------------------------
# Important teaching point:
#   The model was trained on polygon summaries such as ndre_mean or gndvi_q25.
#   For mapping, we create a raster layer with the same name and statistic.
#   We use a simple aggregation window close to the ground plot width.
#
#   This is an approximation for teaching. A research workflow should also map
#   prediction uncertainty and flag extrapolation beyond the training range.
parse_model_term <- function(term) {
  stat_suffixes <- c(
    mean = "_mean",
    median = "_median",
    min = "_min",
    max = "_max",
    q25 = "_q25",
    q75 = "_q75",
    sd = "_sd"
  )

  matched_stat <- names(stat_suffixes)[
    vapply(stat_suffixes, function(suffix) grepl(paste0(suffix, "$"), term), logical(1))
  ]

  if (length(matched_stat) != 1) {
    stop("Could not identify statistic suffix for model term: ", term, call. = FALSE)
  }

  suffix <- stat_suffixes[[matched_stat]]
  base_layer <- sub(paste0(suffix, "$"), "", term)

  list(layer = base_layer, stat = matched_stat)
}

make_stat_raster <- function(base_raster, stat, fact) {
  stat_function <- switch(
    stat,
    mean = mean,
    median = median,
    min = min,
    max = max,
    q25 = q25_fun,
    q75 = q75_fun,
    sd = sd,
    stop("Unknown statistic: ", stat, call. = FALSE)
  )

  terra::aggregate(base_raster, fact = fact, fun = stat_function, na.rm = TRUE)
}

build_prediction_stack <- function(model_terms, base_stack, plot_width_m = 1) {
  # This assumes the raster CRS uses metres. If the raster uses degrees, reproject
  # to a suitable projected CRS before using metre-based plot widths.
  pixel_size <- min(terra::res(base_stack))
  aggregation_factor <- max(1, ceiling(plot_width_m / pixel_size))

  raster_list <- lapply(model_terms, function(term) {
    parsed <- parse_model_term(term)

    if (!parsed$layer %in% names(base_stack)) {
      stop("Base raster layer not found: ", parsed$layer, call. = FALSE)
    }

    out <- make_stat_raster(base_stack[[parsed$layer]], parsed$stat, aggregation_factor)
    names(out) <- term
    out
  })

  prediction_stack <- do.call(c, raster_list)
  names(prediction_stack) <- model_terms
  prediction_stack
}

plot_width_m <- 1
prediction_stack <- build_prediction_stack(
  model_terms = best_model_terms,
  base_stack = predictor_stack,
  plot_width_m = plot_width_m
)

print(prediction_stack)

spad_raster <- terra::predict(prediction_stack, final_lm, na.rm = TRUE)
names(spad_raster) <- "predicted_spad"

# SPAD cannot be negative. This does not solve all extrapolation problems, but it
# avoids impossible negative map values from a linear model.
spad_raster <- terra::clamp(spad_raster, lower = 0, upper = Inf, values = TRUE)

print(spad_raster)

plot(
  spad_raster,
  col = viridis::viridis(100),
  main = "Predicted SPAD"
)
plot(field_vect, add = TRUE, border = "black", lwd = 1)


# 16. Export results ----------------------------------------------------------
# GeoTIFF outputs.


terra::writeRaster(
  spad_raster,
  "outputs/rasters/predicted_spad.tif",
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW")
)


sessionInfo()

normalize_scissor_family <- function(family) {
  match.arg(family, c("gaussian", "binomial", "cox"))
}

coerce_numeric_phenotype <- function(phenotype, family) {
  if (is.matrix(phenotype) || is.data.frame(phenotype)) {
    stop(sprintf("For family = '%s', phenotype must be a vector.", family))
  }
  if (anyNA(phenotype)) {
    stop(sprintf("For family = '%s', phenotype contains missing values.", family))
  }

  phenotype_chr <- if (is.factor(phenotype)) as.character(phenotype) else phenotype
  phenotype_num <- suppressWarnings(as.numeric(phenotype_chr))

  if (length(phenotype_num) != length(phenotype_chr) || anyNA(phenotype_num)) {
    stop(sprintf("For family = '%s', phenotype must be numeric or coercible to numeric without loss.", family))
  }

  phenotype_num
}

prepare_binomial_response <- function(phenotype, tag = NULL) {
  if (is.matrix(phenotype) || is.data.frame(phenotype)) {
    stop("For family = 'binomial', phenotype must be a vector.")
  }
  if (anyNA(phenotype)) {
    stop("For family = 'binomial', phenotype contains missing values.")
  }

  if (is.logical(phenotype)) {
    level_labels <- c("FALSE", "TRUE")
    y <- as.integer(phenotype)
  } else if (is.numeric(phenotype)) {
    level_values <- sort(unique(as.numeric(phenotype)))
    if (length(level_values) != 2) {
      stop("For family = 'binomial', phenotype must contain exactly two distinct values.")
    }
    level_labels <- as.character(level_values)
    y <- as.integer(factor(phenotype, levels = level_values)) - 1L
  } else {
    phenotype_fac <- factor(phenotype)
    if (nlevels(phenotype_fac) != 2) {
      stop("For family = 'binomial', phenotype must contain exactly two levels.")
    }
    level_labels <- levels(phenotype_fac)
    y <- as.integer(phenotype_fac) - 1L
  }

  if (is.null(tag)) {
    tag <- level_labels
  } else {
    tag <- as.character(tag)
    if (length(tag) != 2) {
      stop("For family = 'binomial', tag must have length 2.")
    }
  }

  z <- table(factor(y, levels = 0:1, labels = tag))

  print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
  print("Perform logistic regression on the given phenotypes:")

  as.numeric(y)
}

prepare_gaussian_response <- function(phenotype, tag = NULL) {
  y <- coerce_numeric_phenotype(phenotype, "gaussian")

  if (length(unique(y)) < 2) {
    stop("For family = 'gaussian', phenotype must contain at least two distinct values.")
  }

  if (!is.null(tag)) {
    tag <- as.character(tag)
    unique_values <- sort(unique(y))
    if (length(tag) == length(unique_values)) {
      z <- table(factor(y, levels = unique_values, labels = tag))
      tmp <- paste(z, tag)
      print(paste0("Current phenotype contains ", paste(tmp[1:(length(z) - 1)], collapse = ", "), ", and ", tmp[length(z)], " samples."))
    } else {
      warning("Ignoring `tag` for family = 'gaussian' because its length does not match the number of unique phenotype values.")
    }
  }

  print(sprintf(
    "Current phenotype summary: min = %s, median = %s, max = %s.",
    formatC(min(y), format = "fg", digits = 6),
    formatC(stats::median(y), format = "fg", digits = 6),
    formatC(max(y), format = "fg", digits = 6)
  ))
  print("Perform linear regression on the given phenotypes:")

  y
}

prepare_cox_response <- function(phenotype) {
  y <- as.matrix(phenotype)
  if (ncol(y) != 2) {
    stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
  }

  print("Perform cox regression on the given clinical outcomes:")
  y
}

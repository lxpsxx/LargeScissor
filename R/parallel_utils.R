normalize_parallel_flag <- function(Mthread) {
  isTRUE(Mthread)
}

normalize_parallel_cores <- function(Mcore) {
  cores <- suppressWarnings(as.integer(Mcore)[1])
  if (is.na(cores) || cores < 1L) {
    return(1L)
  }
  cores
}

parallel_backend_enabled <- function(Mthread, Mcore, ntask = NULL) {
  if (!normalize_parallel_flag(Mthread)) {
    return(FALSE)
  }

  cores <- normalize_parallel_cores(Mcore)
  if (cores < 2L || .Platform$OS.type != "unix") {
    return(FALSE)
  }

  if (!is.null(ntask) && ntask < 2L) {
    return(FALSE)
  }

  TRUE
}

resolve_parallel_cores <- function(Mcore, ntask) {
  cores <- normalize_parallel_cores(Mcore)
  detected <- parallel::detectCores(logical = TRUE)
  if (!is.na(detected)) {
    cores <- min(cores, detected)
  }
  min(cores, ntask)
}

parallel_task_apply <- function(tasks, FUN, Mthread = FALSE, Mcore = 1L) {
  ntask <- length(tasks)
  if (!parallel_backend_enabled(Mthread, Mcore, ntask)) {
    return(lapply(tasks, FUN))
  }

  parallel::mclapply(
    tasks,
    FUN,
    mc.cores = resolve_parallel_cores(Mcore, ntask),
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
}

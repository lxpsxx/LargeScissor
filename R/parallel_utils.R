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

estimate_parallel_payload_bytes <- function(...) {
  objects <- list(...)
  if (length(objects) == 0L) {
    return(NA_real_)
  }

  total <- sum(vapply(objects, function(obj) as.numeric(utils::object.size(obj)), numeric(1)), na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(NA_real_)
  }

  total
}

available_parallel_memory_bytes <- function() {
  if (.Platform$OS.type != "unix" || !file.exists("/proc/meminfo")) {
    return(NA_real_)
  }

  meminfo <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(...) character())
  if (!length(meminfo)) {
    return(NA_real_)
  }

  memavailable <- grep("^MemAvailable:", meminfo, value = TRUE)
  if (!length(memavailable)) {
    return(NA_real_)
  }

  fields <- strsplit(trimws(memavailable[1]), "\\s+")[[1]]
  kb <- suppressWarnings(as.numeric(fields[2]))
  if (!is.finite(kb) || kb <= 0) {
    return(NA_real_)
  }

  kb * 1024
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

cap_parallel_cores_by_payload <- function(cores, payload_bytes) {
  if (is.null(payload_bytes) || !is.finite(payload_bytes) || payload_bytes <= 0) {
    return(cores)
  }

  available <- available_parallel_memory_bytes()
  if (is.finite(available) && available > 0) {
    reserve_bytes <- 1024^3
    usable_bytes <- max(available - reserve_bytes, payload_bytes)
    max_cores <- floor(usable_bytes / max(payload_bytes, 1))
    max_cores <- suppressWarnings(as.integer(max_cores))
    if (is.finite(max_cores) && !is.na(max_cores)) {
      return(min(cores, max(1L, max_cores)))
    }
  }

  if (payload_bytes >= 2 * 1024^3) {
    return(min(cores, 2L))
  }
  if (payload_bytes >= 1024^3) {
    return(min(cores, 4L))
  }

  cores
}

resolve_parallel_cores <- function(Mcore, ntask, payload_bytes = NULL) {
  cores <- normalize_parallel_cores(Mcore)
  detected <- parallel::detectCores(logical = TRUE)
  if (!is.na(detected)) {
    cores <- min(cores, detected)
  }
  cores <- min(cores, ntask)
  cores <- cap_parallel_cores_by_payload(cores, payload_bytes)
  max(1L, min(cores, ntask))
}

parallel_task_apply <- function(tasks, FUN, Mthread = FALSE, Mcore = 1L, payload_bytes = NULL) {
  ntask <- length(tasks)
  if (!parallel_backend_enabled(Mthread, Mcore, ntask)) {
    return(lapply(tasks, FUN))
  }

  cores <- resolve_parallel_cores(Mcore, ntask, payload_bytes = payload_bytes)
  if (cores < 2L) {
    return(lapply(tasks, FUN))
  }

  parallel::mclapply(
    tasks,
    FUN,
    mc.cores = cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE,
    mc.allow.recursive = FALSE
  )
}

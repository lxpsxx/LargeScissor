resolve_dense_block_ncol <- function(nrow_matrix, ncol_matrix, target_bytes = 64 * 1024^2) {
  if (is.null(ncol_matrix) || ncol_matrix < 1L) {
    return(1L)
  }

  bytes_per_col <- max(8, as.double(nrow_matrix) * 8)
  block_ncol <- floor(target_bytes / bytes_per_col)
  block_ncol <- suppressWarnings(as.integer(block_ncol))

  if (is.na(block_ncol) || block_ncol < 1L) {
    return(1L)
  }

  min(block_ncol, as.integer(ncol_matrix))
}

build_scissor_quantile_input <- function(bulk_expr, sc_expr) {
  bulk_mat <- as.matrix(bulk_expr)
  n_gene <- nrow(bulk_mat)
  n_bulk <- ncol(bulk_mat)
  n_cell <- ncol(sc_expr)

  dataset0 <- matrix(
    0,
    nrow = n_gene,
    ncol = n_bulk + n_cell,
    dimnames = list(rownames(bulk_mat), c(colnames(bulk_mat), colnames(sc_expr)))
  )
  storage.mode(dataset0) <- "double"
  dataset0[, seq_len(n_bulk)] <- bulk_mat

  if (n_cell > 0L) {
    if (inherits(sc_expr, "sparseMatrix")) {
      block_ncol <- resolve_dense_block_ncol(n_gene, n_cell)
      for (start in seq.int(1L, n_cell, by = block_ncol)) {
        stop_idx <- min(start + block_ncol - 1L, n_cell)
        src_idx <- start:stop_idx
        dst_idx <- n_bulk + src_idx
        dataset0[, dst_idx] <- as.matrix(sc_expr[, src_idx, drop = FALSE])
      }
    } else {
      dataset0[, n_bulk + seq_len(n_cell)] <- as.matrix(sc_expr)
    }
  }

  dataset0
}

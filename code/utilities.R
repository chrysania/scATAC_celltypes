PEAKS_COL="#E85C0D"
REMO_COL="#48CFCB"
CCRE_COL="#CE6DFF"
RUST_COL="#DEA584"
R_COL="#198CE7"

set.seed(1234)

get_knn_purity <- function(emb, idents, method) {
    knp <- data.frame(
        "Score" = knn_purity(emb, idents),
        "Celltype" = idents
    )
    knp <- knp |>
        dplyr::group_by(Celltype) |>
        dplyr::mutate(mn = mean(Score)) |>
        dplyr::ungroup() |>
        dplyr::select(Celltype, mn) |>
        unique()
    knp$Method <- method
    return(knp)
}

get_mean_sil <- function(emb, idents, method) {
    sil <- data.frame(
        "Score" = get_silhouette(emb, idents),
        "Celltype" = idents
    )
    sil <- sil |>
        dplyr::group_by(Celltype) |>
        dplyr::mutate(mn = mean(Score)) |>
        dplyr::ungroup() |>
        dplyr::select(Celltype, mn) |>
        unique()
    sil$Method <- method
    return(sil)
}

get_metrics <- function(remo_obj, peaks_obj, ccre_obj, 
                        remo_emb, peaks_emb, ccre_emb,
                        tissue_name) {

    # KNN purity
    knp_remo <- get_knn_purity(remo_emb, remo_obj$cluster, "REMO")
    knp_peaks <- get_knn_purity(peaks_emb, peaks_obj$cluster, "Peaks")
    knp_ccre <- get_knn_purity(ccre_emb, ccre_obj$cluster, "cCRE")
    knp <- rbind(knp_remo, knp_peaks, knp_ccre)
    knp$Dataset <- tissue_name
    knp$Method <- factor(knp$Method, levels = c("REMO", "Peaks", "cCRE"))

    # silhouette
    sil_remo <- get_mean_sil(remo_emb, remo_obj$cluster, "REMO")
    sil_peaks <- get_mean_sil(peaks_emb, peaks_obj$cluster, "Peaks")
    sil_ccre <- get_mean_sil(ccre_emb, ccre_obj$cluster, "cCRE")
    sil <- rbind(sil_remo, sil_peaks, sil_ccre)
    sil$Dataset <- tissue_name
    sil$Method <- factor(sil$Method, levels = c("REMO", "Peaks", "cCRE"))
    
    # ARI
    ari <- data.frame(
        "Dataset" = tissue_name,
        "Score" = c(
            mclust::adjustedRandIndex(remo_obj$seurat_clusters, remo_obj$cluster),
            mclust::adjustedRandIndex(peaks_obj$seurat_clusters, peaks_obj$cluster),
            mclust::adjustedRandIndex(ccre_obj$seurat_clusters, ccre_obj$cluster)),
        "Method" = c("REMO", "Peaks", "cCRE"),
        "Metric" = "ARI"
    )

    # NMI
    nmi <- data.frame(
        "Dataset" = tissue_name,
        "Score" = c(
            aricode:::NMI(remo_obj$seurat_clusters, remo_obj$cluster),
            aricode:::NMI(peaks_obj$seurat_clusters, peaks_obj$cluster),
            aricode:::NMI(ccre_obj$seurat_clusters, ccre_obj$cluster)),
        "Method" = c("REMO", "Peaks", "cCRE"),
        "Metric" = "NMI"
    )

    return(list("AN" = rbind(ari, nmi), "KNN" = knp, "SIL" = sil))
}

knn_purity <- function(embeddings, clusters, k = 30) {
  nn <- RANN::nn2(data = embeddings, k = k + 1)$nn.idx[, 2:k]
  nn_purity <- vector(mode = "numeric", length = length(x = clusters))
  for (i in seq_len(length.out = nrow(x = nn))) {
    nn_purity[i] <- sum(clusters[nn[i, ]] == clusters[i]) / k
  }
  return(nn_purity)
}

# annotation_matrix = binary matrix, module x annotations
# feature_list = list of module names to test
find_annotations <- function(
    annotation_matrix,
    feature_list,
    p.adjust.method="BH",
    background = NULL
) {
    # subset features matrix & background matrix
    features_idx <- which(rownames(annotation_matrix) %in% feature_list)

    features_matrix <- annotation_matrix[features_idx,]
    # background_matrix <- annotation_matrix[-features_idx,] # background does not include feature set
    if (is.null(x = background)) {
        background_matrix <- annotation_matrix
    } else {
        background_matrix <- annotation_matrix[background, ]
    }

    # calculate foldchange & p-value
    features.counts <- colSums(x = features_matrix)
    background.counts <- colSums(x = background_matrix)

    percent.observed <- features.counts / nrow(features_matrix) * 100
    percent.background <- background.counts / nrow(background_matrix) * 100

    fold.enrichment <- percent.observed / percent.background

    p.list <- vector(mode = "numeric")
    for (i in seq_along(along.with = features.counts)) {
        p.list[[i]] <- phyper(
            q = features.counts[[i]] - 1,
            m = background.counts[[i]],
            n = nrow(x = background_matrix) - background.counts[[i]],
            k = nrow(features_matrix),
            lower.tail = FALSE
        )
    }

    # make result dataframe
    results <- data.frame(
        motif = names(x = features.counts),
        observed = features.counts,
        background = background.counts,
        percent.observed = percent.observed,
        percent.background = percent.background,
        fold.enrichment = fold.enrichment,
        pvalue = p.list
      )
    if (p.adjust.method == 'qvalue') {
        results$p.adjust = qvalue::qvalue(p.list)$qvalue
    } else {
        results$p.adjust = p.adjust(p.list, method = p.adjust.method)
    }
    # results <- results[order(results$fold.enrichment, decreasing=TRUE), ]
    results <- results[order(results$p.adjust, results$fold.enrichment, decreasing=c(FALSE, TRUE)), ]
    return(results)
}

top_experiments <- function(
    experiment_metadata,
    encode_normalized,
    remo,
    all_ccre,
    n=5
) {

    # get CCREs for REMO module
    idents <- ccre[all_ccre$REMO == remo, "id"]

    # ensure experiment metadata and encode data are ordered the same
    exp.keep <- intersect(rownames(experiment_metadata), rownames(encode_normalized))
    encode_normalized <- encode_normalized[exp.keep, ]
    experiment_metadata <- experiment_metadata[exp.keep, ]

    # find experiments targetting the same thing
    all_targets <- experiment_metadata$Biosample.term.name
    rank_sum <- rep(0, nrow(experiment_metadata))
    for (i in idents) {

        ranks <- rank(encode_normalized[, i]) # ranks lowest to highest
        for (j in unique(all_targets)) {
            target_ranks <- ranks[all_targets == j]
            ranks[names(target_ranks)] <- mean(target_ranks)
        }
        rank_sum <- rank_sum + ranks
    }
    
    experiment_names <- rownames(encode_normalized)[head(order(rank_sum, decreasing = TRUE), n)]

    # get biosample target for experiments
    targets <- experiment_metadata[experiment_names, "Biosample.term.name"]
    return(unique(targets))
}

top_frequency_experiments <- function(
    experiment_metadata,
    encode_normalized,
    remo,
    all_ccre,
    n=5
) {

    # get CCREs for REMO module
    idents <- ccre[all_ccre$REMO == remo, "id"]

    # ensure experiment metadata and encode data are ordered the same
    exp.keep <- intersect(rownames(experiment_metadata), rownames(encode_normalized))
    encode_normalized <- encode_normalized[exp.keep, ]
    experiment_metadata <- experiment_metadata[exp.keep, ]

    # find experiments targetting the same thing
    all_targets <- experiment_metadata$Biosample.term.name

    experiments <- vector(mode = "character")

    for (i in idents) {
        # get top experiments and append to vector
        top_experiments <- order(encode_normalized[, i], decreasing = TRUE)[1:10]

        # get biosample target for experiments
        targets <- experiment_metadata[top_experiments, "Biosample.term.name"]
        experiments <- c(experiments, targets)
    }

    # frequency of each experiment in vector
    exp_freq <- sort(table(experiments), decreasing = TRUE)
    return(head(names(exp_freq), n = n))
}

bin_variance <- function(mat, bins=1000, n=40, min.cutoff='q80', verbose = FALSE) {
    df <- data.frame(
        variance = sparseMatrixStats::rowVars(mat),
        mean = rowMeans(mat),
        total_counts = rowSums(mat),
        module = rownames(mat)
    )
    if (is.character(min.cutoff)) {
        percentile.use <- as.numeric(x = sub(pattern = "q", replacement = "", x = as.character(x = min.cutoff)))/100
        count.thresh <- quantile(df$total_counts, probs = percentile.use)[[1]]
    } else {
        count.thresh <- min.cutoff
    }
  df <- df[df$total_counts >= count.thresh, ]
  if (verbose) {
      message("Retained ", nrow(df), " features after count filtering")
  }
  if (nrow(df) == 0) {
    stop("No modules remain after filtering by min_counts.")
  }
  breaks <- unique(quantile(df$mean, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
  df$bin <- findInterval(df$mean, vec = breaks, rightmost.closed = TRUE)
  df <- na.omit(df)
  selected_modules <- vector(mode = 'character')
  for (b in unique(df$bin)) {
    bin_data <- df[df$bin == b, ]
    bin_data <- bin_data[order(-bin_data$variance), ]
    selected_modules <- c(selected_modules, head(bin_data$module, n))
  }
  df$variable <- df$module %in% selected_modules
  if (verbose) {
      message("Selected ", sum(df$variable), " features")
  }
  return(df)
}

lm_variance <- function(
    mat,
    n = 20000,
    min.cutoff = 100,
    weight_mean = 0.5,
    method = "loess",
    span = 0.1,
    bins = 1000,
    sample_per_bin = 50,
    verbose = FALSE
) {
    set.seed(1234)
    rs <- rowSums(mat)
    if (is.character(min.cutoff)) {
        percentile.use <- as.numeric(sub(pattern = "q", replacement = "", x = as.character(min.cutoff))) / 100
        count.thresh <- quantile(rs, probs = percentile.use, na.rm = TRUE)[[1]]
    } else {
        count.thresh <- min.cutoff
    }
    mat <- mat[rs >= count.thresh, ]
    if (verbose) {
        message("Retained ", nrow(mat), " features after count filtering")
    }
    if (nrow(mat) == 0) {
        stop("No modules remain after filtering by min_counts.")
    }
    df <- data.frame(
        variance = sparseMatrixStats::rowVars(mat),
        mean = rowMeans(mat),
        total_counts = rs[rs >= count.thresh],
        module = rownames(mat)
    )
    df$log_mean <- log1p(df$mean)
    df$log_variance <- log1p(df$variance)
    breaks <- seq(min(df$log_mean, na.rm = TRUE), max(df$log_mean, na.rm = TRUE), length.out = bins + 1)
    df$bin <- findInterval(df$log_mean, breaks, rightmost.closed = TRUE)
    sampled_df <- do.call(rbind, lapply(split(df, df$bin), function(subset) {
        if (nrow(subset) > sample_per_bin) {
            subset <- subset[sample(nrow(subset), sample_per_bin), ]
        }
        return(subset)
    }))
    if (method == 'loess') {
        loess_fit <- loess(log_variance ~ log_mean, data = sampled_df, span = span)
        predicted <- predict(loess_fit, newdata = df$log_mean)
        df$residuals <- df$log_variance - predicted
    } else if (method == 'lm') {
        lm_fit <- lm(log_variance ~ log_mean, data = df)
        df$residuals <- resid(lm_fit)
    }
    df$residual_rank <- rank(-df$residuals, ties.method = "average")
    df$mean_rank <- rank(-df$log_mean, ties.method = "average")
    df$combined_rank <- (weight_mean * df$mean_rank) + ((1 - weight_mean) * df$residual_rank)
    selected_modules <- df[order(df$combined_rank), "module"][1:n]
    df$variable <- df$module %in% selected_modules
    if (verbose) {
        message("Selected ", sum(df$variable), " features")
    }
    return(df)
}

select_features <- function(x, min_var = 0.1, max_nz = 0.95, exclude = NULL) {
    x <- NormalizeData(x, scale.factor = 5000)
    row_sd <- sqrt(sparseMatrixStats::rowVars(x))
    high_var <- rownames(x)[row_sd > min_var]
    fraction_nz <- rowSums(x>0) / ncol(x)
    high_nz <- names(fraction_nz)[fraction_nz > max_nz]
    feat_use <- setdiff(high_var, c(high_nz, exclude))
    return(feat_use)
}

cumulative_top_features <- function(mat, threshold = 0.8, plot = FALSE) {
    # rank features
    # return the features that contain the specified fraction of total counts in the matrix
    row_contr <- sort(rowSums(mat) / sum(mat), decreasing = TRUE)
    cumulative_counts <- cumsum(row_contr)
    idx <- which(cumulative_counts > threshold)[1]
    if (plot) {
        plot(cumulative_counts)
        abline(v=idx)
    }
    feat.use <- names(row_contr[1:idx])
    return(feat.use)
}

# Function to extract runtime and max memory from the profiling output
extract_profile_info <- function(profile_file) {
  # Read the profiling file
  lines <- readLines(profile_file)
  
  # Initialize variables to store the results
  elapsed_time <- NA
  max_memory <- NA
  
  # Loop through each line and extract relevant information
  for (line in lines) {
    if (grepl("Command exited with non-zero status 1", line)) {
      df <- data.frame("Time" = NA, "Memory" = NA)
      return(df)
    }
    if (grepl("Elapsed \\(wall clock\\) time", line)) {
      # Extract the elapsed time
      elapsed_time <- convert_time_to_seconds(sub(".*: ", "", line))
    } else if (grepl("Maximum resident set size", line)) {
      # Extract the max memory usage
      max_memory <- as.numeric(sub(".*: ", "", line))
    }
  }

  df <- data.frame("Time" = elapsed_time, "Memory" = max_memory/1024)
  return(df)
}

convert_time_to_seconds <- function(time_str) {
  # Check if the time format is h:mm:ss or m:ss
  if (grepl(":", time_str)) {
    time_parts <- as.numeric(unlist(strsplit(time_str, ":")))
    if (length(time_parts) == 3) {
      # Format is h:mm:ss
      return(time_parts[1] * 3600 + time_parts[2] * 60 + time_parts[3])
    } else if (length(time_parts) == 2) {
      # Format is m:ss
      return(time_parts[1] * 60 + time_parts[2])
    }
  }
  return(as.numeric(time_str))
}


revcomp <- function(x) {
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

ReadCounts <- function(dir) {
  counts <- readMM(paste0(dir, "/matrix.mtx.gz"))
  if (file.exists(paste0(dir, "/barcodes.tsv.gz"))) {
    colnames(counts) <- readLines(paste0(dir, "/barcodes.tsv.gz"))
  } else {
    colnames(counts) <- readLines(paste0(dir, "/barcodes.tsv"))
  }
  rownames(counts) <- readLines(paste0(dir, "/features.tsv.gz"))
  counts <- as(counts, "CsparseMatrix")
  return(counts)
}

downsampleCounts <- function(x, prob) {

    if (!inherits(x, "sparseMatrix")) {
        stop("Input matrix must be a sparse matrix.")
    }

    # convert to triplet form
    x <- as(x, "TsparseMatrix")
    orig_dims <- dim(x)
    orig_dimnames <- dimnames(x)
    ival <- x@i
    jval <- x@j
    counts <- x@x

    # downsample counts
    ds <- rbinom(n = length(counts), size = counts, prob = prob)

    # remove zero positions
    zero_positions <- ds == 0
    i_nz <- ival[!zero_positions]
    j_nz <- jval[!zero_positions]
    ds <- ds[!zero_positions]
    
    # reconstruct sparse matrix
    x <- sparseMatrix(i = i_nz+1, j = j_nz+1, x = ds, dims = orig_dims)
    dimnames(x) <- orig_dimnames
    return(x)
}

run_pca_sparse <- function(x, dim = 50, weight.by.var=TRUE, assay="RNA") {
    d_rowmeans <- rowMeans(x)
    d_sd <- sqrt(sparseMatrixStats::rowVars(x))
    nz_var <- d_sd > 0
    message("Retaining ", sum(nz_var), " features with non-zero variance")
    x <- x[nz_var, ]
    d_rowmeans <- d_rowmeans[nz_var]
    d_sd <- d_sd[nz_var]
    pcs <- irlba::irlba(A = t(x), scale = d_sd, center = d_rowmeans, nv = dim)
    if (weight.by.var) {
        emb <- pcs$u %*% diag(pcs$d)
    }
    loadings <- pcs$v
    rownames(loadings) <- rownames(x)
    colnames(loadings) <- paste0("PC", 1:dim)
    rownames(emb) <- colnames(x)
    colnames(emb) <- colnames(loadings)
    sdev <- pcs$d/sqrt(max(1, ncol(x) - 1))
    dr <- Seurat::CreateDimReducObject(
        embeddings = emb,
        loadings = loadings,
        assay = assay,
        stdev = sdev,
        key = "PC_"
    )
    return(dr)
}

process_obj <- function(
    counts,
    dims = 1:30,
    var_features = NULL,
    exclude = NULL,
    nfeatures = 2000,
    scale.factor=10000,
    normalization="LogNorm",
    resolution = 0.8
) {
  obj <- CreateSeuratObject(counts = counts)
  if (normalization == "LogNorm") {
      obj <- NormalizeData(obj, scale.factor = scale.factor)
  } else if (normalization == "TMM") {
      dge <- edgeR::DGEList(counts = counts)
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      norm_data <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
      LayerData(obj, assay = "RNA", layer = "data") <- norm_data
  }
  if (is.null(var_features)) {
      obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
      vf <- VariableFeatures(obj)
      VariableFeatures(obj) <- setdiff(vf, exclude)
  } else {
      VariableFeatures(obj) <- var_features
  }
  norm_data <- LayerData(obj[['RNA']], layer = 'data')[VariableFeatures(obj),]
  obj[['pca']] <- run_pca_sparse(norm_data, dim=max(50, max(dims)))
  obj <- RunUMAP(obj, reduction = 'pca', dims = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims)
  obj <- FindClusters(obj, algorithm=3, resolution=resolution)
  return(obj)
}

process_atac_obj <- function(counts, dims = 2:30, features = NULL, resolution=0.8) {
  obj <- CreateSeuratObject(counts = counts, assay = "ATAC")
  obj <- RunTFIDF(obj)
  if (is.null(x = features)) {
      obj <- FindTopFeatures(obj)
  } else {
      VariableFeatures(obj) <- features
  }
  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = 'lsi', dims = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = dims)
  obj <- FindClusters(obj, algorithm=3, resolution=resolution)
  return(obj)
}

process_remo <- function (counts, dims = 2:30, features = 1:15000, resolution = 0.8, remo_freq) {
    # normalize by REMO size factor
    counts <- counts / log1p(as.numeric(remo_freq[rownames(counts)]))
    rs <- rowSums(counts)
    feature_rank <- order(rs, decreasing = TRUE)
    counts <- counts[names(rs)[feature_rank[features]], ]
    obj <- CreateSeuratObject(counts = counts, assay = "ATAC")
    obj <- RunTFIDF(obj)
    VariableFeatures(obj) <- rownames(counts)
    obj <- RunSVD(obj)
    obj <- RunUMAP(obj, reduction = "lsi", dims = dims, verbose = FALSE)
    obj <- FindNeighbors(obj, reduction = "lsi", dims = dims)
    obj <- FindClusters(obj, algorithm = 3, resolution = resolution)
    return(obj)
}

get_silhouette <- function(x, idents) {
    dist.matrix <- dist(x)
    sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = idents)), dist = dist.matrix)
    return(sil[, 3])
}

get_procrustes <- function(x, y) {
    procrust <- vegan::procrustes(x, y)
    procrustes_score <- summary(procrust)$rmse
    return(procrustes_score)
}

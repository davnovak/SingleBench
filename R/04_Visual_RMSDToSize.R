
PlotRMSDToSize <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = 1
) {

  if (is.null(idx.run)) {
    stop('idx.run must be specified (or left as 1)')
  }

  .PlotClustering.ValidityChecks(environment())

  rmsd_real <- GetRMSDPerPopulation(benchmark, idx.subpipeline, idx.n_param, idx.run, match_type = 'real')
  rmsd_fixed_cluster <- GetRMSDPerPopulation(benchmark, idx.subpipeline, idx.n_param, idx.run, match_type = 'fixed_cluster')

  sizes_real <- GetPopulationSizes(benchmark)
  sizes_fixed_cluster <- GetMatchedClusterSizes(benchmark, idx.subpipeline, idx.n_param, idx.run, match_type = 'fixed_cluster')

  df.real <-
    data.frame(
      Population = names(rmsd_real),
      ClassType = 'Labelled Population',
      ClusterIndex = NA,
      RMSD = unlist(rmsd_real),
      Size = unlist(sizes_real),
      row.names = NULL
    )
  
  pops <- unlist(
    purrr::map(names(rmsd_fixed_cluster),
               function(n) rep(n, times = sum(!is.na(x <- rmsd_fixed_cluster[[n]])))))
  rmsd_fc <- unlist(rmsd_fixed_cluster)
  rmsd_fc <- rmsd_fc[!is.na(rmsd_fc)]
  sizes_fc <- unlist(sizes_fixed_cluster)
  sizes_fc <- sizes_fc[sizes_fc != 0]
  cl_idcs <- unlist(purrr::map(rmsd_fixed_cluster, names))
  cl_idcs <- cl_idcs[cl_idcs != 'NA']
  cl_idcs <- as.numeric(cl_idcs)
  
  df.fixed_cluster <-
    data.frame(
      Population = pops,
      ClassType = 'Cluster Matched to Population',
      ClusterIndex = cl_idcs,
      RMSD = rmsd_fc,
      Size = sizes_fc,
      row.names = NULL
    )
  data <- rbind(df.real, df.fixed_cluster)
  rownames(data) <- NULL

  nn <- GetNParameterIterationName(benchmark, idx.subpipeline, idx.n_param)
  
  ggplot(data, aes(x = Size, y = RMSD, col = ClassType)) +
    geom_point(size = 4) + geom_text(aes(label = ClusterIndex), size = 2.7, col = 'black') +
    facet_wrap(~ Population) + theme_grey() +
    ggtitle('RMSD/size ratios of labelled populations and matched clusters', nn)
}
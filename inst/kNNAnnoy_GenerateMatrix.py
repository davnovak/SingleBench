
def knn_annoy_python(data, nrow, ncol, tmpfile, ntrees=50, k=100, distf="euclidean"):
  import annoy

  Ann = annoy.AnnoyIndex(ncol, distf)
  Ann.build(ntrees)
  Ann.save(tmpfile)
  
  kNN = AnnoyKnnMatrix(index=Ann, shape=(nrow,ncol), index_path=tmpfile, k=k, include_distances=True, verbose=True)
  kNN.build(data, tmpfile, k=k, metric=distf)
  res = kNN.get_neighbour_indices()
  return res


sister_pairs = function(phylo_dist){
  sister_taxa_list = vector('list', nrow(phylo_dist))
  all_spp_names = row.names(phylo_dist)
  names(sister_taxa_list) = all_spp_names
  for(focal_sp in all_spp_names){
    focal_sp_dists= round(phylo_distance[focal_sp,],5)
    min_phylo_dist = round( min(phylo_distance[focal_sp,][phylo_distance[focal_sp,] != 0]), 5)
    sister_index = which(focal_sp_dists == min_phylo_dist)
    sister_taxa_list[[focal_sp]] = names(focal_sp_dists)[sister_index]
  }
  return(sister_taxa_list)
}
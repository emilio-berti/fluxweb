# script to find all cycles in the metaweb and species involved ---------
# it's not standalone, but you have to run flux-per-pixel until the parallel
# loop to run this.
#
# It's mostly for development purpose in order to understand why `fluxing()` 
# failed for large squared areas and how to account for it.
d <- read_delim(file.path(path, "TetraEU_pairwise_interactions.csv"),
                delim = ";",
                show_col_types = FALSE) %>% 
  filter(grepl("adults", sourceLifestageName, ignore.case = TRUE) | 
           grepl("adults", targetLifestageName, ignore.case = TRUE)) %>% 
  filter(sourceTaxonName != targetTaxonName) %>% 
  dplyr::select(conID = sourceTaxonId, con = sourceTaxonName,
                res = targetTaxonName, resID = targetTaxonId) %>% 
  distinct_all()

library(future)
plan(multisession(workers = 6))
d$isCycle <- FALSE
for (i in seq_len(nrow(d))) {
  cons <- d[i, ][["con"]]
  ress <- d[i, ][["res"]]
  if( (d %>% filter(con == ress, res == cons) %>% nrow()) == 1) d[i, "isCycle"] <- TRUE
}

d %>% 
  filter(isCycle)

d %>% 
  transmute(res = sourceTaxonName, con = targetTaxonName, what = "pass I") %>% 
  distinct_all() %>% 
  left_join(
    d %>% 
      transmute(res = targetTaxonName, con = sourceTaxonName, what = "pass II"),
    by = c("res" = "con", "con" = "res")
  ) %>% 
  distinct_all()
  

tmp_metaweb <- create_metaweb() #load all necessary files and create metaweb
g.metaweb <- tmp_metaweb[[1]]
g.metaweb <- simplify(g.metaweb)
dist <- distances(g.metaweb, mode = "out")
dist[is.infinite(dist)] <- 0
dist[!is.infinite(dist)] <- 1
dist <- dist + t(dist)

cycles <- which(dist > 1, arr.ind = TRUE) %>% 
  as_tibble() %>% 
  filter(row != col) 

cycles <- d %>% 
  left_join(tibble(sourceTaxonId = rownames(dist)[cycles[["row"]]],
                   targetTaxonId = colnames(dist)[cycles[["col"]]],
                   isCycle = TRUE),
            by = c("sourceTaxonId", "targetTaxonId")) %>% 
  filter(isCycle) %>% 
  dplyr::select(sourceTaxonId, sourceTaxonName,
                targetTaxonName, targetTaxonId) %>% 
  distinct_all() %>% 
  arrange(sourceTaxonName, targetTaxonName)

g <- graph_from_data_frame(transmute(cycles, from = sourceTaxonId, to = targetTaxonId))
lay <- layout_components(g)
g$layout <- lay
plot(g, edge.width = 1, edge.arrow.width = .51, edge.arrow.size = .3, curved = 1)

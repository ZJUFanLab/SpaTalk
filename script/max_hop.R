library(SpaTalk)
# human
lr <- lrpairs[lrpairs$species == "Human", ]
ggi_tf <- pathways[pathways$species == "Human", c("src", "dest", "src_tf", "dest_tf")]
ggi_tf <- unique(ggi_tf)
# generate ggi_res
receptor_name <- unique(lr$receptor)
n_hop <- rep(0, length(receptor_name))
for (i in 1:length(receptor_name)) {
    ggi_tf1 <- ggi_tf[ggi_tf$src == receptor_name[i], ]
    if (nrow(ggi_tf1) > 0) {
        ggi_tf1_yes <- ggi_tf1[ggi_tf1$dest_tf == "YES", ]
        n <- 0
        while (nrow(ggi_tf1_yes) == 0) {
            ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest, ]
            if (nrow(ggi_tf1) == 0) {
                break
            }
            ggi_tf1_yes <- ggi_tf1[ggi_tf1$dest_tf == "YES", ]
            n <- n + 1
        }
        n_hop[i] <- n
    }
}
max_hop <- max(n_hop)
### human max_hop = 3

# mouse
lr <- lrpairs[lrpairs$species == "Mouse", ]
ggi_tf <- pathways[pathways$species == "Mouse", c("src", "dest", "src_tf", "dest_tf")]
ggi_tf <- unique(ggi_tf)
# generate ggi_res
receptor_name <- unique(lr$receptor)
n_hop <- rep(0, length(receptor_name))
for (i in 1:length(receptor_name)) {
    ggi_tf1 <- ggi_tf[ggi_tf$src == receptor_name[i], ]
    if (nrow(ggi_tf1) > 0) {
        ggi_tf1_yes <- ggi_tf1[ggi_tf1$dest_tf == "YES", ]
        n <- 0
        while (nrow(ggi_tf1_yes) == 0) {
            ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest, ]
            if (nrow(ggi_tf1) == 0) {
                break
            }
            ggi_tf1_yes <- ggi_tf1[ggi_tf1$dest_tf == "YES", ]
            n <- n + 1
        }
        n_hop[i] <- n
    }
}
max_hop <- max(n_hop)
### mouse max_hop = 4

source("man_functions.R")
source("man_objects.R")


#### RUIZ #####
df_ruiz <- data.frame(
  x1 = c(1,1,0,1,0,0,1,0,1,0),
  x2 = c(0,1,0,1,0,1,1,0,1,0),
  x3 = c(0,1,1,1,0,1,1,0,1,0)
)
df_ruiz[] <- df_ruiz + 1

ruiz_data <- BNDataset(
  data = df_ruiz,
  discreteness = c("d","d","d"),
  variables = c("x1","x2","x3"),
  node.sizes = c(2,2,2)
)

layering_mi_ruiz <- node_order

ruiz_net_bs <- learn.network(
  ruiz_data,
  algo = "mmhc",
  scoring.func = "BDeu",
  layering = layering_mi_ruiz,
  max.parents = num.variables(ruiz_data)-1
)


#### ASIA #####
asia_net_bs <- learn.network(
  asia_bs,
  algo = "mmhc",
  scoring.func = "BDeu",
  layering = asia_order,
  max.parents = num.variables(asia_bs)-1
)


#### CHILD ####
child_net_bs <- learn.network(
  child_bs,
  algo = "mmhc",
  scoring.func = "BDeu",
  layering = child_order,
  max.parents = num.variables(child_bs)-1,
  use.imputed.data = TRUE
)


#### SACHS ####
sachs_vars   <- colnames(sachs)
sachs_nodes   <- sapply(sachs, function(col) length(unique(col)))
sachs_disc_flags   <- rep("d", length(sachs_vars))
sachs_data <- BNDataset(
    data=as.data.frame(sachs),
    discreteness=sachs_disc_flags,
    variables=sachs_vars,
    node.sizes=sachs_nodes
)

sachs_net_bs <- learn.network(
  sachs_data,
  algo         = "mmhc",
  scoring.func = "BDeu",
  layering = sachs_order
)


#### UTILITIES ####
bnstruct_to_bnlearn <- function(bnstruct_net) {
  vars <- variables(bnstruct_net)
  adj  <- dag(bnstruct_net)
  bn   <- empty.graph(nodes=vars)
  for(i in seq_along(vars)) for(j in seq_along(vars))
    if(adj[i,j]==1) 
      bn <- set.arc(bn, from=vars[i], to=vars[j])
  bn
}

# Convert a bnlearn 'bn' object back to a bnstruct 'BN' object using the original BNDataset
# bnlearn_to_bnstruct <- function(bn_learn, bndataset) {
#   # Create an empty BNstruct network with the same dataset metadata
#   bnstruct_net <- BN(dataset = bndataset)
#   # Set its adjacency matrix from the bnlearn object
#   dag(bnstruct_net) <- amat(bn_learn)
#   bnstruct_net
# }



# #### COMPARISON METRICS ####
###### PREPARE DATASETS ######
ruiz_df_fac  <- as.data.frame(lapply(df_ruiz, function(col) as.factor(col)))
asia_df_fac  <- as.data.frame(lapply(asia, function(col) as.factor(col)))
child_df_fac <- as.data.frame(lapply(child, function(col) as.factor(col)))
sachs_df_fac <- as.data.frame(lapply(sachs, function(col) as.factor(col)))
colnames(asia_df_fac) <- gsub("\\.", "-", colnames(asia_df_fac))

###### CONVERT 'BN' OBJECTS TO 'bn' ######
ruiz_bs_bn <- bnstruct_to_bnlearn(ruiz_net_bs)
asia_bs_bn <- bnstruct_to_bnlearn(asia_net_bs)
child_bs_bn <- bnstruct_to_bnlearn(child_net_bs)
sachs_bs_bn <- bnstruct_to_bnlearn(sachs_net_bs)


###### METRICS ######
##### COMPARE STATS #####
compare_stats <- function(net_true, net_learn) {
    # Extract adjacency matrices from either bnstruct::BN or bnlearn::bn objects
    true_mat   <- if (inherits(net_true, "BN")) dag(net_true)
                else if (inherits(net_true, "bn")) amat(net_true)
                else stop("Unsupported network class for reference")
    learn_mat <- if (inherits(net_learn, "BN")) dag(net_learn)
                else if (inherits(net_learn, "bn")) amat(net_learn)
                else stop("Unsupported network class for learned")

    # Ensure learned matrix rows/cols match the true matrix ordering
    common <- intersect(rownames(true_mat), rownames(learn_mat))
    learn_mat <- learn_mat[common, common]
    true_mat  <- true_mat[common, common]

    # Build character vectors of directed edges "from->to"
    true_edges <- apply(which(true_mat == 1, arr.ind = TRUE), 1,
                        function(rc) paste(rownames(true_mat)[rc[1]],
                                           colnames(true_mat)[rc[2]],
                                           sep = "->"))
    learn_edges <- apply(which(learn_mat == 1, arr.ind = TRUE), 1,
                         function(rc) paste(rownames(learn_mat)[rc[1]],
                                            colnames(learn_mat)[rc[2]],
                                            sep = "->"))
    # Identify missing and extra edges
    missing_edges <- setdiff(true_edges, learn_edges)
    extra_edges   <- setdiff(learn_edges, true_edges)

    # Compute true positives, false positives, false negatives
    TP <- sum(true_mat == 1 & learn_mat == 1)
    FP <- sum(true_mat == 0 & learn_mat == 1)
    FN <- sum(true_mat == 1 & learn_mat == 0)

    # Compute precision, recall, F1 (guard against zero denominators)
    Precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
    Recall    <- if ((TP + FN) == 0) NA else TP / (TP + FN)
    F1        <- if (is.na(Precision) || is.na(Recall) || (Precision + Recall) == 0) {
                    NA
                } else {
                    2 * Precision * Recall / (Precision + Recall)
                }

    data.frame(
        TP = TP, 
        FP = FP, 
        FN = FN,
        #Missing = I(list(missing_edges)),
        #Extra   = I(list(extra_edges)),
        Precision = Precision,
        Recall = Recall,
        F1 = F1)
}

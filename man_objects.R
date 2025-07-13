library(bnstruct)

############ RUIZ ############
dummy <- data.frame(
  x1 = c(1,1,0,1,0,0,1,0,1,0),
  x2 = c(0,1,0,1,0,1,1,0,1,0),
  x3 = c(0,1,1,1,0,1,1,0,1,0)
)
node_order <- get_node_ordering_for_K2(dummy, alpha = 0.7)
ruiz_net <- k2(2, dummy, ordering = node_order )


############ ASIA ############
asia_bs <- asia()
asia <- as.data.frame( raw.data(asia_bs) )
colnames(asia) <- variables(asia_bs)
asia_order <- get_node_ordering_for_K2(asia, alpha = 0.9)
asia_net <- k2(ncol(asia) - 1, asia, asia_order)

############ CHILD ############
child_bs <- bnstruct::child()
child_bs <- bnstruct::impute(child_bs,k.impute = 10)
child <- as.data.frame( bnstruct::imputed.data(child_bs) )
colnames(child) <- variables(child_bs)
child_order <- get_node_ordering_for_K2(child, alpha = 0.9)
child_net <- k2(ncol(child) - 1, child, child_order)


############ SACHS ############

sachs <- read_csv('https://www.ccd.pitt.edu/wiki/images/SACHS10k.csv')
sachs_vars <- c("PKC", "Plcg", "PKA", "PIP3", "Raf", "Jnk", "P38", "PIP2", "Mek", "Erk", "Akt")
sachs <- sachs |> select(PKC, Plcg, PKA, PIP3, Raf, Jnk, P38, PIP2, Mek, Erk, Akt)
sachs_order <- get_node_ordering_for_K2(sachs, alpha = 0.9)
sachs_net <- k2(ncol(sachs) - 1, sachs, sachs_order)

## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(warning = FALSE,
                      message = TRUE)

library(yulab.utils)
pload(igraph)
pload(ggplot2)
pload(ggtangle)
pload(aplot)
pload(ggfun)

## -----------------------------------------------------------------------------
library(yulab.utils)

pload(igraph)
pload(ggplot2)
pload(ggtangle)
pload(aplot)
pload(ggfun)

net <- erdos.renyi.game(10, .5)

V(net)$name = letters[1:10]
E(net)$weight = abs(rnorm(length(E(net))))
E(net)$type = sample(LETTERS[1:3], length(E(net)), replace=TRUE)

p <- ggplot(net) + geom_edge()
p

## ----fig.width=15, fig.height=10----------------------------------------------
p1 <- p + geom_point(size=6, color='steelblue')
p2 <- p + geom_point(aes(color = label %in% letters[1:5]), size=6)
p3 <- p + geom_point(aes(size = igraph::degree(net)), color='steelblue') 
p4 <- p + geom_point(aes(shape = label %in% letters[1:5]), color='steelblue', size=8) 

plot_list(p1, p2, p3, p4, ncol=2)

## ----fig.width=12, fig.height=5-----------------------------------------------
p5 <- p + geom_label(aes(label=label, color=label %in% letters[1:5]), size=5)

p6 <- p + geom_label(aes(label=label), size=5, data=ggtree::td_filter(label %in% letters[1:5]))

plot_list(p5, p6)

## ----fig.width=18, fig.height=6.5---------------------------------------------
p7 <- ggplot(net) + geom_edge(aes(linewidth=weight)) 
p8 <- ggplot(net) + geom_edge(aes(color=type)) 
p9 <- ggplot(net) + geom_edge(aes(linetype=type)) 

plot_list(p7, p8, p9, ncol=3)

## ----fig.width=8.8, fig.height=6.5--------------------------------------------
set.seed(123)
expr <- abs(rnorm(10))
d <- data.frame(label=V(net)$name, expression=expr)
p %<+% d + geom_point(aes(color = expression), size=6) +
    scale_color_viridis_c()

## ----fig.width=5.8, fig.height=8----------------------------------------------
flow_info <- data.frame(from = LETTERS[c(1,2,3,3,4,5,6)],
                        to = LETTERS[c(5,5,5,6,7,6,7)])

dd <- data.frame(
    label = LETTERS[1:7],
    v1 = abs(rnorm(7)),
    v2 = abs(rnorm(7)),
    v3 = abs(rnorm(7))
)

g <- igraph::graph_from_data_frame(flow_info)

p <- ggplot(g)  + geom_edge()

pload(scatterpie)

p %<+% dd + 
    geom_scatterpie(cols = c("v1", "v2", "v3")) +
    geom_text(aes(label=label), nudge_y = .2) + 
    coord_fixed()


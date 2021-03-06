library(HiveR)
library(grid)

network.dt <- pair.dt1[, c("mark1", "mark2", "interaction")]

#--- NODE
node1 <- data.frame(id = 1:10,
  lab = c("k4me1", "k4me3", "k9ac", "k9me3", "k27ac", "k27me3", "k36me3", "ctcf", "pol2", "p300"),
  radius = c(4, 4, 9, 9, 27, 27, 36, 4, 8, 12),
  axis = c(1, 1, 2, 1, 2, 1, 1, 3, 3, 3),
  size = c(1, 3, 1, 3, 1, 3, 3, 1, 1, 1))
node1$lab <- as.character(node1$lab)
node1$axis <- as.integer(node1$axis)
node1$color = "#ffffff"

#--- EDGE
edge1 <- network.dt
edge1$id1 <- node1$id[match(edge1$mark1, node1$lab)]
edge1$id2 <- node1$id[match(edge1$mark2, node1$lab)]
edge1$weight <- abs(edge1$interaction)
edge1$color <- sign(edge1$interaction)
edge1$color <- gsub("^1", "darkred", as.character(edge1$color))
edge1$color <- gsub("^-1", "darkgreen", as.character(edge1$color))

edge1 <- edge1[1:10, ]
edge1 <- edge1[-c(1, 5, 9), ]

hpd = list()
hpd$nodes = node1
hpd$edges = edge1
hpd$type = "2D"
hpd$desc = "1st Hive Plot"
hpd$axLabs = c("ME","AC","TF")
class(hpd) = "HivePlotData"

# Check data correctly formatted
chkHPD(hpd, confirm = TRUE)

# plot hive!
pdf('hive.pdf', width=8, height=8)
plotHive(hpd, axLabs = hpd$axLabs, ch = 0.1)
dev.off()
browseURL('hive.pdf')

test <- ranHiveData(nx = 3)
sumHPD(test)
plotHive(test, ch = 10, bkgnd = "white", axLabs = c("TF", "ME", "AC"), axLab.gpar = gpar(col = "pink", fontsize = 14, lwd = 2))
grid.text("males", x = 0, y = 2.3, default.units = "native")
grid.text("females", x = 0, y = -2.3, default.units = "native")
grid.text("Pairing of Eye Color with Hair Color", x = 0, y = 3.75,
          default.units = "native", gp = gpar(fontsize = 18))
grid.text("A test of plotHive annotation options", x = 0, y = 3.25,
          default.units = "native", gp = gpar(fontsize = 12))
grid.text("Images from Wikipedia Commons", x = 0, y = -3.5,
          default.units = "native", gp = gpar(fontsize = 9))

require(plyr)
require(colorspace)
require(classInt)

d = ggplot2::diamonds
d = d[,c(1:4,7)]
head(d); dim(d)

# separate carat-size data into equal interval groups
brks = classIntervals(d$carat, n=11, style="quantile")$brks[1:11] # also try 'equal' style
d$carat = findInterval(d$carat, brks)

## NODES DATA

nodegroups = list()
for(i in 1:4){
  vals = as.numeric(unique(d[[i]]))
  nodegroup = data.frame(id = 1:length(vals), lab = unique(d[[i]]), vals = vals,
                         radius = 100 * vals/max(vals), axis = i)
  sizes = table(d[[i]])
  nodegroup$size = as.numeric(sizes[ match(nodegroup$lab, names(sizes)) ])
  nodegroup$size = 2 * nodegroup$size / max(nodegroup$size)
  if(i>1) nodegroup$id = nodegroup$id + max(nodegroups[[i-1]]$id)
  nodegroups[[ names(d)[i] ]] = nodegroup
}
nodegroups

nodes = rbind(nodegroups[[1]], nodegroups[[2]], nodegroups[[3]], nodegroups[[4]])
nodes$lab = as.character(nodes$lab)
nodes$axis = as.integer(nodes$axis)
nodes$radius = as.numeric(nodes$radius)
nodes$color = "#ffffff"
head(nodes)

## EDGES DATA

# first update edge data with new node IDs
head(d)
for(i in 1:4) {
  header = paste0(names(nodegroups)[i], 'id')
  d[[header]] = nodegroups[[i]]$id[ match(as.numeric(d[[i]]), nodegroups[[i]]$vals) ]
}
head(d)

# edges between the 4 axes in terms of node IDs
for(i in 6:8){
  edgegroup = data.frame(id1 = d[[i]], id2 = d[[i+1]], price = d[[5]])
  if(i==6) all_edges = edgegroup else all_edges = rbind(all_edges, edgegroup)
}
head(all_edges); dim(all_edges)

# summarise edge data
edges = aggregate(all_edges$price, by=list(all_edges$id1, all_edges$id2), FUN='mean')
names(edges) = c('id1','id2','price')
edges = edges[with(edges, order(id1,id2)),]             # reorder

# set edge weights (stroke thickness)
weights = count(all_edges, vars = c('id1', 'id2'))      # summary data
weights = weights[with(weights, order(id1,id2)),]       # reorder to match egdes
all(weights$id1 == edges$id1, weights$id2 == edges$id2) # check all IDs match up
edges$weight = weights$freq * 0.004
edges$weight = pmax(edges$weight, 0.2)  # set min edge weight to still visible
range(weights$freq)
range(edges$weight)

# normalise prices for each group of edges (to utilise full colour range)
p = edges$price
edges$colorvals = 0

for(i in nodegroups[1:3]){
  sel = edges$id1 %in% range(i$id)[1] : range(i$id)[2]
  edges$colorvals[sel] = (p[sel] - min(p[sel])) / (max(p[sel]) - min(p[sel]))
}

edges$color = paste0(hex(HSV(edges$colorvals * 300, 1, 1)), '60')  # set alpha
edges = edges[order(edges$weight, decreasing=T),]   # draw thin edges last

head(edges)

hpd = list()
hpd$nodes = nodes
hpd$edges = edges
hpd$type = "2D"
hpd$desc = "Diamonds"
hpd$axis.cols = rep('#00000000', 4) # make invisible
hpd$axLabs = c("carats","cut","colour","clarity")
class(hpd) = "HivePlotData"

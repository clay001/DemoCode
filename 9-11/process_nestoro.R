library(destiny)
library(rgl)
origin <- read.table('normalisedCountsVariableGenes.txt', header = T, row.names = 1)
# breif show
origin[1:5,1:5]

# variance-stabilizing transformation
logOrigin <- t(log2(origin + 1))

# create a diffusionMap object
dm <- DiffusionMap(logOrigin,distance = "cosine",sigma = 0.16)
plot(dm, c(3,2,1), pch=20)
# length is number of list 1656 (row number)
# 列是4773列基因

#for prettier figure (according to the plot document R)
palette(cube_helix(6))

dpt <- DPT(dm)
plot(dpt, root = 2, paths_to = c(1,3), col_by = 'branch')

transition_nestoro <- as.matrix(dm@transitions)
write.csv(transition_nestoro,file = "transition_nestoro.csv")

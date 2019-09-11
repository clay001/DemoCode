library(readxl)
library(rgl)
library(destiny) 
library(Biobase)
raw_ct <- read_xls('mmc4.xls', 'Sheet1')
raw_ct[1:9, 1:9] #preview of a few rows and columns

# turn tibbles to Expression set
ct <- as.ExpressionSet(as.data.frame(raw_ct))  # 48*442 opposite
ct

###########################################
#cleaning
#remove embryo stage 1 
#see higher than background NA
# 取出列名用正则，再变为整数
num_cells <- gsub('^(\\d+)C.*$', '\\1', ct$Cell) 
ct$num_cells <- as.integer(num_cells)

# cells from 2+ cell embryos
# boolean table
have_duplications <- ct$num_cells > 1
# cells with values ≤ 28
# 2表示对列进行处理
normal_vals <- apply(exprs(ct), 2, function(smp) all(smp <= 28))
# logical judgement,横行基因全部取，纵列细胞做清洗
cleaned_ct <- ct[, have_duplications & normal_vals]
# S4

####################################
#normalization
#endogenous controls Actb and Gapdh
housekeepers <- c('Actb', 'Gapdh')
# 求house基因的列平均
normalizations <- colMeans(exprs(cleaned_ct)[housekeepers, ])
# 改名字
guo_norm <- cleaned_ct
# 减去house的表达量
exprs(guo_norm) <- exprs(guo_norm) - normalizations

#####################################
# dimension reduction
dm <- DiffusionMap(guo_norm)
plot(dm)

####### R mode document #############
palette(cube_helix(6)) 
plot(dm, pch = 20, col_by = 'num_cells', legend_main = 'Cell stage')

################################
# turn to 2D
plot(dm, 1:2,col_by = 'num_cells',legend_main = 'Cell stage')

##################need to use X11############################
# more fancy use
save.image()
plot3d(eigenvectors(dm)[, 1:3],col = log2(guo_norm$num_cells), type = 's',radius = .01) 
view3d(theta = 10, phi = 30, zoom = 0.8)
rgl.close()

##############################can also use ggplot
###hist
hist(exprs(cleaned_ct)['Aqp3', ], breaks = 20,
     xlab = 'Ct of Aqp3', main = 'Histogram of Aqp3 Ct',
     col = palette()[[4]], border = 'white')

#for scRNA-Seq, it is may be necessary to first perform a 
#Principal Component Analysis (PCA) on the data (e.g. using prcomp or princomp) 
#and to calculate the Diffusion Components from the Principal Components
#(typically using the top 50 components yields good results).

###DPT
# 可以将dm转换为dpt对象，有branch和tips分量
# branch递归7层的分支标签，tips指示是否为当前标签的端点
dpt <- DPT(dm)
plot(dpt)
#library(viridis)
# 画图应该是默认取branch的第一列的
plot(dpt, root = 1, paths_to = c(2,3), col_by = 'branch')

transition_p <- as.matrix(dm@transitions)
write.csv(transition_p,file = "transition_p.csv")


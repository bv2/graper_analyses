# color and method defintions
library(RColorBrewer)

get_legend<-function(gg){
  tmp <- ggplot_gtable(ggplot_build(gg))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

methodnm <- "graper"

make_nicenames <- function(nm){
  ifelse(nm=="graper_SS", paste0(methodnm, " (sparse)"),
         ifelse(nm=="graper", paste0(methodnm, " (dense, multiv.)"),
                ifelse(nm=="graper_SS_ungrouped", paste0(methodnm, " (sparse, fact., no groups)"),
                       ifelse(nm=="graper_FF", paste0(methodnm, " (dense)"),
                              ifelse(nm=="adaptiveLasso", "adaptive Lasso",
                                     ifelse(nm=="GroupLasso", "group Lasso",                                   
                                            ifelse(nm=="SparseGroupLasso", "sparse group Lasso",
                                                   ifelse(nm=="IPFLasso", "IPF-Lasso",
                                                          ifelse(nm=="ElasticNet", "elastic net",
                                                                 ifelse(nm=="Ridge", "ridge regression", as.character(nm)))))))))))
}



methods2compare_sparse <- c("graper_SS", "Lasso", "SparseGroupLasso",
"IPFLasso","adaptiveLasso", "ElasticNet",  "varbvs")
methods2compare_dense <- c("Ridge", "graper_FF", "graper", "GRridge", "GroupLasso")

methods2compare_sparse <- sapply(methods2compare_sparse, make_nicenames)
methods2compare_dense <- sapply(methods2compare_dense, make_nicenames)


# cols4methods <- c(c(wes_palette("Darjeeling1"),
#                     wes_palette("Darjeeling2")[-c(1,3,4)])[1:length(methods2compare_sparse)],
#                   brewer.pal(8,"Set1")[c(1,2,4,5,7)])

cols4methods <- rep("gray", 12)
cols4methods[c(1,9,10)] <- c("cornflowerblue","cyan4",  "navy")
cols4methods[c(2,4,5,6)] <- c("darkgoldenrod1", "darkorange3", "coral2", "brown")
cols4methods[c(8,11)] <- c("deeppink2", "darkmagenta")
cols4methods[c(3, 12)] <- c("burlywood3", "chocolate4")
cols4methods[c(7)] <- "darkgrey"

names(cols4methods) <- c(methods2compare_sparse, methods2compare_dense)
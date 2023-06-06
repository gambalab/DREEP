#' @import Matrix
#' @import snow
#' @import progress
make.gene.sets = function(xx,n.genes=250,nc=4)
{
  message("Extracting top scoring genes from each cell..")
  sq = seq(1,ncol(xx),ifelse(ncol(xx)<1000,floor(ncol(xx)/2),1000))

  pb <- progress_bar$new(
    format = "  processing [:bar] :current/:total (:percent) [eta::eta]",
    total = length(sq), clear = FALSE, width= 60)
  pb$tick(0)

  cl = snow::makeCluster(nc)
  l = NULL
  for (i in 2:length(sq))
  {
    x = as.matrix(xx[,sq[i-1]:(sq[i]-1)])
    #print("Exporting data...")
    snow::clusterExport(cl,'x',envir = environment())
    snow::clusterExport(cl,'n.genes',envir = environment())
    #print("Started...")
    l = c(l,snow::parLapply(cl,seq_len(ncol(x)), function(idx) {names(x[order(x[,idx],decreasing = T)[1:n.genes],idx])}))
    pb$tick()
  }

  if (sq[length(sq)]<=ncol(xx))
  {
    x = as.matrix(xx[,sq[i]:ncol(xx)])
    #print("Exporting data...")
    snow::clusterExport(cl,'x',envir = environment())
    snow::clusterExport(cl,'n.genes',envir = environment())
    #print("Started...")
    l = c(l,snow::parLapply(cl,seq_len(ncol(x)), function(idx) {names(x[order(x[,idx],decreasing = T)[1:n.genes],idx])}))
  }
  snow::stopCluster(cl)
  pb$tick()
  message("Completed!!")
  names(l) = colnames(xx)
  return(l)
}


#' Drug Response Estimation from single-cell Expression Profiles (DREEP)
#'
#' Drug Response Estimation from single-cell Expression Profiles
#'
#' @import fastmatch
#' @import snow
#' @import fgsea
#' @export
runDREEP = function(M,n.markers=250,cores=0,gsea="multilevel",gpds.signatures=c("CTRP2","GDSC")) {
  gpds.signatures <- toupper(gpds.signatures)
  gpds.signatures = base::match.arg(arg = gpds.signatures,choices = c("CTRP2","GDSC","PRISM"),several.ok = TRUE)
  gsea = base::match.arg(arg = gsea,choices = c("simple","multilevel"),several.ok = FALSE)

  if (cores==0) {cores = detectCores()}

  message("Loading GPDS signatures..")
  meta.drug = NULL
  M.drug.mkrs = NULL

  pb <- progress_bar$new(
    format = "  processing [:bar] :current/:total (:percent) [eta::eta]",
    total = length(gpds.signatures), clear = FALSE, width= 60)
  pb$tick(0)

  for(i in gpds.signatures) {
    M.drug.mkrs = cbind(M.drug.mkrs,readRDS(paste0(find.package("DREEP"),"/data/",i,".gpds.rds")))
    meta.drug = rbind(meta.drug,readRDS(paste0(find.package("DREEP", quiet = FALSE),"/data/",i,".drug.metadata.rds")))
    pb$tick()
  }
  message("DONE!")

  M = M[rownames(M)%in%rownames(M.drug.mkrs),]
  M.drug.mkrs = M.drug.mkrs[rownames(M.drug.mkrs) %in% rownames(M),]

  l.mkrs = make.gene.sets(M,n.markers,1)

  message("Running DREEP...")
  cl = snow::makeCluster(cores)
  genes = rownames(M.drug.mkrs)
  snow::clusterExport(cl,c("l.mkrs","genes"),envir = environment())
  r = snow::parApply(cl,M.drug.mkrs, 2, function(x){
    names(x) = genes
    if (gsea=="simple"){
      res = fgsea::fgseaSimple(pathways = l.mkrs,stats = x,nproc = 1,gseaParam = 0, nperm = 1000)
    } else {
      res = fgsea::fgseaMultilevel(pathways = l.mkrs,stats = x,nproc = 1,gseaParam = 0, eps = 0)
    }
    res = as.data.frame(res[,c("pval","ES")])
  })
  snow::stopCluster(cl)
  message("FINISHED!!")

  df = do.call("rbind",lapply(r, function(x) data.frame(sens=sum(x[,"pval"] < 0.05 & x[,"ES"]<0)/nrow(x),res=sum(x[,"pval"] < 0.05 & x[,"ES"]>0)/nrow(x),med=median(x[,"ES"]))))
  df$conpound.id = sapply(strsplit(x = names(r),split = "_",fixed = T),function(x) x[[2]])
  df$conpound = sapply(strsplit(x = df$conpound.id,split = ":",fixed = T),function(x) x[[1]])
  df$dataset = sapply(strsplit(x = names(r),split = "_",fixed = T),function(x) x[[1]])
  rownames(df) = df$conpound.id
  df = df[order(df$med,decreasing = F),]
  df$target = meta.drug$target[fastmatch::fmatch(df$conpound.id,meta.drug$id)]
  df$moa = meta.drug$moa[fastmatch::fmatch(df$conpound.id,meta.drug$id)]
  df$smiles = meta.drug$smiles[fastmatch::fmatch(df$conpound.id,meta.drug$id)]
  df$dataset = meta.drug$dataset[fastmatch::fmatch(df$conpound.id,meta.drug$id)]
  M.es = sapply(r, function(x) x[,"ES"])
  M.pval = sapply(r, function(x) x[,"pval"])
  rownames(M.es) = rownames(M.pval) = colnames(M)

  return(list("df"=df,"es.mtx"=M.es,"es.pval"=M.pval))
}

#' Differential Drug Analysis
#'
#' Identify a set of compounds specific for two cell populations
#'
#' @import fastmatch
#' @import ggplot2
#' @export
runDiffDrugAnalisys = function(dreep.data,cell.set1,cell.set2,fdr.th=0.1,show.plot=T) {
  tmp <- data.frame(conpound.id=sapply(strsplit(x = colnames(dreep.data$es.mtx),split = "_",fixed = T),function(x) x[[2]]))
  rownames(tmp) = tmp$conpound.id
  tmp$conpound <- sapply(strsplit(x = tmp$conpound.id,split = ":",fixed = T),function(x) x[[1]])
  tmp$dataset = sapply(strsplit(x = colnames(dreep.data$es.mtx),split = "_",fixed = T),function(x) x[[1]])
  tmp$med1 <- apply(dreep.data$es.mtx[cell.set1,],2,median)
  tmp$med2 <- apply(dreep.data$es.mtx[cell.set2,],2,median)

  tmp$specific <- sign(tmp$med1) != sign(tmp$med2)
  oldw <- getOption("warn")
  options(warn = -1)
  tmp$p.value = apply(dreep.data$es.mtx, 2, function(x,c1=cell.set1,c2=cell.set2) ks.test(x[c1],x[c2])$p.value)
  options(warn = oldw)
  tmp$fdr <- p.adjust(tmp$p.value,method = "fdr")
  tmp$significant <- factor(as.character(tmp$fdr<fdr.th & tmp$specific),levels = c("TRUE","FALSE"))

  if(show.plot) {
    p <- ggplot(data = tmp,aes(x=med1,y=med2,color=significant,size=significant)) +
          geom_point() + geom_abline(color="red") +
          theme_bw() + coord_equal() +
          xlab("Drug median ES (cell set 1)") + ylab("Drug median ES (cell set 2)") +
          scale_color_manual(values = c("red","gray")) +
          scale_size_manual(values = c(.5,.1))
    print(p)
  }

  tmp$target = dreep.data$df$target[fastmatch::fmatch(tmp$conpound.id,dreep.data$df$conpound.id)]
  tmp$moa = dreep.data$df$moa[fastmatch::fmatch(tmp$conpound.id,dreep.data$df$conpound.id)]
  tmp$smiles = dreep.data$df$smiles[fastmatch::fmatch(tmp$conpound.id,dreep.data$df$conpound.id)]
  tmp <- tmp[order(tmp$fdr),]
  tmp$population <- "cell.set1"
  tmp$population[tmp$med2<0] <- "cell.set2"
  dreep.data$DiffDrugRes <- tmp[tmp$fdr<fdr.th & tmp$specific,-which(colnames(tmp)%in%c("specific","significant"))]
  return(dreep.data)
}

#' Cell Dimensionality Reduction
#'
#' Cell Dimensionality Reduction in the drug space
#'
#' @import ggplot2
#' @importFrom uwot umap tumap
#' @import Rtsne
#' @export
runDrugReduction = function(dreep.data,pval.th=0.05,drug.subset=NULL,cores=0,seed=180582,verbose=T,reduction="umap",cellDistAbsolute=T,storeCellDist=F, ...) {
  if (cores==0) {cores = ifelse(detectCores()>1,detectCores()-1,1)}
  reduction = base::match.arg(arg = reduction,choices = c("umap","tsne","tumap"),several.ok = FALSE)

  if(is.null(drug.subset)) {
    tmp.es <- dreep.data$es.mtx
    } else {
      tmp.es <- dreep.data$es.mtx[,drug.subset]
  }
  tmp.es[dreep.data$es.pval[,colnames(tmp.es)]<pval.th] <- 0
  tmp.es = tmp.es[rowSums(tmp.es!=0)>0,]

  d = cellDistWrapper(t(tmp.es),cores,cellDistAbsolute,verbose)
  base::set.seed(seed)
  if (reduction=="umap") {
    dreep.data$embedding = as.data.frame(uwot::umap(X = d,
                                              scale = F,
                                              n_threads = cores,
                                              verbose = verbose,
                                              ret_model = F,
                                              ...)
                                         )
  }

  if (reduction=="tumap") {
    dreep.data$embedding = as.data.frame(uwot::tumap(X = d,
                                                    scale = F,
                                                    n_threads = cores,
                                                    verbose = verbose,
                                                    ret_model = F,
                                                    ...)
    )
  }

  if (reduction=="tsne") {
    dreep.data$embedding = base::as.data.frame(Rtsne::Rtsne(X = d,
                                                            dims = 2,
                                                            pca = F,
                                                            verbose = verbose,
                                                            max_iter=1000,
                                                            num_threads=cores,
                                                            is_distance = TRUE,
                                                            ...)$Y)
  }

  rownames(dreep.data$embedding) = base::rownames(tmp.es)
  colnames(dreep.data$embedding) = base::c("X","Y")

  if(storeCellDist) {
    dreep.data$cellDist=d
    rm(d,tmp.es);gc()
  } else {
    rm(d,tmp.es);gc()
  }

  return(dreep.data)
}


#' Detect the Number of CPU Cores
#'
#' Detects the number of (logical) CPU cores.
#'
detectCores <- function() {
  .Call("detectCoresCpp")
}


cellDistWrapper <- function(m,cores,cellDistAbsolute,verbose) {
  message("Computing cell distances in the drug space..")
  as.dist(cellDist(m = Matrix(data = m,sparse = T),ncores = cores,verbose = verbose,full = F,diag = F,absolute=cellDistAbsolute))
}

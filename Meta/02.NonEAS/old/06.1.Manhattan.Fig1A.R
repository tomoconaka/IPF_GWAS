library(lattice)
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}
library(data.table)
library(dplyr)
library(tidyr)
d <- fread("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/Meta_scratch/ALL_meta_IPF.b38_filtered.txt.gz")
#d <- fread(args[1])
#d <- d %>% filter(Pvalue < 0.05)
d$Pvalue <- as.numeric(d$Pvalue)
d <- d %>% mutate(Pvalue = ifelse(Pvalue < 1e-65, 1e-65, Pvalue))
ann<-rep(1, length(d$Pvalue))
ann[with(d, CHR==1 & POS>=9109660 - 500000 & POS<=9109660 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==1 & POS>=150579566 - 500000 & POS<=150579566 + 500000 & Pvalue < 5e-8)] <-3
ann[with(d, CHR==1 & POS>=155199564 - 500000 & POS<=155199564 + 500000 & Pvalue < 5e-8)] <-3
ann[with(d, CHR==1 & POS>=163373599 - 500000 & POS<=163373599 + 500000 & Pvalue < 5e-8)] <-3
ann[with(d, CHR==1 & POS>=214487380 - 500000 & POS<=214487380 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==2 & POS>=64652534 - 500000 & POS<=64652534 + 500000 & Pvalue < 5e-8)] <-3
ann[with(d, CHR==3 & POS>=44805230 - 500000 & POS<=44805230 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==3 & POS>=169768720 - 500000 & POS<=169768720 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==4 & POS>=88915668 - 500000 & POS<=88915668 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==5 & POS>=1282299 - 500000 & POS<=1282299 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==5 & POS>=169477001 - 500000 & POS<=169477001 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==6 & POS>=7562999 - 500000 & POS<=7562999 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==6 & POS>=27698141 - 500000 & POS<=27698141 + 500000 & Pvalue < 5e-8)] <-2
ann[with(d, CHR==6 & POS>=32436600 - 500000 & POS<=32436600 + 500000 & Pvalue < 5e-8)] <-3
ann[with(d, CHR==6 & POS>=43385242 - 500000 & POS<=43385242 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==7 & POS>=1936821 - 500000 & POS<=1936821 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==7 & POS>=100032719 - 500000 & POS<=100032719 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==8 & POS>=119931204 - 500000 & POS<=119931204 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==9 & POS>=106717987 - 500000 & POS<=106717987 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==10 & POS>=103920828 - 500000 & POS<=103920828 + 500000 & Pvalue < 5e-8)] <-4
ann[with(d, CHR==11 & POS>=1219991 - 500000 & POS<=1219991 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==12 & POS>=101892202 - 500000 & POS<=101892202 + 500000 & Pvalue < 5e-8)] <-2
ann[with(d, CHR==13 & POS>=112886335 - 500000 & POS<=112886335 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==15 & POS>=40426335 - 500000 & POS<=40426335 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==15 & POS>=85744679 - 500000 & POS<=85744679 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==16 & POS>=112241 - 500000 & POS<=112241 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==16 & POS>=67895674 - 500000 & POS<=67895674 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==17 & POS>=46253848 - 500000 & POS<=46253848 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==17 & POS>=76029575 - 500000 & POS<=76029575 + 500000 & Pvalue < 5e-8)] <-3
ann[with(d, CHR==18 & POS>=666625 - 500000 & POS<=666625 + 500000 & Pvalue < 5e-8)] <-3
ann[with(d, CHR==19 & POS>=4717660 - 500000 & POS<=4717660 + 500000 & Pvalue < 5e-8)]<-4
ann[with(d, CHR==20 & POS>=63638397 - 500000 & POS<=63638397 + 500000 & Pvalue < 5e-8)]<-4
ann<-factor(ann, levels=1:4)
png(paste0(paste0("/mnt/gpfs1/project/Interstitial_Pneumonia/ANALYSIS/tnakanishi/repo/IPF_GWAS/Figures/Fig1a.Manhattan.png")), width=1200, height=500)
manhattan.plot(d$CHR, d$POS, d$Pvalue, sig.level=5e-8, annotate = ann)
dev.off()






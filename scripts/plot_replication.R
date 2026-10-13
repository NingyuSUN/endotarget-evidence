#!/usr/bin/env Rscript
# Standard base-R figures from report-builder tables. No inference or packages.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4) stop("Usage: plot_replication.R data_dir output_dir prefix label")
data_dir<-args[1]; out<-args[2]; prefix<-args[3]; label<-args[4]
read_table<-function(name) read.delim(file.path(data_dir,name), stringsAsFactors=FALSE, check.names=FALSE)
is_true<-function(x) tolower(as.character(x))=="true"
render<-function(name,fun,w=8,h=5.5) {
  for (ext in c("pdf","png")) {
    dest<-file.path(out,paste0(prefix,"_",name,".",ext))
    if(file.exists(dest)) stop("Refuse overwrite: ",dest)
    if(ext=="pdf") pdf(dest,width=w,height=h,useDingbats=FALSE) else png(dest,width=w*160,height=h*160,res=160)
    par(mar=c(4.5,4.8,3.6,1),las=1)
    fun()
    mtext(label,side=3,line=2,cex=.7)
    dev.off()
  }
}
counts<-read_table("counts.tsv")
render("counts",function() {
  par(mfrow=c(1,2),mar=c(4.5,4.2,3.6,1))
  for(k in c(20,50)) {
    s<-counts[counts$k==k,]; s<-s[match(c("B0","B1","G1"),s$method),]
    barplot(t(as.matrix(s[,c("agree","opposite","unevaluable")])),names.arg=s$method,
      col=c("#2A788E","#D95F02","#BDBDBD"),border=NA,ylim=c(0,k*1.3),
      ylab=paste0("Genes / fixed ",k," selected"),main=paste0("Top ",k))
    legend("top",legend=c("Same direction","Opposite","Unevaluable"),
      fill=c("#2A788E","#D95F02","#BDBDBD"),bty="n",cex=.7)
  }
})
scatter<-read_table("scatter.tsv")
render("scatter",function() {
  if(nrow(scatter)==0) {plot.new(); text(.5,.5,"No finite common-gene effects"); return()}
  sel<-is_true(scatter$selected)
  plot(scatter$discovery_logFC,scatter$external_logFC,pch=16,cex=.45,
    col=adjustcolor("#9E9E9E",alpha.f=.3),xlab="Discovery ENDO-minus-non-ENDO log2 fold change",
    ylab="HECA ENDO-minus-non-ENDO log2 fold change",
    main=paste0("Common fitted network genes: n=",nrow(scatter)))
  abline(h=0,v=0,col="#666666",lty=2)
  points(scatter$discovery_logFC[sel],scatter$external_logFC[sel],pch=21,bg="#2A788E",col="#173D45",cex=.8)
  legend("topleft",legend=c("All common fitted B0-universe genes","Frozen B0/B1/G1 top20 union"),
    pch=c(16,21),col=c("#9E9E9E","#173D45"),pt.bg=c(NA,"#2A788E"),bty="n",cex=.8)
  mtext(paste0("Finite selected points: ",sum(sel),"; missing outcomes not plotted; associations, not intervention effects"),side=1,line=3.4,cex=.65)
},w=8,h=6)
if(file.exists(file.path(data_dir,"lodo.tsv"))) {
  folds<-read_table("lodo.tsv")
  render("lodo",function() {
    par(mfrow=c(1,2),mar=c(4.5,4.8,3.6,1))
    for(metric in c("retention","jaccard")) {
      plot(NA,xlim=c(.5,3.5),ylim=c(-.09,1.05),xaxt="n",
        xlab="Discovery ranking method",ylab=if(metric=="retention") "Intersection / original 20" else "Intersection / union",
        main=if(metric=="retention") "Top20 retention" else "Top20 Jaccard")
      axis(1,at=1:3,labels=c("B0","B1","G1"))
      abline(h=c(0,.5,1),col="#E0E0E0",lty=2)
      for(i in 1:3) {
        s<-folds[folds$method==c("B0","B1","G1")[i],]
        good<-is_true(s$valid); offset<-(as.numeric(s$fold)-mean(as.numeric(s$fold)))*.013
        points(i+offset[good],s[[metric]][good],pch=16,col="#2A788E",cex=.65)
        if(any(!good)) points(i+offset[!good],rep(-.05,sum(!good)),pch=4,col="#D95F02",lwd=2)
        text(i,1.04,paste0(sum(good),"/",nrow(s)," valid"),cex=.7)
      }
      mtext("Crosses below zero mark failed folds; failures retained",side=1,line=3.4,cex=.65)
    }
  },w=9,h=5.5)
}
if(file.exists(file.path(data_dir,"bootstrap.tsv"))) {
  boot<-read_table("bootstrap.tsv"); limits<-read_table("bootstrap_interval.tsv")
  valid<-is_true(boot$valid); v<-boot$G1_minus_B0_top20[valid]
  if(length(v)<450 || limits$failed>50) stop("Bootstrap interval suppression required")
  render("bootstrap",function() {
    if(diff(range(v))==0) {
      plot(v,rep(1,length(v)),pch=16,col=adjustcolor("#2A788E",alpha.f=.25),
        xlim=range(v)+c(-.025,.025),xlab="Paired G1-minus-B0 top20 concordance difference",
        ylab="Identical valid draws (one displayed level)",main=paste0(length(v)," valid / 500 total draws"))
    } else hist(v,breaks=30,col="#2A788E",border="white",
      xlab="Paired G1-minus-B0 top20 concordance difference",ylab="Valid draws",
      main=paste0(length(v)," valid / 500 total draws"))
    abline(v=0,col="#666666",lty=2)
    abline(v=c(limits$low,limits$high),col="#D95F02",lwd=2)
    mtext(paste0("Failed=",limits$failed,"/500; type7 conditional interval; fixed discovery lists/directions"),side=1,line=3.4,cex=.65)
  },w=8,h=5.5)
}
cat("R=",R.version.string,"; base graphics; no model fitting\n",sep="")


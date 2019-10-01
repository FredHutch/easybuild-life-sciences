PEER_plotModel <- function(model){
    par(mfrow=c(2,1))
    #bounds = PEER_getBounds(model)
    bounds = (1:10) + rnorm(10)
    #vars = PEER_getResidualVars(model)
    vars = -(1:10)/10. + rnorm(10,0,0.1)
    par(mar=c(5,4,4,5)+.1)
    plot(bounds, type="l", col="red", lwd=2, xlab="Iterations", ylab="Lower bound")
    par(new=TRUE)
    plot(vars,,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
    axis(4)
    mtext("Residual variance",side=4,line=3)
    legend("left",col=c("red","blue"),lty=1,legend=c("Lower bound","Residual variance"))
    alpha = PEER_getAlpha(model)
    plot(alpha,xlab="Factors",ylab="Inverse variance of factor weights", type="b", col="blue", lwd=4, xaxp=c(1,length(alpha), length(alpha)-1))
}

# LOD scores (base 10) [n_snps x n_genes], snp locations => plot of gene discoveries based on LOD cutoff, density plot
PEER_plotQtl <- function(lods, snploc, filterfun, cutoff=5){
    par(mfrow=c(2,1))
    lod_cutoffs = as.matrix(1:30)
    max_lods = apply(lods, 2, max) # largest LODs per gene    
    plot(lod_cutoffs, apply(lod_cutoffs, 1, function(x){sum(x < max_lods)}), "l", lwd=6, col="blue", xlab="LOD cutoff", ylab="# genes with eQTL", cex.lab=1.5, cex.axis=1.2)
    plot(snploc, apply(lods>cutoff, 1, sum), "l", lwd=6, col="red", xlab="Genome", ylab="# eQTLs")    
}


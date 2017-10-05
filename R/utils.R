#=============================================================================#
#=============================================================================#
#' Plot 2 x 2 confusion matrix
#'
#' This function creates an illustrated confusion matrix for two-class problems
#'
#' @param confMatrix a 2 x 2 confusion matrix (classA | classB | classError)
#' predicted vs observed
#' @param plotLabels a character vector of length two with the ordered name of the
#' two classes (default: colnames(confMatrix)[1:2])
#' @param plotTitle main title of plot (default: "")
#' @param colDiag colour of the diagonal segments (correctly classified) 
#' (default: "darkblue")
#' @param colOffDiag colour of the off diagnonal segments (misclassified)
#' (default: "grey")
#' @param colAxes colour of the x-y axes (default: "blue")
#' @param colTicks colour of the x-y ticks (default: "red")
#' @return None, a plot is created
#' @export
#' 
plotconfusionmatrix <- function(confMatrix, plotLabels=colnames(confMatrix)[1:2], plotTitle="", 
                                colDiag="darkblue", colOffDiag="grey",
                                colAxes="blue", colTicks="red")
{
    # Parse confusion matrix (a numeric matrix)
    proptab = prop.table(confMatrix, margin=1)
    r = as.vector(proptab)
    par(pty="s")
    
    # Creat plot and axes
    plot(NULL, type="n", xlab="", ylab="", xaxt='n',yaxt='n',main=plotTitle,
         xlim=c(0,4),ylim=c(0,4),xaxs = "i", yaxs = "i",asp=1,cex.lab=1.2,font.lab=2)
    axis(1, at=c(1, 3), labels=plotLabels,cex.axis=1.2)
    axis(2, at = c(1, 3), labels=rev(plotLabels),cex.axis=1.2)
    mtext("Predicted", side=1, line=3, cex=1.1, font=2)
    mtext("Observed", side=2, line=3, cex=1.1, font=2)
    abline(h=2, lwd=0.8, col=colAxes)
    abline(v=2, lwd=0.8, col=colAxes)
    
    # Create ticks
    segments(x0=seq(0.05,1.95,by=0.38), y0=1.98, y1=2.02, col=colTicks, lwd=0.8)
    segments(x0=seq(2.05,3.95,by=0.38), y0=1.98, y1=2.02, col=colTicks, lwd=0.8)
    segments(y0=seq(0.05,1.95,by=0.38), x0=1.98, x1=2.02, col=colTicks, lwd=0.8)
    segments(y0=seq(2.05,3.95,by=0.38), x0=1.98, x1=2.02, col=colTicks, lwd=0.8)
    
    # Plot circles
    circlize::draw.sector(start.degree=90, end.degree=180, center=c(1.95,2.05), 
                          rou1=1.9*r[1], col=colDiag, clock.wise=F)
    circlize::draw.sector(start.degree=180, end.degree=270, center=c(1.95,1.95), 
                          rou1=1.9*r[2], col=colOffDiag, clock.wise=F)
    circlize::draw.sector(start.degree=270, end.degree=360, center=c(2.05,1.95), 
                          rou1=1.9*r[4], col=colDiag, clock.wise=F)
    circlize::draw.sector(start.degree=0, end.degree=90, center=c(2.05,2.05), 
                          rou1=1.9*r[3], col=colOffDiag, clock.wise=F)
    
    # Add misclassification rates
    text(0.1, 3.9, paste(format(100*r[1],digits = 1), "% (N=",confMatrix[1, 1],")",sep=''), 
         adj=c(0,1),font=2)
    text(0.1, 0.1, paste(format(100*r[2],digits = 1), "% (N=",confMatrix[2, 1],")",sep=''), 
         adj=c(0,0),font=2)
    text(3.9, 3.9, paste(format(100*r[3],digits = 1), "% (N=",confMatrix[1, 2],")",sep=''), 
         adj=c(1,1),font=2)
    text(3.9, 0.1, paste(format(100*r[4],digits = 1), "% (N=",confMatrix[2, 2],")",sep=''), 
         adj=c(1,0),font=2)
    
    par(pty="m")
} 
#=============================================================================#
#=============================================================================#
#' Plot ROC curves
#'
#' This function plots ROC curves for each fold and their average
#'
#' @param ROC a list of ROC curves produced by ROCR and named ("Fold1, Fold2,...,Average")
#' @param AUC a list of area under the ROC curve (same size as ROC)
#' @param colAvg colour for the average ROC curve
#' @param colFold colour for each fold ROC curve
#' @param colRandom colour for randomly guession (y=x) line
#' @return None, a plot is created
#' @export
#' 
plotROC <- function(ROC, AUC, colAvg="black", colFold="grey", colRandom="red")
{
    if (!all(names(ROC)==names(AUC)))
    {
        print("ROC and AUC should have the same named columns!")
    } else {
        # Plot average ROC
        par(pty="s", cex.axis=1.5)
        plot(ROC[["Average"]], col=colAvg, lwd=4, font.lab=1.5, cex.lab=1.5)
        legStr <- c()
        
        # Plot all other folds
        folds <- setdiff(names(ROC), "Average")
        for (fold in folds)
        {
            k <- substr(fold, nchar(fold), nchar(fold)) # extract fold number
            plot(ROC[[fold]], col=colFold, lwd=2, lty=2, add=T)
            legStr <- c(legStr, paste("AUC", k, "=", format(AUC[[fold]], digits=2)))
        }
        # Add mean ROC legend
        legStr <- c(legStr, paste("AUC (avg)", "=", format(AUC[["Average"]], digits=2))) 
        
        # Add randomly guessing line
        abline(a=0, b=1, col=colRandom, lty=3, lwd=2)
        
        # Add legend
        legend("bottomright", legStr, col=c(rep(colFold, length(folds)), colAvg), 
               lwd=c(rep(2, length(folds)), 4), lty=1, cex=1.2)
    }
    par(pty="m")
}
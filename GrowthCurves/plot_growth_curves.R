#-----------------------------------------------------------------
# Plot growth curves
#-----------------------------------------------------------------
library(ggplot2)
library(grid)

asfIDs = c(356,360,361,492,500,502,519)

averageGrowthCurve <- function(fileName)
{  
  gc = read.table(fileName,sep='\t',header=F,row.names=NULL)
  gcdf = data.frame(gc)
  names(gcdf) = c("time","OD","group")
  
  # Determine average growth curve and plot a standard error around it    
  t1 = gcdf$time[gcdf$group == 1]
  od1 = gcdf$OD[gcdf$group == 1]
  od1[od1<0.001] = 0.001
  t2 = gcdf$time[gcdf$group == 2]
  od2 = gcdf$OD[gcdf$group == 2]
  t3 = gcdf$time[gcdf$group == 3]
  od3 = gcdf$OD[gcdf$group == 3]
  t4 = gcdf$time[gcdf$group == 4]
  od4 = gcdf$OD[gcdf$group == 4]
  if (length(t2) != length(t1))
  {
    fn2 = splinefun(t2, y = od2, method = "natural", ties = mean)
    od2 = fn2(t1)
    od2[od2<0.001] = 0.001
  }
  if (length(t3) != length(t1))
  {
    fn3 = splinefun(t3, y = od3, method = "natural", ties = mean)
    od3 = fn3(t1)
    od3[od3<0.001] = 0.001
  }
  if (length(t4) != length(t1))
  {
    fn4 = splinefun(t4, y = od4, method = "natural", ties = mean)
    od4 = fn4(t1)
    od4[od4<0.001] = 0.001
  }
  ODs = matrix(0,ncol=4,nrow=length(t1))
  ODs[,1] = od1
  ODs[,2] = od2
  ODs[,3] = od3
  ODs[,4] = od4
  
  means = c()
  stdError = c()
  for (k in 1:nrow(ODs))
  {
    curRow = ODs[k,]
    means[k] = mean(as.numeric(curRow))
    stdError[k] = sd(as.numeric(curRow))/2
  }
  curvedf = data.frame(time=t1,ODave=means,ODstdError=stdError)
  return(curvedf)
}

tiff(filename = "ASF_SpentMedia_GrowthCurves.tiff", width = 1000, height = 600, units = "px")

layout = matrix(seq(1, 7 * 8), ncol = 8, nrow = 7)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
lengths = c()
for(i in asfIDs)
{
  # Get the maximum OD for each row
  maxOD = 0.1
  for(j in c(0,asfIDs))
  {
    fileName = paste0("SpentMediaGrowthCurveData\\","ASF",i,"in",j,".tsv")
    gc = read.table(fileName,sep='\t',header=F,row.names=NULL)
    print(max(gc[,2]))
    if(max(gc[,2]) > maxOD)
    {
      maxOD = max(gc[,2])
    }
  }
  
  # Plot
  for(j in c(0,asfIDs))
  {
    fileName = paste0("SpentMediaGrowthCurveData\\ASF",i,"in",j,".tsv")
    curvedf= averageGrowthCurve(fileName)
    
    p = ggplot(curvedf, aes(time,ODave)) +
      geom_ribbon(fill="gray50", aes(ymin=ODave-ODstdError, ymax=ODave+ODstdError)) + 
      geom_line(color="black") +
      theme_bw() +
      ylim(0,maxOD) +
      xlim(0,72) +
      ylab(paste("ASF",i)) +
      theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank())
    print(p, vp = viewport(layout.pos.row = match(i,asfIDs), layout.pos.col = match(j,c(0,asfIDs))))
    
  }  
}

dev.off()

#-----------------------------------------------------------------
# Calculate area between growth curves (as a measure of growth inhibition)
#-----------------------------------------------------------------
library(zoo)

auc_inhibition_mat = matrix(0,ncol=7,nrow=7)
for(i in asfIDs)
{
  fileName = paste0("SpentMediaGrowthCurveData\\ASF",i,"in",0,".tsv")
  curve0_df = averageGrowthCurve(fileName)
  auc0 = sum(diff(curve0_df$time)*rollmean(abs(curve0_df$ODave-0.001),2))
    
  for(j in asfIDs)
  {
    fileName = paste0("SpentMediaGrowthCurveData\\ASF",i,"in",j,".tsv")
    curve_df= averageGrowthCurve(fileName)
    aucj = sum(diff(curve_df$time)*rollmean(abs(curve_df$ODave-0.001),2))
    inhibition = -(auc0-aucj)/auc0
    auc_inhibition_mat[match(i,asfIDs),match(j,asfIDs)] = inhibition
    
    # Visualize, just for QC
    print(paste(fileName,": ",inhibition))
    plot(curve0_df$time,curve0_df$ODave,col='red',type='l')
    lines(curve_df$time,curve_df$ODave,col='blue',type='l')

  }
}
colnames(auc_inhibition_mat) = paste0("Spent",asfIDs)
rownames(auc_inhibition_mat) = paste0("ASF",asfIDs)
write.table(auc_inhibition_mat, file = "auc_inhibitions_all_pairs.tsv", sep = "\t", quote=FALSE, row.names = TRUE, col.names = TRUE)


#-----------------------------------------------------------------
# Write AUC inhibition table in edge-wise format for plotting
#-----------------------------------------------------------------
fileConn = file("ASF_AUC_interaction_net.txt")
lineList = c()
for(i in asfIDs)
{
  for(j in asfIDs)
  {
    # Species i grown in spent j
    type = -1
    if(auc_inhibition_mat[match(i,asfIDs),match(j,asfIDs)] > 0)
    {
      type = 1
    }
    
    lineList = c(lineList, paste(j,abs(auc_inhibition_mat[match(i,asfIDs),match(j,asfIDs)]),type,'TRUE',i) )
  }
}

writeLines(lineList, fileConn)
close(fileConn)




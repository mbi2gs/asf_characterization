# Plotting NMR results
#
# Written by Matt Biggs, November 2015
library(ggplot2)
library(reshape2)

#------------------------------------------------
# Read in data
#------------------------------------------------
relativeIntensities = read.table('NMR_Spectra\\Integrals_UVA.csv',sep=',',header=T,row.names=NULL)
relativeIntensities = relativeIntensities[,-c(2,3,4,5,6)]
relativeIntensities[,1] = as.character(relativeIntensities[,1])

peakInfo = read.table('NMR_Spectra\\Peak_name_associations.txt',sep='\t',header=T,row.names=NULL)

sampleInfo = read.table('NMR_Spectra\\sample_info_7Oct15.tsv',sep='\t',header=F,row.names=NULL)
names(sampleInfo) = c("Grower","Medium","Rep","Round")
sampleInfo = sampleInfo[-253,]  # Last sample in "sample_info_7Oct15.tsv" is a blank buffer, so no need to keep it

# Read in 492 re-do data
relativeIntensities_492 = read.table('NMR_Spectra\\Integrals_492_redo.csv',sep=',',header=T,row.names=NULL)
relativeIntensities_492 = relativeIntensities_492[,-c(2,3,4,5,6)]
relativeIntensities_492[,1] = as.character(relativeIntensities_492[,1])

sampleInfo_492 = read.table('NMR_Spectra\\sample_info_492redo.tsv',sep='\t',header=F,row.names=NULL)
names(sampleInfo_492) = c("Grower","Medium","Rep","Round")

#------------------------------------------------
# Rename rows and calculate z-scores
# Keep one copy as individual data points ("RI_zs")
# and one copy as means (for the heatmap plot) ("RI_zm")
#------------------------------------------------
RI_zs = data.frame(matrix(ncol = 86, nrow = 0))
RI_zm = data.frame(matrix(ncol = 86, nrow = 0))
media = c(0,356,360,361,492,500,502,519)
species = c(356,360,361,492,500,502,519)
blankMedia = relativeIntensities[(sampleInfo$Grower==0 & sampleInfo$Medium==0),-1]
spentMedia = relativeIntensities[(sampleInfo$Grower>0 | sampleInfo$Medium>0),-1]
spentMediaInfo = sampleInfo[(sampleInfo$Grower>0 | sampleInfo$Medium>0),]
FMmeans = colMeans(blankMedia)
FMstdevs = apply(blankMedia,2,sd)

# Calculate z-scores
for(r in 1:nrow(spentMedia))
{
  # Rename 
  tmpName = paste0(spentMediaInfo$Grower[r],'in',spentMediaInfo$Medium[r])
  
  # Calculate z-scores
  tmpZscores = (spentMedia[r,] - FMmeans) / FMstdevs
  tmpdf = data.frame(matrix(ncol = 86, nrow = 1))
  tmpdf[1,2:86] = tmpZscores
  tmpdf[1,1] = tmpName
  RI_zs = rbind(RI_zs,tmpdf)
}
names(RI_zs) = c('Name',as.character(peakInfo$Name))
rowNames_RI_zs = RI_zs[,1]

# Calculate average z-scores
for(s in species)
{
  for(m in media)
  {   
    # Rename 
    tmpName = paste0(s,'in',m)
    
    # Calculate average z-scores for heatmap
    tmp = RI_zs[(spentMediaInfo$Grower==s & spentMediaInfo$Medium==m),2:86]
    if(nrow(tmp) > 4)
    {
      tmp = RI_zs[(spentMediaInfo$Grower==s & spentMediaInfo$Medium==m & spentMediaInfo$Round==2),2:86]
    }
    tmpMeans = colMeans(tmp)
    tmpdf = data.frame(matrix(ncol = 86, nrow = 1))
    tmpdf[1,2:86] = tmpMeans
    tmpdf[1,1] = tmpName
    RI_zm = rbind(RI_zm,tmpdf)
  }
}
names(RI_zm) = c('Name',as.character(peakInfo$Name))
tmpRnames = RI_zm[,1]
row.names(RI_zm) = tmpRnames
RI_zm = RI_zm[,-1]


# Rename rows and calculate z-scores for 492 re-do samples
# Keep one copy as individual data points ("RI_zs_492")
# and one copy as means (for the heatmap plot) ("RI_zm")
RI_zs_492 = data.frame(matrix(ncol = 86, nrow = 0))
RI_zm_492 = data.frame(matrix(ncol = 86, nrow = 0))
blankMedia_492 = relativeIntensities_492[(sampleInfo_492$Grower==0 & sampleInfo_492$Medium==0),-1]
spentMedia_492 = relativeIntensities_492[(sampleInfo_492$Grower>0 | sampleInfo_492$Medium>0),-1]
spentMediaInfo_492 = sampleInfo_492[(sampleInfo_492$Grower>0 | sampleInfo_492$Medium>0),]
FMmeans_492 = colMeans(blankMedia_492)
FMstdevs_492 = apply(blankMedia_492,2,sd)

# Calculate z-scores
for(r in 1:nrow(spentMedia_492))
{
  # Rename 
  tmpName = paste0(spentMediaInfo_492$Grower[r],'in',spentMediaInfo_492$Medium[r])
  
  # Calculate z-scores
  tmpZscores = (spentMedia_492[r,] - FMmeans_492) / FMstdevs # Use the stdevs from previous data so that results are comparable
  tmpdf = data.frame(matrix(ncol = 86, nrow = 1))
  tmpdf[1,2:86] = tmpZscores
  tmpdf[1,1] = tmpName
  RI_zs_492 = rbind(RI_zs_492,tmpdf)
}
names(RI_zs_492) = c('Name',as.character(peakInfo$Name))
rowNames_RI_zs_492 = RI_zs_492[,1]

# Calculate average z-scores
for(r in seq(1,nrow(spentMediaInfo_492),4))
{
  s = spentMediaInfo_492[r,1]
  m = spentMediaInfo_492[r,2]
  
  # Rename 
  tmpName = paste0(s,'in',m)
  
  # Calculate average z-scores for heatmap
  tmp = RI_zs_492[(spentMediaInfo_492$Grower==s & spentMediaInfo_492$Medium==m),2:86]
  tmpMeans = colMeans(tmp)
  tmpdf = data.frame(matrix(ncol = 86, nrow = 1))
  tmpdf[1,2:86] = tmpMeans
  tmpdf[1,1] = tmpName
  RI_zm_492 = rbind(RI_zm_492,tmpdf)
}
names(RI_zm_492) = c('Name',as.character(peakInfo$Name))
tmpRnames_492 = RI_zm_492[,1]
row.names(RI_zm_492) = tmpRnames_492
RI_zm_492 = RI_zm_492[,-1]

#------------------------------------------------
# Add the 492 re-do data to the previous data
#------------------------------------------------
RI_zs_v2 = RI_zs
RI_zm_v2 = RI_zm

for(m in c(0,360,361,492,500,502))
{
  tmpName = paste0(492,'in',m)
  
  tmpReplacement = RI_zm_492[row.names(RI_zm_492) == tmpName,1:85]
  RI_zm_v2[row.names(RI_zm_v2) == tmpName,1:85] = tmpReplacement
  
  tmpReplacement = RI_zs_492[(spentMediaInfo_492$Grower==492 & spentMediaInfo_492$Medium==m),1:85]
  RI_zs_v2[(spentMediaInfo$Grower==492 & spentMediaInfo$Medium==m),1:85] = tmpReplacement
}

#-----------------------------------------------------------------
# Plotting Metabolomics Z-Scores (only known metabolites) by species
# Compare all duplicated conditions to see 
# what, if anything, changed.
#-----------------------------------------------------------------
library(gridExtra)

species = c(356,360,361,492,500,502,519)
RI_zs_known = RI_zs[,-grep('Unknown',names(RI_zs))]
RI_zs_known = RI_zs_known[,-1]

RI_zs_492_known = RI_zs_492[,-grep('Unknown',names(RI_zs_492))]
RI_zs_492_known = RI_zs_492_known[,-1]

plotList = list()
for(s in species)
{
  plotList[[match(s,species)]] = list()
  for(m in media)
  {   
    tmpName = paste0(s,'in',m)
    tmpTitle = paste0('ASF',s,' in Spent',m)
    if(m == 0)
    {
      tmpTitle = paste0('ASF',s,' in Fresh Media')
    }
    
    # Old data
    tmp = RI_zs_known[rowNames_RI_zs == tmpName,]
    tmp[tmp > 20] = 20
    tmpdf = melt(tmp,id=c())
    
    p = ggplot(tmpdf, aes(x=factor(variable), y=value)) + 
        geom_hline(yintercept=0,color="red") +
        geom_boxplot(outlier.size = 0) + 
        geom_jitter(color='gray',alpha=0.7) +
        ylim(-6,20) +
        ylab("Z-score") +
        ggtitle(tmpTitle) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1,size=7), axis.title.x = element_blank())
    
    # 492 re-do data
    if(sum(rowNames_RI_zs_492 == tmpName) > 0)
    {
      tmp_492 = RI_zs_492_known[rowNames_RI_zs_492 == tmpName,]
      tmp_492[tmp_492 > 20] = 20
      tmpdf_492 = melt(tmp_492,id=c())
      
      p = p + geom_boxplot(data = tmpdf_492,outlier.size = 0,color="red",alpha=0.5) + 
          geom_jitter(data = tmpdf_492,color='red',alpha=0.5)
    }  
    #print(p)
    #ggsave(paste(tmpName,'.tiff'),width=6,height=4) 
    plotList[[match(s,species)]][[match(m,media)]] = p
  }
}

pdf("ASF_metabolite_zscores_plots.pdf", onefile = TRUE)
for(s in species)
{
  grid.arrange(plotList[[match(s,species)]][[1]],plotList[[match(s,species)]][[2]],plotList[[match(s,species)]][[3]],plotList[[match(s,species)]][[4]],ncol=2)
  grid.arrange(plotList[[match(s,species)]][[5]],plotList[[match(s,species)]][[6]],plotList[[match(s,species)]][[7]],plotList[[match(s,species)]][[8]],ncol=2)
}
dev.off()

#-----------------------------------------------------------------
# Plotting HeatMap of Metabolomics Z-Scores (all, including unknown, metabolites)
#-----------------------------------------------------------------
library(pheatmap)
library(RColorBrewer)

col.pal.fn = colorRampPalette(c('blue','red'),bias=1,space='rgb',interpolate ='linear')
col.pal = col.pal.fn(10)
fontsize = 15

RI_zm_thresholded = RI_zm_v2
RI_zm_thresholded[RI_zm_thresholded > 6] = 6
RI_zm_thresholded[RI_zm_thresholded < -6] = -6
RI_zm_thresholded[RI_zm_thresholded < 2 & RI_zm_thresholded > -2] = 0

hm.parameters <- list(RI_zm_thresholded, 
                      color = col.pal,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 30,
                      treeheight_col = 30,
                      fontsize = fontsize, fontsize_row = fontsize,
                      fontsize_col = fontsize,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      main = "",
                      clustering_method = "average",
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_rows = "euclidean", 
                      clustering_distance_cols = "euclidean",
                      filename = "relative_intensities_all_mets_heat_map.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)

#-----------------------------------------------------------------
# Plotting HeatMap of Metabolomics Z-Scores (only known metabolites)
#-----------------------------------------------------------------
RI_zm_thresh_known = RI_zm_thresholded
RI_zm_thresh_known = RI_zm_thresh_known[,-grep('Unknown',names(RI_zm_thresh_known))]

hm.parameters <- list(RI_zm_thresh_known, 
                      color = col.pal,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 30,
                      treeheight_col = 30,
                      fontsize = fontsize, fontsize_row = fontsize,
                      fontsize_col = fontsize,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      main = "",
                      clustering_method = "average",
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_rows = "euclidean", 
                      clustering_distance_cols = "euclidean",
                      filename = "relative_intensities_known_mets_heat_map.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)

#-----------------------------------------------------------------
# Plotting HeatMap of Metabolomics Z-Scores 
# (only known metabolites and single spent media)
#-----------------------------------------------------------------
RI_zm_tk_ss = RI_zm_thresh_known
RI_zm_tk_ss = RI_zm_tk_ss[grep('in0',row.names(RI_zm_thresh_known)),]

hm.parameters <- list(RI_zm_tk_ss, 
                      color = col.pal,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 30,
                      treeheight_col = 30,
                      fontsize = fontsize, fontsize_row = fontsize,
                      fontsize_col = fontsize,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      main = "",
                      clustering_method = "average",
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_rows = "euclidean", 
                      clustering_distance_cols = "euclidean",
                      filename = "relative_intensities_single_spent_known_mets_heat_map.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)


#-----------------------------------------------------------------
# Put all interaction data (growth and metabolomics) in a single plot
#-----------------------------------------------------------------

# Plot the underlying heatmap
growth_rate_inhibition = read.table('..\\GrowthCurves\\auc_inhibitions_all_pairs.tsv',sep='\t',header=T,row.names=1)
growth_rate_inhibition_df = melt(growth_rate_inhibition,id=c())
names(growth_rate_inhibition_df) = c("media","growth_rate_inhibition")
growth_rate_inhibition_df$growth_rate_inhibition = growth_rate_inhibition_df$growth_rate_inhibition + 1
growth_rate_inhibition_df$species = rep(row.names(growth_rate_inhibition),times=7)

p = ggplot(growth_rate_inhibition_df, aes(x=factor(media),y=factor(species))) +
  geom_tile(aes(fill=growth_rate_inhibition)) +
  scale_fill_gradient(low="black", high="white",na.value = "blue") +
  ylim(rev(levels(factor(growth_rate_inhibition_df$species)))) +
  theme(axis.text.y = element_text(size=8,color="black"), axis.text.x = element_text(angle = 0,vjust=0.5,hjust=0.5,size=8,color="black"),
        axis.title = element_blank(), axis.ticks = element_blank(), legend.position="none", panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA))

tiff("metabolomics_growth_compound_figure.tif", width=13, height=12, units="cm", res=600)
print(p)

# Create smaller metabolite plots
for(s in species)
{
  for(m in species)
  {
    # Name of SpentA
    tmpSpentA_name = paste0(m,'in',0)
    # Name of BinSpentA
    tmpBinSpentA_name = paste0(s,'in',m)
    # Make empty data frame
    xc = match(m,species)
    yc = match(s,species)
    tmp_df = data.frame(x=rep(1:36,times=3),y=rep(1:3,each=36),z=rep(NA,108))
    
    # Populate data frame with z-scores from spentA and BinSpentA
    tmpSpentA = RI_zm_thresh_known[row.names(RI_zm_thresh_known) == tmpSpentA_name,]
    tmpBinSpentA = RI_zm_thresh_known[row.names(RI_zm_thresh_known) == tmpBinSpentA_name,]
    tmp_df$z[tmp_df$y == 2] = as.numeric(tmpSpentA)
    tmp_df$z[tmp_df$y == 3] = as.numeric(tmpBinSpentA)
    
    p1 = ggplot(tmp_df, aes(x=x, y=y)) + 
      geom_tile(aes(fill=z)) + 
      scale_fill_gradient(limits=c(min(RI_zm_thresh_known),max(RI_zm_thresh_known)),low="blue", high="red", na.value = NA) +
      coord_polar() +
      theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), legend.position="none",
            panel.background = element_rect(fill = "transparent",colour = NA), 
            plot.background = element_rect(fill = "transparent",colour = NA))
    
    # Plot
    #print(p1)
    print(p1, vp=viewport(width=.24,height=.24,x=.316+(xc-1)*(0.1145),y=0.025+(7-yc)*(0.1217), just=c("right","bottom")) )
  }
}
dev.off()


#-----------------------------------------------------------------
# Count the various categories of metabolite profiles
#-----------------------------------------------------------------
# Convert thresholded relative intensity z-scores to simple up/down/no-change indicators
RI_zm_thresh2_known = RI_zm_thresh_known
RI_zm_thresh2_known = sign(RI_zm_thresh2_known)

# Growth is anything over 10%
gri = (growth_rate_inhibition + 1) > 0.1

# Counting
noChange = c(0,0,0,0,0,0)
lowerInDoubleSpent = c(0,0,0,0,0,0)
higherInDoubleSpent = c(0,0,0,0,0,0)

for(s in species)
{
  for(m in species)
  {
    if(m != s)
    {
      # Growth is anything over 10%
      grew = gri[match(s,species),match(m,species)]
        
      # Name of SpentA
      tmpSpentA_name = paste0(m,'in',0)
      # Name of BinSpentA
      tmpBinSpentA_name = paste0(s,'in',m)
      
      # Find the correct conditions
      tmpSpentA = RI_zm_thresh2_known[row.names(RI_zm_thresh2_known) == tmpSpentA_name,]
      tmpBinSpentA = RI_zm_thresh2_known[row.names(RI_zm_thresh2_known) == tmpBinSpentA_name,]
      
      # Infer profile classification
      add = 3
      if(grew)
      {
        add = 0
      }    
      
      # No change
      noChangeHigh = sum((tmpSpentA == tmpBinSpentA) & (tmpBinSpentA == 1))
      noChangeMed = sum((tmpSpentA == tmpBinSpentA) & (tmpBinSpentA == 0))
      noChangeLow = sum((tmpSpentA == tmpBinSpentA) & (tmpBinSpentA == -1))
      noChange[1+add] = noChange[1+add] + noChangeHigh
      noChange[2+add] = noChange[2+add] + noChangeMed
      noChange[3+add] = noChange[3+add] + noChangeLow
      
      # Lower in double spent media
      lowerInDoubleSpent_high2low = sum((tmpSpentA > tmpBinSpentA) & (tmpSpentA == 1) & (tmpBinSpentA == -1))
      lowerInDoubleSpent_high2med = sum((tmpSpentA > tmpBinSpentA) & (tmpSpentA == 1) & (tmpBinSpentA == 0))
      lowerInDoubleSpent_med2low = sum((tmpSpentA > tmpBinSpentA) & (tmpSpentA == 0) & (tmpBinSpentA == -1))
      lowerInDoubleSpent[1+add] = lowerInDoubleSpent[1+add] + lowerInDoubleSpent_high2low
      lowerInDoubleSpent[2+add] = lowerInDoubleSpent[2+add] + lowerInDoubleSpent_high2med
      lowerInDoubleSpent[3+add] = lowerInDoubleSpent[3+add] + lowerInDoubleSpent_med2low
      
      # Higher in double spent media
      higherInDoubleSpent_low2high = sum((tmpSpentA < tmpBinSpentA) & (tmpSpentA == -1) & (tmpBinSpentA == 1))
      higherInDoubleSpent_low2med = sum((tmpSpentA < tmpBinSpentA) & (tmpSpentA == -1) & (tmpBinSpentA == 0))
      higherInDoubleSpent_med2high = sum((tmpSpentA < tmpBinSpentA) & (tmpSpentA == 0) & (tmpBinSpentA == 1))
      higherInDoubleSpent[1+add] = higherInDoubleSpent[1+add] + higherInDoubleSpent_low2high
      higherInDoubleSpent[2+add] = higherInDoubleSpent[2+add] + higherInDoubleSpent_low2med
      higherInDoubleSpent[3+add] = higherInDoubleSpent[3+add] + higherInDoubleSpent_med2high
    }
  }
}

d = 7*6*36 # 7 species x 6 media (excluding own) x 36 known metabolites
print(paste("noChange total: ",sum(noChange)/d))
print(paste("lowerInDoubleSpent total: ",sum(lowerInDoubleSpent)/d))
print(paste("higherInDoubleSpent total: ",sum(higherInDoubleSpent)/d))

print(rbind(noChange,lowerInDoubleSpent,higherInDoubleSpent)/d)
print(rbind(noChange,lowerInDoubleSpent,higherInDoubleSpent))

#-----------------------------------------------------------------
# Look for emergent metabolites, where emergent metabolites 
# are those which are produced (or consumed) only when species A
# is growing in the spent media of another species
# Interesting (unique) emergent metabolites are those which only
# occur in one particular interaction.
#-----------------------------------------------------------------
# Convert thresholded relative intensity z-scores to simple up/down/no-change indicators
RI_zm_thresh2_known = RI_zm_thresh_known
RI_zm_thresh2_known = sign(RI_zm_thresh2_known)
howManyEmergent = 0
totalPossibleCases = 0

outputEmergentMetsFileName = "asfEmergentMets.txt"
write("Emergent Metabolites", file = outputEmergentMetsFileName, append = FALSE, sep = "\t")
write("Key: 2 = Didn't produce metabolite in spent (produced only in fresh)", file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
write("     3 = Newly consumed metabolite in spent (consumed only in spent)", file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
write("     5 = Not only didn't produce, but consumed metabolite in spent", file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
write("     7 = Newly produced in spent (produced only in spent)", file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
write("     11 = Didn't consume in spent (consumed only in fresh)", file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
write("     13 = Not only didn't consume, but produced metabolite in spent \n", file = outputEmergentMetsFileName, append = TRUE, sep = "\t")

for(s in species)
{
  # All samples generated by species "s"
  tmpMedia = RI_zm_thresh2_known[substring(row.names(RI_zm_thresh2_known),1,3) == as.character(s),]
  inFresh = tmpMedia[1,]
  inSpent = tmpMedia[-1,]
  emergents = inSpent
  
  # Normalize "double spent" media to originating spent media
  # Exclude own spent media
  for(r in 1:nrow(inSpent))
  {
    derivedFromSpentMediaName = paste0(substring(row.names(inSpent[r,]),6,8),'in0')
    originatingSpent = RI_zm_thresh2_known[row.names(RI_zm_thresh2_known) == derivedFromSpentMediaName,]
    ownSpentMediaName = paste0(s,'in0')
    ownSpent = RI_zm_thresh2_known[row.names(RI_zm_thresh2_known) == ownSpentMediaName,]
    doubleSpent = inSpent[r,]
    tmpDiff = sign(inSpent[r,] - originatingSpent)
    row.names(tmpDiff) = c("diff")    
    
    # Identify metabolites produced or consumed only in spent media
    emergent = tmpDiff
    row.names(emergent) = c("emergent")    
    for(m in 1:ncol(inSpent))
    {
       tmpEmergent = 0
       compare2ownSpent = tmpDiff[,m] - ownSpent[,m]
       if(compare2ownSpent == -1 & tmpDiff[,m] == 0 & ownSpent[,m] == 1)
       {
         if(doubleSpent[,m] < ownSpent[,m])
         {
           # Didn't produce m in spent (produced only in fresh)
           tmpEmergent = 2
         }
       }
       if(compare2ownSpent == -1 & tmpDiff[,m] == -1 & ownSpent[,m] == 0)
       {
         # Newly consumed m in spent
         tmpEmergent = 3
       }
       if(compare2ownSpent == -2 & tmpDiff[,m] == -1 & ownSpent[,m] == 1)
       {
         # Not only didn't produce, but consumed m in spent
         tmpEmergent = 5
       }
       if(compare2ownSpent == 1 & tmpDiff[,m] == 1 & ownSpent[,m] == 0)
       {
         # Newly produced in spent
         tmpEmergent = 7
       }
       if(compare2ownSpent == 1 & tmpDiff[,m] == 0 & ownSpent[,m] == -1)
       {
         if(doubleSpent[,m] > ownSpent[,m])
         {
           # Didn't consume in spent (consumed only in fresh)
           tmpEmergent = 11
         }
       }
       if(compare2ownSpent == 2 & tmpDiff[,m] == 1 & ownSpent[,m] == -1)
       {
         # Not only didn't consume, but produced m in spent
         tmpEmergent = 13
       }
       emergent[,m] = tmpEmergent
    }
    #print(rbind(originatingSpent,doubleSpent,tmpDiff,ownSpent,emergent))
    row.names(emergent) = row.names(inSpent[r,])
    emergents[r,] = emergent    
  }
  ownDoubleSpentMediaName = paste0(s,'in',s)
  emergents = emergents[row.names(emergents) != ownDoubleSpentMediaName,]
  howManyEmergent = howManyEmergent + sum(sign(emergents))
  totalPossibleCases = totalPossibleCases + 216
  write(paste0("\nASF",s), file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
  tmpColNames = names(emergents)
  tmpColNames[1] = paste0('\t',tmpColNames[1])
  write.table(emergents, file = outputEmergentMetsFileName, append = TRUE, sep = "\t", quote = FALSE, col.names = tmpColNames)
    
  # Identify metabolites only produced/consumed in one spent media condition
  write("Unique Emergent Metabolites:", file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
  columnSignedSums = colSums(sign(emergents))
  for(em in 1:ncol(emergents))
  {
    if(columnSignedSums[em] == 1)
    {
      write(paste0(names(emergents)[em],"\t",colSums(emergents)[em],"\t",row.names(emergents)[emergents[,em] > 0]), file = outputEmergentMetsFileName, append = TRUE, sep = "\t")
    }
  }  
  
}





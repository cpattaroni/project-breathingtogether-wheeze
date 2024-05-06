###############################################################################################
########################################### MICROBIOME ########################################
###############################################################################################

# Transform phyloseq object into dataframes
SampleDataframe <- function(phylo) {return(as(phyloseq::sample_data(phylo), 'data.frame'))}
OTUDataframe <- function(phylo) {return(as.matrix(data.frame(otu_table(phylo))))}

# Add ASV name based on the highest taxonomic class 
highest_taxClass <- function (x) {
  y <- ifelse(is.na(x$Species), # Check if the Species column is empty
              ifelse(is.na(x$Genus), # If yes, check if the Genus column is empty
                     ifelse(is.na(x$Family), # If yes, check if the Family column is empty
                            ifelse(is.na(x$Order), # If yes, check if the Order column is empty
                                   ifelse(is.na(x$Class), # If yes, check if the Class column is empty
                                          ifelse(is.na(x$Phylum), # If yes, check if the Phylum column is empty
                                                 paste0(x$Kingdom, '_', c(1:dim(x)[1])), # If yes, assign the Kingdom to the row
                                                 paste0(x$Phylum, '_', c(1:dim(x)[1])) ), # If not, assign the Phylum to the row
                                          paste0(x$Class, '_', c(1:dim(x)[1])) ), # If not, assign the Class to the row
                                   paste0(x$Order, '_', c(1:dim(x)[1])) ), # If not, assign the Order to the row
                            paste0(x$Family, '_', c(1:dim(x)[1])) ), # If not, assign the Family to the row
                     paste0(x$Genus, '_', c(1:dim(x)[1])) ), # If not, assign the Genus to the row
              paste0(x$Genus, '_', x$Species, '_', c(1:dim(x)[1])) ) # If Species is not empty, assign Genus_Species to the row
  return(y)
}

# Add taxonomy after tip_glom
addtaxonomyafteraggregation <- function(phylo) {
  # Convert tax_table to data frame and add ASV column
  taxtab <- as.data.frame(tax_table(phylo))
  taxtab$ASV <- rownames(taxtab)
  
  # Loop through each row and extract genus and species information from ASV column
  for (i in 1:nrow(taxtab)) {
    if (str_count(taxtab$ASV[i], "_") == 2) {
      taxtab$Genus[i] <- str_split(taxtab$ASV, "_")[[i]][1]
      taxtab$Species[i] <- str_split(taxtab$ASV, "_")[[i]][2]
    }
  }
  
  # Update the tax_table in the phylo object with the new taxonomy information
  tax_table(phylo) <- tax_table(as.matrix(taxtab))
  return(phylo)
}

# Drop samples with incomplete metadata
ps_drop_incomplete <- function(ps, vars = NA, verbose = FALSE) {
  # Convert phyloseq object to data.frame
  df <- SampleDataframe(ps)
  
  # If no specific variables are specified, use all sample variables in the phyloseq object
  if (identical(vars, NA)) vars <- phyloseq::sample_variables(ps)
  
  # Subset the data.frame to only include the specified sample variables
  df_sub <- df[, vars, drop = FALSE]
  
  # Remove rows (samples) with missing data
  df_sub <- df_sub[stats::complete.cases(df_sub), , drop = FALSE]
  
  # If verbose is not set to FALSE, print a message indicating the number of samples with missing data
  if (!isFALSE(verbose)) {
    incomplete <- nrow(df) - nrow(df_sub)
    if (incomplete > 0 || identical(verbose, "max")) {
      message("Dropping samples with missings: ", incomplete)
    }
  }
  
  # If verbose is set to "max", print a message indicating the number of missing values for each variable
  if (identical(verbose, "max")) {
    for (v in vars) {
      n_missings <- sum(is.na(df[[v]]))
      if (n_missings > 0) message(v, " has NAs: ", n_missings)
    }
  }
  
  # Get the names of the samples to keep (those without missing data)
  keepers <- rownames(df_sub)
  
  # Subset the phyloseq object to only include the samples without missing data
  phyloseq::prune_samples(samples = keepers, x = ps)
}

# Transform date into "hot/cold" months
TimeOfYear <- function(datevector) {
  months <- month(as.Date(datevector))
  season <- ifelse(months >= 5 & months <= 10, 'Warm', 'Cold')
  return(season)
}

# Plot contaminants
PlotContaminants <- function(phylo, contam, title) {
  # Transform the sample counts by setting any abundance greater than 0 to 1
  pa <- transform_sample_counts(phylo, function(abund) 1 * (abund > 0))
  
  # Subset the resulting object into two phyloseq objects: pa_neg containing negative control samples and pa_pos containing all other samples.
  pa_neg <- prune_samples(phyloseq::sample_data(pa)$NegCtrl == TRUE, pa)
  pa_pos <- prune_samples(phyloseq::sample_data(pa)$NegCtrl == FALSE, pa)
  
  # Create a data frame df_pa that contains the sum of taxonomic abundances for each sample type (positive and negative controls) and contaminant.
  df_pa <- data.frame(pa_pos = taxa_sums(pa_pos), pa_neg = taxa_sums(pa_neg), contaminant = contam$contaminant)
  
  # Create a scatter plot using the ggplot2 package.
  contam_plot <- ggplot(df_pa, aes(x = pa_neg, y = pa_pos, col = contaminant, fill = contaminant)) +
    # Add points to the plot
    geom_point(size = 2, shape = 21, position = position_dodge2(0.2), stroke = 0.25) +
    # Set the title and axis labels
    labs(title = title,
         x = 'Prevalence in negative controls (Swab, extraction and PCR)',
         y = 'Prevalence in samples') +
    # Set the fill colors for the contaminants
    scale_fill_manual(values = adjustcolor(c('darkseagreen', 'indianred'), 0.5)) +
    # Set the point colors for the contaminants
    scale_color_manual(values = adjustcolor(c('darkseagreen', 'indianred'), 1)) +
    # Remove the legend
    theme(legend.position = c(1, 0), legend.box.background=element_rect(), text=element_text(size=12))
  
  # Return the ggplot object containing the contaminant plot
  return(contam_plot)
}

# Boxplots for metadata
BoxplotMetadataGroup <- function(phylo, parameter, title, yaxis, palette) {
  # Drop samples with missing data for the specified parameter
  phylo <- ps_drop_incomplete(phylo, parameter)
  # Convert to a SampleDataframe
  df <- SampleDataframe(phylo)
  # Find the maximum value of the parameter for setting the y-axis limit
  max <- max(as.numeric(df[, parameter]))
  
  # Create the plot using ggplot2
  plot <- ggplot(df, aes(x=Group, y=as.numeric(df[, parameter]))) +
    # Add the boxplot (white fill with colored outline)
    geom_boxplot(aes(group=Group), width=0.5, outlier.shape=NA, col=adjustcolor(palette, alpha=1), fill=adjustcolor('white', alpha=1)) +
    # Add the boxplot with transparent fill
    geom_boxplot(aes(group=Group), width=0.5, outlier.shape=NA, col=adjustcolor(palette, alpha=1), fill=adjustcolor(palette, alpha=0.25)) +
    # Add individual points (colored by group)
    geom_point(aes(y=as.numeric(df[, parameter]), fill=Group, color=Group),
               size=2, shape=21, position=position_dodge2(0.2), stroke=0.25) +
    # Set the fill and color scales
    scale_fill_manual(values=adjustcolor(palette, alpha=0.5)) +
    scale_color_manual(values=adjustcolor(palette, alpha=1)) +
    # Add significance bars between groups
    geom_signif(comparisons=list(c('Yes', 'No'), step_increase=0.2), map_signif_level=c('***'=0.001, '**'=0.01, '*'=0.05)) +
    # Set the x-axis labels
    scale_x_discrete(labels=c(paste0('Controls', '\n', 'n = ', table(df$Group)[1]), paste0('Wheezers', '\n', 'n = ', table(df$Group)[2]))) +
    # Set the plot title and axis labels
    ggtitle(title) + xlab('') + ylab(title) + ylab(yaxis) +
    # Remove the legend and set the y-axis limit
    theme(legend.position = c(1, 0), legend.box.background=element_rect(), text=element_text(size=12)) + ylim(0, (max + max * 0.1))
  
  return(plot)
}

# Boxplots for metadata
BoxplotMetadataParameters <- function(phylo, parameter.num, parameter.bin, title, yaxis, palette, paletteboxplot, labels){
  # Drop incomplete rows for binary parameter
  phylo <- ps_drop_incomplete(phylo, parameter.bin)
  
  # Create data frame from phyloseq object
  df <- SampleDataframe(phylo)
  
  # Convert numeric parameter to numeric type and binary parameter to factor type
  df$parameter.num <- as.numeric(df[, parameter.num])
  df$parameter.bin <- as.factor(df[, parameter.bin])
  
  # Get maximum value for numeric parameter for ylim
  max <- max(df$parameter.num)
  
  # Create plot with ggplot
  plot <- ggplot(df, aes(x = parameter.bin, y = parameter.num)) +
    # Add boxplots
    geom_boxplot(width = 0.5, outlier.shape = NA, col = adjustcolor(paletteboxplot, alpha = 1),
                 fill = adjustcolor('white', alpha = 1)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, col = adjustcolor(paletteboxplot, alpha = 0.25),
                 fill = adjustcolor(paletteboxplot, alpha = 0.1)) +
    # Add points with jitter and fill/color by group
    geom_point(aes(x = parameter.bin, y = parameter.num, color = parameter.bin, fill = parameter.bin),
               size = 2, shape = 21, position = position_jitter(0.1), stroke = 0.25) +
    # Add color and fill scales
    scale_fill_manual(values = adjustcolor(palette, alpha = 0.5)) +
    scale_color_manual(values = adjustcolor(palette, alpha = 1)) +
    # Add significance bars
    geom_signif(comparisons = list(c(levels(df$parameter.bin)[1], levels(df$parameter.bin)[2]), step_increase = 0.2),
                map_signif_level = c('***' = 0.001, '**' = 0.01, '*' = 0.05)) +
    # Add plot labels and formatting
    ggtitle(title) + xlab('') + ylab(title) + ylab(yaxis) + theme(legend.position = 'none') +
    scale_x_discrete(labels = c(paste0(levels(df$parameter.bin)[1], '\n', 'n = ', table(df$parameter.bin)[1]),
                                paste0(levels(df$parameter.bin)[2], '\n', 'n = ', table(df$parameter.bin)[2]))) +
    ylim(0, (max + max * 0.1)) + theme(legend.position = c(1, 0), legend.box.background=element_rect(), text=element_text(size=12))
  
  return(plot)
}

# Age plots for metadata
BoxplotMetadataAge <- function(phylo, parameter, title, xaxis, yaxis, palette) {
  # Drop incomplete rows based on the given parameter
  phylo <- ps_drop_incomplete(phylo, parameter)
  
  # Convert data to a sample dataframe
  df <- SampleDataframe(phylo)
  
  # Remove unused factor levels from the dataframe
  df <- droplevels(df)
  Medians <- data.frame(Median=c(median(df[df$TimepointSimplified == 0, parameter]), median(df[df$TimepointSimplified == 3, parameter]),
                                 median(df[df$TimepointSimplified == 6, parameter]), median(df[df$TimepointSimplified == 12, parameter])),
                        TimePoint=c(0, 3, 6, 12))
  # Create the plot with ggplot2
  plot <- ggplot(df, aes(x = TimepointSimplified, y = as.numeric(df[, parameter]))) +
    geom_line(data=Medians, aes(x = TimePoint, y = Median), color = adjustcolor(palette, alpha = 0.5), size = 5) +
    geom_boxplot(aes(group = as.factor(TimepointSimplified)), fill = "white", color = "white", width = 2, outlier.shape = NA) +
    geom_boxplot(aes(group = as.factor(TimepointSimplified), fill = as.factor(TimepointSimplified), color = as.factor(TimepointSimplified)), width = 2, outlier.shape = NA) +
    # Add points for each group with jitter
    geom_point(aes(fill = as.factor(TimepointSimplified), color = as.factor(TimepointSimplified)), size = 2, shape = 21, stroke = 0.25,
               position = position_jitterdodge(jitter.width = 2, jitter.height = 0.02)) +
    # Set the fill and color scales with the given color palette
    scale_fill_manual(values = adjustcolor(palette, alpha = 0.5))  +
    scale_color_manual(values = adjustcolor(palette, alpha = 1)) +
    # Set the x-axis breaks to show the time points
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15)) +
    ylim(c(0, max(as.numeric(df[, parameter])+as.numeric(df[, parameter])*0.1))) +
    # Add the title, x-axis label, and y-axis label
    ggtitle(title) + xlab(xaxis) + ylab(yaxis) +
    # Add a text annotation with the sample sizes for each group
    annotate('text', x = max(df$TimepointMonth) - 2, y = min(as.numeric(df[, parameter])),
             label = paste0('1 week n = ', table(df$TimepointSimplified)[1], '\n', 
                            '3 months n = ', table(df$TimepointSimplified)[2], '\n',
                            '6 months n = ', table(df$TimepointSimplified)[3], '\n',
                            '12 months n = ', table(df$TimepointSimplified)[4])) +
    # Remove the legend
    theme(legend.position = "none", text=element_text(size=12))
  
  # Return the plot
  return(plot)
}

# Ordination age
PlotOrdinationAge <- function(phylo, title, palette){
  
  # Perform ordination analysis and add the PCoA scores to the sample metadata
  ordi <- ordinate(phylo, 'PCoA', 'wunifrac')
  phyloseq::sample_data(phylo)$PCoA1 <- ordi$vectors[,1]
  phyloseq::sample_data(phylo)$PCoA2 <- ordi$vectors[,2]
  phylo <- SampleDataframe(phylo)
  
  # Define axis labels using the relative eigenvalues of the PCoA
  axisx <- paste0('PCoA1 ',round(ordi$values$Relative_eig*100)[1],'%')
  axisy <- paste0('PCoA1 ',round(ordi$values$Relative_eig*100)[2],'%')
  
  # Create the plot using ggplot2
  plot <- ggplot(phylo, aes(x=PCoA1, y=PCoA2)) +
    
    # Add x and y side densities, filled by timepoint and colored by timepoint
    geom_xsidedensity(aes(y=after_stat(density), fill=Timepoint, color=Timepoint), size=0.25) +
    geom_ysidedensity(aes(x=after_stat(density), fill=Timepoint, color=Timepoint), size=0.25) +
    
    # Add confidence ellipses, filled by timepoint and alpha=0.2
    stat_ellipse(aes(fill=Timepoint), geom='polygon', alpha=0.2, type='t', level=0.75) + scale_shape_identity() +
    
    # Add points filled and colored by timepoint, with size=2 and stroke=0.25
    #geom_point(aes(fill=Timepoint, color=Timepoint), shape=ifelse(phylo$Group=='Yes' & phylo$TimepointFactor=='T4-1year', 25, 21), size=2, stroke=0.25) + 
    geom_point(aes(fill=Timepoint, color=Timepoint), size=2, shape=21, stroke=0.25) +
    
    # Define color and fill scales using the palette argument
    scale_color_manual(values=adjustcolor(palette, alpha=1)) +
    scale_fill_manual(values=adjustcolor(palette, alpha=0.5)) + 
    
    # Define the legend position and axis labels
    theme(legend.position='bottom') +
    xlab(axisx) + ylab(axisy) + ggtitle(title) +
    
    # Define theme elements such as panel scale and axis ticks
    theme(ggside.panel.scale = 1/4, axis.text.x=element_blank(), 
          axis.text.y=element_blank(), axis.ticks=element_blank(), legend.position = c(1, 0), legend.box.background=element_rect(), text=element_text(size=12))
  
  # Return the plot object
  return(plot)
}

# Calculate distance transitions
UnifracTimepoint <- function(phylo, title) {
  
  # Extract sample data and keep longitudinal samples with more than 2 timepoints
  samples.data <- SampleDataframe(phylo)
  keep.longitudinal <- names(which(rowSums(table(samples.data$Patient, samples.data$TimepointFactor)) > 2))
  phylo <- prune_samples(phyloseq::sample_data(phylo)$Patient %in% keep.longitudinal, phylo)
  
  # Calculate distance matrix and prepare sample data for plotting
  distance.matrix <- as.matrix(phyloseq::distance(phylo, 'wunifrac'))
  samples.data <- SampleDataframe(phylo)
  samples.data$barcode <- rownames(samples.data)
  
  # Create a data frame to store the pairwise distances between time points
  transition.df <- data.frame(
    'T1:1week-3months' = rep(NA, length(keep.longitudinal)),
    'T2:3months-6months' = rep(NA, length(keep.longitudinal)),
    'T3:6months-1year' = rep(NA, length(keep.longitudinal))
  )
  
  # Add patient IDs and group information to the transition data frame
  transition.df$Patient <- keep.longitudinal
  transition.df$Group <- ifelse(
    keep.longitudinal %in% samples.data[samples.data$Group == 'Yes',]$Patient == TRUE,
    'Yes',
    'No'
  )
  
  # Loop through all longitudinal samples and calculate pairwise distances
  for (i in 1:length(keep.longitudinal)) {
    T1 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T1-1week')]
    T2 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T2-3months')]
    T3 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T3-6months')]
    T4 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T4-1year')]
    transition.df$T1.1week.3months[i] <- ifelse(
      c(T1,T2) %in% colnames(distance.matrix),
      distance.matrix[T1,T2],
      NA
    )
    transition.df$T2.3months.6months[i] <- ifelse(
      c(T2,T3) %in% colnames(distance.matrix),
      distance.matrix[T2,T3],
      NA
    )
    transition.df$T3.6months.1year[i] <- ifelse(
      c(T3,T4) %in% colnames(distance.matrix),
      distance.matrix[T3,T4],
      NA
    )
  }
  
  # Melt the transition data frame to long format for plotting
  transition.melt <- melt(transition.df)
  
  # Remove missing values from the data frame
  transition.melt <- transition.melt[!is.na(transition.melt$value),]

  # Return plot
  return(transition.melt)}

# Plot distance transitions
PlotUnifracTimepoint <- function(phylo, title) {
  
  # Extract sample data and keep longitudinal samples with more than 2 timepoints
  samples.data <- SampleDataframe(phylo)
  keep.longitudinal <- names(which(rowSums(table(samples.data$Patient, samples.data$TimepointFactor)) > 2))
  phylo <- prune_samples(phyloseq::sample_data(phylo)$Patient %in% keep.longitudinal, phylo)
  
  # Calculate distance matrix and prepare sample data for plotting
  distance.matrix <- as.matrix(phyloseq::distance(phylo, 'wunifrac'))
  samples.data <- SampleDataframe(phylo)
  samples.data$barcode <- rownames(samples.data)
  
  # Create a data frame to store the pairwise distances between time points
  transition.df <- data.frame(
    'T1:1week-3months' = rep(NA, length(keep.longitudinal)),
    'T2:3months-6months' = rep(NA, length(keep.longitudinal)),
    'T3:6months-1year' = rep(NA, length(keep.longitudinal))
  )
  
  # Add patient IDs and group information to the transition data frame
  transition.df$Patient <- keep.longitudinal
  transition.df$Group <- ifelse(
    keep.longitudinal %in% samples.data[samples.data$Group == 'Yes',]$Patient == TRUE,
    'Yes',
    'No'
  )
  
  # Loop through all longitudinal samples and calculate pairwise distances
  for (i in 1:length(keep.longitudinal)) {
    T1 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T1-1week')]
    T2 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T2-3months')]
    T3 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T3-6months')]
    T4 <- samples.data$barcode[which(samples.data$Patient == keep.longitudinal[i] & samples.data$TimepointFactor == 'T4-1year')]
    transition.df$T1.1week.3months[i] <- ifelse(
      c(T1,T2) %in% colnames(distance.matrix),
      distance.matrix[T1,T2],
      NA
    )
    transition.df$T2.3months.6months[i] <- ifelse(
      c(T2,T3) %in% colnames(distance.matrix),
      distance.matrix[T2,T3],
      NA
    )
    transition.df$T3.6months.1year[i] <- ifelse(
      c(T3,T4) %in% colnames(distance.matrix),
      distance.matrix[T3,T4],
      NA
    )
  }
  
  # Melt the transition data frame to long format for plotting
  transition.melt <- melt(transition.df)
  
  # Remove missing values from the data frame
  transition.melt <- transition.melt[!is.na(transition.melt$value),]
  
  # Plot
  plot <- ggplot(transition.melt, aes(x=variable, y=value, fill=variable, color=variable)) +
    geom_boxplot(aes(group = variable), fill = "white", color = "white", width = 0.5, outlier.shape = NA) +
    geom_boxplot(aes(group = variable), width = 0.5, outlier.shape = NA) +
    geom_point(size = 2, shape = 21, stroke = 0.25, position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +
    scale_fill_manual(values=adjustcolor(custom_palette_age, alpha=0.5)) +
    scale_color_manual(values=adjustcolor(custom_palette_age, alpha=1)) +
    geom_signif(comparisons=list(c('T1.1week.3months','T2.3months.6months')),
                map_signif_level=c('***'=0.001, '**'=0.01, '*'=0.05), 
                col='black') +
    scale_x_discrete(labels=c(paste0('1 week - 3 months', '\n', 'n = ', table(transition.melt$variable)[1]),
                              paste0('3 months - 6 months', '\n', 'n = ', table(transition.melt$variable)[2]),
                              paste0('6 months - 1 year', '\n', 'n = ', table(transition.melt$variable)[3]))) +
    ylim(0,(max(transition.melt$value)+max(transition.melt$value)*0.1)) +
    ggtitle(title) + xlab('') + ylab('Δ distance between pairs (wunifrac)') + theme(legend.position='none', text=element_text(size=12))
  
  # Return plot
  return(plot)}

# ANOSIM
ANOSIMunifrac <- function(phylo, design) {
  
  # Set the random seed to 2
  set.seed(2)
  
  # Set the distance metric to 'wunifrac'
  distance <- 'wunifrac'
  
  # Perform an ANOSIM test using the specified design and distance metric
  test <- adonis2(formula = as.formula(paste('phyloseq::distance(phylo, method = distance) ~', design)),
                  data = as.data.frame(as.matrix(phyloseq::sample_data(phylo))), 
                  distance = distance, dfun = vegdist, na.action = na.omit)
  
  # Return the test results
  return(test)
}

# ANOSIM
ANOSIMunifrac <- function(phylo, design) {
  
  # Set the random seed to 2
  set.seed(2)
  
  # Set the distance metric to 'wunifrac'
  distance <- 'wunifrac'
  
  # Perform an ANOSIM test using the specified design and distance metric
  test <- adonis2(formula = as.formula(paste('phyloseq::distance(phylo, method = distance) ~', design)),
                  data = as.data.frame(as.matrix(phyloseq::sample_data(phylo))), 
                  distance = distance, dfun = vegdist, na.action = na.omit)
  
  # Return the test results
  return(test)
}

# DA testing with LINDA
DAtesting <- function(phylo, variables, model, level) {
  require(phyloseq)
  # Remove samples with missing data in the specified parameter
  phylo <- ps_drop_incomplete(phylo, variables)
  
  # Convert phyloseq object to meco object
  meco <- phyloseq2meco(phylo)
  
  # Tidy the taxonomic table
  meco$tax_table <- tidy_taxonomy(meco$tax_table)

  # Perform classification using the LINDA method
  meco <- trans_diff$new(dataset = meco, method = 'linda', group = model, taxa_level = level)
  
  # Return the results of classification as a data frame
  return(meco)
}

# DA testing: LINDA barplot
DAtestingPlotBar <- function(meco, comparison, palette, number_features, title) {
  
  # Filter the dataframe based on the predefined comparison and significance criteria
  df <- as.data.frame(meco$res_diff)
  filtered_df <- df[df$Comparison == comparison & df$Significance != "ns",]
  filtered_df <- filtered_df[order(filtered_df$P.adj), ]  # Sort by smallest P.adj
  
  # Select the top number_features
  filtered_df <- filtered_df[1:number_features, ]
  
  # Reorder based on log2FoldChange
  filtered_df <- filtered_df[order(filtered_df$log2FoldChange),]
  
  # Add direction
  filtered_df$Direction <- ifelse(filtered_df$log2FoldChange<0, strsplit(comparison, " - ")[[1]][2], strsplit(comparison, " - ")[[1]][1])
  
  # Calculate transparency based on log10 of P.adj
  filtered_df$Transparency <- scales::rescale((1/log10(filtered_df$P.adj)), c(0.9, 1))
  
  # Rename taxa
  filtered_df$Taxa <- gsub(".*\\|a__", "", filtered_df$Taxa)
  
  # Create the barplot
  p <- ggplot(filtered_df, aes(x = factor(Taxa, levels = unique(filtered_df$Taxa)), y = log2FoldChange, 
                               color = Direction, fill = Direction, alpha = Transparency)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    coord_flip() +
    labs(x = "ASV", y = "log2 Fold Change") +
    ggtitle(title) +
    theme(legend.position = c(1,0),
          legend.box.background = element_rect(), 
          text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm")) +
    guides(alpha = "none")
  return(p)
}

# Record DA testing results for some covariates
DAfeaturesCounts <- function(phylo, paramlist, level = 'ASV') {

  # Create a data frame to store the number of DA features
  df <- data.frame(nbfeatures = rep(NA, length(paramlist)))
  
  # Loop through each parameter in the paramlist
  for (i in 1:length(paramlist)) {
    da <- DAtesting(phylo, variables = paramlist[i], model = paste0("~", paramlist[i]),  level = level)
    df[i, ] <- length(which(da$res_diff$Significance!="ns"))}
  # Add the parameter list as a column to the data frame
  df$Covariate <- paramlist
    # Return the data frame
  return(df)}

# Redundancy analysis for variables selection
dbRDAselectionUnifrac <- function(phylo, variables, distance='wunifrac', direction='both') {
  set.seed(2)
  # Remove samples with missing values for the specified variables
  phylo <- ps_drop_incomplete(phylo, var=variables, verbose=TRUE)
  # Extract the metadata table from the phyloseq object
  meta <- SampleDataframe(phylo)
  # Perform a dbRDA with a constant term only
  dbrda.0 <- dbrda(phyloseq::distance(phylo, method=distance) ~ 1, data=meta)
  # Perform a dbRDA with all specified variables
  dbrda.all <- dbrda(as.formula(paste('phyloseq::distance(phylo, method=distance) ~', paste(variables, collapse=" + "))), data=meta)
  # Perform stepwise variable selection using the ordiR2step function
  select <- ordiR2step(dbrda.0, scope=formula(dbrda.all), adjR2thresh=RsquareAdj(dbrda.all)$adj.r.squared)
  # Make a copy of the selected model and adjust the p-values for multiple comparisons
  select.adj <- select
  select.adj$anova$`Pr(>F)` <- p.adjust(select$anova$`Pr(>F)`, method='fdr', n=length(variables))
  # Return the results of the adjusted model and the unadjusted model
  return(select.adj$anova)
}

# Plot presence of viruses with given metadata
plotanyViruses <- function(phylo, parameter, title){
  phylo <- ps_drop_incomplete(phylo, var=parameter, verbose=TRUE)
  vir.data <- SampleDataframe(phylo)
  vir.data$parameter <- vir.data[parameter][,1]
  p <- ggplot(vir.data %>% dplyr::count(AnyVirus, parameter) %>% mutate(pct = n/sum(n)), 
              aes(parameter, n, fill=AnyVirus)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 5) +
    labs(title = title, x = parameter, y = "% of samples within group") +
    scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
    scale_x_discrete(labels = c(paste0(names(table(vir.data[parameter]))[1], '\n', 'n = ', table(vir.data[parameter])[1]),
                                paste0(names(table(vir.data[parameter]))[2], '\n', 'n = ', table(vir.data[parameter])[2]))) +
    scale_fill_manual(values = c("#BF812D", "#DFC27D", "#F5F5F5")) +
    theme(legend.position = c(1,0),
          legend.box.background = element_rect(),
          text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm")) +
    guides(alpha = "none")
  return(p)
}

# Plot number of DE/DA features
PlotfeaturesCounts <- function(countmatrix, color, title){
  p <- ggplot(countmatrix, aes(x=Newname, y=nbfeatures, color=color, fill=color)) +
    geom_bar(stat="identity") + coord_flip() +
    scale_fill_manual(values=scales::alpha(color, 0.75)) +
    scale_color_manual(values=color) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
    labs(y="Number of features", x="Covariate") +
    ggtitle(title) +
    theme(legend.position = "none",
          text = element_text(size = 12)) +
    guides(alpha = "none")
  return(p)
}

###############################################################################################
###################################### GENE EXPRESSION ########################################
###############################################################################################

# Longitudinal boxplots ############ DELETE ############
LongitudinalBoxplots <- function(dge, Gene, title) {
  # Prepare dataframe
  data <- data.frame(Abundance=dge$E[rownames(dge$E) %in% Gene,], Group=dge$targets$Group, AgeFactor=dge$targets$TimepointSimplified)
  
  # Statistics for each timepoint
  p1 <- wilcox.test(data[data$AgeFactor==0,]$Abundance ~ data[data$AgeFactor==0,]$Group)
  p2 <- wilcox.test(data[data$AgeFactor==3,]$Abundance ~ data[data$AgeFactor==3,]$Group)
  p3 <- wilcox.test(data[data$AgeFactor==6,]$Abundance ~ data[data$AgeFactor==6,]$Group)
  p4 <- wilcox.test(data[data$AgeFactor==12,]$Abundance ~ data[data$AgeFactor==12,]$Group)
  
  # Add stats as labels
  labels <- c(round(p1$p.value, digits=2), round(p2$p.value, digits=2), round(p3$p.value, digits=2), round(p4$p.value, digits=2))
  
  # Create the ggplot object with basic layers
  p <- ggplot(data, aes(x=AgeFactor, y=Abundance, fill=Group)) +
    geom_smooth(aes(color=Group, fill=Group), method='lm', level=0.75) +
    geom_point(aes(color=Group), size=2, shape=21, position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0.1), stroke=0.25) +
    geom_boxplot(aes(group=interaction(Group, AgeFactor), color=Group, fill=Group), outlier.shape=NA, width=0.5) +
    geom_smooth(aes(group = Group, color = Group, fill = Group), method = 'lm', level = 0.8) +
    # Add labels at specific x and y coordinates for each label
    geom_label(label=labels[1], x=0, y=max(data$Abundance), color="black", fill="white", label.size=0) +
    geom_label(label=labels[2], x=3, y=max(data$Abundance), color="black", fill="white", label.size=0) +
    geom_label(label=labels[3], x=6, y=max(data$Abundance), color="black", fill="white", label.size=0) +
    geom_label(label=labels[4], x=12, y=max(data$Abundance), color="black", fill="white", label.size=0) +
    # Set the color and fill scales
    scale_fill_manual(values=adjustcolor(custom_palette_2groups, alpha=0.5)) +
    scale_color_manual(values=adjustcolor(custom_palette_2groups, alpha=1)) +
    # Set the x-axis breaks
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18)) + 
    # Set the plot title and axis labels
    ggtitle(title) + xlab('Age (months)') + ylab(paste(Gene, "(voom normalised)")) +
    # Remove the legend
    theme(legend.position='none')
  
  # Return the ggplot object
  return(p)
}

# Count number of DE genes for some parameters
DEfeaturesCounts <- function(dge, paramlist) {
  
  # Create a data frame to store the number of DA features
  df <- data.frame(nbfeatures = rep(NA, length(paramlist)))
  
  # Loop through each parameter in the paramlist
  for (i in 1:length(paramlist)) {
    # Remove samples with NA data
    dge <- dge[,which(!is.na(dge$targets[paramlist[i]]))]
    de.results <- topTable(eBayes(lmFit(dge$E, model.matrix(formula(paste0("~", paramlist[i])), data=dge$targets))), coef=2, number="Inf")
    df[i, ] <- length(which(de.results$adj.P.Val<FDR_param))}
  
  # Add the parameter list as a column to the data frame
  df$Covariate <- paramlist
  
  # Return the data frame
  return(df)}

# Plot Volcano 
PlotVolcanoGenes <- function(DE_matrix, Direction, palette, title) {
  
  plot <- ggplot(data=DE_matrix, aes(x=logFC, y=-log10(adj.P.Val), label=label)) + # Create a ggplot object with log2FC on x-axis and -log10(p-value) on y-axis, with label as the label of each point
    
    # Add horizontal line at -log10(FDR_param) and vertical lines at -FC_param and FC_param
    geom_hline(yintercept=-log10(FDR_param), col='grey', linetype="dashed") + 
    geom_vline(xintercept=c(-FC_param, FC_param), col='grey', linetype="dashed") + 
    
    # Add points and labels for significant differentially expressed genes
    geom_point(aes(fill=Direction, color=Direction), shape=ifelse(DE_matrix$Immune=="Yes", 22, 21), size=ifelse(DE_matrix$Immune=="Yes", 3, 2), stroke=0.25) + 
    geom_label_repel(aes(fill=Direction), color='black', label.size=0, label.r=0.5, size=4, label.padding=0.25, max.overlaps=30, 
                     segment.size=0.25, segment.alpha=0.95, segment.color="grey") + 
    
    # Add manual color scales for fill and color
    scale_fill_manual(values=adjustcolor(palette, alpha=0.5)) +
    scale_color_manual(values=adjustcolor(palette, alpha=1)) +
    
    # Add title and move legend to bottom
    ggtitle(title) + theme(legend.position = c(1, 0), legend.box.background=element_rect(), text=element_text(size=12)) +
    
    # Set x-axis limits to be symmetric around 0
    xlim(c(-max(abs(DE_matrix$logFC)), max(abs(DE_matrix$logFC))))
  
  return(plot)
}

# Transform non-numeric data into numeric
SelectIntoNumeric <- function(metadata, variables){
  # Subset the data frame to include only the selected variables
  metadata <- metadata[,colnames(metadata) %in% variables]
  
  # Convert any character columns to factor columns
  cols <- which(sapply(metadata, is.character)) 
  if (length(cols>=1)) {
    metadata[cols] <- lapply(metadata[cols], as.factor)
  }
  
  # Convert any non-numeric or non-integer columns to numeric
  cols <- which(!sapply(metadata,class)%in% c('numeric', 'integer'))
  if (length(cols>=1)) {
    cols.factor <- which(sapply(metadata,class)=='factor')
    metadata[cols] <- lapply(metadata[cols], as.numeric)
  }
  
  # Reverse the order of the selected variables in the data frame
  metadata <- metadata[,rev(variables)]
  
  # Return the modified data frame
  return(metadata)
}

# Calculate distance between timepoints
DistanceTimepoint <- function(se, title){
  
  # Get the sample metadata
  samples.data <- se$targets
  
  # Keep only samples from patients with more than two timepoints
  keep.longitudinal <- names(which(rowSums(table(samples.data$Patient, samples.data$TimepointFactor))>2))
  se <- se[,se$targets$Patient %in% keep.longitudinal] 
  
  # Calculate the distance matrix between samples
  distance.matrix <- as.matrix(dist(t(se$E), method='maximum'))
  
  # Add a "barcode" column to the metadata
  samples.data$barcode <- rownames(samples.data)
  
  # Initialize a data frame to store the distance transitions
  transition.df <- data.frame('T1:1week-3months'=rep(NA, length(keep.longitudinal)),
                              'T2:3months-6months'=rep(NA, length(keep.longitudinal)),
                              'T3:6months-1year'=rep(NA, length(keep.longitudinal)))
  transition.df$Patient <- keep.longitudinal
  
  # Assign each patient to a "Group" based on whether they belong to a certain subgroup
  transition.df$Group <- ifelse(keep.longitudinal %in% samples.data[samples.data$Group=='Yes',]$Patient==TRUE, 'Yes', 'No')
  
  # Loop over each patient and calculate the distances between their timepoints
  for (i in 1:length(keep.longitudinal)) {
    T1 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T1-1week')]
    T2 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T2-3months')]
    T3 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T3-6months')]
    T4 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T4-1year')]
    transition.df$T1.1week.3months[i] <- ifelse(c(T1,T2) %in% colnames(distance.matrix), distance.matrix[T1,T2], NA)
    transition.df$T2.3months.6months[i] <- ifelse(c(T2,T3) %in% colnames(distance.matrix), distance.matrix[T2,T3], NA)
    transition.df$T3.6months.1year[i] <- ifelse(c(T3,T4) %in% colnames(distance.matrix), distance.matrix[T3,T4], NA)
  }
  
  # Melt the transition data frame to long format
  transition.melt <- melt(transition.df)
  
  # Remove rows with missing values
  transition.melt <- transition.melt[!is.na(transition.melt$value),]
  return(transition.melt)}

# Plot distance transitions
PlotDistanceTimepoint <- function(se, title){
  
  # Get the sample metadata
  samples.data <- se$targets
  
  # Keep only samples from patients with more than two timepoints
  keep.longitudinal <- names(which(rowSums(table(samples.data$Patient, samples.data$TimepointFactor))>2))
  se <- se[,se$targets$Patient %in% keep.longitudinal] 
 
  # Calculate the distance matrix between samples
  distance.matrix <- as.matrix(dist(t(se$E), method='maximum'))
  
  # Add a "barcode" column to the metadata
  samples.data$barcode <- rownames(samples.data)
  
  # Initialize a data frame to store the distance transitions
  transition.df <- data.frame('T1:1week-3months'=rep(NA, length(keep.longitudinal)),
                              'T2:3months-6months'=rep(NA, length(keep.longitudinal)),
                              'T3:6months-1year'=rep(NA, length(keep.longitudinal)))
  transition.df$Patient <- keep.longitudinal
  
  # Assign each patient to a "Group" based on whether they belong to a certain subgroup
  transition.df$Group <- ifelse(keep.longitudinal %in% samples.data[samples.data$Group=='Yes',]$Patient==TRUE, 'Yes', 'No')
  
  # Loop over each patient and calculate the distances between their timepoints
  for (i in 1:length(keep.longitudinal)) {
    T1 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T1-1week')]
    T2 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T2-3months')]
    T3 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T3-6months')]
    T4 <- samples.data$barcode[which(samples.data$Patient==keep.longitudinal[i] & samples.data$TimepointFactor=='T4-1year')]
    transition.df$T1.1week.3months[i] <- ifelse(c(T1,T2) %in% colnames(distance.matrix), distance.matrix[T1,T2], NA)
    transition.df$T2.3months.6months[i] <- ifelse(c(T2,T3) %in% colnames(distance.matrix), distance.matrix[T2,T3], NA)
    transition.df$T3.6months.1year[i] <- ifelse(c(T3,T4) %in% colnames(distance.matrix), distance.matrix[T3,T4], NA)
  }
  
  # Melt the transition data frame to long format
  transition.melt <- melt(transition.df)
  
  # Remove rows with missing values
  transition.melt <- transition.melt[!is.na(transition.melt$value),]
  
  # Create a boxplot with overlaid points and significance bars
  plot <- ggplot(transition.melt, aes(x=variable, y=value, fill=variable, color=variable)) +
    geom_boxplot(aes(group = variable), fill = "white", color = "white", width = 0.5, outlier.shape = NA) +
    geom_boxplot(aes(group = variable), width = 0.5, outlier.shape = NA) +
    geom_point(size = 2, shape = 21, stroke = 0.25, position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +
    scale_fill_manual(values=adjustcolor(custom_palette_age, alpha=0.5)) +
    scale_color_manual(values=adjustcolor(custom_palette_age, alpha=1)) +
    geom_signif(comparisons=list(c('T1.1week.3months','T2.3months.6months')),
                map_signif_level=c('***'=0.001, '**'=0.01, '*'=0.05), 
                col='black') +
    scale_x_discrete(labels=c(paste0('1 week - 3 months', '\n', 'n = ', table(transition.melt$variable)[1]),
                              paste0('3 months - 6 months', '\n', 'n = ', table(transition.melt$variable)[2]),
                              paste0('6 months - 1 year', '\n', 'n = ', table(transition.melt$variable)[3]))) +
    ylim(0,(max(transition.melt$value)+max(transition.melt$value)*0.1)) +
    ggtitle(title) + xlab('') + ylab('Δ distance between pairs (maximum)') + theme(legend.position='none', text=element_text(size=12))
  return(plot)}

# Redundancy analysis for variables selection
dbRDAselectionMaximum <- function(dge, variables, distance='maximum', direction='both') {
  set.seed(2)
  # Remove samples with missing values for the specified variables
  meta <- dge$targets
  meta <- meta[,colnames(meta) %in% variables]
  meta <- na.omit(meta)  
  dge <- dge[,which(colnames(dge) %in% rownames(meta))]
  # Calculate the distance matrix between samples
  distance.matrix <- as.matrix(dist(t(dge$E), method=distance))
  # Perform a dbRDA with a constant term only
  dbrda.0 <- dbrda(distance.matrix ~ 1, data=meta)
  # Perform a dbRDA with all specified variables
  dbrda.all <- dbrda(as.formula(paste('distance.matrix ~', paste(variables, collapse=" + "))), data=meta)
  # Perform stepwise variable selection using the ordiR2step function
  select <- ordiR2step(dbrda.0, scope=formula(dbrda.all), adjR2thresh=RsquareAdj(dbrda.all)$adj.r.squared, direction=direction)
  # Make a copy of the selected model and adjust the p-values for multiple comparisons
  select.adj <- select
  select.adj$anova$`Pr(>F)` <- p.adjust(select$anova$`Pr(>F)`, method='fdr', n=length(variables))
  # Return the results of the adjusted model and the unadjusted model
  return(select.adj$anova)
}

# PathfindR on gene list
ImmunePathwayAnalysis <- function(DE_table){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  genes.upregulated <- DE_table[!is.na(DE_table$label) & DE_table$logFC>0,]$label
  genes.downregulated <- DE_table[!is.na(DE_table$label) & DE_table$logFC<0,]$label
  genes.upregulated <- bitr(genes.upregulated, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')
  genes.downregulated <- bitr(genes.downregulated, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')
  convert_entrez_to_symbol <- function(entrez_string) {
    entrez_ids <- strsplit(entrez_string, '/')[[1]]
    tryCatch({
      gene_symbols <- bitr(entrez_ids, fromType='ENTREZID', toType='SYMBOL', OrgDb='org.Hs.eg.db', drop=FALSE)
      return(paste(gene_symbols$SYMBOL, collapse="/"))
    }, error = function(e) {
      return(NA)
    })
  }
  
  # Run pathway analysis
  pathway.upregulated.KEGG <- as.data.frame(enrichKEGG(genes.upregulated$ENTREZID, pvalueCutoff = 0.05))
  if (nrow(pathway.upregulated.KEGG)>0) {
    pathway.upregulated.KEGG$Description <- paste("", pathway.upregulated.KEGG$Description)
    pathway.upregulated.KEGG$Direction <- "Upregulated"
    pathway.upregulated.KEGG$Gene_ratio <- sapply(strsplit(pathway.upregulated.KEGG$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  }
  pathway.downregulated.KEGG <- as.data.frame(enrichKEGG(genes.downregulated$ENTREZID, pvalueCutoff = 0.05))
  if (nrow(pathway.downregulated.KEGG)>0) {
    pathway.downregulated.KEGG$Description <- paste("", pathway.downregulated.KEGG$Description)
    pathway.downregulated.KEGG$Direction <- "Downregulated"
    pathway.downregulated.KEGG$Gene_ratio <- sapply(strsplit(pathway.downregulated.KEGG$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  }
  pathway.upregulated.WP <- as.data.frame(enrichWP(genes.upregulated$ENTREZID, pvalueCutoff = 0.05, organism = "Homo sapiens"))
  if (nrow(pathway.upregulated.WP)>0) {
    pathway.upregulated.WP$Description <- paste("WP:", pathway.upregulated.WP$Description)
    pathway.upregulated.WP$Direction <- "Upregulated" 
    pathway.upregulated.WP$Gene_ratio <- sapply(strsplit(pathway.upregulated.WP$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  }
  pathway.downregulated.WP <- as.data.frame(enrichWP(genes.downregulated$ENTREZID, pvalueCutoff = 0.05, organism = "Homo sapiens"))
  if (nrow(pathway.downregulated.WP)>0) {
    pathway.downregulated.WP$Description <- paste("WP:", pathway.downregulated.WP$Description)
    pathway.downregulated.WP$Direction <- "Downregulated"
    pathway.downregulated.WP$Gene_ratio <- sapply(strsplit(pathway.downregulated.WP$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  }
  pathway.upregulated.GO <- as.data.frame(enrichGO(genes.upregulated$ENTREZID, ont = "BP", pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db))
  if (nrow(pathway.upregulated.GO)>0) {
    pathway.upregulated.GO$Description <- paste("GO:", pathway.upregulated.GO$Description)
    pathway.upregulated.GO$Direction <- "Upregulated"
    pathway.upregulated.GO$Gene_ratio <- sapply(strsplit(pathway.upregulated.GO$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  }
  pathway.downregulated.GO <- as.data.frame(enrichGO(genes.downregulated$ENTREZID, ont = "BP", pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db))
  if (nrow(pathway.downregulated.GO)>0) {
    pathway.downregulated.GO$Description <- paste("GO:", pathway.downregulated.GO$Description)
    pathway.downregulated.GO$Direction <- "Downregulated"
    pathway.downregulated.GO$Gene_ratio <- sapply(strsplit(pathway.downregulated.GO$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    }
  
  # Bind all results together
  #path.output.all <- rbind(pathway.upregulated.KEGG, pathway.downregulated.KEGG,
   #                        pathway.upregulated.GO, pathway.downregulated.GO,
    #                       pathway.upregulated.WP, pathway.downregulated.WP)
  
  # Just KEGG
  pathway.upregulated.KEGG$Direction <- "Up"
  pathway.downregulated.KEGG$Direction <- "Down"
  path.output.all <- rbind(pathway.upregulated.KEGG, pathway.downregulated.KEGG)
  path.output.all$geneSymbol <- sapply(path.output.all$geneID, convert_entrez_to_symbol)
  
  return(path.output.all)
}

# Plot immune pathways
PlotImmunePathways <- function(pathways_table, palette, number_pathways, min_genecounts, title){
  
  # Filter top pathways
  filtered_df <- pathways_table[order(pathways_table$p.adjust), ]
  filtered_df <- filtered_df[filtered_df$Count >= min_genecounts, ]
  filtered_df <- filtered_df[1:number_pathways, ]
  
  # Calculate transparency based on log10 of P.adj
  filtered_df$Transparency <- scales::rescale((1/log10(filtered_df$p.adjust)), c(0.8, 0.9))
  
  # Reorder based on log2FoldChange
  filtered_df <- filtered_df[order(filtered_df$Gene_ratio),]
  
  # Create the barplot
  p <- ggplot(filtered_df, aes(x = factor(Description, levels = unique(filtered_df$Description)), y = Gene_ratio, 
                               color = Direction, fill = Direction, alpha = Transparency)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    coord_flip() +
    labs(x = "Pathway", y = "Fold enrichment") +
    ggtitle(title) +
    theme(legend.position = c(1,0),
          legend.box.background = element_rect(), 
          text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm")) +
    guides(alpha = "none")
  return(p)
}

###############################################################################################
###################################### MULTI-OMICS ############################################
###############################################################################################

# Plot variance explained
MOFAvarianceexplained <- function(MOFArun, palette, title) {
  r2f <- as.data.frame(melt(MOFArun@cache$variance_explained$r2_per_factor[[1]]))
  # Create the plot
  p <- ggplot(r2f, aes(x = Var1, y = value, color = Var2, fill = Var2)) +
    geom_bar(position='stack', stat='identity') +
    scale_fill_manual(values=adjustcolor(palette, alpha=0.75)) +
    scale_color_manual(values=adjustcolor(palette, alpha=1)) +
    labs(x="", y='% of variance explained') +
    ggtitle(title) +
    theme(legend.position = c(1,0),
          legend.box.background = element_rect(), 
          text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm")) +
    guides(alpha = "none")
  return(p)}

# Plot top lodings for a given factor of a given view
MOFAfactorloadings <- function(MOFArun, factor, omics, palette, number_features, title){
  
  # Get factor values
  df <- get_weights(MOFArun, scale = FALSE)[[omics]]
  filtered_df <- as.data.frame(df[ ,factor])
  filtered_df$Feature <- rownames(filtered_df)
  filtered_df <- filtered_df[order(abs(filtered_df[,1]), decreasing = TRUE), ]
  colnames(filtered_df) <- c("Factor", "Feature")
  
  # Add direction
  filtered_df$Direction <- ifelse(filtered_df$Factor < 0, "-", "+")

  # Select the top number_features and sort by value
  filtered_df <- filtered_df[1:number_features, ]
  filtered_df <- filtered_df[order(filtered_df$Factor), ]
  
  # Create the barplot
  p <- ggplot(filtered_df, aes(x = factor(Feature, levels = unique(filtered_df$Feature)), y = Factor, 
                               color = Direction, fill = Direction, alpha = 0.5)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values=adjustcolor(palette, alpha=0.5)) +
    scale_color_manual(values=adjustcolor(palette, alpha=1)) +
    coord_flip() +
    labs(x = "Feature", y = factor) +
    ggtitle(title) +
    theme(legend.position = c(1,0),
          legend.box.background = element_rect(), 
          text = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm")) +
    guides(alpha = "none")
  return(p)}

# Boxplots for metadata
BoxplotBinaryMOFAfactor <- function(MOFArun, factor, parameter, title, palette, paletteboxplot) {
  
  # Extract factor values and metadata information from the MOFA object
  MOFA <- MOFArun
  factorvalues <- get_factors(MOFA, factors = factor)
  metadata <- MOFA@samples_metadata[parameter]
  group <- MOFA@samples_metadata['Group']
  
  # Combine factor values, metadata information, and group information into a data frame
  df <- cbind(factorvalues, metadata, group)
  colnames(df) <- c('Factor', 'Parameter', 'Group')
  
  # Remove rows with missing values in the metadata column
  df.complete <- df[df$Parameter != 'NA', ]
  
  # Create the plot using ggplot2
  plot <- ggplot(df.complete, aes(x = Parameter, y = Factor)) +
    
    # Add the boxplot with solid outlines
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      col = adjustcolor(paletteboxplot, alpha = 1),
      fill = adjustcolor('white', alpha = 1)
    ) +
    
    # Add the boxplot with dotted outlines
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      col = adjustcolor(paletteboxplot, alpha = 0.25),
      fill = adjustcolor(paletteboxplot, alpha = 0.1)
    ) +
    
    # Add points representing the data, colored by group
    geom_point(
      aes(color = Group, fill = Group),
      size = 2,
      shape = 21,
      position = position_jitter(0.1),
      stroke = 0.25
    ) +
    
    # Set the fill and outline colors for the points
    scale_fill_manual(values = adjustcolor(palette, alpha = 0.5)) +
    scale_color_manual(values = adjustcolor(palette, alpha = 1)) +
    
    # Add significance labels using a Wilcoxon test and adjust the font size based on the significance level
    geom_signif(
      test = 'wilcox.test',
      comparisons = list(
        c(levels(as.factor(df.complete$Parameter))[1],levels(as.factor(df.complete$Parameter))[2]),
        step_increase = 0.2
      ),
      map_signif_level = c('***' = 0.001, '**' = 0.01, '*' = 0.05)
    ) +
    
    # Set the plot title and axis labels
    ggtitle(title) +
    xlab('') +
    ylab(title) +
    ylab(factor) +
    
    # Remove the legend
    theme(legend.position = c(1, 0), legend.box.background=element_rect(), text=element_text(size=12))
  
  # Return the plot
  return(plot)
}

# Return MOFA factors for a given parameter
ReturnMOFAfactor2 <- function(MOFArun, factor, parameter, title, palette, paletteboxplot) {
  
  # Extract factor values and metadata information from the MOFA object
  MOFA <- MOFArun
  factorvalues <- get_factors(MOFA, factors = factor)
  metadata <- MOFA@samples_metadata[parameter]
  group <- MOFA@samples_metadata['Group']
  
  # Combine factor values, metadata information, and group information into a data frame
  df <- cbind(factorvalues, metadata, group)
  colnames(df) <- c('Factor', 'Parameter', 'Group')
  
  # Remove rows with missing values in the metadata column
  df.complete <- df[df$Parameter != 'NA', ]
  
  # Return the plot
  return(df.complete)
}

# MOFA boxplot 4 groups
BoxplotBinaryMOFAfactor4groups <- function(MOFArun, factor, parameter, title, palette, paletteboxplot) {
  
  # Extract factor values and metadata information from the MOFA object
  MOFA <- MOFArun
  factorvalues <- get_factors(MOFA, factors = factor)
  metadata <- MOFA@samples_metadata[parameter]
  group <- MOFA@samples_metadata['Group']
  
  # Combine factor values, metadata information, and group information into a data frame
  df <- cbind(factorvalues, metadata, group)
  colnames(df) <- c('Factor', 'Parameter', 'Group')
  
  # Remove rows with missing values in the metadata column
  df.complete <- df[df$Parameter != 'NA', ]
  
  # Create the plot using ggplot2
  plot <- ggplot(df.complete, aes(x = Parameter, y = Factor)) +
    
    # Add the boxplot with solid outlines
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      col = adjustcolor(paletteboxplot, alpha = 1),
      fill = adjustcolor('white', alpha = 1)
    ) +
    
    # Add the boxplot with dotted outlines
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      col = adjustcolor(paletteboxplot, alpha = 0.25),
      fill = adjustcolor(paletteboxplot, alpha = 0.1)
    ) +
    
    # Add points representing the data, colored by group
    geom_point(
      aes(color = Group, fill = Group),
      size = 2,
      shape = 21,
      position = position_jitter(0.1),
      stroke = 0.25
    ) +
    
    # Set the fill and outline colors for the points
    scale_fill_manual(values = adjustcolor(palette, alpha = 0.5)) +
    scale_color_manual(values = adjustcolor(palette, alpha = 1)) +
    
    # Add significance labels using a Wilcoxon test and adjust the font size based on the significance level
    geom_signif(
      test = 'wilcox.test',
      comparisons = list(
        c(levels(as.factor(df.complete$Parameter))[1],levels(as.factor(df.complete$Parameter))[2]),
        #c(levels(as.factor(df.complete$Parameter))[1],levels(as.factor(df.complete$Parameter))[3]),
        #c(levels(as.factor(df.complete$Parameter))[1],levels(as.factor(df.complete$Parameter))[4]),
        #c(levels(as.factor(df.complete$Parameter))[2],levels(as.factor(df.complete$Parameter))[3]),
        #c(levels(as.factor(df.complete$Parameter))[2],levels(as.factor(df.complete$Parameter))[4]),
        c(levels(as.factor(df.complete$Parameter))[3],levels(as.factor(df.complete$Parameter))[4])), step_increase = 0.1,
      map_signif_level = c('***' = 0.001, '**' = 0.01, '*' = 0.05)
      ) +
    
    # Set the plot title and axis labels
    ggtitle(title) +
    xlab('') +
    ylab(title) +
    ylab(factor) +
    
    # Remove the legend
    theme(legend.position = c(1, 0), legend.box.background=element_rect(), text=element_text(size=12))
  
  # Return the plot
  return(plot)
}




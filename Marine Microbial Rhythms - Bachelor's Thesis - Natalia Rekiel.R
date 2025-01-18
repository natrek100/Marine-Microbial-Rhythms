#Bachelor's Thesis - Marine Microbial Rhythms: Yearly Patterns in the Northwestern Mediterranean (2007-2015)


#Part 1: Basic data analysis - the SOLA data

#Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(fossil)
library(gridExtra)
library(grid)
library(pheatmap)
library(igraph)
library(tibble)
library(ggrepel)
library(circlize)
library(zoo)

#Import data:
#abundances and taxonomic information
abundance <- read.csv("./abundance.csv")
taxa_info <- read.csv("./taxa_info.csv")


# set the date column as date
abundance$X <- as.Date(abundance$X, format="%Y-%m-%d")

#add a "year" column to the abundance data
abundance$year <- format(as.Date(abundance$X, "%Y-%m-%d"), "%Y") 

# Summarise all columns by year - to achieve the yearly abundance 
yearly_abundance <- abundance %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

yearly_abundance <- yearly_abundance %>%
  mutate(across(-1, ~ round(.x, digits = 2)))

#Part 2: The ratio based on domain - per count and percentage


### The count ratio ###

#now we create a matrix that tells us which ASVs are present per year -> 0.0 abundance = not present = 0, any abundance higher than 0 = present = 1

ASV_presence <- yearly_abundance %>%
  mutate(across(-year, ~ ifelse(. > 0, 1, 0)))

#divide ASV_abundance per domain, and sum the ASVs presence

euk_count <- ASV_presence %>%
  dplyr::select(year, starts_with("Euk")) %>% # Explicitly use dplyr::select
  group_by(year) %>%
  summarise(euk_presence = rowSums(across(where(is.numeric)))) %>%
  mutate(across(everything(), as.numeric))


bac_count <- ASV_presence %>%
  dplyr::select(year, starts_with("Bac")) %>%
  group_by(year) %>%
  summarise(bac_presence = rowSums(across(where(is.numeric)))) %>%
  mutate(across(everything(), as.numeric))

arc_count <- ASV_presence %>%
  dplyr::select(year, starts_with("Arc")) %>%
  group_by(year) %>%
  summarise(arc_presence = rowSums(across(where(is.numeric)))) %>%
  mutate(across(everything(), as.numeric))

#join the counts together again
ASV_count <- euk_count %>%
  left_join(bac_count, by = "year") %>%
  left_join(arc_count, by = "year")

long_format_asv_counts <- ASV_count %>%
  pivot_longer(cols = c(euk_presence, bac_presence, arc_presence),
               names_to = "Domain", values_to = "ASV_count" )




# Create the bar plot. Black outline, domain bars stacked for each year, coloured, with counts inside of them. The grid in gray.
ggplot(long_format_asv_counts, aes(x = factor(year), y = ASV_count, fill = Domain)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.6) +  #black outline
  geom_text(aes(label = ASV_count),  
            position = position_stack(vjust = 0.5),  #stacking
            color = "white", size = 5, fontface = "bold") + #white counts inside bars
  labs(title = "Number of ASVs Present Per Year by Domain", 
       x = "Year", y = "ASV Count") +
  scale_fill_manual(values = c("euk_presence" = "#1f77b4",  
                               "bac_presence" = "#2ca02c", 
                               "arc_presence" = "#d62728"),
                    labels = c("Archaea", "Bacteria", "Eukaryota")) +  #bar colors based on domains
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),  
    legend.title = element_text(face = "bold"),  
    legend.position = "top",  
    panel.grid.major = element_line(color = "gray"),  #gray grid
    panel.grid.minor = element_blank()  
  )






### The percentage ratio ###

# Ensure presence counts are numeric
ASV_count <- ASV_count %>%
  mutate(
    euk_presence = as.numeric(euk_presence),
    bac_presence = as.numeric(bac_presence),
    arc_presence = as.numeric(arc_presence)
  )

# Calculate total presence and percentages
ASV_count_percent <- ASV_count %>%
  mutate(
    total = euk_presence + bac_presence + arc_presence,  # Total counts
    perc_euk = ifelse(total > 0, (euk_presence / total) * 100, 0),
    perc_bac = ifelse(total > 0, (bac_presence / total) * 100, 0),
    perc_arc = ifelse(total > 0, (arc_presence / total) * 100, 0)
  )
mean(ASV_count_percent$total)
# View the result
print(ASV_count_percent)


#"reshape" to long format
percentage_long <- ASV_count_percent %>%
  pivot_longer(cols = c(perc_euk, perc_bac, perc_arc),
               names_to = "Domain", values_to = "percentage" ) %>%
  select(year, Domain, percentage, total)

#visualize using a barplot, same design as previously.
ggplot(percentage_long, aes(x = factor(year), y = percentage, fill = Domain)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.6) +  
  geom_text(aes(label = round(percentage, 1)),  
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") + 
  labs(title = "Percentage of ASVs Present Per Year by Domain", 
       x = "Year", y = "Percentage (%)") +  
  scale_fill_manual(values = c("perc_euk" = "#1f77b4", 
                               "perc_bac" = "#2ca02c", 
                               "perc_arc" = "#d62728"),
                    labels = c("Archaea", "Bacteria", "Eukaryota")) + 
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),  
    legend.title = element_text(face = "bold"), 
    legend.position = "top",  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank()  
  )






#Part 3:  Taxa composition plots 

#firstly let's transform the yearly_abundance into long format
long_format <- yearly_abundance %>%
  pivot_longer(-year, names_to = "ASV", values_to = "Abundance")

#add taxa_info to long format data

long_format_taxa <- merge(long_format, taxa_info, by=("ASV"))

#sum the abundances based on the phylum and year
summed_data <- aggregate(Abundance ~ year + Phylum, data = long_format_taxa, sum)


# select top 5 phyla based on each year. For the remaining phyla, sum their abundances to create an "other" bar.
top_5_and_remaining <- summed_data %>%
  group_by(year) %>%
  
  mutate(rank = rank(-Abundance, ties.method = "first")) %>%
  
  mutate(Phylum = ifelse(rank > 5, "Other", Phylum)) %>%
  
  group_by(year, Phylum) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()




# Calculate the relative abundance (percentage) for each phylum per year
top_5_and_remaining <- top_5_and_remaining %>%
  group_by(year) %>%
  mutate(Percentage = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# put other phyla at the end of the list
top_5_and_remaining$Phylum <- factor(top_5_and_remaining$Phylum, 
                                    levels = c(setdiff(unique(top_5_and_remaining$Phylum), "Other"), 
                                               "Other"))



# Define colors for classes with distinct shades for better contrast
phyla_colors = c(
  "Ochrophyta" = "chartreuse3",    
  "Chlorophyta" = "forestgreen",                
  "Thaumarchaeota" = "darkgreen",         
  "Dinoflagellata" = "orange",           
  "Proteobacteria" = "orange3",             
  "Euryarchaeota" = "mediumpurple",             
  "Other" = "#A9A9A9"               
)
ordered_classes <- c("Thaumarchaeota", "Chlorophyta", "Ochrophyta", "Proteobacteria", "Dinoflagellata", "Euryarchaeota", "Other")
top_5_and_remaining$Phylum <- factor(top_5_and_remaining$Phylum, levels = ordered_classes)

# stacked barplot representing top phyla - minimalistic, without percentage labelled
ggplot(top_5_and_remaining, aes(x = factor(year), y = Percentage, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) +  
  scale_fill_manual(values = phyla_colors) +  
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") + 
  labs(x = "Year", y = "Percentage of Total Abundance", 
       title = "Top 5 Phyla by Relative Abundance for Each Year", 
       fill = "Phylum") +  
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18, face = "bold"),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank(),  
    legend.position = "bottom",  
    legend.direction = "horizontal",  
    legend.justification = "center",  
    legend.text = element_text(size = 16)
   )








#Part 4 : Alpha diversity

#merge the data in long format and the taxa information datasets by ASVs, and sort by year
long_format_div <- merge(long_format, taxa_info, by=("ASV"))
long_format_div <- long_format_div %>%
  arrange(year)
#long_format_div <- long_format_div %>%
#  mutate(Abundance = round(Abundance, digits = 0))
#calculate the alpha diversity indices for the dataset - ASV count, shannon entropy, simpson index and chao1 richness.
alpha_div <- long_format_div %>%
  group_by(year) %>%
  summarise(
    # Number of unique ASVs (species richness)
    Sprichness = sum(Abundance > 0),  # Count of ASVs with non-zero abundance
    
    # Calculate total abundance for the year to use for diversity calculations
    total_abundance = sum(Abundance), 
    
    # Shannon index
    Shannon = -sum((Abundance / total_abundance) * log(Abundance / total_abundance + 1e-10), na.rm = TRUE), 
    
    # Simpson index
    Simpson = 1 - sum((Abundance / total_abundance) ^ 2, na.rm = TRUE), 
    
    # Chao1
    Chao = chao1(Abundance),
    
    # Maximum Shannon index (log(S), where S is the number of unique ASVs)
    max_Shannon = log(Sprichness)
  ) %>%
  # Add a normalized Shannon value
  mutate(shannon_normalized = Shannon / max_Shannon)


# Display the results
print(alpha_div)

# defining colors with quite good contrast for the plots.
div_colors <- c("2007" = "#FF6F61",  # Year 2007
                  "2008" = "#6B5B93",  # Year 2008
                  "2009" = "#88B04B",  # Year 2009
                  "2010" = "#7986CB",  # Year 2010
                  "2011" = "#A0522D",  # Year 2011
                  "2012" = "#64B5F6",  # Year 2012
                  "2013" = "#FF8C00",  # Year 2013
                  "2014" = "#A1887F",  # Year 2014
                  "2015" = "#D32F2F")  # Year 2015

#Shannon bar plot
shannon <- ggplot(alpha_div, aes(x = year, y = shannon_normalized, fill = year)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) +  
  geom_text(aes(label = round(shannon_normalized, 2)),  
            position = position_stack(vjust = 0.5),  
            color = "white", size = 6, fontface = "bold") + 
  labs(y = "Shannon Entropy") +
  scale_fill_manual(values = div_colors) +  
  theme_minimal(base_size = 18) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.title.x = element_blank(), 
    legend.position = "none",  
  )

# ASV count (Sprichness) plot
sprichness <- ggplot(alpha_div, aes(x = year, y = Sprichness, fill = year)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) + 
  geom_text(aes(label = Sprichness),  
            position = position_stack(vjust = 0.5),  
            color = "white", size = 6, fontface = "bold") +  
  labs( y = "Number of ASVs") +
  scale_fill_manual(values = div_colors) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none",  
  )

#Simpson plot
simpson <- ggplot(alpha_div, aes(x = year, y = Simpson, fill = year)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) + 
  geom_text(aes(label = round(Simpson, 2)),  
            position = position_stack(vjust = 0.5),  
            color = "white", size = 6, fontface = "bold") +  
  labs( x = "Year", y = "Simpson Index") +
  scale_fill_manual(values = div_colors) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )

# Chao1 Richness plot
chao1_richness <- ggplot(alpha_div, aes(x = year, y = Chao, fill = year)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) + 
  geom_text(aes(label = round(Chao, 0)),  
            position = position_stack(vjust = 0.5),  
            color = "white", size = 6, fontface = "bold") +  
  labs(x = "Year", y = "Chao1 Richness") +
  scale_fill_manual(values = div_colors) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none",  
  )

# make two plots - ASV count and chao1 richness aligned together, and a separate plot with Shannon entropy and Simpson index
grid.arrange(sprichness, chao1_richness, 
             ncol = 1, nrow = 2)

grid.arrange(shannon, simpson, 
             ncol = 1, nrow = 2)




#Part 5: Beta Diversity
#use our yearly_abundance but without the year column
abundances_without_years <- yearly_abundance %>%
  select(-year)

#then save as a matrix
abundances_without_years <- as.matrix(abundances_without_years)

#also make sure that values are numerical
as.data.frame(lapply(abundances_without_years,as.numeric))

#calculate bray-curtis distances using the vegan package
bray_distances <- vegdist(abundances_without_years, method = "bray")

#save as a matrix
bray_distances <- as.matrix(bray_distances)

#set row/col names
rownames(bray_distances) <- yearly_abundance$year
colnames(bray_distances) <- yearly_abundance$year

# Define the desired order of years 
desired_order <- c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015)

#the bray-curtis dissimilarity matrix as a heatmap 
pheatmap(bray_distances, 
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         labels_row = desired_order,  
         labels_col = desired_order,
         display_numbers = TRUE,
         number_color = "black",  
         fontsize_number = 15,    
         fontsize_row = 18,       
         fontsize_col = 18,       
         xlab = "Year", 
         ylab = "Year",
         cellwidth = 68,  
         cellheight = 68,
)

#the same heatmap but including clustering
pheatmap(bray_distances,
         display_numbers = TRUE,
         number_color = "black",  
         fontsize_number = 15,   
         fontsize_row = 18,       
         fontsize_col = 18,        
         main = "Beta (Bio-)diversity - Bray-Curtis",
         xlab = "Year", ylab = "Year",
         cellwidth = 50,  
         cellheight = 50,
)

#NMDS plot
# Perform NMDS analysis on the Bray-Curtis distance matrix
nmds_result <- metaMDS(bray_distances, distance = "bray", k = 2, trymax = 100)

# Extract NMDS coordinates and add year information for plotting
nmds_points <- as.data.frame(nmds_result$points)
nmds_points$Year <- as.factor(desired_order)  # Add year as a factor for coloring


# Create the NMDS plot 
ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = Year)) +
  geom_point(size = 6, shape = 21, fill = div_colors, stroke = 1.2) +  
  geom_text_repel(aes(label = Year), color = "black", size = 7, fontface = "bold") +  
  scale_color_manual(values = div_colors) + 
  labs(
    title = "NMDS Plot - Beta Diversity (Bray-Curtis)",
    subtitle = paste("Stress:", round(nmds_result$stress, 3)),
    x = "NMDS Dimension 1",
    y = "NMDS Dimension 2",
    color = "Year"
  ) +
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    plot.subtitle = element_text(hjust = 0.5, size = 18, face = "italic"),
    legend.position = "right",  
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 18, face = "bold")  
  ) +
  theme(axis.line = element_line(size = 0.5, color = "gray"),  
        axis.ticks = element_line(color = "gray"))










#Part 6: Analysis of the cytoscape files

#cytoscape tables
cytoscape_taxa <- read.csv("./cytoscape_taxa.csv")
cytoscape_data <- read.csv("./cytoscape_yearly.csv")
cytoscape_ferret_CCM <- read.csv("./Edge_Table.csv")

#Max Abundance Year: a loop that goes through each ASV and their abundances for all years, then selects the largest number (abundance) and the corresponding year. Saved as: ASV name, the maximal count and the corresponding year.
max_abundance_year <- data.frame(Nodes = as.character(), max_abundance = as.numeric(), max_abundance_year = character(), stringsAsFactors = F)
max_abu_year_data <- data.frame(from = as.character(),  max_abundance_year = character(), stringsAsFactors = F)


for(asv in colnames(yearly_abundance)[-1]) {
  max_abu <- max(yearly_abundance[[asv]])
  max_year <- as.numeric(rownames(yearly_abundance)[which.max(yearly_abundance[[asv]])])
  max_abundance_year <- rbind(max_abundance_year, data.frame(Nodes = asv, max_abundance = max_abu, max_abundance_year = max_year ))
}

max_year <- max_abundance_year[,c("Nodes", "max_abundance_year")]
merge_cytoscape_taxa <- merge(cytoscape_taxa, max_year, by = "Nodes", all.x = T)
max_data_year <- max_abu_year_data[,c("from", "max_abundance_year")]
merge_cytoscape_data <- merge(cytoscape_data, max_year, by.x = "from", by.y = "Nodes", all.x = T)


max_year <- max_year %>%
  mutate(max_abundance_year = as.character(2006 + max_abundance_year)) 


#Part 7: Chord Diagrams



CCM_split <- cytoscape_ferret_CCM %>% select(name, corr, from_clu, to_clu) %>%
  separate(name, into = c("from", "to"), sep = " \\(interacts with\\) ")

CCM_split <- merge(CCM_split, max_year, by.x = "to", by.y = "Nodes", all.x = TRUE) #for the to nodes
names(CCM_split)[names(CCM_split) == "max_abundance_year"] <- "Year_to"

CCM_split <- merge(CCM_split, max_year, by.x = "from", by.y = "Nodes", all.x = TRUE) #for the from nodes
names(CCM_split)[names(CCM_split) == "max_abundance_year"] <- "Year_from"

par(cex = 2.5)

#chord diagram including interactions in the same cluster/year
chord_data <- with(CCM_split, table(Year_from, Year_to))
chordDiagram(chord_data, grid.col = div_colors, annotationTrack = c("name", "grid"))

#chord diagram only with interactions from different clusters/years
diff_int <- CCM_split[CCM_split$Year_from != CCM_split$Year_to,] #keep only interactions that are not from the same year
chord_diff <- with(diff_int, table(Year_from, Year_to))
chordDiagram(chord_diff, grid.col = div_colors, annotationTrack = c("name", "grid"))






#Part 8: In and Out degree

#all counts
#select the from or to column, group by year and count rows for each year to calculate the degree.
year_counts_out <- CCM_split %>% select(Year_from) %>%
  group_by(Year_from) %>%
  summarise(out_degree = n())

year_counts_in <- CCM_split %>% select(Year_to) %>%
  group_by(Year_to) %>%
  summarise(in_degree = n())

#merge the counts
in_out_degree_all <- merge(year_counts_out, year_counts_in, by.x = "Year_from", by.y = "Year_to") 
names(in_out_degree_all)[names(in_out_degree_all) == "Year_from"] <- "Year"

#calculate percentage
degree_percent_all <- in_out_degree_all %>%
  group_by(Year) %>%
  mutate(
    total = out_degree + in_degree,  # Total counts
    perc_in = ifelse(total > 0, (out_degree / total) * 100, 0),
    perc_out = ifelse(total > 0, (in_degree / total) * 100, 0)
  )



#count based stacked  barplot
degree_long_count_all <- degree_percent_all %>%
  pivot_longer(cols = c(out_degree, in_degree),
               names_to = "degree", values_to = "count") %>%
  select(Year, degree, count)


all_count <- ggplot(degree_long_count_all, aes(x = factor(Year), y = count, fill = degree)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.6) +  
  geom_text(aes(label = round(count, 1)),  
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") + 
  labs(title = "In and Out Degree (all connections)", y = "Count") +  
  scale_fill_manual(values = c("in_degree" = "#ff7f0e", 
                               "out_degree" = "#17becf"), 
                    labels = c("In Degree", "Out Degree")) +
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16), 
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16),  
    legend.position = "top",  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank()  
  )


#percentage based plot
degree_percent_long_all <- degree_percent_all %>%
  pivot_longer(cols = c(perc_in, perc_out),
               names_to = "degree", values_to = "percentage") %>%
  select(Year, degree, percentage)



perc_all_count <- ggplot(degree_percent_long_all, aes(x = factor(Year), y = percentage, fill = degree)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.6) +  
  geom_text(aes(label = round(percentage, 1)),  
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") + 
  labs( x = "Year", y = "Percentage (%)") +  
  scale_fill_manual(values = c("perc_in" = "#ff7f0e", 
                               "perc_out" = "#17becf"), 
                    labels = c("In Degree", "Out Degree")) +
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16),  
    legend.position = "none",  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank()  
  )


grid.arrange(all_count, perc_all_count, 
             ncol = 1, nrow = 2)




#the same thing but with only counts between different years
#select the from or to column, group by year and count rows for each year to calculate the degree.
year_counts_out_ex <- diff_int %>% select(Year_from) %>%
  group_by(Year_from) %>%
  summarise(out_degree = n())

year_counts_in_ex <- diff_int %>% select(Year_to) %>%
  group_by(Year_to) %>%
  summarise(in_degree = n())

#merge the counts
in_out_degree_ex <- merge(year_counts_out_ex, year_counts_in_ex, by.x = "Year_from", by.y = "Year_to") 
names(in_out_degree_ex)[names(in_out_degree_ex) == "Year_from"] <- "Year"

#calculate percentage
degree_percent_ex <- in_out_degree_ex %>%
  group_by(Year) %>%
  mutate(
    total = out_degree + in_degree,  # Total counts
    perc_in = ifelse(total > 0, (out_degree / total) * 100, 0),
    perc_out = ifelse(total > 0, (in_degree / total) * 100, 0)
  )



#count based plot
degree_long_count_ex <- degree_percent_ex %>%
  pivot_longer(cols = c(out_degree, in_degree),
               names_to = "degree", values_to = "count") %>%
  select(Year, degree, count)


count_excluded <- ggplot(degree_long_count_ex, aes(x = factor(Year), y = count, fill = degree)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.6) +  
  geom_text(aes(label = round(count, 1)),  
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") + 
  labs(title = "In and Out Degree (between different years)", y = "Count") +  
  scale_fill_manual(values = c("in_degree" = "#ff7f0e", 
                               "out_degree" = "#17becf"), 
                    labels = c("In Degree", "Out Degree")) +
theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16), 
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16),  
    legend.position = "top",  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank()  
  )


#percentage based plot
degree_percent_long_ex <- degree_percent_ex %>%
  pivot_longer(cols = c(perc_in, perc_out),
               names_to = "degree", values_to = "percentage") %>%
  select(Year, degree, percentage)



perc_count_excluded <- ggplot(degree_percent_long_ex, aes(x = factor(Year), y = percentage, fill = degree)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.6) +  
  geom_text(aes(label = round(percentage, 1)),  
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") + 
  labs( x = "Year", y = "Percentage (%)") +  
  scale_fill_manual(values = c("perc_in" = "#ff7f0e", 
                                  "perc_out" = "#17becf"), 
                       labels = c("In Degree", "Out Degree")) +
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16),  
    legend.position = "none",  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank()  
  )


grid.arrange(count_excluded, perc_count_excluded, 
             ncol = 1, nrow = 2)



  
#Part 9: Closeness and Betweenness
CCM_split_years_and_asvs_only <- CCM_split %>% select("from", "Year_from") #select the connections
CCM_split_igraph <- CCM_split %>% select("from", "to")


graph_ccm_all <- graph_from_data_frame(CCM_split_igraph, directed = TRUE) #recreate the network

#closeness
closeness_ccm_all <- data.frame(closeness(graph_ccm_all))
closeness_ccm_all <- tibble::rownames_to_column(closeness_ccm_all)
names(closeness_ccm_all) = c("ASV", "closeness centrality")
#betweenness
betweenness_ccm_all <- data.frame(betweenness(graph_ccm_all))
betweenness_ccm_all <- tibble::rownames_to_column(betweenness_ccm_all)
names(betweenness_ccm_all) = c("ASV", "betweenness centrality")

closeness_ccm_all <- merge(closeness_ccm_all, CCM_split_years_and_asvs_only, by.x = "ASV", by.y = "from")
betweenness_ccm_all <- merge(betweenness_ccm_all, CCM_split_years_and_asvs_only, by.x = "ASV", by.y = "from")



# Reshape to long format for plots
closeness_long_all <- closeness_ccm_all %>%
  mutate(metric = "closeness_centrality") %>%
  select(ASV, Year_from, metric, closeness_value = `closeness centrality`)


betweenness_long_all <- betweenness_ccm_all %>%
  mutate(metric = "betweenness_centrality") %>%
  select(ASV, Year_from, metric, betweenness_value = `betweenness centrality`)




#closeness plot

ggplot(closeness_long_all, aes(x = factor(Year_from), 
                               y = closeness_value, 
                               fill = factor(Year_from))) +
  geom_violin(trim = FALSE) + 
  ylim(0,0.00004) +
  labs(title = "Closeness Centrality Violin Plot",
       x = "Year", 
       y = "Centrality Value",
       fill = "Year") +  
  scale_fill_manual(values = div_colors) +
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16, angle = 45, hjust = 2), 
    axis.text.y = element_text(size = 16),  
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.position = "none",  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank()  
  )



#betweenness plot

ggplot(betweenness_long_all, aes(x = factor(Year_from), 
                               y = betweenness_value, 
                               fill = factor(Year_from))) +
  geom_violin(trim = FALSE) + 
  ylim(0,100000) +
  labs(title = "Betweenness Centrality Violin Plot",
       x = "Year", 
       y = "Centrality Value",
       fill = "Year") +  
  scale_fill_manual(values = div_colors) +
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
    axis.text.x = element_text(size = 16, angle = 45, hjust = 2), 
    axis.text.y = element_text(size = 16),  
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(size = 18), 
    legend.position = "none",  
    panel.grid.major = element_line(color = "grey"),  
    panel.grid.minor = element_blank()  
  )


#Yearly graph
CCM_split_yearly_igraph <- CCM_split %>% select("Year_from", "Year_to")
par(cex = 1)
graph_ccm_yearly <- graph_from_data_frame(CCM_split_yearly_igraph, directed = TRUE)
vcount(graph_ccm_yearly)
ecount(graph_ccm_yearly)

#calculate closeness centrality for the yearly network
closeness_values <- closeness(graph_ccm_yearly)
closeness_df <- data.frame(
  Year = as.numeric(V(graph_ccm_yearly)$name), 
  Closeness = closeness_values
)

# Calculate betweenness centrality for each node (year)
betweenness_values <- betweenness(graph_ccm_yearly)
betweenness_df <- data.frame(
  Year = as.numeric(V(graph_ccm_yearly)$name), 
  Betweenness = betweenness_values
)
#closeness centrality plot
ggplot(closeness_df, aes(x = as.factor(Year), y = Closeness, fill = as.factor(Year))) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) +
  geom_text(aes(label = round(Closeness, 3)), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 6, fontface = "bold") +
  labs(title = "Closeness Centrality", x = "Year", y = "Closeness Centrality") +
  scale_fill_manual(values = div_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title.x = element_text(size = 18, face = "bold"),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 16),  
    legend.position = "none", )

# Plot betweenness centrality as a barplot with white text inside bars
ggplot(betweenness_df, aes(x = as.factor(Year), y = Betweenness, fill = as.factor(Year))) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) +
  geom_text(aes(label = round(Betweenness, 2)), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 6, fontface = "bold") +
  labs(title = "Betweenness Centrality", x = "Year", y = "Betweenness Centrality") +
  scale_fill_manual(values = div_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title.x = element_text(size = 18, face = "bold"),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 16),  
    legend.position = "none",  

  )



#Part 10: NMI

# Create the violin plot with a boxplot inside
ggplot(CCM_split, aes(x = "", y = corr, fill = as.factor(Year_from))) +  
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, alpha = 0.7) + 
  geom_boxplot(aes(x = "", y = corr), 
               width = 0.1, 
               color = "black", 
               fill = "white", 
               alpha = 0.5, 
               outlier.shape = NA) +  
  labs(
    title = "Normalized Mutual Information - NMI (out) per year",
    x = NULL, y = "NMI (out)",
    fill = "Year"
  ) +
  ylim(0, 1) +
  scale_fill_manual(values = div_colors) +
  facet_wrap(~ Year_from, nrow = 1, strip.position = "bottom") +  
  theme_minimal(base_size = 18) + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(size = 16),  
    legend.position = "none", 
    strip.text = element_text(size = 16),  
    strip.background = element_blank()  
  )




#Part 11: Environmental Data analysis

#read the data, clean the missing values and set the date as date and add a year column
biogeochemical <- read.csv("./SOLA biogeochemical variables filtered.csv")
biogeochemical[biogeochemical == 999999] <- NA
biogeochemical[biogeochemical == 999999.000] <- NA
biogeochemical$DATE <- as.Date(biogeochemical$DATE, format = "%m/%d/%Y")
biogeochemical$Year <- format(biogeochemical$DATE, "%Y")


# fill missing values by interpolation
biogeochemical$T <- na.approx(biogeochemical$T)
biogeochemical$S <- na.approx(biogeochemical$S)
biogeochemical$O <- na.approx(biogeochemical$O)
biogeochemical$PH <- na.approx(biogeochemical$PH)
biogeochemical$NH4 <- na.approx(biogeochemical$NH4)
biogeochemical$NO3 <- na.approx(biogeochemical$NO3)
biogeochemical$NO2 <- na.approx(biogeochemical$NO2)
biogeochemical$CHLA <- na.approx(biogeochemical$CHLA)
biogeochemical$PO4 <- na.approx(biogeochemical$PO4)
biogeochemical$SIOH4<- na.approx(biogeochemical$SIOH4)
biogeochemical$COP <- na.approx(biogeochemical$COP)
biogeochemical$NOP <- na.approx(biogeochemical$NOP)

#calculate yearly averages
biogeochemical_avg <- biogeochemical_long %>%
  group_by(Year, Variable) %>%
  summarize(Average_Value = mean(Value, na.rm = TRUE))


# Reshape the data into long format
biogeochemical_long <- biogeochemical %>%
  pivot_longer(cols = c(T, S, O, PH, CHLA, NH4, NO3, NO2, PO4, SIOH4, COP, NOP ), names_to = "Variable", values_to = "Value")

# Add Year to the data
biogeochemical_long$Year <- format(biogeochemical_long$DATE, "%Y")

# Assign specific colors to variables
color_palette <- c(
  "T" = "coral",            
  "S" = "cadetblue1",           
  "O" = "cadetblue",       
  "PH" = "darkolivegreen2",       
  "CHLA" = "darkgreen",       
  "NH4" = "darkmagenta",      
  "NO3" = "darkorange",   
  "NO2" = "darkorchid4",    
  "PO4" = "aquamarine4",     
  "SIOH4" = "cornsilk3", 
  "COP" = "cadetblue3",      
  "NOP" = "darkolivegreen"          
)

#plot all variables in one plot
ggplot(biogeochemical_long, aes(x = DATE, y = Value, color = Variable)) +
  geom_line(size = 1.1) + 
  scale_color_manual(values = color_palette) + 
  geom_smooth(aes(group = Variable), color = "gray50", size = 1.2, linetype = "solid", se = FALSE) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  facet_grid(Variable ~ ., scales = "free_y") +  
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20),  
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),  
    panel.grid.major = element_line(color = "gray80", size = 0.5, linetype = "solid"), 
    panel.grid.minor = element_blank(), 
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 20, color = "black", face = "bold")  
  ) +
  labs(
    title = "Environmental Variables Across Years",
    y = "Value",
    x = "Year"
  )



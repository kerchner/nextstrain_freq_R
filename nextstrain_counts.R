library(dplyr)
library(ggplot2)
library(R.utils) # Needed for gunzip()
library(readr)

dir.create('data', showWarnings = FALSE)

# Control parameters
data_source <- 'gisaid' # 'open'
redownload <- FALSE # switch to TRUE to get fresh data

if (redownload) {
  # Download data and gunzip
  download.file(url = paste0('https://data.nextstrain.org/files/workflows/forecasts-ncov/',
                             data_source, '/nextstrain_clades/global.tsv.gz'),
                destfile = 'data/global.tsv.gz')
  
  gunzip('data/global.tsv.gz', remove = FALSE, overwrite = TRUE)
}

# Use readr::read_delim instead of read.table to avoid problems with the UTF-8 encoding
df <- read_delim('data/global.tsv',
                 delim = '\t') %>%
  rename(country = location)

# Extract 2-digit year-month (20-12)
df <- df %>%
  mutate(date_ym = substr(date, 3, 7))

# Group into variants
df <- df %>% 
  mutate(variant = case_when(clade %in% 
                               c("19A", "19B", "20A", "20B", "20C", "20D", "20F", "20G", "20E") ~ "Pre_alpha",
                             clade %in% 
                               c("20I") ~ "Alpha",
                             clade %in%
                               c("20H") ~ "Beta", 
                             clade %in%
                               c("21B") ~ "Kappa", 
                             clade %in%
                               c("21D") ~ "Eta",
                             clade %in% 
                               c("21E") ~ "Theta",
                             clade %in% 
                               c("20J", "20J") ~ "Gamma",
                             clade %in% 
                               c("21A", "21I", "21J") ~ "Delta",
                             clade %in% 
                               c("21C") ~ "Epsilon",
                             clade %in% 
                               c("21F") ~ "Iota",
                             clade %in% 
                               c("21G") ~ "Lambda",
                             clade %in% 
                               c("21H") ~ "Mu",
                             clade %in% 
                               c("21K", # BA.1
                                 "21L", # BA.2
                                 "21M", # Omicron
                                 "22A", # BA.4
                                 "22B", # BA.5
                                 "22C", # BA.2.12.1
                                 "22D", # BA.2.75
                                 "22E", # BQ.1
                                 "23C", # CH.1.1
                                 "23I") # BA.2.86 - highly mutated vs. Omicron but emerged late 2023
                             ~ "Omicron",
                             clade %in%
                               c("22F", # XBB
                                 "23A", # XBB.1.5
                                 "23B", # XBB.1.16
                                 "23D", # XBB.1.9
                                 "23E", # XBB.2.3
                                 "23F", # EG.5.1
                                 "23G", # XBB.1.5.70
                                 "23H") # HK.3
                             ~ "XBB",
                             clade %in%
                               c("recombinant") ~ "recombinant"))

df <- df %>%
  filter(country %in% c("Spain", "Mexico", "USA", "Italy", "Kenya", "Hong Kong",
                         "France", "Malaysia", "Switzerland", "Turkey", "Saudi Arabia", "India")) %>%
  arrange(country, date, clade)

# Check: Which clades don't yet have a variant mapped?
# Should be empty if we've chosen a variant for all clades.  Otherwise, need to fix.
df_unknown <- df %>% filter(is.na(variant))
table(df_unknown$clade)

# Compute monthly totals by country and variant
df_monthly <- df %>%
  mutate(date_ym = substr(date, 3, 7)) %>%
  group_by(date_ym, country, variant) %>%
  summarize(count = sum(sequences, na.rm = TRUE))

# Compute relative ("normalized") variant frequencies per month, location
df_monthly <- df_monthly %>%
  group_by(date_ym, country) %>%
  mutate(nor = round(count/sum(count, na.rm = TRUE), digits = 3))

# Create plot (without studies)
p1 <- ggplot() +
  geom_bar(data = df_monthly, stat="identity",
           mapping = aes(x = date_ym, fill = variant, y=nor)) +
  theme_light() + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1, size = 6),
        axis.text.y = element_text(size = 6)) + 
  facet_wrap(~ country,  nrow = 12, ncol = 1)#, scale = "free_x")

ggsave(filename = paste0("plot1", data_source, ".pdf"), plot = p1,
       width = 8, height = 20, units = "in")  



library(dplyr)
library(tidyr)
library(ggplot2)

minfactor <- function(number) {
  for (i in 2:10) {
    if (number %% i == 0) {
      return(i)
    }
  }
  return(NA)
}


# chrom_lengths <- read.delim("data/chromosome_lengths.tsv")
# blocks <- read.delim('../SE_SRWW_PHG/haploblocks/test2.prune.in', header = F) %>%
blocks <- read.delim('../haploblocks/recombination_thin.prune.in', header = F) %>%
  separate(V1, c('Chromosome', 'Position', 'ref'))  %>%
  mutate(Position = as.numeric(Position))
# blocks$color <- as.factor(rep(1:7))

blocks_segments <- blocks %>%
  arrange(Chromosome, Position) %>%
  group_by(Chromosome) %>%
  mutate(BlockStart = Position+1,
         BlockEnd   = lead(Position)) %>%
  filter(!is.na(BlockEnd)) 
coln=3
blocks_segments$arbitrary_color <- as.factor(rep(1:coln, length = nrow(blocks_segments)))

ggplot() +
  geom_segment(data = blocks_segments,
               aes(x = BlockStart, xend = BlockEnd,
                   y = Chromosome, yend = Chromosome,
                   color = arbitrary_color),
               linewidth = 4, lineend = "butt") + 
  theme_minimal() +
  labs(x = "Position (bp)", y = "Chromosome", title="Haploblock sizes")
# ggsave('../SE_SRWW_PHG/figures/plink_haploblocks.png', width=14, height=6,)
ggsave('../figures/plink_haploblocks_2.png', width=14, height=6,)


first_blocks <- blocks %>%
  arrange(Chromosome, Position) %>%
  group_by(Chromosome) %>%
  slice_head(n = 1) %>%
  mutate(BlockStart = 0,
         BlockEnd   = Position) %>%
  ungroup()

# combine with the rest
blocks_segments <- bind_rows(blocks_segments, first_blocks) %>%
  select(Chromosome, BlockStart, BlockEnd) %>%
  mutate(ID = paste0("TraesCS", Chromosome, "03", BlockStart, "_", BlockEnd),
         X = 0,
         Strand = ".") %>%
  arrange(Chromosome, BlockStart)

# write.table(blocks_segments, '../SE_SRWW_PHG/data/haploblocks.txt',
write.table(blocks_segments, '../data/haploblocks.bed',
            quote = F, sep = "\t", 
            row.names = F, col.names = F)


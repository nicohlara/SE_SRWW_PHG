
library(dplyr)
library(tidyr)
library(ggplot)

minfactor <- function(number) {
  for (i in 2:10) {
    if (number %% i == 0) {
      return(i)
    }
  }
  return(NA)
}


# chrom_lengths <- read.delim("data/chromosome_lengths.tsv")
blocks <- read.delim('../SE_SRWW_PHG/haploblocks/test2.prune.in', header = F) %>%
  separate(V1, c('Chromosome', 'Position', 'ref'))  %>%
  mutate(Position = as.numeric(Position))
# blocks$color <- as.factor(rep(1:7))

blocks_segments <- blocks %>%
  arrange(Chromosome, Position) %>%
  group_by(Chromosome) %>%
  mutate(BlockStart = Position,
         BlockEnd   = lead(Position)) %>%
  filter(!is.na(BlockEnd)) 
coln=minfactor(nrow(blocks_segments))
blocks_segments$color <- as.factor(rep(1:coln, times = nrow(blocks_segments)/coln))

ggplot() +
  geom_segment(data = blocks_segments,
               aes(x = BlockStart, xend = BlockEnd,
                   y = Chromosome, yend = Chromosome,
                   color = color),
               linewidth = 4, lineend = "butt")
ggsave('../SE_SRWW_PHG/figures/plink_haploblocks.png')


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
         Strand = ".")

write.table(blocks_segments, '../SE_SRWW_PHG/data/haploblocks', 
            quote = F, sep = "\t", 
            row.names = F, col.names = F)


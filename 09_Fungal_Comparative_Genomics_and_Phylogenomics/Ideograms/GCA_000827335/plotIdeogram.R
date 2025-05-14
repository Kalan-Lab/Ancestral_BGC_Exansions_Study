library(ggideogram)
library(dplyr)
library(ggplot2)
library(tidyverse)

hi.dat <- as_tibble(read.table('HI_Positions.txt', header=T, sep='\t'))
bgc.dat <- as_tibble(read.table('Ideogram_Coordinates.txt', header=T, sep='\t'))

# Chrom	Start	End	Stain	Arm
# Protein	Chrom	Position	Type

pdf("Ideogram.pdf", height=3, width=6)

mutate(bgc.dat, Chr = factor(Chrom)) %>%
    ggplot() +
    geom_ideogram(aes(x = Chrom, ymin = Start, ymax = End, 
                      chrom = Chrom, fill = Stain), 
                  radius = unit(4, 'pt'), width = 0.5, 
                  linewidth = 1, chrom.col = 'grey80',
                  show.legend = FALSE) + scale_fill_manual(values=c('grey', 'gold')) + 
    geom_point(data = hi.dat,
               aes(x = Chrom, y = Position, 
                   shape = Type, colour = Type), show.legend=F,
               position = position_nudge(x = 0.0)) +
    scale_colour_manual(values = c("#7A3E56", "#AA5ED6")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    theme(
        axis.title = element_blank(),
	#axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
    )


dev.off()

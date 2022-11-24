library(dplyr)
library(tidyverse)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


source("RIO.R")

input <- function(infile) {
         pfix = prefix()
   parameters <<- readParameters(infile)
   kraken <<- paste(pfix, parameters['kraken', 2], sep="/")
   sckmer <<- paste(pfix, parameters['sckmer', 2], sep="/")
}

run <- function() {}

output <- function(outfile) {

# kraken report
report = read.delim(kraken, header = F)
report$V8 = trimws(report$V8)
report[report$V8 %in% c('Homo sapiens', 'Bacteria', 'Fungi', 'Viruses'), ]

# sckmer data
kmer_data = read.table(sckmer, header = T)
head(kmer_data)

length(unique(report$V8[report$V6 %in% c('G', 'S')]))

c = kmer_data %>%
  subset(kmer > 1) %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn > 3) %>%
  group_by(taxid) %>%
  summarize(r = cor.test(kmer, uniq, method = 'spearman')$estimate,
            p = cor.test(kmer, uniq, method = 'spearman')$p.value,
            .groups='keep') %>%
  mutate(p = p.adjust(p))

c$name = report$V8[match(c$taxid, report$V7)] # add taxa names
print(c)

}

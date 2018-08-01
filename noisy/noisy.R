library(tidyverse)
library(lme4)

rm(list = ls())

d <- read.csv('~/Desktop/infotheory/noisy/out.txt', sep = '\t')

ggplot(data = d, aes(mutual.info, group = lex.type, fill = lex.type)) + 
  geom_histogram(position = 'stack', bins = 40) +
  xlab('Mutual Information') + ylab('Count') +
  scale_fill_discrete(name="Lexicon Type", labels = c('Original', 'Scrambled Frequencies')) +
  geom_segment(aes(x = d$mutual.info[1], y = 25, xend = d$mutual.info[1], yend = 15),
                                                             arrow = arrow(length = unit(0.5, "cm"))) +
  theme_set(theme_gray(base_size = 16))
ggsave('~/Documents/alc2018/noisy.png')
d$real.greater <- d$mutual.info < d$mutual.info[1]
sum(d$real.greater)


d <- read.csv('~/Desktop/infotheory/noisy/dutch_out.txt', sep = '\t')

ggplot(data = d, aes(mutual.info, group = lex.type, fill = lex.type)) + 
  geom_histogram(position = 'stack', bins = 40) +
  xlab('Mutual Information') + ylab('Count') +
  scale_fill_discrete(name="Lexicon Type", labels = c('Original', 'Scrambled Frequencies')) +
  geom_segment(aes(x = d$mutual.info[1], y = 50, xend = d$mutual.info[1], yend = 35),
               arrow = arrow(length = unit(0.5, "cm"))) +
  theme_set(theme_gray(base_size = 16))
ggsave('~/Documents/alc2018/dutch_noisy.png')
d$real.greater <- d$mutual.info < d$mutual.info[1]
sum(d$real.greater)


# lmer for effect of word prob on prob of noiselessness (accurate transmission)
e <- read.csv('~/Desktop/infotheory/noisy/by_word.output.txt', sep = '\t')
summary(lmer(data = e, scale(p.as.self) ~ scale(freq) + (1|skeleton), REML = F))

e_dut <- read.csv('~/Desktop/infotheory/noisy/dutch_by_word.output.txt', sep = '\t')
summary(lmer(data = e_dut, scale(p.as.self) ~ scale(freq) + (1|skeleton), REML = F))

for (x in unique(e$skeleton)) {
  xd <- subset(e, skeleton == x)
  if (nrow(xd) > 5 & coefficients(lm(data = xd, p.as.self ~ frequency))[2] > .0001) {
    print(x)
    print(coefficients(lm(data = xd, p.as.self ~ frequency)))
    ggplot(data = sample_n(xd, min(15, nrow(xd))), aes(x = frequency, y = p.as.self, label = word)) + geom_text() + geom_smooth(method = 'lm') + ggtitle(x)
    ggsave(paste0('~/Desktop/infotheory/noisy/graphs/', x, '.png'))
  }
}


# CVFCF
ggplot(data = subset(e, skeleton == 'CVFCF'), aes(x = frequency, y = p.as.self, label = word)) + geom_text() + geom_smooth(method = 'lm')
# CVFW
ggplot(data = subset(e, skeleton == 'CVFW'), aes(x = frequency, y = p.as.self, label = word)) + geom_text() + geom_smooth(method = 'lm')
# VCWCFV
ggplot(data = subset(e, skeleton == 'VCWCFV'), aes(x = frequency, y = p.as.self, label = word)) + geom_text() + geom_smooth(method = 'lm')

options(echo=FALSE)
options(error=function()traceback(2))

library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(reshape2)
library(Hmisc)
library(directlabels)
library(plot3D)

Sys.setlocale("LC_ALL", "en_US.UTF-8")
args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    stop('No measurements file specified.')
} else if (length(args) == 1) {
    stop('Missing plot name.')
}


plot.median <- function(y)
{
    m <- median(y)
    s <- ifelse(m >= 1e5, format(m, digits=2, scientific=TRUE), paste(round(m, 2)))
    return (data.frame(y = m, label = s))
}

prettify <- function(df)
{
    if ('Algorithm' %in% colnames(df)) {
        df$Algorithm <- ordered(df$Algorithm,
                                levels=c('memcpy', 'branching',                 'pirk',                     'vectorized',               'register',                     'rewired (smallpages)',                 'rewired (hugepages)',              'rewired simd (smallpages)',                'rewired simd (hugepages)',                 'mixed cracking',     'pirk (MT)',                  'rewired (hugepages+MT)',       'hybrid (hugepages+MT)'),
                                labels=c('memcpy', 'branching\ncrack-in-two',    'predicated\ncrack-in-two',  'vectorized\ncrack-in-two',  'predicated++\ncrack-in-two',    'rewired\ncrack-in-two (smallpages)',    'rewired\ncrack-in-two\n(hugepages)', 'SIMD rewired\ncrack-in-two (smallpages)',   'SIMD rewired\ncrack-in-two\n(hugepages)',    'mixed crack-in-two', 'refined partition & merge',  'rewired partition & merge',    'rewired partiton & merge*'))
    }
    if ('Element' %in% colnames(df)) {
        df$Element <- ordered(df$Element,
                              levels=c(4, 8),
                              labels=c('4+4 Byte', '8+8 Byte'))
    }

    df$Value[df$Attribute=='mem_read'] <- df$Value[df$Attribute=='mem_read'] / (1024 * 1024)
    df$Value[df$Attribute=='mem_write'] <- df$Value[df$Attribute=='mem_write'] / (1024 * 1024)

    if ('Attribute' %in% colnames(df)) {
        # TODO
    }
    return(df)
}


m <- read.csv(args[1], header=TRUE, sep=',', fileEncoding='UTF-8', strip.white=TRUE, na.strings=c('NA', 'nan'))
print(levels(m$Algorithm))
m$Element <- ordered(m$Element, levels=c('4', '8'))
m$Selectivity <- as.numeric(as.character(m$Selectivity))
m$Value <- as.numeric(as.character(m$Value))
m$Threads <- factor(m$Threads)

{
    branch <- m[m$Attribute=='branch',]
    brmiss <- m[m$Attribute=='brmiss',]
    brmissratio <- branch
    brmissratio$Attribute <- revalue(brmissratio$Attribute, c('branch'='brmissratio'))
    brmissratio$Value <- brmiss$Value / branch$Value
    m <- rbind(m, brmissratio)
}

{
    mem_read <- m[m$Attribute=='mem_read',]
    mem_write <- m[m$Attribute=='mem_write',]
    mem_IO <- mem_read
    mem_IO$Attribute <- revalue(mem_IO$Attribute, c('mem_read'='mem_IO'))
    mem_IO$Value <- (mem_read$Value + mem_write$Value)
    m <- rbind(m, mem_IO)
}

str(m)

sub <- subset(m, Attribute == 'time' & Threads %in% c(1,2,4,6,8), select=c('Algorithm', 'Threads', 'Value'))

sub$Value = sub$Value / 1000 # to seconds

ggplot(prettify(sub), aes(Algorithm, Value, group=Threads, fill=Threads)) +
    geom_bar(show.legend=TRUE, stat='summary', fun.y='median', position=position_dodge(width=.8), width=.6) +
    ylab('Time [s]') +
    theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(panel.grid.major.x=element_blank()) +
    theme(legend.position=c(1,1), legend.direction='horizontal', legend.justification=c(1,1), legend.background=element_blank()) +
    guides(fill=guide_legend(nrow=1))
ggsave('scaleup.pdf', device='pdf', width=6, height=2)

str(sub)

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
                                levels=c('memcpy', 'branching',                 'pirk',                 'vectorized',           'register',                     'rewired (smallpages)',                 'rewired (hugepages)',              'rewired simd (smallpages)',                'rewired simd (hugepages)',                 'mixed cracking', 'pirk (MT)',                  'rewired (hugepages+MT)',       'hybrid (hugepages+MT)'),
                                labels=c('memcpy', 'branching crack-in-two',    'predicated crack-in-two',  'vectorized crack-in-two',  'predicated++ crack-in-two',    'rewired crack-in-two (smallpages)',    'rewired crack-in-two', 'SIMD rewired crack-in-two (smallpages)',   'SIMD rewired crack-in-two (hugepages)',    'mixed cracking', 'refined partition & merge',  'rewired partition & merge',    'rewired partiton & merge*'))
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

# Get the median of memcpy
memcpy <- subset(m, Experiment=='misc' & Algorithm=='memcpy' & Element=='4' & Attribute=='time' & Threads==1,
                 select=c(Algorithm, Value, Attribute))
medians <- ddply(memcpy, c('Attribute'), plyr::summarise, Time=median(Value))
memcpy_median <- as.numeric(medians[1,][2])
mem_bandwidth <- .426 # sec per 4GiB on deeprig01
print(paste('median (memcpy) is', memcpy_median))


# Show the behaviour of branching and predicated for different selectivities and in relation to memcpy.
{
    print('comparison of branching and Pirk by selectivity to memcpy')

    sub <- subset(m,
                  (Layout=='AoS' | Layout=='SoA') &
                      Experiment=='partition' & Distribution=='uniform dense' & Attribute=='time' &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='vectorized') &
                      (Element=='4' | Element=='8'),
                  select=c(Layout, Algorithm, Element, Selectivity, Value, Attribute))

    summarized_sub <- ddply(sub, c('Layout', 'Algorithm', 'Selectivity', 'Element'), plyr::summarise, time=median(Value))
    summarized_sub$time <- summarized_sub$time / 1000 # to seconds

    ggplot(prettify(summarized_sub), aes(Selectivity, time, colour=Algorithm, linetype=Element)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_wrap(~ Layout) +
        geom_hline(aes(yintercept=mem_bandwidth), size=.5, linetype='F1', show.legend=NA) +
        annotate('text', x=0, y=mem_bandwidth, label='single core bandwidth', hjust=-.05, vjust=1.4, family='serif', colour='black', size=3) +
        labs(y='Time [s]') +
        scale_linetype_discrete(name='Element size') +
        scale_y_continuous(sec.axis=dup_axis(name='')) +
        theme(legend.position='bottom', legend.box='horizontal', legend.direction='vertical') +
        guides(colour=guide_legend(order=1)) +
        expand_limits(y=0)
    ggsave(paste(args[1], '_selectivity.pdf', sep=''), device='pdf', width=5, height=3.8)
}

# Show the behaviour of branching and predicated for different selectivities and in relation to memcpy.
{
    print('comparison of more algorithms by running time over selectivity')

    # AoS + SoA
    sub <- subset(m,
                  Experiment=='partition' & Distribution=='uniform dense' & Attribute=='time' &
                      (Layout=='AoS' | Layout=='SoA') &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='register' | grepl('^rewired.*pages)', Algorithm) | Algorithm=='vectorized') &
                      (Element=='4' | Element=='8'),
                  select=c(Layout, Algorithm, Element, Selectivity, Value, Attribute))

    summarized_sub <- ddply(sub, c('Layout', 'Algorithm', 'Selectivity', 'Element'), plyr::summarise, time=median(Value))
    summarized_sub$time <- summarized_sub$time / 1000 # ms to s

    ggplot(prettify(summarized_sub), aes(Selectivity, time, colour=Algorithm)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_grid(Element ~ Layout) +
        geom_hline(aes(yintercept=mem_bandwidth), size=.5, linetype='F1', show.legend=NA) +
        annotate('text', x=0, y=mem_bandwidth, label='single core bandwidth', hjust=-.05, vjust=1.4, family='serif', colour='black', size=3) +
        labs(y='Time [s]') +
        scale_y_continuous(breaks=seq(0, 20, 1), minor_breaks=seq(0, 20, .5), sec.axis=dup_axis(name='')) +
        theme(legend.title=element_blank(), legend.position='bottom', legend.box='horizontal', legend.direction='vertical', legend.margin=margin(-12,0,0,0,'pt')) +
        guides(colour=guide_legend(order=1, ncol=2, title.position='top')) +
        expand_limits(y=0)
    ggsave(paste(args[1], '_selectivity_all.pdf', sep=''), device='pdf', width=5, height=6)
}

# Compare branching vs. Pirk in terms of time, branch misses, and stalls.
{
    print('comparison of branching vs. naive')

    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Distribution=='uniform dense' &
                      (Algorithm=='branching' | Algorithm=='pirk') &
                      (Element=='4' | Element=='8') &
                      (Attribute=='time' | Attribute=='brmiss' | Attribute=='stalled'),
                  select=c(Algorithm, Element, Selectivity, Value, Attribute))
    sub[sub$Attribute=='brmiss', 'Value'] <- sub[sub$Attribute=='brmiss', 'Value'] / 1e6
    sub[sub$Attribute=='stalled', 'Value'] <- sub[sub$Attribute=='stalled', 'Value'] / 1e6
    sub[sub$Attribute=='time', 'Value'] <- sub[sub$Attribute=='time', 'Value'] / 1000
    sub$Attribute <- revalue(sub$Attribute, c(brmiss='Branch Misses (Mil.)', stalled='Stalled Cycles (Mil.)', time='Running Time [s]'))

    ggplot(prettify(sub), aes(Selectivity, Value, colour=Algorithm, linetype=Element)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_wrap(~ Attribute, scales='free_y', nrow=1) +
        expand_limits(y=0) +
        scale_linetype_discrete(name='Element size') +
        theme(legend.position='right', legend.direction='vertical') +
        theme(panel.spacing.x=unit(1, 'lines')) +
        guides(colour=guide_legend(order=1)) +
        scale_x_continuous(breaks=seq(0, 1, .5)) +
        #theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
        theme(legend.position='bottom', legend.box='horizontal', legend.direction='vertical', legend.margin=margin(t = -3, unit='mm'))
    ggsave(paste(args[1], '_comparison_branching.pdf', sep=''), device='pdf', width=5, height=2.5)

    # Compute average IPC for all algorithms, and print them to the console.
    ipc_sub <- subset(m,
                      Experiment=='partition' & Distribution=='uniform dense' &
                          (Element=='4' | Element=='8') &
                          Attribute=='ipc',
                      select=c(Algorithm, Value))
    ipcs <- ddply(ipc_sub, c('Algorithm'), plyr::summarise, IPC=mean(Value))
    print(ipcs)
}

# Plot a comparison of Pirk and register crack-in-two for all element sizes
{
    print('comparison of various predicated algorithms as boxplots for all element sizes')

    # AoS
    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Distribution=='uniform dense' & Element==4 &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='register'),
                  select=c(Algorithm, Attribute, Element, Selectivity, Value))

    ggplot(prettify(sub), aes(Selectivity, Value, fill=Algorithm, linetype=Element)) +
        geom_boxplot(position='dodge') +
        facet_wrap(~ Attribute, scales='free_y') +
        expand_limits(y = 0) +
        xlab('Selectivity') + ylab('Value')
    ggsave(paste(args[1], '_comparison_predication_AoS.pdf', sep=''), device='pdf', width=12, height=8)
}

# AoS vs. SoA
{
    print('comparing algorithms AoS vs. SoA')

    sub <- subset(m,
                  (Layout=='AoS' | Layout== 'SoA') &
                      Experiment=='partition' &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='register' | grepl('^rewired.*pages)', Algorithm)) &
                  (Selectivity >= .2 & Selectivity <= .8) &
                  Element==4 & Attribute=='time',
            select=c(Layout, Algorithm, Value))

    ggplot(prettify(sub), aes(Algorithm, Value, fill=Layout)) +
        geom_bar(show.legend=TRUE, stat='summary', fun.y='median', position=position_dodge(width=.8), width=.6) +
        labs(y='Time [ms]') +
        theme(axis.text.x=element_text(angle=15, hjust=1))

    ggsave(paste(args[1], '_comparison_layout.pdf', sep=''), device='pdf', width=9, height=9)
}

# Register vs. Rewired threshold
if (F) {
    print('threshold between register and rewired')

    sub <- subset(m,
                  Experiment=='increasing' & Attribute=='time' & Element==4 &
                      (Algorithm=='register' | grepl('^rewired.*pages)', Algorithm)),
                  select=c(Layout, Algorithm, Size, Value))
    sub$Size <- sub$Size / (1024 * 1024)

    ggplot(prettify(sub), aes(Size, Value, colour=Algorithm, fill=Algorithm)) +
        facet_wrap(~ Layout, ncol=2) +
        stat_summary(geom='line', fun.y=mean, size=.3) +
        stat_summary(geom='ribbon', fun.data=mean_cl_normal, fun.args=list(conf.int=.95), colour=NA, alpha=.3) +
        labs(x='Size [MiB]', y='Time [ms]') +
        scale_x_continuous(breaks=seq(0, 128, 16), minor_breaks=seq(0, 128, 4)) +
        scale_y_continuous(minor_breaks=seq(0, 120, 5)) +
        theme(legend.title=element_blank(), legend.position='bottom', legend.direction='vertical', legend.box='horizontal', legend.margin=margin(-12,0,0,0,'pt')) +
        guides(colour=guide_legend(ncol=2))
    ggsave(paste(args[1], '_threshold.pdf', sep=''), device='pdf', width=5.3, height=4)
}

# Cracking
{
    print('comparing cracking')

    # AoS + SoA
    sub <- subset(m,
                  Experiment=='crack' &
                      (Layout=='AoS' | Layout=='SoA') &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='register' | grepl('^rewired.*hugepages)', Algorithm)) &
                      Attribute=='time',
            select=c(Layout, Algorithm, Element, Value))

    sub$Value <- sub$Value / 1000 # to seconds

    ggplot(prettify(sub), aes(Element, Value, fill=Algorithm)) +
        geom_bar(show.legend=TRUE, stat='summary', fun.y='median', position=position_dodge(width=.8), width=.6) +
        facet_wrap(~ Layout) +
        labs(x='Element size', y='Time [s]') +
        scale_y_continuous(breaks=seq(0, 50, 10), minor_break=seq(0, 50, 5), sec.axis=dup_axis(name='')) +
        theme(legend.title=element_blank(), legend.position='bottom', legend.direction='vertical', legend.box='horizontal', legend.background=element_blank(), legend.margin=margin(-12,0,0,0,'pt')) +
        theme(axis.ticks.x=element_blank()) +
        theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
        guides(fill=guide_legend(ncol=2), title.position='top')
    ggsave(paste(args[1], 'cracking_all.pdf', sep='_'), device='pdf', width=5, height=4)
}

# Cracking with mixed
{
    print('comparing cracking with mixed')

    # AoS + SoA
    sub <- subset(m,
                  Experiment=='crack' &
                      (Layout=='AoS' | Layout=='SoA') &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='register' | grepl('^rewired.*pages)', Algorithm) | Algorithm=='mixed cracking') &
                      Attribute=='time',
            select=c(Layout, Algorithm, Element, Value))

    sub$Value <- sub$Value / 1000 # to seconds

    ggplot(prettify(sub), aes(Element, Value, fill=Algorithm)) +
        geom_bar(show.legend=TRUE, stat='summary', fun.y='median', position=position_dodge(width=.8), width=.6) +
        facet_wrap(~ Layout) +
        labs(x='Element size', y='Time [s]') +
        scale_y_continuous(breaks=seq(0, 50, 10), minor_break=seq(0, 50, 5), sec.axis=dup_axis(name='')) +
        guides(fill=guide_legend(ncol=2)) +
        theme(legend.position='bottom', legend.box='horizontal', legend.background=element_blank())
    ggsave(paste(args[1], 'cracking_mixed.pdf', sep='_'), device='pdf', width=7, height=5)
}

if (F) {
    print('Time per query when cracking')

    sub <- subset(m, Experiment=='per_query' &
                  (Algorithm=='branching' | Algorithm=='register' | Algorithm=='register' | grepl('^rewired.*hugepages)', Algorithm)),
                  select=c(Layout, Element, Algorithm, Num, Value))

    ggplot(prettify(sub), aes(Num, Value, colour=Algorithm, fill=Algorithm)) +
        stat_summary(fun.y=median, geom='point', shape=20, size=.01, alpha=.5, show.legend=FALSE) +
        geom_smooth(method='loess', size=.4, n=500, span=.3) +
        facet_grid(Element ~ Layout) +
        scale_y_log10(breaks=10^seq(0,3), minor_breaks=c(tcrossprod(1:10, 10^(-1:4)))) +
        labs(x='# of query', y='Time [ms]') +
        theme(legend.title=element_blank(), legend.position='bottom', legend.direction='vertical', legend.box='horizontal', legend.margin=margin(-12,0,0,0,'pt')) +
        guides(colour=guide_legend(ncol=2))
    ggsave(paste(args[1], '_per_query.pdf', sep=''), device='pdf', width=5.1, height=5)
}

# Mixed cracking
{
    sub <- subset(m,
                  Experiment=='crack' & Layout=='SoA' & Element==4 &
                      (Algorithm=='mixed cracking' | Algorithm=='rewired simd (hugepages)'),
                  select=c(Algorithm, Attribute, Value))

    ggplot(prettify(sub), aes(Algorithm, Value)) +
        facet_wrap(~ Attribute, scales='free_y') +
        geom_boxplot() +
        scale_y_continuous(limits=c(0,NA)) +
        theme(axis.text.x=element_text(angle=30,hjust=1))
    ggsave(paste(args[1], '_comparison_cracking.pdf', sep=''), device='pdf')
}

# Parallel crack-in-two
{
    print('Parallel crack-in-two')

    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Element==4 & Attribute=='time' & Threads %in% c(2,4,6,8) &
                      (Algorithm=='pirk (MT)' | Algorithm=='rewired (hugepages+MT)'),
                  select=c(Experiment, Algorithm, Threads, Selectivity, Value))
    partitioning <- ddply(sub, c('Experiment', 'Algorithm', 'Threads', 'Selectivity'), plyr::summarise, Value=median(Value))

    merging <- subset(m,
                  Layout=='AoS' & Experiment=='merge' & Element==4 & Attribute=='time' & Threads %in% c(2,4,6,8) &
                      (Algorithm=='pirk (MT)' | Algorithm=='rewired (hugepages+MT)'),
                  select=c(Experiment, Algorithm, Threads, Selectivity, Value))

    partitioning <- within(merge(partitioning, merging, by=c('Algorithm', 'Threads', 'Selectivity')),
                           { Value.x <- Value.x - Value.y })[,c('Experiment.x', 'Algorithm', 'Threads', 'Selectivity', 'Value.x')]
    partitioning <- rename(partitioning, c('Experiment.x'='Experiment', 'Value.x'='Value'))
    all <- rbind(partitioning, merging)
    all$Experiment <- revalue(all$Experiment, c('partition'='parallel partitioning'))

    p_meds <- ddply(prettify(partitioning), c('Experiment', 'Algorithm', 'Threads'), plyr::summarise, med=round(median(Value), 2))
    p_meds$Experiment <- revalue(p_meds$Experiment, c('partition'='parallel partitioning'))

    ggplot(prettify(all), aes(Selectivity, Value, group=Experiment, fill=Experiment)) +
        geom_area(alpha=.7) +
        facet_grid(Algorithm ~ Threads) +
        theme(panel.spacing=unit(1,'lines')) +
        geom_hline(aes(yintercept=mem_bandwidth*1000), size=.5, linetype='F1', show.legend=NA) +
        annotate('text', x=0, y=mem_bandwidth*1000, label='single core bandwidth', hjust=-.05, vjust=1.4, family='serif', colour='black', size=3) +
        labs(y='Time [ms]') +
        theme(legend.position='bottom', legend.direction='horizontal', legend.margin=margin(-11,0,0,0,'pt')) +
        scale_fill_discrete(name='Phase') +
        geom_text(data=p_meds, aes(x=.5, y=0, label=paste(med, 'ms')), vjust=-.5, show.legend=FALSE)
    ggsave(paste(args[1], '_partition_parallel.pdf', sep=''), device='pdf', width=7, height=4.1)

    merging <- subset(merging, Algorithm!='hybrid (hugepages+MT)', select=c(Experiment:Value))
    m_meds <- ddply(prettify(merging), c('Experiment', 'Algorithm', 'Threads'), plyr::summarise, Max=round(max(Value), 2))

    ggplot(prettify(merging), aes(Selectivity, Value, group=Experiment, fill=Experiment)) +
        geom_area(alpha=.7, show.legend=FALSE) +
        facet_grid(Algorithm ~ Threads) +
        theme(panel.spacing=unit(1,'lines')) +
        labs(y='Time [ms]') +
        theme(legend.position='bottom', legend.direction='horizontal') +
        geom_text(data=m_meds, aes(x=.5, y=Max, label=paste('max =', Max, 'ms')), vjust=-.5, show.legend=FALSE) +
        expand_limits(y=160)
    ggsave(paste(args[1], '_partition_parallel_merge.pdf', sep=''), device='pdf', width=7, height=4.1)
}

# Comparison of parallel cracking on all attributes
{
    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Element==4 &
                      (Algorithm=='pirk (MT)' | Algorithm=='rewired (hugepages+MT)'),
                  select=c(Algorithm, Threads, Attribute, Value))

    sub2 <- subset(m,
                   Experiment=='misc' & Algorithm=='memcpy' & Element==4,
                   select=c(Algorithm, Threads, Attribute, Value))
    str(sub2)

    sub <- rbind(sub, sub2)

    ggplot(sub, aes(Threads, Value, colour=Algorithm)) +
        geom_boxplot() +
        facet_wrap(~ Attribute, scales='free_y') +
        expand_limits(y = 0)
    ggsave(paste(args[1], '_comparison_parallel.pdf', sep=''), device='pdf', width=15, height=6)
}

# Comapre Pirk vs. register on all attributes
{
    sub <- subset(m,
                  Layout=='SoA' & Experiment=='partition' & Algorithm=='pirk' & Element==4,
                  select=c(Selectivity, Attribute, Value))

    ggplot(prettify(sub), aes(Selectivity, Value)) +
        facet_wrap(~ Attribute, scales='free') +
        geom_smooth(method='loess', fill=NA) +
        expand_limits(x=0, y=0)
    ggsave(paste(args[1], '_pirk.pdf', sep=''), device='pdf')

    sub <- subset(m,
                  Layout=='SoA' & Experiment=='partition' & Algorithm=='register' & Element==4,
                  select=c(Selectivity, Attribute, Value))

    ggplot(prettify(sub), aes(Selectivity, Value)) +
        facet_wrap(~ Attribute, scales='free') +
        geom_smooth(method='loess', fill=NA) +
        expand_limits(x=0, y=0)
    ggsave(paste(args[1], '_register.pdf', sep=''), device='pdf')
}

# for the presentation, emit an evaluation of branching on 4B AoS alone
{
    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Distribution=='uniform dense' & Algorithm=='branching' &
                      Element=='4' & (Attribute=='time' | Attribute=='brmiss' | Attribute=='stalled'),
                  select=c(Algorithm, Element, Selectivity, Value, Attribute))
    sub[sub$Attribute=='brmiss', 'Value'] <- sub[sub$Attribute=='brmiss', 'Value'] / 1e6
    sub[sub$Attribute=='stalled', 'Value'] <- sub[sub$Attribute=='stalled', 'Value'] / 1e6
    sub[sub$Attribute=='time', 'Value'] <- sub[sub$Attribute=='time', 'Value'] / 1000
    sub$Attribute <- revalue(sub$Attribute, c(brmiss='Branch Misses (Mil.)', stalled='Stalled Cycles (Mil.)', time='Running Time [s]'))

    ggplot(prettify(sub), aes(Selectivity, Value, colour=Algorithm, linetype=Element)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_wrap(~ Attribute, scales='free_y', nrow=1) +
        expand_limits(y=0) +
        scale_linetype_discrete(name='Element size') +
        theme(panel.spacing.x=unit(1, 'lines')) +
        guides(colour=guide_legend(order=1)) +
        scale_x_continuous(breaks=seq(0, 1, .5)) +
        theme(legend.position='bottom', legend.direction='vertical', legend.margin=margin(t = -3, unit='mm')) +
        guides(colour=guide_legend(order=1,title=NULL),linetype=FALSE) +
        theme(plot.background=element_rect(fill="transparent",colour=NA)) +
        theme(legend.background=element_rect(fill="transparent",colour=NA))
    ggsave(paste(args[1], '_eval_branching.pdf', sep=''), device='pdf', width=5, height=2.3)
}
{
    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Distribution=='uniform dense' &
                      (Algorithm=='branching' | Algorithm=='pirk') &
                      Element=='4' & (Attribute=='time' | Attribute=='brmiss' | Attribute=='stalled'),
                  select=c(Algorithm, Element, Selectivity, Value, Attribute))
    sub[sub$Attribute=='brmiss', 'Value'] <- sub[sub$Attribute=='brmiss', 'Value'] / 1e6
    sub[sub$Attribute=='stalled', 'Value'] <- sub[sub$Attribute=='stalled', 'Value'] / 1e6
    sub[sub$Attribute=='time', 'Value'] <- sub[sub$Attribute=='time', 'Value'] / 1000
    sub$Attribute <- revalue(sub$Attribute, c(brmiss='Branch Misses (Mil.)', stalled='Stalled Cycles (Mil.)', time='Running Time [s]'))

    ggplot(prettify(sub), aes(Selectivity, Value, colour=Algorithm, linetype=Element)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_wrap(~ Attribute, scales='free_y', nrow=1) +
        expand_limits(y=0) +
        scale_linetype_discrete(name='Element size') +
        theme(panel.spacing.x=unit(1, 'lines')) +
        scale_x_continuous(breaks=seq(0, 1, .5)) +
        theme(legend.position='bottom', legend.direction='vertical', legend.margin=margin(t = -3, unit='mm')) +
        guides(colour=guide_legend(order=1,title=NULL),linetype=FALSE) +
        theme(plot.background=element_rect(fill="transparent",colour=NA)) +
        theme(legend.background=element_rect(fill="transparent",colour=NA))
    ggsave(paste(args[1], '_eval_pirk.pdf', sep=''), device='pdf', width=5, height=2.5)
}
{
    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Distribution=='uniform dense' & Attribute=='time' &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='vectorized') &
                      Element=='4',
                  select=c(Layout, Algorithm, Element, Selectivity, Value, Attribute))

    summarized_sub <- ddply(sub, c('Layout', 'Algorithm', 'Selectivity', 'Element'), plyr::summarise, time=median(Value))
    summarized_sub$time <- summarized_sub$time / 1000 # to seconds

    ggplot(prettify(summarized_sub), aes(Selectivity, time, colour=Algorithm)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_wrap(~ Layout) +
        geom_hline(aes(yintercept=mem_bandwidth), size=.5, linetype='F1', show.legend=NA) +
        annotate('text', x=0, y=mem_bandwidth, label='single core bandwidth', hjust=-.05, vjust=1.4, family='serif', colour='black', size=3) +
        labs(y='Time [s]') +
        theme(legend.position='bottom', legend.direction='vertical', legend.margin=margin(t = -3, unit='mm')) +
        guides(colour=guide_legend(order=1,title=NULL),linetype=FALSE) +
        expand_limits(y=0) +
        theme(plot.background=element_rect(fill="transparent",colour=NA)) +
        theme(legend.background=element_rect(fill="transparent",colour=NA))
    ggsave(paste(args[1], '_eval_vectorized.pdf', sep=''), device='pdf', width=2.2, height=3.4)
}
{
    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Distribution=='uniform dense' & Attribute=='time' &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='vectorized' | Algorithm=='register') &
                      Element=='4',
                  select=c(Layout, Algorithm, Element, Selectivity, Value, Attribute))

    summarized_sub <- ddply(sub, c('Layout', 'Algorithm', 'Selectivity', 'Element'), plyr::summarise, time=median(Value))
    summarized_sub$time <- summarized_sub$time / 1000 # to seconds

    ggplot(prettify(summarized_sub), aes(Selectivity, time, colour=Algorithm)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_wrap(~ Layout) +
        geom_hline(aes(yintercept=mem_bandwidth), size=.5, linetype='F1', show.legend=NA) +
        annotate('text', x=0, y=mem_bandwidth, label='single core bandwidth', hjust=-.05, vjust=1.4, family='serif', colour='black', size=3) +
        labs(y='Time [s]') +
        theme(legend.position='bottom', legend.direction='vertical', legend.margin=margin(t = -3, unit='mm')) +
        guides(colour=guide_legend(order=1,title=NULL),linetype=FALSE) +
        expand_limits(y=0) +
        theme(plot.background=element_rect(fill="transparent",colour=NA)) +
        theme(legend.background=element_rect(fill="transparent",colour=NA))
    ggsave(paste(args[1], '_eval_register.pdf', sep=''), device='pdf', width=2.2, height=3.7)
}
{
    sub <- subset(m,
                  Layout=='AoS' & Experiment=='partition' & Distribution=='uniform dense' & Attribute=='time' &
                      (Algorithm=='branching' | Algorithm=='pirk' | Algorithm=='vectorized' | Algorithm=='register' |
                       Algorithm=='rewired (hugepages)') &
                      Element=='4',
                  select=c(Layout, Algorithm, Element, Selectivity, Value, Attribute))

    summarized_sub <- ddply(sub, c('Layout', 'Algorithm', 'Selectivity', 'Element'), plyr::summarise, time=median(Value))
    summarized_sub$time <- summarized_sub$time / 1000 # to seconds

    ggplot(prettify(summarized_sub), aes(Selectivity, time, colour=Algorithm)) +
        geom_smooth(method='loess', fill=NA, span=.5) +
        facet_wrap(~ Layout) +
        geom_hline(aes(yintercept=mem_bandwidth), size=.5, linetype='F1', show.legend=NA) +
        annotate('text', x=0, y=mem_bandwidth, label='single core bandwidth', hjust=-.05, vjust=1.4, family='serif', colour='black', size=3) +
        labs(y='Time [s]') +
        theme(legend.position='bottom', legend.direction='vertical', legend.margin=margin(t = -3, unit='mm')) +
        guides(colour=guide_legend(order=1,title=NULL),linetype=FALSE) +
        expand_limits(y=0) +
        theme(plot.background=element_rect(fill="transparent",colour=NA)) +
        theme(legend.background=element_rect(fill="transparent",colour=NA))
    ggsave(paste(args[1], '_eval_rewired.pdf', sep=''), device='pdf', width=2.8, height=4.3)
}

quit()

# Generate a set of plots that are useful for internal analysis.  These are not meant for publication, as they provide
# an unfiltered view on data.
for (layout in levels(m$Layout))
{
    if (layout == 'None')
        next

    for (exp in levels(m$Experiment))
    {
        for (dist in levels(m$Distribution))
        {
            for (sel in levels(factor(m$Selectivity)))
            {
                # Grouped boxplot
                sub <- subset(m, Layout==layout & Experiment==exp & Distribution==dist & Selectivity==sel,
                              select=c(Element, Algorithm, Attribute, Value))

                if (nrow(sub) == 0)
                    next

                print(paste(layout, exp, dist, sel, 'grouped'))
                p <- ggplot(sub, aes(x=Element, y=Value)) +
                    geom_boxplot(show.legend=TRUE, aes(fill=factor(Element)), position=position_dodge(width=.5), width=.6, outlier.shape=NA) +
                    stat_summary(fun.y=median, geom='line', aes(group=1), color='#BBBBBBA0') +
                    stat_summary(fun.data=plot.median, geom='text', size=2, vjust=-1) +
                    facet_grid(Attribute ~ Algorithm, scales='free_y') +
                    geom_blank() +
                    labs(x='Algorithm', y='Property') +
                    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
                    theme(plot.margin=unit(c(1,1,1,1), 'cm')) +
                    scale_fill_discrete(name='Element size (Byte)') +
                    ggtitle(paste('Machine:', args[2], ', Set up:', layout, exp, dist))
                ggsave(paste(args[1], layout, exp, dist, sel, 'grouped.pdf', sep='_'), device='pdf', width=15, height=15)
            }
        }
    }
}

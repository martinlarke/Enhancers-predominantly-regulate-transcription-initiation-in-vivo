# heatmap tss coverage for martin

library(tidyverse)
library(RColorBrewer)
library(cowplot)

# set output dir ----------------------------------------------------------
out.dir <- "~/fusessh/other_projects/heatmaps_for_martin/"
setwd(out.dir)

# load data ---------------------------------------------------------------
df5 <- as.tibble(read.table("local_links/homer.tss.5primecov.txt"))
names(df5) <- c("chr", "start", "end", "tag", "homer.summary", "strand", "pos", "cov")
df5 <- unique(df5)

df3 <- as.tibble(read.table("local_links/homer.tss.3primecov.txt"))
names(df3) <- c("chr", "start", "end", "tag", "homer.summary", "strand", "pos", "cov")
df3 <- unique(df3)

dfM <- as.tibble(read.table("local_links/mnase.cov.txt"))
names(dfM) <- c("chr", "start", "end", "tag", "homer.summary", "strand", "pos", "cov")
dfM <- unique(dfM)

# readjust positional window (+/- 200) ----------------------------------
df5 <- df5 %>%
  mutate(pos = pos - 500) %>%
  filter(pos >= -500 & pos <= 500)

df3 <- df3 %>%
  mutate(pos = pos - 500) %>%
  filter(pos >= -500 & pos <= 500)

dfM <- dfM %>%
  mutate(pos = pos - 500) %>%
  filter(pos >= -500 & pos <= 500)

# # split plus minus ------------------------------------------------------
# p.df5 <- df5 %>% filter(strand == '+')
# m.df5 <- df5 %>% filter(strand == '-')
# 
# p.df3 <- df3 %>% filter(strand == '+')
# m.df3 <- df3 %>% filter(strand == '-')
# 
# p.dfM <- dfM %>% filter(strand == '+')
# m.dfM <- dfM %>% filter(strand == '-')

# invert minus strand
idf5 <- df5 %>% 
  mutate(pos = if_else(strand == "-", -pos, pos))

idf3 <- df3 %>% 
  mutate(pos = if_else(strand == "-", -pos, pos))

idfM <- dfM %>% 
  mutate(pos = if_else(strand == "-", -pos, pos))

# add a ranks (abs cov on 5  prime -----------------------------------
# get ranks from 5' coverage
ranked.idf <- idf5 %>%
  group_by(tag, strand) %>%
  summarise(abs.cov = sum(cov)) %>%
  ungroup() %>%
  mutate(rank.abs.cov = row_number(desc(abs.cov)))

idf5 <- left_join(idf5, ranked.idf, by = c("tag", "strand"))

idf3 <- left_join(idf3, ranked.idf, by = c("tag", "strand"))

idfM <- left_join(idfM, ranked.idf, by = c("tag", "strand"))

# combine 5 and  3 prime ---------------------------------------------
idf <- rbind(idf5 %>% mutate(loc = "5prime"), 
              idf3 %>% mutate(loc = "3prime")) %>% 
  mutate(loc = factor(loc, levels = c("5prime", "3prime")))


# for plotting 5 and 3 prime in the same plot -------------------------
idf.combined <- idf %>%
  mutate(log.cov = log(cov+1)) %>%
  mutate(log.cov = if_else(loc == "3prime", -log.cov, log.cov)) %>% # invert 3prime
  group_by(tag, strand, pos, rank.abs.cov) %>%
  summarise(log.cov.combined = sum(log.cov)) %>%
  ungroup()

# plot ------------------------------------------------------------------
cpal <- brewer.pal(3, "Set1")

p.prime <- idf.combined %>%
  # filter(rank.abs.cov < 100) %>%
  mutate(log.cov.combined = if_else(log.cov.combined > 4, 4, log.cov.combined)) %>%
  mutate(log.cov.combined = if_else(log.cov.combined < -4, -4, log.cov.combined)) %>%
  ggplot(aes(x = pos, y = rank.abs.cov, fill = log.cov.combined)) +
  geom_tile() +
  ggtitle("5' - 3' coverage") +
  labs(fill = "log(5p+1)-log(3p+1)") +
  scale_fill_gradient2(low = cpal[1], mid = "white", high = cpal[2], midpoint = 0) +
  scale_y_reverse() +
  coord_cartesian(xlim = c(-100, 100)) +
  theme(strip.background = element_blank())

# p.prime
# Save ------------------------------------------------------------------  
ggsave(p.prime, filename = "plot_5prime_3prime_merged_limit.png", width = 7.5, height = 15)


# MNase 3prime plot ---------------------------------------------------
# median of 3 prime max
max.3prime.signal <- pull(idf3 %>%
  group_by(pos) %>%
  summarise(cov.sum = sum(cov)) %>%
  filter(cov.sum == max(cov.sum)) %>%
  select(pos))

# select max per tss  
idf3.max <- idf3 %>% 
  group_by(tag) %>%
  filter(cov == max(cov)) %>%
  arrange(rank.abs.cov) %>%
  ungroup()

p.mnase.median <- idfM %>%
  # filter(rank.abs.cov < 100) %>%
  ggplot(aes(x = pos, y = rank.abs.cov, fill = log(cov+1))) +
  geom_tile() +
  # geom_vline(xintercept = max.3prime.signal, linetype = "dotted") +
  ggtitle("MNase vs Median of 3' Max") +
  scale_fill_gradient2(low = "white", mid = "black", high = "black", midpoint = 3) +
  annotate("segment", x = 100, xend = 100, y = 0, yend = -100, colour = "black") +
  annotate("segment", x = -100, xend = -100, y = 0, yend = -100, colour = "black") +
  annotate("segment", x = 100, xend = 100, y = 8000, yend = 8100, colour = "black") +
  annotate("segment", x = -100, xend = -100, y = 8000, yend = 8100, colour = "black") +
  scale_y_reverse() +
  theme(strip.background = element_blank()) + 
  coord_cartesian(xlim=c(-500,500))

# pm.plus
ggsave(p.mnase.median, filename = "plot_3prime_max_vs_mnase_median_of_max2.png", width = 7.5, height = 15)

# p.mnase.dots <- idfM %>%
#   # filter(rank.abs.cov < 100) %>%
#   ggplot(aes(x = pos, y = rank.abs.cov, fill = log(cov+1))) +
#   geom_tile() +
#   geom_point(data = idf3.max, col = cpal[1], size = .1) +
#   ggtitle("MNase vs Max 3' per TSS") +
#   # scale_fill_gradientn(colours = brewer.pal(9, "Greys")) +
#   scale_fill_gradient2(low = "white", mid = "black", high = "black", midpoint = 2) +
#   scale_y_reverse() +
#   theme(strip.background = element_blank()) + 
#   coord_cartesian(xlim=c(-500,500))
# 
# ggsave(p.mnase.dots, filename = "plot_3prime_max_vs_mnase.png", width = 7.5, height = 15)


# combinations -----------------------------------------------
p <- plot_grid(p.prime + theme(
                  legend.pos = "bottom",
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.y = element_blank()
                ) + coord_cartesian(ylim = c(-100, 8100), xlim = c(-100, 100)), 
               p.mnase.median + theme(
                 legend.pos = "bottom",
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.line.y = element_blank()
                 ) + coord_cartesian(ylim = c(-100, 8100), xlim = c(-500, 500)),
               align = "h")
ggsave(p, filename = "plot_combined_5_3_merged_3prime_max_vs_mnase_median_of_max.png", width = 7.5, height = 15)


# Separate 5' 3' MNase pltos for Supplementary ------------------------
idf <- idf %>% filter(cov != 0)

prime5.col <- c(153,0,204)
prime5.col <- rgb(prime5.col[1], prime5.col[2], prime5.col[3], maxColorValue=255)

p.5prime <- idf %>%
  # filter(rank.abs.cov < 100) %>%
  filter(loc == "5prime") %>%
  mutate(log.cov = log(cov+1)) %>%
  mutate(log.cov = if_else(log.cov > 4, 4, log.cov)) %>%
  ggplot(aes(x = pos, y = rank.abs.cov, fill = log.cov)) +
  geom_tile() +
  ggtitle("5' - 3' coverage") +
  labs(fill = "log(5prime + 1)") +
  # scale_fill_gradientn(colours = brewer.pal(9, "Blues")) +
  scale_fill_gradient(low = "white", high = prime5.col) +
  scale_y_reverse() +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 8100)) +
  theme(strip.background = element_blank())

# prime3.col <- c(230,159,0)
prime3.col <- c(213,94,0)
prime3.col <- rgb(prime3.col[1], prime3.col[2], prime3.col[3], maxColorValue=255)

p.3prime <- idf %>%
  # filter(rank.abs.cov < 100) %>%
  filter(loc =="3prime") %>%
  mutate(log.cov = log(cov+1)) %>%
  mutate(log.cov = if_else(log.cov > 4, 4, log.cov)) %>%
  ggplot(aes(x = pos, y = rank.abs.cov, fill = log.cov)) +
  geom_tile() +
  ggtitle("5' - 3' coverage") +
  labs(fill = "log(3prime + 1)") +
  # scale_fill_gradientn(colours = brewer.pal(9, "Reds")) +
  scale_fill_gradient(low = "white", high = prime3.col) +
  scale_y_reverse() +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 8100)) +
  theme(strip.background = element_blank())


psub <- plot_grid(
  p.5prime + theme(
    legend.pos = "left",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    text = element_text(size = 10, family="Arial"),
    legend.text=element_text(size=10)
  ) + coord_cartesian(ylim = c(-100, 8100), xlim = c(-100, 100)),
  p.3prime + theme(
    legend.pos = "left",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    text = element_text(size = 10, family="Arial"),
    legend.text=element_text(size=10)
  ) + coord_cartesian(ylim = c(-100, 8100), xlim = c(-100, 100)),
  p.mnase.median + theme(
    legend.pos = "left",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    text = element_text(size = 10, family="Arial"),
    legend.text=element_text(size=10)
  ) + coord_cartesian(ylim = c(-100, 8100), xlim = c(-500, 500)),
align = "v", nrow = 3)


ggsave(psub, filename = "plot_supplementary_5_3_3prime_max_vs_mnase_median_of_max_5.png", width = 5, height = 20)



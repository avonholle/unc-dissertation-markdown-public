# plot.sim.lgmm.R
# Plot the simulated lgmm
library(ggplot2)
library(reshape)
require(gridExtra)

# read in the data from the C:\Users\vonholle\Documents\dissertation\unc-dissertation-markdown\includes\scripts\power\aim2
# power-aim2.inp Mplus program

setwd("C:/Users/vonholle/Documents/dissertation/unc-dissertation-markdown/includes/scripts/power/aim1.2")

df.1 = read.table("sim2_rev1-or1-5.dat", header=F)
head(df.1)
colnames(df.1)

colnames(df.1) = c("y", "m01", "m02", "m03", "m04", "m05", "m06",
                   "m07", "c")
df.1$id = rownames(df.1)
head(df.1)

# plot data according to class

# melt the data from wide to long
df.1.long = melt(df.1, id=c("y", "c", "id"))
head(df.1.long)
df.1.long$time = as.numeric(with(df.1.long, substr(variable,2,3)))-1
head(df.1.long)
table(df.1.long$time)

table(df.1.long$c)
df.1.long$c.f = factor(df.1.long$c, labels=c("Group 1", "Group 2"))

# plot
sub.1 = df.1.long[df.1.long$id %in% c(seq(1,20,1)),]
sub.1[1:20,]

p1 = ggplot(data=sub.1,
       aes(x=time, y=value, group=id, colour=c.f)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks=unique(sub.1$time)) +
    scale_colour_discrete(name="Latent Class") +
    theme_bw() +
    labs(x="Time (months)", y="Weight (lbs)") 
p1

# add the population average line to the plot
# plot a cubic curve to get an idea of proper coefficients for power analsyis
# (adapted from plot.sim.R)
# ..........................................................
i = 6
s = 2.6
q2 = -0.12
c = -0.005

i2 = 6
s2 = 2
q2.2 = -0.05

age = seq(0, 6, by=1)

r.quad = i + s*age + q2*age^2
r.quad.2 = i2 + s2*age + q2.2*age^2

df.2 = data.frame(age=age, value=r.quad, c.f="Group 1", id=600)
df.2

df.3 = data.frame(age=age, value=r.quad.2, c.f="Group 2", id=601)
df.3

df.4 = rbind.data.frame(df.2, df.3)
df.4

p1 + geom_line(aes(x=age, y=value, group=id, colour=c.f), 
               data=df.4, size=2) +
  theme(legend.position="bottom")

# No groups
# ..............................................................
p2 = ggplot(data=sub.1,
            aes(x=time, y=value, group=id)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks=unique(sub.1$time)) +
  theme_bw() +
  labs(x="Time (months)", y="Weight (lbs)") 
p2

p1.extra =  p1 + geom_line(aes(x=age, y=value, group=id, colour=c.f), 
                           data=df.4, size=2) +
  theme_bw(base_size=18) +
  theme(legend.position="bottom", legend.text=element_text(size=18))


  theme(legend.position="bottom", base_size=18,
        legend.text=element_text(size=18))

grid.arrange(p2, p1.extra, ncol=2)

# output plot

png("sim.lgmm-2.png", width=2400, height=2400, res=300)
#grid.arrange(p2, p1.extra, ncol=2) 
p1.extra
dev.off()





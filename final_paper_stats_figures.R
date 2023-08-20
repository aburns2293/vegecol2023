# code written by Anna Burns (contact: anna.burns@stud.uni-greifswald.de)
# last edited on 18.08.2023

library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(MASS)
library(EnvStats)

# box-cox transformation function from 
# https://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/Transformations.html

BCTransform <- function(y, lambda=0) {
  if (lambda == 0L) { log(y) }
  else { (y^lambda - 1) / lambda }
}

# data available upon request from anna.burns@stud.uni-greifswald.de

veg.allperiods<-read.table("veg_allperiods.csv", header=T, sep=";",dec=".")
head.allperiods<-read.table("head_allperiods.csv", header=T, sep=";",dec=".")
data23<-read.table("all_data_2023.csv", header=T, sep=";",dec=".")
head23<-read.table("header_data2023.csv", header=T, sep=";",dec=".")

# question: How is the relationship between leaf nutrient concentration (N, P) 
# and site soil nutrient concentration?

# data cleaning and tidying

data23 <- data23 %>% mutate(P.vol = bulk_mean*P_cal_soil,
                            K.vol = bulk_mean*K_cal_soil)

## Phosphorous

p.data <- data23 %>% select(plot_id, location, P.vol, P_Holcus, P_Plantago)

p.long <- p.data %>% pivot_longer(!c(plot_id, location, P.vol), 
                                  names_to = "Species",
                                  values_to = "P")

p.long$Species <- str_sub(p.long$Species, 3)
p.long$ID <- str_c(p.long$plot_id,"_",p.long$location,"_",p.long$Species)
  
## Nitrogen

n.data <- data23 %>% select(plot_id, location, Cn_soil, N_Holcus, N_Plantago)

n.long <- n.data %>% pivot_longer(!c(plot_id, location, Cn_soil), 
                                  names_to = "Species",
                                  values_to = "N")

n.long$Species <- str_sub(n.long$Species, 3)
n.long$ID <- str_c(n.long$plot_id,"_",n.long$location,"_",n.long$Species)

## Full dataset

nutrient.combined <- left_join(n.long, p.long, by = join_by(ID))

nutrient.combined <- nutrient.combined %>% select(-c(plot_id.x, location.x, Species.x))
nutrient.combined <- nutrient.combined %>% rename(location = location.y)
nutrient.combined <- nutrient.combined %>% rename(plot_id = plot_id.y)
nutrient.combined <- nutrient.combined %>% rename(Species = Species.y)

nutrient.combined.long <- pivot_longer(nutrient.combined, 
                                       !c(ID, plot_id, Cn_soil, P.vol, location, Species),
                                       names_to = "Nutrient",
                                       values_to = "% in leaves")

# initial values

nutrient.combined.long %>% group_by(Species, Nutrient) %>% 
  summarise(mean = mean(`% in leaves`), sd = sd(`% in leaves`), n=n())

# initial visualizations

ggplot(p.long, aes(x=P.vol, y=P, col = Species, group = Species)) + 
  geom_point() +
  geom_smooth(se=FALSE) +
  theme_classic() +
  labs(title = "Fig. 3: Phosphorous Concentrations in Leaves and Soil", 
       x="Soil P mg/100cm³",
       y="% P in leaf biomass")

ggplot(n.long, aes(x=Cn_soil, y=N, col = Species, group = Species)) + 
  geom_point() +
  geom_smooth(se=FALSE) +
  theme_classic() +
  labs(title = "Fig. 2: Nitrogen Concentrations in Leaves and Soil", 
       x="Soil C/N",
       y="% N in leaf biomass")

ggplot(nutrient.combined.long) +
  geom_boxplot(aes(x = Nutrient, y = `% in leaves`, col = Species)) +
  theme_classic() +
  labs(title = "Fig. 1: Leaf Nutrient Concentrations", y = "% of leaf biomass")

# phosphorous statistics

## holcus

### testing conditions

p.hol <- p.long %>% filter(Species == "Holcus")

p.hol.fit1 <- lm(P~P.vol, data = p.hol)
summary(p.hol.fit1)
plot(p.hol.fit1)

### trying boxcox transformation
p.hol.bc <- boxcox(p.hol.fit1)
p.hol.bc.power <- p.hol.bc$x[which.max(p.hol.bc$y)]

p.hol$p.bc <- BCTransform(p.hol$P, p.hol.bc.power)

p.hol.fit2 <- lm(p.bc~P.vol, data=p.hol)
summary(p.hol.fit2)
plot(p.hol.fit2)

### trying paired logs - happiest with this model

p.hol.fit3 <- lm(log(P)~log(P.vol), data = p.hol)
summary(p.hol.fit3)

plot(p.hol.fit3)

ggplot(p.hol, aes(x=log(P.vol), y=log(P))) + geom_point() + geom_smooth(se=FALSE)
ggplot(p.hol, aes(x=P.vol, y=P)) + geom_point() + geom_smooth(se=FALSE)

### trying gaussian fit

p.hol.fit4 <- glm(P~P.vol, family = "gaussian", data=p.hol)
summary(p.hol.fit4)

plot(p.hol.fit4)

## plantago

### testing conditions

p.plantago <- p.long %>% filter(Species == "Plantago")

p.pla.fit1 <- lm(P~P.vol, data = p.plantago)
summary(p.pla.fit1)
plot(p.pla.fit1)

### trying paired logs - happiest with this model

p.pla.fit2 <- lm(log(P)~log(P.vol), data = p.plantago)
summary(p.pla.fit2)

plot(p.pla.fit2)

ggplot(p.plantago, aes(x=log(P.vol), y=log(P))) + geom_point() + geom_smooth(se=FALSE)
ggplot(p.plantago, aes(x=P.vol, y=P)) + geom_point() + geom_smooth(se=FALSE)

## relationship between species

hol.p.vec <- p.long[p.long$Species == "Holcus",]$P
pla.p.vec <- p.long[p.long$Species == "Plantago",]$P

### correlation test

cor.test(hol.p.vec, pla.p.vec)

### permutation test

p.perm <- twoSamplePermutationTestLocation(hol.p.vec, pla.p.vec, paired = TRUE)

p.perm$p.value

p.perm$estimate

# nitrogen statistics

## holcus

### testing conditions

n.hol <- n.long %>% filter(Species == "Holcus")

ggplot(n.hol, aes(x=Cn_soil, y=N)) + geom_point()

n.hol.fit1 <- lm(N~Cn_soil, data = n.hol)
summary(n.hol.fit1)
plot(n.hol.fit1)
n.hol.bc <- boxcox(n.hol.fit1)

n.hol.bc.power <- n.hol.bc$x[which.max(n.hol.bc$y)]

n.hol$n.bc <- BCTransform(n.hol$N, n.hol.bc.power)

n.hol.fit2 <- lm(n.bc~Cn_soil, data=n.hol)
summary(n.hol.fit2)
plot(n.hol.fit2)

## plantago

n.plantago <- n.long %>% filter(Species == "Plantago")

ggplot(n.plantago, aes(x=Cn_soil, y=N)) + geom_point()

n.pla.fit1 <- lm(N~Cn_soil, data = n.plantago)
summary(n.pla.fit1)
plot(n.pla.fit1)

n.pla.bc <- boxcox(n.pla.fit1)
n.pla.bc.power <- n.pla.bc$x[which.max(n.pla.bc$y)]

n.plantago$n.bc <- BCTransform(n.plantago$N, n.pla.bc.power)

n.pla.fit2 <- lm(n.bc~Cn_soil, data=n.plantago)
summary(n.pla.fit2)
plot(n.pla.fit2)

## relationship between species

hol.n.vec <- n.long[n.long$Species == "Holcus",]$N
pla.n.vec <- n.long[n.long$Species == "Plantago",]$N

## correlation test

cor.test(hol.n.vec, pla.n.vec)

## permutation test

n.perm <- twoSamplePermutationTestLocation(hol.n.vec, pla.n.vec, paired = TRUE)

n.perm$p.value

n.perm$estimate

# relationship between N/P

## leaf N - soil P

ggplot(nutrient.combined.long[nutrient.combined.long$Nutrient == "N",], 
       aes(x = P.vol, y = `% in leaves`, col = Species, group = Species)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_classic() +
  xlab("Soil P mg/100 cm³") +
  ylab("Leaf N (% of biomass)") +
  ggtitle("Relationship between leaf N and soil P concentrations")

hol.leafn.soilp <- nutrient.combined %>% filter(Species == "Holcus")

hol.leafn.soilp$N.bc <- BCTransform(hol.leafn.soilp$N, -1.4343434)

hol.leafn.soilp.fit1 <- lm(N.bc ~ Cn_soil + P.vol, 
                           data = hol.leafn.soilp)
summary(hol.leafn.soilp.fit1)
plot(hol.leafn.soilp.fit1)

pla.leafn.soilp <- nutrient.combined %>% filter(Species == "Plantago")
pla.leafn.soilp.fit1 <- lm(N ~ Cn_soil + P.vol, data = pla.leafn.soilp)
summary(pla.leafn.soilp.fit1)
plot(pla.leafn.soilp.fit1)

## leaf P - soil C/N

ggplot(nutrient.combined.long[nutrient.combined.long$Nutrient == "P",], 
       aes(x = Cn_soil, y = `% in leaves`, col = Species, group = Species)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_classic() +
  xlab("Soil C/N ratio") +
  ylab("Leaf P (% of biomass)") +
  ggtitle("Relationship between leaf P and soil C/N concentrations")

hol.leafp.soilcn <- nutrient.combined %>% filter(Species == "Holcus")

hol.leafp.soilcn.fit1 <- lm(log(P) ~ log(Cn_soil), 
                           data = hol.leafp.soilcn)
summary(hol.leafp.soilcn.fit1)
plot(hol.leafp.soilcn.fit1)

pla.leafp.soilcn <- nutrient.combined %>% filter(Species == "Plantago")
pla.leafp.soilcn.fit1 <- lm(P ~ Cn_soil, data = pla.leafp.soilcn)
summary(pla.leafp.soilcn.fit1)

## soil P - soil C/N

ggplot(nutrient.combined.long, 
       aes(x = Cn_soil, y = P.vol)) +
  geom_point() +
  geom_smooth(se = FALSE, col = "black") +
  theme_classic() +
  xlab("Soil C/N ratio") +
  ylab("Soil P mg/100 cm³") +
  ggtitle("Fig. 5: Relationship between soil P and C/N")

## plant P - plant N

ggplot(nutrient.combined, 
       aes(x = N, y = P, group = Species, col = Species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  xlab("Leaf N (% of biomass)") +
  ylab("Leaf P (% of biomass)") +
  ggtitle("Fig. 6: Relationship between leaf P and N")

pla.pn <- nutrient.combined %>% filter(Species == "Plantago")

pla.pn.fit1 <- lm(P~N, data = pla.pn)
plot(pla.pn.fit1)
summary(pla.pn.fit1)

hol.pn <- nutrient.combined %>% filter(Species == "Holcus")

hol.pn.fit1 <- lm(P~N, data = hol.pn)
plot(hol.pn.fit1)
summary(hol.pn.fit1)

# c/n:p ratio

nutrient.combined2 <- nutrient.combined %>% mutate(soil_ratio = Cn_soil/P)

nutrient.combined.long2 <- pivot_longer(nutrient.combined, 
                                        !c(ID, plot_id, soil_ratio, Cn_soil, P.vol, location, Species),
                                                                  names_to = "Nutrient",
                                                                  values_to = "% in leaves")

ggplot(nutrient.combined.long2, aes(x = soil_ratio, y = `% in leaves`, col = Species, group = Species)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(vars(Nutrient)) +
  theme_classic() +
  xlab("soil C/N : soil P volume (mg/100 cm³)") +
  ylab("% of leaf biomass") +
  ggtitle("Fig. 7: Relationship between soil C/N : P ratio and leaf nutrients")

## holcus

### phosphorous

hol.ratio.p <- nutrient.combined.long2 %>% filter(Species =="Holcus" &  Nutrient == "P")

hol.ratio.p.fit1 <- lm(`% in leaves`~ soil_ratio, data = hol.ratio.p)
summary(hol.ratio.p.fit1)
plot(hol.ratio.p.fit1)

hol.ratio.p.bc <- MASS::boxcox(hol.ratio.p.fit1)
hol.ratio.p.bc.power <- hol.ratio.p.bc$x[which.max(hol.ratio.p.bc$y)]

hol.ratio.p$P.bc <- BCTransform(hol.ratio.p$`% in leaves`, hol.ratio.p.bc.power)

hol.ratio.p.fit2 <- lm(P.bc~soil_ratio, data = hol.ratio.p)
summary(hol.ratio.p.fit2)
plot(hol.ratio.p.fit2)

### nitrogen

hol.ratio.n <- nutrient.combined.long2 %>% filter(Species =="Holcus" &  Nutrient == "N")

hol.ratio.n.fit1 <- lm(`% in leaves`~ soil_ratio, data = hol.ratio.n)
summary(hol.ratio.n.fit1)
plot(hol.ratio.n.fit1)

## plantago

### phosphorous

pla.ratio.p <- nutrient.combined.long2 %>% filter(Species =="Plantago" &  Nutrient == "P")

pla.ratio.p.fit1 <- lm(`% in leaves`~ soil_ratio, data = pla.ratio.p)
summary(pla.ratio.p.fit1)
plot(pla.ratio.p.fit1)

pla.ratio.p.bc <- MASS::boxcox(pla.ratio.p.fit1)
pla.ratio.p.bc.power <- pla.ratio.p.bc$x[which.max(pla.ratio.p.bc$y)]

pla.ratio.p$P.bc <- BCTransform(pla.ratio.p$`% in leaves`, pla.ratio.p.bc.power)

pla.ratio.p.fit2 <- lm(P.bc~soil_ratio, data = pla.ratio.p)
summary(pla.ratio.p.fit2)

### nitrogen

pla.ratio.n <- nutrient.combined.long2 %>% filter(Species =="Plantago" &  Nutrient == "N")

pla.ratio.n.fit1 <- lm(`% in leaves`~ soil_ratio, data = pla.ratio.n)
summary(pla.ratio.n.fit1)
plot(pla.ratio.n.fit1)



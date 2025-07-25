library(ggplot2)
library(MuMIn)
library(lme4)
library(DHARMa)
library(effects)
library(dplyr)
library(mgcv)
library(tidyr)
library(dplyr)
library(writexl)
library(ggplot2)
library(devtools)
library(ggtext)
library(gridExtra)
library(grid)
library(corrplot)
library(car)
library(gratia)
library(MASS)
library(emmeans)
library(margins)
library(lmerTest)
library(ggeffects)
library(parameters)

#############
###fig s1
#plot of spat through time to demonstrate July 2nd settlement date

ScallopSettlement <- H0

ScallopSettlement$Year <- as.character(ScallopSettlement$Year)
ScallopSettlement$Date <- as.Date(ScallopSettlement$Date)

#per chat gpt, these colors are color blind friendly
ScallopSettlementPlot<-ggplot(ScallopSettlement, aes(x=Date, y=AverageSpatAbundance)) +
  geom_line(aes(color = Year), size = .5, alpha=0.75) +
  scale_color_manual(values = c("gold", "#66c2a5", "#fc8d62", "#8da0cb")) +
  geom_point(aes(color = Year), size = 1, alpha=0.9)+
  xlab("Date")+
  ylab("Average Spat Abundance (count per bag)")+
  stat_smooth(aes(x = Date, y = AverageSpatAbundance), method = "lm",
             formula = y ~ poly(x, 4), se = FALSE, linetype = "dashed", color="gray25") +
  theme_bw()+
  theme(legend.position = c(.8, 0.8))
ScallopSettlementPlot
  
#ESA rules: 8.5 cm, max = 18 by 24 tall
ggsave("ScallopSettlementPlot.tiff", ScallopSettlementPlot, dpi = 600, bg = "white",
       width = 18,
       height = 14,
       units = "cm")


#############
###fig 2 (H1)

#mixed nb glm with year as a catgorical random factor, model did not converge with site as a factor
H1$year <- as.character(H1$year)

SAVscallopsurveyWeight <- H1

#site would not converge 
m<-glmer.nb(scallop_freq ~ species_grp * spun_wt +(1|year), data=SAVscallopsurveyWeight)

res1 <- simulateResiduals(m)
plot(res1)

summary(m)
Anova(m, type = "III")
print(m, corr = FALSE)

multip<-emmeans(m, ~  species_grp)
pairs(multip)

species_colors <- c("Codium" = "darkgreen", 
                    "Thick Fleshy Red" = "magenta4", 
                    "Filamentous Red" = "lightcoral")

pr <- ggpredict(m, c("spun_wt", "species_grp"))

pr <- pr %>%
  rename(species_grp = group)

pr_filtered <- pr %>%
  filter((species_grp == "Filamentous Red" & x <= 34.5) |
           (species_grp == "Thick Fleshy Red" & x <= 15) |
           species_grp == "Codium")

scallopfreqSAVlineplot <- ggplot(pr_filtered, aes(x = x, y = predicted, color = species_grp, fill = species_grp, linetype = species_grp)) +
  geom_line() +
  geom_ribbon(aes(ymin = pmin(conf.low, 15), ymax = pmin(conf.high, 15), fill = species_grp), alpha = 0.3, color = NA) +
  geom_point(data = SAVscallopsurveyWeight, aes(x = spun_wt, y = scallop_freq, color = species_grp), show.legend = FALSE) +
  labs(x = "Spun Weight (mg)", y = "SAV Scallop Frequency (counts per plant)", color = "SAV Group", fill = "SAV Group", linetype = "SAV Group") +
  theme_bw() +
  scale_color_manual(values = species_colors, 
                     labels = c("Spongy Green", "Filamentous Red", "Thick Fleshy Red")) +
  scale_fill_manual(values = species_colors, 
                    labels = c("Spongy Green", "Filamentous Red", "Thick Fleshy Red")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("Spongy Green", "Filamentous Red", "Thick Fleshy Red")) +
  guides(color = guide_legend(title = "SAV Group"), 
         fill = guide_legend(title = "SAV Group"),
         linetype = guide_legend(title = "SAV Group")) +
  theme(legend.position = c(.7, 0.8)) +
  coord_cartesian(ylim = c(0, 15), xlim=c(0,60)) 
scallopfreqSAVlineplot

# repeat above but for surface area, not weight
m_sa<-glm.nb(scallop_freq ~ species_grp * surface_area, data=SAVscallopsurveyWeight)

res1 <- simulateResiduals(m_sa)
plot(res1)

summary(m_sa)
Anova(m_sa, type = "III")
print(m_sa, corr = FALSE)

multip_sa<-emmeans(m_sa, ~  species_grp)
pairs(multip_sa)

pr_sa <- ggpredict(m_sa, c("surface_area", "species_grp"))

pr_sa <- pr_sa %>%
  rename(species_grp = group)

# making plot to have on reserve but we wont inculde in final figure bc its ugly and shows next to nothong

m<-glm.nb(scallop_freq ~ species_grp, data=SAVscallopsurveyWeight)
summary(m)
Anova(m, type = "II")
print(m_sa, corr = FALSE)

multip<-emmeans(m, ~  species_grp)
pairs(multip)


scallopfreqSAV_SAplot<-ggplot(pr_sa, aes(x = x, y = predicted, color = species_grp, fill = species_grp, linetype = species_grp)) +
  geom_line() +
  geom_ribbon(aes(ymin = pmin(conf.low, 15), ymax = pmin(conf.high, 15), fill = species_grp), alpha = 0.3, color = NA) +
  geom_point(data = SAVscallopsurveyWeight, aes(x = surface_area, y = scallop_freq, color = species_grp), show.legend = FALSE) +
  labs(x = "Surface Area (cm^2)", y = "SAV Scallop Frequency (counts per plant)", color = "SAV Group", fill = "SAV Group", linetype = "SAV Group") +
  theme_bw() +
  scale_color_manual(values = species_colors, 
                     labels = c("Spongy Green", "Filamentous Red", "Thick Fleshy Red")) +
  scale_fill_manual(values = species_colors, 
                    labels = c("Spongy Green", "Filamentous Red", "Thick Fleshy Red")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("Spongy Green", "Filamentous Red", "Thick Fleshy Red")) +
  guides(color = guide_legend(title = "SAV Group"), 
         fill = guide_legend(title = "SAV Group"),
         linetype = guide_legend(title = "SAV Group")) +
  theme(legend.position = "none")+
coord_cartesian(ylim = c(0, 15), xlim=c(0,10000))  


# box plot for sav group 
scallopfreqSAVboxplot <- ggplot(SAVscallopsurveyWeight, aes(x = species_grp, y = scallop_freq, color = species_grp)) +
  geom_boxplot(outlier.colour = "white") + 
  geom_point(position = position_jitter(width = 0.15), alpha = 0.5, size = 2) +  
  theme_bw() +
  guides(color = "none") +
  scale_color_manual(values = species_colors, 
                    labels = c("Spongy Green", "Filamentous Red", "Thick Fleshy Red")) +
  labs(x = "SAV Group", y = "SAV Scallop Frequency (counts per plant)")+
  scale_x_discrete(labels = c("Codium" = "Spongy Green", 
                              "Filamentous Red" = "Filamentous Red", 
                              "Thick Fleshy Red" = "Thick Fleshy Red"))
scallopfreqSAVboxplot

scallopfreqSAV<-grid.arrange(scallopfreqSAVlineplot, scallopfreqSAVboxplot, 
                                ncol = 1, nrow = 2)

ggsave("scallopfreqSAV.tiff",scallopfreqSAV, dpi = 600, bg = "white",
       width = 18,
       height = 22,
       units = "cm")


#############
###fig 3 (H2 and H3)

#hog neck and SB only
#remove unneeded columns

SAVscallopsurveysHNSB <- H2

SAVscallopsurveysHNSB2 <- SAVscallopsurveysHNSB[, c("Year", "date", "dayspost", "Site", "sample", "Species", "Species_Group", "attachment_ht", "canopy_ht")]

#remove NAs for attachment height
SAVscallopsurveysHNSB3 <- SAVscallopsurveysHNSB2[!is.na(SAVscallopsurveysHNSB2$attachment_ht), ]
SAVscallopsurveysHNSB3$Year <- as.character(SAVscallopsurveysHNSB3$Year)

#residuals are bad when using year as a random factor
m<-lmer(attachment_ht ~ Species_Group+(1|Site),
      data=SAVscallopsurveysHNSB3)
summary(m)

#residuals look good
res1 <- simulateResiduals(m)
plot(res1)

m<-aov(attachment_ht ~ Error(Site)+ Species_Group,
       data=SAVscallopsurveysHNSB3)
summary(m)

#non parametric test (since data are not balanced) give the same outcomes as the anova
kruskal.test(attachment_ht ~ Site, data = SAVscallopsurveysHNSB3)
kruskal.test(attachment_ht ~ Species_Group, data = SAVscallopsurveysHNSB3)

#not color blind friendly, but colors are redundant with the x axis labels
species_colors2 <- c("Codium" = "darkgreen", 
                     "Grinnellia" = "violetred2",
                     "Thick Fleshy Red" = "magenta4", 
                     "Filamentous Red" = "lightcoral")

speciesattachheightboxplot<-ggplot(SAVscallopsurveysHNSB3, aes(x = Species_Group, y = attachment_ht, color=Species_Group)) +
  geom_boxplot(outlier.colour = "white") +  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 8.4 - 1.02, ymax = 8.4 + 1.02,
           fill = "darkgreen", alpha = 0.15) + 
  geom_point(position = position_jitter(width = 0.15), alpha = 0.5, size = 0.8) +  
  theme_bw() +
  guides(color = "none")+
  scale_color_manual(values = species_colors2) +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank()) + 
  labs(y = "Scallop Attachment Height (cm)") +
  ylim(0, 13.5) +
  geom_hline(yintercept = 8.4, linetype = "dashed", color = "darkgreen", size = 1)+
  annotate("text", x = 2.1, y = 8.5, label = "'Eelgrass mean attachment height'", 
           parse = TRUE, hjust = -0.1, vjust = -0.1, color = "darkgreen", size = 3.6)
speciesattachheightboxplot

#working with original file with all sites (ie sites other than HN and SB)
#remove NAs
SAVscallopsurveys <- H3
SAVscallopsurveys2 <- SAVscallopsurveys[!is.na(SAVscallopsurveys$canopy_ht), ]
SAVscallopsurveysHNSB3$Year <- as.character(SAVscallopsurveysHNSB3$Year)

#remove species for which there are two of fewer observstions
SAVscallopsurveys2 <- SAVscallopsurveys2 %>%
  group_by(species) %>%
  filter(n() > 2) %>%
  ungroup()

#log transform bc residuals looked bad
SAVscallopsurveys2 <- SAVscallopsurveys2 %>%
  mutate(log_canopy_ht = log(canopy_ht))

#get the mean zostera height to print on plot as a reference point
mean_canopy_ht_Zostera <- SAVscallopsurveys2 %>%
  filter(species == "Zostera") %>%
  summarise(mean_canopy_ht = mean(canopy_ht, na.rm = TRUE)) %>%
  pull(mean_canopy_ht)
print(mean_canopy_ht_Zostera) # it is 36.3

#now get only the sites of interest
SAVscallopsurveys3 <- SAVscallopsurveys2 %>%
  filter(site %in% c("SB", "HN"))

#bingo, bootiful residuals
m<-lmer(log_canopy_ht ~ species_grp+(1|site),
      data=SAVscallopsurveys3)
res1 <- simulateResiduals(m)
plot(res1)
summary(m)

m<-aov(log_canopy_ht ~ Error(site)+ species_grp,
       data=SAVscallopsurveys3)
summary(m)

multip<-emmeans(m, ~  species_grp)
pairs(multip)

#non parametric test (since data are not balanced) give the same outcomes as the anova
kruskal.test(canopy_ht ~ site, data = SAVscallopsurveys3)
kruskal.test(canopy_ht ~ species_grp, data = SAVscallopsurveys3)

#create custom labels 
custom_labels <- levels(SAVscallopsurveys3$species_grp)
names(custom_labels) <- custom_labels
custom_labels["Codium"] <- expression("Spongy Green")
custom_labels["Grinnellia"] <- expression("Folioise Red")

#not color blind friendly, but colors are redundant with the x axis labels
species_colors2 <- c("Codium" = "darkgreen", 
                    "Grinnellia" = "violetred2",
                    "Thick Fleshy Red" = "magenta4", 
                    "Filamentous Red" = "lightcoral")

Canopy_SpeciesGroupPlot<-ggplot(SAVscallopsurveys3, aes(x = species_grp, y = canopy_ht, color=species_grp)) +
  geom_boxplot(outlier.colour = "white") +  # Add boxplot
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 36.3 - 2.2, ymax = 36.3 + 2.2,
           fill = "darkgreen", alpha = 0.15) + 
  geom_point(position = position_jitter(width = 0.15), alpha = 0.5,size=.8) +  
  theme_bw()+
  theme(axis.title = element_text(size = 11))+
  labs(x = "SAV Species Group", y = "SAV Canopy Height (cm)")+
  scale_y_continuous(limits = c(0,40))+
  scale_x_discrete(labels = custom_labels) +
  guides(color = "none")+
  scale_color_manual(values = species_colors2) +
  geom_hline(yintercept = 36.3, linetype = "dashed", color = "darkgreen", size = 1)+
  annotate("text", x = 2.11, y = 36.6, label = "'Eelgrass mean canopy height'", 
           parse = TRUE, hjust = -0.1, vjust = -0.1, color = "darkgreen", size = 3.6)
Canopy_SpeciesGroupPlot

SpeciesHeightPlot_error<-grid.arrange(speciesattachheightboxplot, Canopy_SpeciesGroupPlot,
                                ncol = 1, nrow = 2)

ggsave("SpeciesHeightPlot_error.tiff", SpeciesHeightPlot_error, dpi = 600, bg = "white",
       width = 18,
       height = 22,
       units = "cm") 


############
#### H4
#There were no differences in position of scallop attachment to SAVs (periphery vs middle)

SAVscallopposition<- H4
chisq.test(table(SAVscallopposition$position))


#############
###fig 4 (H5 and H6)

#hog neck and SB only
SAVscallopfreqHNSB  <- H5

# avg of all sites 
SAVscallopfreqHNSB2 <- SAVscallopfreqHNSB %>%
  group_by(dayspost, Year) %>%
  summarize(scallop_freq = mean(scallop_freq, na.rm = TRUE))

SAVscallopfreqHNSB2$Year <- as.factor(SAVscallopfreqHNSB2$Year)

#make model
m<-gam(scallop_freq ~ s(dayspost), tw(), data=SAVscallopfreqHNSB2)
res1 <- simulateResiduals(m)
plot(res1)
summary(m)
DS <- data_slice(m, dayspost = evenly(dayspost, n = 100))
FV <- fitted_values(m, data = DS, scale = "response")

#make plot 
savscallopfreqplot<-ggplot(FV, aes(x = dayspost, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3, fill = "gray33") +
  geom_line(size = 0.5,colour="black") +
  geom_point(data = SAVscallopfreqHNSB2, 
             aes(x = dayspost, y = scallop_freq,color = Year), 
             size = 3) +   xlab("Days Post July 2")+
  ylab("SAV Canopy Scallop Frequency (counts per plant)")+
  theme_bw()+
  theme(legend.position = c(0.85, 0.75),
        axis.title.x = element_blank()      
  ) +
  scale_color_manual(values = c("gold", "#66c2a5", "#fc8d62", "#8da0cb")) 
savscallopfreqplot

#only includes HN and SB
#this df differs from the SAVscallopsurveysHNSB in that it only includes observations for which 
#there are shell height values

SAVscallopshellheightHNSB2 <-H6 %>%
  group_by(date, site, Year, dayspost) %>%
  summarize(shell_ht = mean(shell_ht))


SAVscallopshellheightHNSB2 <- subset(SAVscallopshellheightHNSB2, dayspost <= 35)

SAVscallopshellheightHNSB2$Year <- as.character(SAVscallopshellheightHNSB2$Year)

SAVscallopshellheightHNSB2 <- SAVscallopshellheightHNSB2[SAVscallopshellheightHNSB2$dayspost >= 0, ]



#doesnt converge with the site factor/estimtes site factor to be zero, so remove
m<-glmer(shell_ht ~ dayspost + (1 | site) + (1 | Year), family = gaussian(link = "log"), data=SAVscallopshellheightHNSB2)
# does converge
m<-glmer(shell_ht ~ dayspost + (1 | Year), family = gaussian(link = "log"), data=SAVscallopshellheightHNSB2)

summary(m)
r.squaredGLMM(m)

res1 <- simulateResiduals(m)
plot(res1)

efct <- effect("dayspost", mod = m, xlevels = list(dayspost = seq(0, 35, length.out = 100)))

efct <- as.data.frame(efct)

savscallopshellheightplot <- ggplot() +
  geom_ribbon(data = efct, aes(x = dayspost, ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(data = efct, aes(x = dayspost, y = fit)) +
  geom_point(data = SAVscallopshellheightHNSB2, aes(x = dayspost, y = shell_ht, colour = Year, shape = site), size = 3) +
  xlab("Days Post July 2") +
  ylab("SAV Canopy Shell Height (mm)") +
  theme_bw() +
  labs(shape = "Site", colour = "Year") +
  scale_color_manual(values = c("gold", "#66c2a5", "#fc8d62", "#8da0cb")) +
  guides(
    shape = guide_legend(order = 1), 
    colour = guide_legend(order = 2) 
  ) +
  theme(
    legend.position = c(.2, 0.8),
    legend.box = "horizontal"  
  )

savscallopshellheightplot

SAVscallopYEARplot<-grid.arrange(savscallopfreqplot, savscallopshellheightplot,
                             ncol = 1, nrow = 2)

ggsave("SAVscallopYEARplot.tiff",SAVscallopYEARplot, dpi = 600, bg = "white",
       width = 18,
       height = 22,
       units = "cm")


#############
###fig 5
ordered_levels <- c("Foliose Red","Filamentous Red","Thick Fleshy Red","Foliose Green", "Spongy Green", "Spongy Green*",
                    "Eelgrass", "Eelgrass†") 

maxshellheights <- F5
#apply  specified order to the savtype factor
maxshellheights <- maxshellheights %>%
  mutate(savtype = factor(savtype, levels = ordered_levels))

maxshellheightsplot <- ggplot(maxshellheights, aes(x = savtype, y = barheight, fill = savgroup)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5, color = "black") +  
  geom_text(aes(label = ninbar), vjust = 1.3, size = 4) +  
  geom_text(aes(label = nabovebar), vjust = -.3, size = 5) +  
  ylim(0, 50) +  
  xlab("SAV Group") + 
  ylab("Maximum Size of Scallops (mm) Attached to Given SAV") + 
  theme_bw() +
  scale_fill_manual(
    values = c(
      "chloro" = "darkgreen",        
      "rhodo" = "lightcoral",        
      "grass" = "darkolivegreen3" )
  ) +
  guides(fill = "none") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)  
  )
maxshellheightsplot

ggsave("maxshellheightplot.tiff", maxshellheightsplot, dpi = 600, bg = "white",
       width = 18,
       height = 14,
       units = "cm")


#############
###fig 6 (H7 and H8)

#number of benthic scallops increased after July 2
#only includes HN and SB
#poisson dist for count data
benthicscallopfreq <- H7
benthicscallopfreq <- subset(benthicscallopfreq, dayspost <= 37)

benthicscallopfreq $Year <- as.character(benthicscallopfreq$Year)

m<-glmer(scallop_freq ~ dayspost + (1 | site) + (1 | Year), family=poisson(), data=benthicscallopfreq)
summary(m)

res1 <- simulateResiduals(m)
plot(res1)

r.squaredGLMM(m)

efct <- effect("dayspost", mod = m, xlevels = list(dayspost = seq(0, 37, length.out = 100)))

efct <- as.data.frame(efct)

benthicscallopfreqplot <- ggplot() +
  geom_ribbon(data = efct, aes(x = dayspost, ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(data = efct, aes(x = dayspost, y = fit)) +
  geom_point(data = benthicscallopfreq, aes(x = dayspost, y = scallop_freq, color = Year, shape = site), size = 3) +
  labs(
    x = "Days Post July 2", 
    y = expression(Benthic~Scallop~Frequency~(counts~per~0.08~m^2)), 
    color = "Year", 
    shape = "Site"
  ) +
  theme_bw() +
  theme(
    legend.position = c(.7, 0.85),
    legend.box = "horizontal",  
    axis.title.x = element_blank(),   
    axis.text.x = element_blank()     
  ) +
  scale_color_manual(values = c("#fc8d62", "#8da0cb")) +
  guides(
    shape = guide_legend(order = 1),  
    color = guide_legend(order = 2)  
  )
benthicscallopfreqplot    

#############
#hog neck only bc there were only four other observations for other sites

benthicscallopSHclean <- H8

benthicscallopSHclean$Year <- as.character(benthicscallopSHclean$Year)

benthicscallopSHclean<- benthicscallopSHclean %>%
  filter(site == "Hog Neck")

m<-lmer(shell_ht ~ dayspost +(1 | Year),data=benthicscallopSHclean)
summary(m)

r.squaredGLMM(m)
plot(m)
car::Anova(m, type=3)

res1 <- simulateResiduals(m)
plot(res1)

efct <- effect("dayspost", mod = m, xlevels = list(dayspost = seq(0, 82, length.out = 100)))

efct <- as.data.frame(efct)

benthicscallopSHplot <- ggplot() +
  geom_ribbon(data = efct, aes(x = dayspost, ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(data = efct, aes(x = dayspost, y = fit)) +
  geom_point(data = benthicscallopSHclean, aes(x = dayspost, y = shell_ht, colour = Year), size = 3) +
  ylab("Benthic Scallop Shell Height (mm)") +
  xlab("Days Post July 2")+
  theme_bw() +
  theme(
    legend.position = c(.2, 0.75)) +
  scale_color_manual(values = c("gold", "#66c2a5", "#fc8d62", "#8da0cb"))

benthicscallopSHplot

benthicplot<-grid.arrange(benthicscallopfreqplot,  benthicscallopSHplot,
                                 ncol = 1, nrow = 2)

ggsave("benthicplot2.tiff",benthicplot, dpi = 600, bg = "white",
       width = 18,
       height = 22,
       units = "cm")



#### LAB EXPERIMENTS


##############
#hypothesis 9
#Lab Survival of juvenile scallops with mud crab predators varied by size and substrate
MudCrabSAVeaten <- H9
m<-glmer(scallop_eaten ~ Substrate * SIZE + (1|TANK), family=binomial, data=MudCrabSAVeaten)
summary(m)
car::Anova(m, type=3)

res1 <- simulateResiduals(m)
plot(res1)

r.squaredGLMM(m)


#create interaction term to more easily get mult comparisons output
library(multcomp)
MudCrabSAVeaten$SubBySize <- interaction(MudCrabSAVeaten$Substrate, MudCrabSAVeaten$SIZE)
m<-glmer(scallop_eaten ~ SubBySize + (1|TANK), family=binomial, data=MudCrabSAVeaten)
summary(glht(m, mcp(SubBySize = "Tukey")))

scallop_summary <- MudCrabSAVeaten %>%
  group_by(SIZE, Substrate) %>%
  summarise(percent_eaten = mean(scallop_eaten) * 100) %>%
  ungroup()

ScallopsEatenPlot <- ggplot(scallop_summary, aes(x = SIZE, y = percent_eaten, fill = Substrate)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("Codium" = "darkgreen", "Zostera" = "darkseagreen3", "Sand" = "tan")) +
  labs(x = "Size Group", y = "Percent Scallops Eaten (%)") +
  ylim(0, 35) +
  scale_x_discrete(labels = c("B" = "5-7 mm", "C" = "9-12 mm")) +
  theme_bw()+
  theme(legend.position = c(.1, 0.5))
ScallopsEatenPlot

ggsave("ScallopsEatenPlot2.tiff", ScallopsEatenPlot, dpi = 600, bg = "white",
       width = 18,
       height = 14,
       units = "cm")


##############
#hypothesis 10
#Lab Height of attachment of scallops did not differ by SAV or with presence/absence of mud crabs 

MudCrabSAVattchht<- H10
m<-lm(ATTACHMENT ~ SAV + CRAB + SAV * CRAB, data=MudCrabSAVattchht)
res1 <- simulateResiduals(m)
plot(res1)
summary(m)
m2<-aov(m)
summary(m2)

m<-lm(ATTACHMENT ~ SAVBYCRAB, data=MudCrabSAVattchht)
m2<-aov(m)
TukeyHSD(m2)


MudCrabSAVattchht$SAV <- factor(MudCrabSAVattchht$SAV, levels = c("Zostera", "Codium"))

MudCrabAttachPlot<-ggplot(MudCrabSAVattchht, aes(x = factor(CRAB), y = ATTACHMENT, color = SAV)) +
  geom_boxplot(outlier.colour = "white", position = position_dodge(width = 0.95)) +  # Add boxplot with dodging
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.95), alpha = 0.7, size = 0.9) +  
  scale_color_manual(values = c("Codium" = "darkgreen", "Zostera" = "darkseagreen3")) +
  labs(x = "Crab Presence (No = Absence, Yes = Presence)", y = "Scallop Attachment Height in SAV (cm)") +
  theme_bw()+
  ylim(0, 24) +
  theme(legend.position = c(.1, 0.89))
MudCrabAttachPlot

ggsave("MudCrabAttachPlot2.png", MudCrabAttachPlot, dpi = 200, bg = "white",
       width = 18,
       height = 14,
       units = "cm")





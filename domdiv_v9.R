library(tidyverse)
library(plyr)
library(ggeffects)
library(ggthemes)
library(MuMIn)
library(nlme)
library(vegan)
library(lmerTest)
library(codyn)
library(emmeans)

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", force = TRUE)
library(pairwiseAdonis)


removal_full <- read.csv("C:/Users/ohler/Desktop/Sevilleta_RemovalPlots_biomass_19Nov2022.csv")
removal <- removal_full%>%
  subset(kartez != "litter")%>%
  subset(kartez != "soil")%>%
  subset(kartez != "total")%>%
  subset( treatment != "R_all" & kartez != "litter" & kartez != "soil")%>%
  subset( site == "removal_blue" | site == "removal_black")


#calculate without compensation
compensation_dummy <- subset(removal, treatment =="control")

#remove dominant species 
compensation_dummy$subordinate_species <- ifelse(compensation_dummy$site == "removal_black" & compensation_dummy$kartez == "BOER4", NA, ifelse(compensation_dummy$site == "removal_blue" &compensation_dummy$kartez == "BOGR2",NA,compensation_dummy$biomass.BIM))

#re-name treatment to w/o compensation
compensation_dummy$treatment <- "subordinate_species"
compensation_dummy <- compensation_dummy[!is.na(compensation_dummy$subordinate_species),]
compensation_dummy <- unite(compensation_dummy, "plot", c("plot", "treatment"),remove = FALSE)

#sum biomass
removal.1 <- merge(compensation_dummy, removal, by = c("site","year","plot","treatment","kartez","genus","sp.epithet","family","LifeHistory","PhotoPath","FunctionalGroup","cover","biomass.BM","biomass.BIM","SiteCluster","MetStation","season.precip","GDD","SPEI"), all = TRUE)

shrub_frame <- removal.1 %>%
  subset(FunctionalGroup == "shrub")%>%
  ddply(c("site", "year", "plot","treatment"),
        function(x)data.frame(
          shrub.mass=sum(x$biomass.BIM)
        ))
shrub_frame$year <- shrub_frame$year + 1

#ddply to reduce coremass to plot mass
removal_plot <- ddply(removal.1, c( "site", "year", "plot", "treatment"),
                      function(x)data.frame(
                        biomass = sum(x$biomass.BIM),
                        dominance = max(x$cover)/sum(x$cover),        
                        EQ = community_structure(
                          x,
                          time.var = NULL,
                          abundance.var = "cover",
                          replicate.var = NULL,
                          metric = c("EQ"))$EQ,
                        boer4_mass = sum(x$biomass.BIM[x$kartez %in% "BOER4"]),
                        bogr2_mass = sum(x$biomass.BIM[x$kartez %in% "BOGR2"]),
                        boer4_cover = sum(x$cover[x$kartez %in% "BOER4"]),
                        bogr2_cover = sum(x$cover[x$kartez %in% "BOGR2"]),
                        veg_cover = sum(x$cover),
                        richness = length(x$kartez),
                        SPEI = mean(x$SPEI, na.rm=TRUE)
                      ))%>%
  merge(shrub_frame, by=c("site", "year", "plot","treatment"),all.x=TRUE)%>%
  subset(year != 1995 & year != 2001)%>% #1995 and 2001 data are not IDd to species so species richness and dominance metrics are not attainable
  subset(year!=2014)#%>% #'14 don't have previous year data so shrub biomass cannot be accounted for
#subset(year!=1996) #no SPEI data for 1996

removal_plot$shrub.mass[is.na(removal_plot$shrub.mass)] <- 0
removal_plot$npp <- removal_plot$biomass - removal_plot$shrub.mass
removal_plot <- subset(removal_plot, npp <300) #remove outlier NPP values on upper tail of distribution


#split the two sites
removal_plot$treatment <- factor(removal_plot$treatment, levels = c("control", "subordinate_species", "R_BOER4", "R_BOGR2"))

blue_removal <- subset(removal_plot, site == "removal_blue")
black_removal <- subset(removal_plot, site == "removal_black")


#####
#SR change at both sites

#blue
mod <- lme(richness~treatment*as.factor(year),data=subset(blue_removal, treatment != "control"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
summary(mod)
write.csv(data.frame(pairs(emmeans(mod, ~ treatment | year))), "C:/Users/ohler/Documents/Misc Grad school/DomDiv/blue_richness.csv")

#black
mod <- lme(richness~treatment*as.factor(year),data=subset(black_removal, treatment != "control"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
summary(mod)
write.csv(data.frame(pairs(emmeans(mod, ~ treatment | year))), "C:/Users/ohler/Documents/Misc Grad school/DomDiv/black_richness.csv")



#NPP change at both sites
#blue
mod <- lme(npp~treatment*as.factor(year),data=subset(blue_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
summary(mod)
write.csv(data.frame(pairs(emmeans(mod, ~ treatment | year))), "C:/Users/ohler/Documents/Misc Grad school/DomDiv/blue_npp.csv")


temp <- ddply(blue_removal, .(year, treatment),
              function(x)data.frame(
                npp = mean(x$npp),
                se.npp = sd(x$npp)/sqrt(length(x$npp))
              ))

ggplot(subset(temp, treatment != "subordinate_species" ), aes(year, npp, shape = treatment))+
  geom_pointrange(aes(ymin = npp-se.npp, ymax = npp+se.npp),size = 1)+
  geom_line()+
  ylim(-5,210)+
  ylab("ANPP (g/m2)")+
  xlab("")+
  scale_shape_manual(values = c(17, 2))+ 
  guides(shape = FALSE)+
  theme_base()
#ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/blue_ANPP_year_v4.pdf",
#       device = "pdf",
#       width = 5,
#       height = 5)




#black
mod <- lme(npp~treatment*as.factor(year),data=subset(black_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
summary(mod)
write.csv(data.frame(pairs(emmeans(mod, ~ treatment | year))), "C:/Users/ohler/Documents/Misc Grad school/DomDiv/black_npp.csv")

temp <- ddply(black_removal, .(year, treatment),
              function(x)data.frame(
                npp = mean(x$npp),
                se.npp = sd(x$npp)/sqrt(length(x$npp))
              ))

ggplot(subset(temp, treatment != "subordinate_species" ), aes(year, npp, shape = treatment))+
  geom_pointrange(aes(ymin = npp-se.npp, ymax = npp+se.npp),size = 1)+
  geom_line()+
  ylim(-6,240)+
  ylab("ANPP (g/m2)")+
  xlab("")+
  scale_shape_manual(values = c(19, 21))+ 
  guides(shape = FALSE)+
  theme_base()
#ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/black_ANPP_year_v4.pdf",
#       device = "pdf",
#       width = 5,
#       height = 5)



#SR-NPP by site and treatment
#first a quick experiment without separating ecosystem
removal_plot <- removal_plot%>%
  dplyr::mutate(treatment_revised = case_when(
  treatment == "control" ~ "control",
  treatment == "R_BOGR2" ~ "removal",
  treatment == "R_BOER4" ~ "removal",
))

removal_plot <- removal_plot%>%
                            unite("site.plot",c("site", "plot"), remove = FALSE)

mod <- lme(npp~richness*treatment_revised*site,data=subset(removal_plot, treatment != "subordinate_species"),random=~1|site.plot,correlation=corAR1(form=~year),method="ML")
summary(mod)
anova(mod)
pairs(emtrends(mod, ~ site | richness | treatment_revised, var = "richness"))
pairs(emtrends(mod, ~ richness | treatment_revised, var = "richness"))


#blue
mod <- lme(npp~richness*treatment,data=subset(blue_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
mod.1 <- lme(npp~poly(richness,2)*treatment,data=subset(blue_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
mod.2 <- lme(npp~poly(richness,3)*treatment,data=subset(blue_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
mod.3 <- lme(npp~poly(richness,4)*treatment,data=subset(blue_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
anova(mod, mod.1, mod.2, mod.3)
AICc(mod, mod.1, mod.2, mod.3)

summary(mod)
r.squaredGLMM(mod)


x <- ggpredict(mod, c("richness","treatment"))
plot(x, ci=FALSE)+
  ylim(0, 200)+
  geom_ribbon(aes(ymin=predicted-std.error,ymax=predicted+std.error, linetype = group, fill = group), alpha =.5)+
  scale_fill_manual(values = c("#00ced1","#ffd700"))+
  scale_linetype_manual(values = c(1, 2))+
  scale_color_manual(values = c("#00ced1", "#ffd700"))+
  ggtitle("")+
  xlab("Species richness")+
  ylab("Net primary production (g/m2)")+
  geom_point(data = subset(blue_removal, treatment != "subordinate_species"), aes( x = richness, y = npp, shape = treatment, color = treatment), inherit.aes = FALSE, alpha = .9, size = 2.5)+
  scale_shape_manual(values = c(17, 2))+ 
  guides(color = FALSE, linetype = FALSE, fill = FALSE, shape = FALSE)+
    theme_base()

#ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/blue_ANPP_sr_v4.pdf",
#       device = "pdf",
#       width = 5,
#      height = 5)


#black
mod <- lme(npp~richness*treatment,data=subset(black_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
mod.1 <- lme(npp~poly(richness,2)*treatment,data=subset(black_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
mod.2 <- lme(npp~poly(richness,3)*treatment,data=subset(black_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
mod.3 <- lme(npp~poly(richness,4)*treatment,data=subset(black_removal, treatment != "subordinate_species"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
anova(mod, mod.1, mod.2, mod.3)
AICc(mod, mod.1, mod.2, mod.3)

summary(mod)
r.squaredGLMM(mod)

x <- ggpredict(mod, c("richness","treatment"))
plot(x, ci=FALSE)+
  #ylim(-50,190)+
  ylim(0, 200)+
  geom_ribbon(aes(ymin=predicted-std.error,ymax=predicted+std.error, linetype = group, fill = group), alpha =.5)+
  scale_fill_manual(values = c("#a67c52","#bf0a30"))+
  scale_color_manual(values = c("#a67c52", "#bf0a30"))+
  scale_linetype_manual(values = c(1, 2))+
  ggtitle("")+
  xlab("Species richness")+
  ylab("Net primary production (g/m2)")+
  geom_point(data = subset(black_removal, treatment != "subordinate_species"), aes( x = richness, y = npp, shape = treatment, color = treatment), inherit.aes = FALSE, alpha = .9, size = 2.5)+
  scale_shape_manual(values = c(19, 21))+
  guides(color = FALSE, fill = FALSE, linetype = FALSE, shape = FALSE)+
  theme_base()
#ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/black_ANPP_sr_v4.pdf",
#       device = "pdf",
#       width = 5,
#       height = 5)



####Richness over time
blue_comb <- blue_removal%>%
  ddply(.(site, year, treatment), function(x)data.frame(
    total_cover = mean(x$veg_cover),
    cover_se = sd(x$veg_cover)/sqrt(length(x$veg_cover)),
    total_npp = mean(x$npp),
    npp_se = sd(x$npp)/sqrt(length(x$npp)),
    total_richness = mean(x$richness),
    richness_se = sd(x$richness)/sqrt(length(x$richness))
  ))

black_comb <- black_removal%>%
  ddply(.(site, year, treatment), function(x)data.frame(
    total_cover = mean(x$veg_cover),
    cover_se = sd(x$veg_cover)/sqrt(length(x$veg_cover)),
    total_npp = mean(x$npp),
    npp_se = sd(x$npp)/sqrt(length(x$npp)),
    total_richness = mean(x$richness),
    richness_se = sd(x$richness)/sqrt(length(x$richness))
  ))


##Blue
ggplot(subset(blue_comb, treatment != "control"), aes(year, total_richness, shape = treatment ))+
  geom_pointrange(aes(ymin = total_richness - richness_se, ymax = total_richness + richness_se),size = 1)+
  geom_line()+
  ylab("Subordinate species richness")+
  xlab("")+
  ylim(0,23)+
  ggtitle("")+
  scale_shape_manual(values = c(17, 2))+
  guides( shape = FALSE)+
  theme_base()
#ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/blue_sr_year_v4.pdf",
#       device = "pdf",
#       width = 5,
#       height = 5)

##Black
ggplot(subset(black_comb, treatment != "control"), aes(year, total_richness, shape = treatment, ))+
  geom_pointrange(aes(ymin = total_richness - richness_se, ymax = total_richness + richness_se),size = 1)+
  scale_shape_manual(values = c(19, 21))+
    geom_line()+
  ylab("Subordinate species richness")+
  xlab("")+
  ylim(0,23)+
  ggtitle("")+
  guides(shape = FALSE)+
  theme_base()
#ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/black_sr_year_v4.pdf",
#       device = "pdf",
#       width = 5,
#       height = 5)


# ####Veg cover over time
# ##Blue
# mod <- lme(veg_cover~treatment*as.factor(year),data=subset(blue_removal, treatment != "control"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
# summary(mod)
# emmeans(mod, ~ treatment | year)
# 
# ggplot(subset(blue_comb, treatment != "control"), aes(year, total_cover, shape = treatment, ))+
#   geom_pointrange(aes(ymin = total_cover - cover_se, ymax = total_cover + cover_se),size = 1)+
#   scale_shape_manual(values = c(19, 21))+
#   geom_line()+
#   ylab("Subordinate species cover")+
#   xlab("")+
#   ylim(0,64)+
#   ggtitle("")+
#   guides(shape = FALSE)+
#   theme_base()
# #ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/blue_cover_year_v3.pdf",
#   #     device = "pdf",
#   #     width = 5,
#  #      height = 5)
# 
# 
# 
# ##Black
# mod <- lme(veg_cover~treatment*as.factor(year),data=subset(black_removal, treatment != "control"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
# summary(mod)
# 
# ggplot(subset(black_comb, treatment != "control"), aes(year, total_cover, shape = treatment ))+
#   geom_pointrange(aes(ymin = total_cover - cover_se, ymax = total_cover + cover_se),size = 1)+
#   scale_shape_manual(values = c(19, 21))+
#   geom_line()+
#   ylab("Subordinate species cover")+
#   xlab("")+
#   ylim(0,64)+
#   ggtitle("")+
#   guides(shape = FALSE)+
#   theme_base()
# #ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/black_cover_year_v3.pdf",
# #       device = "pdf",
# #       width = 5,
#  #      height = 5)

####################### 
#Which functional groups win when dominants are removed

removal_full <- read.csv("C:/Users/ohler/Desktop/Sevilleta_RemovalPlots_biomass_19Nov2022.csv")
removal <- removal_full%>%
  subset( treatment != "R_all" & kartez != "litter" & kartez != "soil" & kartez != "total")%>%
  subset( site == "removal_blue" | site == "removal_black")


removal$FG_revised <- ifelse(removal$kartez == "ASTER", "forb",
                             ifelse(removal$kartez == "SPORO", "grass",
                                    
                                    
                                    
                                    as.character(removal$FunctionalGroup)))


removal$FG_revised <- ifelse(removal$kartez == "ASTER", "forb",
                             ifelse(removal$kartez == "SPORO", "grass",  as.character(removal$FunctionalGroup)))

#x <- subset(removal, site == "removal_black" & treatment == "R_BOER4" & year == "2021" & kartez == "BOER4") #this was just to get some numbers about treatment efficacy

just_functional_group <- removal%>%
  ddply(.(site,year,plot,treatment,FG_revised),
        function(x)data.frame(
          cover = sum(x$cover),
          richness = length(x$kartez)
        ))

functional_group <-  ddply(just_functional_group, .(site,year,plot,treatment),
                           function(x)data.frame(
                             plot_cover = sum(x$cover)
                           ))%>%
  merge(just_functional_group, by = c("site","year","plot","treatment"))%>%
  subset(plot_cover < 100)

functional_group$relative_cover <- functional_group$cover / functional_group$plot_cover






##BLACK GRASS
mod <- lme(cover~treatment*as.factor(year),data=subset(functional_group, site == "removal_black" & FG_revised == "grass"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
summary(mod)
write.csv(data.frame(pairs(emmeans(mod, ~ treatment | year))), "C:/Users/ohler/Documents/Misc Grad school/DomDiv/black_grass.csv")


##BLUE GRASS
mod <- lme(cover~treatment*as.factor(year),data=subset(functional_group, site == "removal_blue" & FG_revised == "grass"),random=~1|plot,correlation=corAR1(form=~year),method="ML")
summary(mod)
write.csv(data.frame(pairs(emmeans(mod, ~ treatment | year))), "C:/Users/ohler/Documents/Misc Grad school/DomDiv/blue_grass.csv")

##BLUE AND BLACK FORB #there's enough data but nearly no differences so ignoring for now
##BLUA AND BLACK SHRUB #there's not enough data to run these models for shrubs. In the last 10 years there's nearly no shrubs in any of the plots


#################################
######dataset information to be included in methods and results 
blue_removal %>%
  subset(treatment == "control")%>%
  ddply(c("year"),
        function(x)data.frame(
          richness = mean(x$richness)
        ))


black_removal %>%
  subset(treatment == "control")%>%
  ddply(c("year"),
        function(x)data.frame(
          richness = mean(x$richness)
        ))



############################
#########


functional_group_avg <- functional_group%>%
  ddply(.(site, year, treatment, FG_revised),function(x)data.frame(
    cover = mean(x$cover),
    se.cover = sd(x$cover)/sqrt(length(x$cover)),
    relative_cover = mean(x$relative_cover),
    se.relative_cover = sd(x$relatve_cover)/sqrt(length(x$relative_cover))
  ))


functional_group_avg$site <- factor(functional_group_avg$site, levels = c("removal_blue", "removal_black"))

ggplot(subset(functional_group_avg, FG_revised == "grass"), aes(year, cover, shape = paste(treatment, site)))+
  facet_wrap(~site)+
  scale_shape_manual(values = c(1, 2, 16, 17))+
  geom_pointrange(aes(ymin = cover-se.cover, ymax = cover+se.cover),size = 1)+
  geom_line()+
  xlab("")+
  ylab("Percent cover")+
  guides(shape = FALSE)+
  theme_base()
#ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/grass_cover_year_v4.pdf",
#       device = "pdf",
#      width = 10,
#       height = 5)


###########################################
###############################Which species are increasing

species <- removal%>%
  subset(year != "1995" & year != "2001")%>% #data not to species these years
  ddply(.(site, year, kartez, treatment),function(x)data.frame(
    cover = mean(x$cover)
  ))%>%
  subset(treatment == "R_BOGR2" | treatment == "R_BOER4")%>%
  subset(cover>10)



###################Fun with NMDS###########

##blue PERMANOVA
removal.temp <- removal[,c("site", "year","plot", "treatment", "kartez", "cover")]

removal.temp <- subset(removal.temp, cover >= 1 & kartez != "annual.forb" & kartez != "shrub" & kartez != "perennial.grass" & kartez != "perennial.forb" & site == "removal_blue"
)%>%
  subset(year != "1995" & year != "2001") #data not to species these years #only single instance of MUTO2 and t unnecessarily dominates the NMDS

removal_spread <- spread(removal.temp, "kartez", "cover", fill = "0")
removal_spread$treatment <- plyr::revalue(removal_spread$treatment, c(R_BOER4 = "removal", R_BOGR2 = "removal"))

#This creates a matrix with only the species columns which is required to run ordination analyses
removal_matrix_1 <- removal_spread[,5:64]
removal_matrix <- sapply( removal_matrix_1, as.numeric )

perm <- how(nperm = 999)
setBlocks(perm) <- with(removal_spread, removal_spread$plot)
removal_permanova_blue <- pairwise.adonis2(removal_matrix~treatment*year, data=removal_spread, method="bray", permutations = perm)
removal_permanova_blue #this could be a table



##black PERMANOVA
removal.temp <- removal[,c("site", "year","plot", "treatment", "kartez", "cover")]

removal.temp <- subset(removal.temp, cover >= 1 & kartez != "annual.forb" & kartez != "shrub" & kartez != "perennial.grass" & kartez != "perennial.forb" & site == "removal_black"
)%>%
  subset(year != "1995" & year != "2001") #data not to species these years #only single instance of MUTO2 and t unnecessarily dominates the NMDS

removal_spread <- spread(removal.temp, "kartez", "cover", fill = "0")
removal_spread$treatment <- plyr::revalue(removal_spread$treatment, c(R_BOER4 = "removal", R_BOGR2 = "removal"))

#This creates a matrix with only the species columns which is required to run ordination analyses
removal_matrix_1 <- removal_spread[,5:54]
removal_matrix <- sapply( removal_matrix_1, as.numeric )

perm <- how(nperm = 999)
setBlocks(perm) <- with(removal_spread, removal_spread$plot)
removal_permanova_black <- pairwise.adonis2(removal_matrix~treatment*year, data=removal_spread, method="bray", permutations = perm)
removal_permanova_black #this could be a table



##############################################
#############try to calculate beta diversity within a site/treatment for each year

#blue beta
removal.reduced <- removal[,c("site", "year","plot", "treatment", "kartez", "cover")]
removal.reduced <- subset(removal.reduced, cover >= 1 & site == "removal_blue")

removal_spread <- spread(removal.reduced, "kartez", "cover", fill = "0")

test_master <- {}

for (i in unique(removal_spread$year)) {
  
  temp.df <- subset(removal_spread, year == i)
  removal_matrix_temp.1 <- temp.df[,5:69]
  removal_matrix_temp <- sapply( removal_matrix_temp.1, as.numeric )
  
  dis <- vegdist(removal_matrix_temp)
  
  groups <- interaction(subset(removal_spread, year == i)$treatment, subset(removal_spread, year == i)$site)
  ## Calculate multivariate dispersions
  mod <- betadisper(dis, groups, type = "centroid")
  
  temp.test <- data.frame(mod$group.distances, i)
  temp.test$treatment <- row.names(temp.test)
  
  
  test_master <- rbind(test_master, temp.test)
  
  rm(temp.df,temp.test)
  
}

beta_div_blue <- test_master%>%
  tidyr::separate( col = treatment, into = c("treatment", "ecosystem"), sep = "[.]")%>%
mutate(treatment = plyr::revalue(treatment, c(R_BOGR2 = "removal", R_BOER4 = "removal")))



##black beta
removal.reduced <- removal[,c("site", "year","plot", "treatment", "kartez", "cover")]
removal.reduced <- subset(removal.reduced, cover >= 1 & site == "removal_black")

removal_spread <- spread(removal.reduced, "kartez", "cover", fill = "0")

test_master <- {}

for (i in unique(removal_spread$year)) {
  
  temp.df <- subset(removal_spread, year == i)
  removal_matrix_temp.1 <- temp.df[,5:58]
  removal_matrix_temp <- sapply( removal_matrix_temp.1, as.numeric )
  
  dis <- vegdist(removal_matrix_temp)
  
  groups <- interaction(subset(removal_spread, year == i)$treatment, subset(removal_spread, year == i)$site)
  ## Calculate multivariate dispersions
  mod <- betadisper(dis, groups, type = "centroid")
  
  temp.test <- data.frame(mod$group.distances, i)
  temp.test$treatment <- row.names(temp.test)
  
  
  test_master <- rbind(test_master, temp.test)
  
  rm(temp.df,temp.test)
  
}

beta_div_black <- test_master%>%
  tidyr::separate( col = treatment, into = c("treatment", "ecosystem"), sep = "[.]")%>%
  mutate(treatment = plyr::revalue(treatment, c(R_BOGR2 = "removal", R_BOER4 = "removal")))



#combine black and blue beta
beta_div_summary <- rbind(beta_div_blue, beta_div_black)%>%
                      ddply( .(treatment, ecosystem),
                          function(x)data.frame(
                            mean.dist = mean(x$mod.group.distances),
                            se.dist = sd(x$mod.group.distances)/sqrt(length(x$mod.group.distances))
                          ))
ggplot(beta_div_summary, aes(treatment, mean.dist, shape = paste(treatment, ecosystem)))+
  facet_wrap(~ecosystem)+
  geom_pointrange(aes(ymin = mean.dist - se.dist, ymax = mean.dist+se.dist), size = 0.75)+
  geom_point(data = rbind(beta_div_blue, beta_div_black), aes(y=mod.group.distances ), alpha = 0.5)+
  ylim(0,0.6)+
  scale_shape_manual(values = c(16, 17,1, 2))+
  ylab("Beta Diversity")+
  xlab("")+
  theme_base()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("C:/Users/ohler/Documents/Misc Grad School/DomDiv/Figures/betadiv_v5.pdf",
       device = "pdf",
       width = 4,
       height = 4)

#beta diversity hypothesis test
mod.blue <- gls(mod.group.distances~treatment,correlation = corAR1(form=~i|treatment),  data = beta_div_blue)
summary(mod.blue)


mod.black <- gls(mod.group.distances~treatment,correlation = corAR1(form=~i|treatment),  data = beta_div_black)
summary(mod.black)










                     
##
###CSF trial

mod <- lme(npp~SPEI, random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "control" & SPEI != "NaN"))
mod1 <- lme(npp~poly(SPEI,2), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "control" & SPEI != "NaN"))
mod2 <- lme(npp~poly(SPEI,3), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "control" & SPEI != "NaN"))
mod3 <- lme(npp~poly(SPEI,4), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "control" & SPEI != "NaN"))
AICc(mod, mod1, mod2, mod3)
summary(mod3)

mod <- lme(npp~SPEI, random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "R_BOGR2" & SPEI != "NaN"))
mod1 <- lme(npp~poly(SPEI,2), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "R_BOGR2" & SPEI != "NaN"))
mod2 <- lme(npp~poly(SPEI,3), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "R_BOGR2" & SPEI != "NaN"))
mod3 <- lme(npp~poly(SPEI,4), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_blue" & treatment == "R_BOGR2" & SPEI != "NaN"))
AICc(mod, mod1, mod2, mod3)
summary(mod3)


mod <- lme(npp~SPEI, random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "control" & SPEI != "NaN"))
mod1 <- lme(npp~poly(SPEI,2), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "control" & SPEI != "NaN"))
mod2 <- lme(npp~poly(SPEI,3), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "control" & SPEI != "NaN"))
mod3 <- lme(npp~poly(SPEI,4), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "control" & SPEI != "NaN"))
AICc(mod, mod1, mod2, mod3)
summary(mod3)

mod <- lme(npp~SPEI, random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "R_BOER4" & SPEI != "NaN"))
mod1 <- lme(npp~poly(SPEI,2), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "R_BOER4" & SPEI != "NaN"))
mod2 <- lme(npp~poly(SPEI,3), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "R_BOER4" & SPEI != "NaN"))
mod3 <- lme(npp~poly(SPEI,4), random = ~1|year, corr = corAR1(form = ~ 1|year),data = subset(removal_plot, site == "removal_black" & treatment == "R_BOER4" & SPEI != "NaN"))
AICc(mod, mod1, mod2, mod3)
summary(mod3)

###INDICATOR SPECIES ANALYSIS TRIAL

library(indicspecies)

#removal.1

wide <- removal.1%>%
        dplyr::select(site, year, plot, treatment, kartez, cover)%>%
        pivot_wider( names_from = "kartez", values_from = "cover", values_fill = 0)

  x <- subset(wide, site == "removal_blue" & treatment != "subordinate_species")[,5:50]

  indval <- multipatt(x
          , subset(wide, site == "removal_blue" & treatment != "subordinate_species")$treatment, 
                      control = how(nperm=999)) 
summary(indval)





x <- subset(wide, site == "removal_black" & treatment != "subordinate_species")[,5:50]

indval <- multipatt(x
                    , subset(wide, site == "removal_black" & treatment != "subordinate_species")$treatment, 
                    control = how(nperm=999)) 
summary(indval)








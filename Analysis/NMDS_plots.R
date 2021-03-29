#NMDS plots
library(ggplot2)
library(vegan)
scores <- read.csv("Processing/1b_out_ChemSummaries_Indiv.csv", header = T)

#################################
#Plots

#colored by latitude, shaped by age
plot = ggplot(scores,aes(x=NMDS1, y=NMDS2, colour = lat, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_gradient2(low="gold",mid = "red", high="blue", midpoint = 26) + scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("nmds1_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
# setEPS()
# postscript("20180517_nmds1_all_legend.eps")
# vv + labs(shape="Leaf age", color = "Latitude")
# dev.off()
plot = ggplot(scores,aes(x=NMDS3, y=NMDS4, colour = lat, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_gradient2(low="gold",mid = "red", high="blue", midpoint = 26) + scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("nmds3_all.eps")
# vv + theme(legend.position = "none")
# dev.off()

#re-order levels for key, also makes listing colors easier
scores$region <- factor(scores$region, levels =c("north temperate", "temperate","subtropical", "tropical") )
levels(scores$region)

#colored by region
plot = ggplot(scores,aes(x=NMDS1, y=NMDS2, colour = region, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))+ scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv
# setEPS()
# postscript("nmds1_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
# setEPS()
# postscript("nmds1_all_legend.eps")
# vv
# dev.off()
#axis 1 is mostly separating by latitude (but also age for tropical); for mature leaves, it's tropical/subtropical vs everything north of FL; for young leaves, it's similar but tropical and subtropical are separating out more. axis 2 separates by leaf age for PHAM but not PHRI. THIS IS SO COOL because it shows that more divergence is happening in young leaves! And it's crazy that FL is grouping with CR!
plot = ggplot(scores,aes(x=NMDS3, y=NMDS4, colour = region, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))+ scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv
# setEPS()
# postscript("20180531_nmds2_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
#for the most part leaf ages still cluster geographically, but there is more geographic overlap in general. Axis 3 is separating young from mature, to differing degrees depending on region. Axis 4 separates PHAM from PHRI (a bit)

####################
#will draw an ellipse for each leaf age in each region
scores$reg_age <- as.factor(paste(scores$region,scores$age,sep="_"))

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
#ordination ellipse object here allows you to change settings
ord<-ordiellipse(ord = scores[12:13], groups = scores$reg_age, display = "sites", 
                 kind = "se", conf = 0.95, label = T)



df_ellipse <- data.frame()
for(g in levels(scores$reg_age)){
  df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(scores[scores$reg_age==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,reg_age=g))
}

df_ellipse <- data.frame()
for(g in levels(scores$reg_age)){
  df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(scores[scores$reg_age==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,reg_age=g))
}

df_ellipse <- cbind(df_ellipse, read.table(text = df_ellipse$reg_age,sep="_",colClasses = "factor",col.names=c("region","age")))

plot = ggplot(scores,aes(x=NMDS1, y=NMDS2, colour = region, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))+ scale_shape_manual(values = c(19,1)) +  geom_path(data=df_ellipse, aes(x=NMDS1, y=NMDS2, linetype=age))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv

ggplot(scores,aes(x=NMDS1,y=NMDS2)) + geom_path(data=df_ellipse,aes(x=NMDS1,y=NMDS2))

#THIS WORKS WITHOUT AGE
ggplot(scores,aes(x=NMDS1, y=NMDS2, colour = region))+  geom_point(size=4, stroke = 2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2")) +  geom_path(data=df_ellipse, aes(x=NMDS1, y=NMDS2)) + geom_path(data=df_ellipse,aes(x=NMDS1,y=NMDS2)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))

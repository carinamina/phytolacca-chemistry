#NMDS plots

scores <- read.csv("Processing/1b_out_ChemSummaries.csv", row.names = T)

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
scores$region <- factor(scores$region, levels =c("temperate", "n/a","subtropical", "tropical") )
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
#for the most part leaf ages still cluster geographically, but there is more geographic overlap in general. Axis 3 is separating young from mature, to differing degrees depending on region. Axis 4 separates PHAM from PHRI (more or less)

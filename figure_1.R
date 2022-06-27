library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(gridGraphics)

source(file = "./import_data.R")

map<-getData('GADM', country='MAYOTTE', level=1)

map2 = fortify(map)
map2$type = "Outer"
map2$type[map2$id %in% c(170797,170813,170814,170820,170828)] = "Central"

p3 = ggplot()
p3 = p3 + geom_polygon(data=map2, aes(long,lat,fill=type,group=group), col="black")
p3 = p3 + scale_fill_manual(name=NULL, values=c("Central"="#1f78b4", "Outer"="#a6cee3"))
p3 = p3 + theme(aspect.ratio=1.3,
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.background = element_rect(fill="white"),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                plot.background = element_blank(),
                legend.position="none",
                axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.title = element_blank())

dat_plot = dat1
dat_plot$annot = ifelse(dat_plot$type=="cas", "A", "B")
dat_plot$type[dat_plot$type=="cas"] = "Number of reported cases"
dat_plot$type[dat_plot$type=="sero"] = "Number of samples collected"
p1 = ggplot(data = dat_plot, aes(x=week, fill=aggr_geo3))
p1 = p1 + geom_bar(col="black", width=0.8)
p1 = p1 + facet_wrap(~type,ncol=1,scales = "free_y", switch = "y")
p1 = p1 + geom_text(aes(label=annot, x=1.5, y=-Inf), vjust=-1, size=7)
p1 = p1 + scale_fill_manual(values=c("#1f78b4", "#a6cee3"), labels=c("Central communes", "Outer communes", "Missing geographical\ninformation"), name=NULL, na.value="white")
p1 = p1 + xlab("Weeks") + ylab(NULL)
p1 = p1 + scale_x_discrete(limits = c(paste0("2018-", 44:52), paste0("2019-0", 1:9), paste0("2019-", 10:32)))
p1 = p1 + theme_bw()
p1 = p1 + theme(strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(size = 13),
                legend.title=element_blank(),
                legend.text=element_text(size=11),
                legend.position=c(0.15, 0.9),
                legend.background = element_rect(fill = "white", color = "black"),
                legend.direction="vertical",
                axis.text.x = element_text(angle = 90, size=c(10,0)),
                axis.title.x = element_text(size=13))

pdf(file="./res/fig_1.pdf", width=8, height=7)
p1
print(p3, vp=viewport(width = 0.3, height = 0.3, x = 0.8, y = 0.83))
dev.off()

# Supplementary Table S1:

tab_suppl_1 <- as.data.frame(dat_plot
                           %>% filter(type == "Number of reported cases", age >= 15)
                           %>% group_by(week, aggr_geo3)
                           %>% summarise(n = n_distinct(id))
                           %>% dcast(week ~ aggr_geo3, value.var="n"))
tab_suppl_1[is.na(tab_suppl_1)] <- 0
tab_suppl_1 <- tab_suppl_1[,1:3]
tab_suppl_1

tab_suppl_2 <- as.data.frame(dat_plot
                           %>% filter(type == "Number of samples collected", age >= 15)
                           %>% group_by(week, aggr_geo3)
                           %>% summarise(n = n_distinct(id))
                           %>% dcast(week ~ aggr_geo3, value.var="n"))
tab_suppl_2[is.na(tab_suppl_2)] <- 0
tab_suppl_2 <- tab_suppl_2[,1:3]
tab_suppl_2

tab_suppl_3 <- as.data.frame(dat_plot
                             %>% filter(type == "Number of samples collected", inf == "anc", age >= 15)
                             %>% group_by(week, aggr_geo3)
                             %>% summarise(n = n_distinct(id))
                             %>% dcast(week ~ aggr_geo3, value.var="n"))
tab_suppl_3[is.na(tab_suppl_3)] <- 0
tab_suppl_3 <- tab_suppl_3[,1:3]
tab_suppl_3

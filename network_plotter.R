library(RTN)
library(RedeR)
library(TreeAndLeaf)
library(data.table)
library(RColorBrewer)
library(classInt)
library(gdata)
library(stringr)
library(igraph)

load("../rdatas/shh_rtn.RData")

interesse_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")

tfs <- rtni@regulatoryElements
tree <- tni.graph(rtni, tnet = "dpi", gtype = "rmap", 
                  regulatoryElements = interesse_regs)
tna <- tna.get(rtna, what = "mra")

#---- Defining network MRs ----

index <- which(tfs %in% tna$Regulon)
temp <- data.frame(gene = tfs[-index])
temp$is_mr <- "no"

temp <- rbind(temp, data.frame(gene = tna$Regulon, is_mr = "yes"))

#---- Initial plot ----

filepath <- system.file(package = "RedeR", "java/reder_v2.jar")
cmd <- "java -jar"
command <- paste(cmd, shQuote(filepath), "openshellDcall", 9091, sep = ' ')
system(command, wait = FALSE)

gg <- hclust2igraph(tree$hcl)

gg <- formatTree(gg, theme = 1)


gg <- att.mapv(gg, dat = temp, refcol = 0)
gg <- att.setv(gg, from = "name", to = "nodeAlias")

V(gg)$nodeAlias[str_detect(V(gg)$nodeAlias, "^N[0-9]+$")] <- ""

rdp <- RedPort()

calld(rdp)

treeAndLeaf(rdp, gg)

addLegend.color(rdp, colvec = c("#58139B", "#D80916"), title = "",
                labvec = c("Master Regulator", "Regulon of Interest"), 
                dxborder = 30, dyborder = 30, vertical = T,
                position = "bottomright")

#---- Relax it ----

for(i in 1:100){
relax(rdp, p1 = 10, p2 = 100, p3 = 280, p4 = 280, 
      p5 = 50, p6 = 10, p7 = 180, p8 = 200)
}

#Let this run for a couple times. It works wonders.

for(x in 1:5){
  
  for(i in 1:1000){relax(rdp, p1 = 10, p2 = 100, p3 = 180, p4 = 180, 
                         p5 = 50, p6 = 10, p7 = 150, p8 = 200)}
  
  print("sleep!")
  if(x < 30){
    Sys.sleep(12) # Initial iterations need more time in between
  }else{
    Sys.sleep(4)
  }
  print("wake up!")
}

for(x in 1:10){
  for(i in 1:100){relax(rdp, p1 = 2, p2 = 100, p3 = 30, p4 = 30, p5 = 10, 
                      p6 = 10, p7 = 30, p8 = 10)}
  print("sleep!")
  Sys.sleep(1)
  print("wake up!")
}
  
#---- Selection of MR regulons, so it is possible to change its color ----

interesse_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")

selectNodes(rdp, tna$Regulon)
selectNodes(rdp, interesse_regs)



#---- Size legend is necessary to read the amount of elements in regulon ----


# scl <- round(as.numeric(gg$legNodeSize$scale))
# 
# leg <- as.character(round(as.numeric(gg$legNodeSize$legend)))
# 
# addLegend.size(rdp, sizevec = scl, labvec = leg, 
#                title = "Node Size (Regulon size)", position = "bottomright", 
#                intersp = 10, ftsize = 10, vertical = F)  



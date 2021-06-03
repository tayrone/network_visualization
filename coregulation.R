library(RTN)
library(RedeR)
library(RColorBrewer)
library(classInt)
library(gdata)
library(stringr)
library(igraph)
library(purrr)

load("../rdatas/g3_rtn.RData")

interesse_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")

tree <- tni.graph(rtni, tnet = "dpi", gtype = "amap", 
                  regulatoryElements = interesse_regs)


filepath <- system.file(package = "RedeR", "java/reder_v2.jar")
cmd <- "java -jar"
command <- paste(cmd, shQuote(filepath), "openshellDcall", 9091, sep = ' ')
system(command, wait = FALSE)

rdp <- RedPort()

calld(rdp)

addGraph(rdp, tree)


#---- Generate network for each pair of regulons, and then plot each one ----

highly_connected <- interesse_regs

for(regulon1 in highly_connected){
  for(regulon2 in highly_connected){
    tree <- tni.graph(rtni, tnet = "dpi", gtype = "rmap", 
                      regulatoryElements = c(regulon1, regulon2))
    
    assign(paste0(regulon1, "_", regulon2), tree)
    
  }
}

addGraph(rdp, ZNF365_ZNF385B)



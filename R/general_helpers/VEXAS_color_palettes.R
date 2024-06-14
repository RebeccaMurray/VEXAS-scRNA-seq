library(ggsci)
library(pals)

############################## Cell Types ################################

cell.type.palette <- c(
  HSC = "#2061AD", # pal
  LMPP = "#2E5050", # pal
  EMP = "#8C1A1C",
  BaEoMa = "#F7F49C",
  `Early Eryth` = "red3", # pal
  `Late Eryth` = "red2", # pal # #CD212A
  MkP = "#F5B936", # pal
  `CD4 T` = "#842D8D",
  `CD8 T` = "darkorchid4",
  `NK` = "darkorchid3",
  CLP = "#91357d",
  B = "#D53D96",
  Plasma = "#EB977A",
  GMP = "#2A9098",
  `CD14 Mono` = "#208D43",
  cDC = "palegreen4",
  pDC = "#8A772D",
  `pre-mDC` = "darkgoldenrod3"
)

# 
# # Based off scales::show_col(pals::stepped2(n=20))
# cell.type.palette <- c(
#   HSC = "#393B79", # pal
#   LMPP = "#5254A3", # pal
#   EMP = "#6B6ECF",
#   BaEoMa = "#E7CB94",
#   `Late Eryth` = "brown4", # pal
#   `Early Eryth` = "#AD494A", # pal
#   MkP = "#D6616B", # pal
#   `CD4 T` = "#6e7939",
#   `CD8 T` = "#8CA252",
#   `NK` = "#B5CF6B",
#   CLP = "#91357d",
#   B = "#b14d8e",
#   Plasma = "#ca699d",
#   GMP = "#8C6D31",
#   `CD14 Mono` = "#813753",
#   cDC = "#a65461",
#   pDC = "gold1",
#   `pre-mDC` = "darkgoldenrod3"
# )
# 
# # Based off scales::show_col(pals::stepped2(n=20))
# cell.type.palette <- c(
#   HSC = "#393B79", # pal
#   LMPP = "#5254A3", # pal
#   EMP = "#6B6ECF",
#   BaEoMa = "#bd925a",
#   `Late Eryth` = "brown4", # pal
#   `Early Eryth` = "#AD494A", # pal
#   MkP = "#D6616B", # pal
#   `CD4 T` = "#6e7939",
#   `CD8 T` = "#8CA252",
#   `NK` = "#B5CF6B",
#   CLP = "#b14d8e",
#   B = "#dd88ac",
#   Plasma = "#ca699d",
#   GMP = "darkgoldenrod4",
#   `CD14 Mono` = "darkgoldenrod",
#   cDC = "darkgoldenrod3",
#   pDC = "goldenrod2",
#   `pre-mDC` = "darkgoldenrod1"
# )
# 
# cell.type.palette <- c(
#   HSC = "#393B79", # pal
#   LMPP = "#5254A3", # pal
#   EMP = "#6B6ECF",
#   BaEoMa = "#bd925a",
#   `Late Eryth` = "brown4", # pal
#   `Early Eryth` = "#AD494A", # pal
#   MkP = "#D6616B", # pal
#   `CD4 T` = "#6e7939",
#   `CD8 T` = "#8CA252",
#   `NK` = "#B5CF6B",
#   CLP = "#b14d8e",
#   B = "#dd88ac",
#   Plasma = "#ca699d",
#   GMP = "darkgoldenrod4",
#   `CD14 Mono` = "#A16928",
#   cDC = "#bd925a",
#   pDC = "#d6bd8d",
#   `pre-mDC` = "#e0c2a2"
# )
# 
# # magentas <- c("#f3cbd3", "#eaa9bd","#dd88ac","#ca699d", "#b14d8e","#91357d", "#6c2167")
# # browns <- c("#ede5cf","#e0c2a2",'#d39c83',"#c1766f","#a65461","#813753","#541f3f")
# # "#A16928","#bd925a",#d6bd8d,#edeac2,#b5c8b8,#79a7ac,#2887a1
# # 
# # "#A16928","#bd925a","#d6bd8d","#edeac2",#b5c8b8,#79a7ac,#2887a1
# 
# 
# cell.type.palette <- c(
#   HSPC = "dodgerblue4", # pal,
#   HSC = "dodgerblue4", # pal
#   GMP = "slategray3", # pal
#   LMPP = "dodgerblue3", # pal
#   `Early Eryth` = "skyblue3", # pal
#   EMP = "deepskyblue3", # pal
#   `Prog Mk` = "cornflowerblue", # pal
#   Other = "grey80",
#   BaEoMa = "lightsalmon3", # pal
#   `Late Eryth` = "dodgerblue", # pal
#   `B cells` = "violet", # pal
#   `CD8 T cells` = "#A55194", # pal
#   `NK` = "#CE6DBD", # pal
#   `CD4 T cells` = "#DE9ED6", # pal,
#   `CD4/CD8/NK cells` = "#DE9ED6", # pal
#   Plasma = "thistle3", # pal
#   pDC = "sienna2", # pal
#   `CD14 Mono` = "darkorange4", # pal
#   cDC2 = "sienna1", # pal
#   Macrophage = "sienna3", # pal
#   `Stromal` = "tan3", # pal,
#   `pre-pDC` = "sienna4"
# )
# 
# cell.type.palette <- c(
#   HSPC = "skyblue4", # pal,
#   HSC = "skyblue4", # pal
#   GMP = "sienna3", # pal
#   LMPP = "dodgerblue3", # pal
#   `Early Eryth` = "skyblue3", # pal
#   EMP = "deepskyblue3", # pal
#   `Prog Mk` = "cornflowerblue", # pal
#   Other = "grey80",
#   BaEoMa = "lightsalmon", # pal
#   `Late Eryth` = "dodgerblue", # pal
#   `B cells` = "violet", # pal
#   `CD8 T cells` = "#A55194", # pal
#   `NK` = "#CE6DBD", # pal
#   `CD4 T cells` = "#DE9ED6", # pal,
#   `CD4/CD8/NK cells` = "#DE9ED6", # pal
#   Plasma = "thistle3", # pal
#   pDC = "lightsalmon1", # pal
#   `CD14 Mono` = "lightsalmon2", # pal
#   cDC2 = "lightsalmon3", # pal
#   Macrophage = "sienna1", # pal
#   `Stromal` = "tan3", # pal,
#   `pre-pDC` = "sienna2"
# )
# 
# cell.type.palette <- c(
#   HSPC = "skyblue4", # pal,
#   HSC = "skyblue4", # pal
#   GMP = "darkseagreen", # pal
#   LMPP = "dodgerblue3", # pal
#   `Early Eryth` = "skyblue3", # pal
#   EMP = "deepskyblue3", # pal
#   `Prog Mk` = "cornflowerblue", # pal
#   Other = "grey80",
#   BaEoMa = "seagreen3", # pal
#   `Late Eryth` = "dodgerblue", # pal
#   `B cells` = "violet", # pal
#   `CD8 T cells` = "#A55194", # pal
#   `CD4/CD8 T cells` = "#A55194", # pal
#   `NK` = "#CE6DBD", # pal
#   `NK cells` = "#CE6DBD", # pal
#   `CD4 T cells` = "#DE9ED6", # pal,
#   `CD4/CD8/NK cells` = "#DE9ED6", # pal
#   Plasma = "thistle3", # pal
#   pDC = "darkseagreen2", # pal
#   `CD14 Mono` = "darkseagreen4", # pal
#   `CD16 Mono` = "palegreen3", # pal
#   cDC2 = "darkseagreen3", # pal
#   cDC = "darkseagreen3", # pal
#   Macrophage = "palegreen3", # pal
#   `Stromal` = "palegreen2", # pal,
#   `pre-pDC` = "seagreen2"
# )
# 
# cell.type.palette <- c(
#   HSC = "#2260AC", # pal
#   GMP = "#2D8E96", # pal
#   `GMP-Neu` = "#2D8E96", # pal
#   `GMP-Mono` = "#2D8E96", # pal
#   LMPP = "darkslategray", # pal
#   `Early Eryth` = "firebrick3", # pal
#   `Late Eryth` = "firebrick2",
#   EMP = "firebrick4", # pal
#   `MkP` = "indianred3", # pal
#   Platelet = "indianred2",
#   BaEoMa = "sienna2", # pal
#   `B` = "palevioletred", # pal
#   Plasma = "palevioletred4", # pal
#   `CLP` = "palevioletred1",
#   `cDC` = "palevioletred2",
#   pDC = "palevioletred3", # pal
#   `CD4 T` = "salmon1", # pal
#   `CD8 T` = "salmon2", # pal
#   `NK` = "salmon3", # pal
#   `CD14 Mono` = "rosybrown2"
#   )
# #
# # cell.type.palette <- scales::hue_pal()(17)
# # 
# cell.type.palette <- c(HSC = "#b76352",
#                                LMPP = "#61b94e",
#                                EMP = "#a859cc",
#                                `Early Eryth` = "#b4b331",
#                                `Late Eryth` = "#6567d7",
#                                MkP = "#d78b38",
#                                BaEoMa = "#6d95de",
#                                CLP = "#d14c34",
#                                GMP = "#46aecc",
#                                `CD14 Mono` = "#d14169",
#                                `B` = "#5fc397",
#                                `Plasma` = "#cf4ba4",
#                                `pre-mDC` = "#618534",
#                                pDC = "#c68cd7",
#                                `CD4 T` = "#378457",
#                                `CD8 T` = "#974c81",
#                                NK = "#c0ad66",
#                                "#5a63aa")
#                                # "#866d2c",
#                                # "#e184a1")
# 
# 
# # #
# cell.type.palette <- c("#855C75",
# "#D9AF6B",
# "#AF6458",
# "#736F4C",
# "#526A83",
# "#625377",
# "#68855C",
# "#9C9C5E",
# "#427551",
# "#A06177",
# "#8C785D",
# "#467378",
# "#7C7C7C",
# 
# # "#88CCEE",
# "#AA4499",
# "#2b314e",
# "#CC6677",
# "#DDCC77",
# "#117733",
# "#332288",
# "#44AA99",
# "#999933",
# "#882255",
# "#661100",
# "#6699CC"
# # "#888888"
# )
# 
# cell.type.palette <- c(
#   "#427551",
#   "#512956",
#   "#345620",
#   "#5c5e96",
#   "#756429",
#   "#38688b",
#   "#6b222d",
#   "#3e6c6c",
#   "#91506f",
#   "#253a1e",
#   "#905646",
#   "#533518",
#   "#1c3d3e",
#   "#66664b",
#   "#4b302e",
#   "#645c6c"
# )


scales::show_col(cell.type.palette)

############################## HSC subcluster colors ########################

cluster.colors <- pal_npg()(7)
names(cluster.colors) <- seq(0, 6)

############################## VEXAS vs control colors ######################
# vexas.color <- "#3C5488FF"
# control.color <- "#F39B7FFF"

vexas.color = "#333c92"
control.color <- "#add9e8"

sample.type.palette <- c(VEXAS = vexas.color, Control = control.color)



############################## Genotyping colors #############################

# ## GoT paper colors
# wt.color <- "#333c92"
# mut.color <- "#add9e8"
# na.color <- "#bfbdbd"

# # NPG colors
# mut.color <- npg.pal[5]
# wt.color <- npg.pal[4]
# 
# # viridis colors
# mut.color <- viridis.pal[1]
# wt.color <- viridis.pal[6]
# 
# simple colors
# mut.color <- "red"
# wt.color <- "#ADD8E6"

## Colors from Franco
# wt.color = "#38459C"
# het.color = "#F89622"
# mut.color = "#CB2026"
# na.color = "#D3D3D3"

## Colors from Gotcha paper  
mut.color = "#BA242A"
wt.color = "#333D84"
na.color = "#D3D3D3"

vexas.genotyping.palette <- c(`WT` = wt.color, `MUT` = mut.color, `NA` = na.color)


############################## Misc ###############################
# vexas.vs.control.palettes <- list(
#   `HSC` = list(VEXAS = npg.pal[4], Control = npg.pal[6]),
#   `EMP` = list(VEXAS = npg.pal[1], Control = npg.pal[5]),
#   `Early Eryth` = list(VEXAS = npg.pal[1], Control = npg.pal[5]),
#   `GMP` = list(VEXAS = npg.pal[3], Control = npg.pal[7]),
#   `CD14 Mono` = list(VEXAS = npg.pal[3], Control = npg.pal[7]),
#   `LMPP` = list(VEXAS = npg.pal[9], Control = npg.pal[10])
# )
# vexas.HSC.color <- npg.pal[4]
# control.HSC.color <- npg.pal[6]
# 
# vexas.mono.color <- npg.pal[1]
# control.mono.color <- npg.pal[5]
# npg.pal <- pal_npg()(15)
# 

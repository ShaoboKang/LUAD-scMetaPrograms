library(scatterpie)
library(gtools)

CARD.visualize.pie2 <- function(proportion,spatial_location,colors = NULL,round_size){
  res_CARD = as.data.frame(proportion)
  res_CARD = res_CARD[,mixedsort(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  location$x = -location$x
  if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorCandidate = c("#1e77b4","#ff7d0b","#ceaaa3","#2c9f2c","#babc22","#d52828","#9267bc",
                     "#8b544c","#e277c1","#d42728","#adc6e8","#97df89","#fe9795","#4381bd","#f2941f","#5aa43a","#cc4d2e","#9f83c8","#91675a",
                     "#da8ec8","#929292","#c3c237","#b4e0ea","#bacceb","#f7c685",
                     "#dcf0d0","#f4a99f","#c8bad8",
                     "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
                     "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
                     "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
                     "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785")
  if(is.null(colors)){
    #colors = brewer.pal(11, "Spectral")
    if(ncol(res_CARD) > length(colorCandidate)){
      colors = colorRampPalette(colorCandidate)(ncol(res_CARD))
    }else{
      colors = colorCandidate[sample(1:length(colorCandidate),ncol(res_CARD))]
    }
  }else{
    colors = colors
  }
  data = cbind(res_CARD,location)
  ct.select = colnames(res_CARD)
  p = suppressMessages(ggplot() + geom_scatterpie(aes(x=y, y=x,r = round_size),data=data,
                                                  cols=ct.select,color=NA) + coord_fixed(ratio = 1) +
                         scale_fill_manual(values =  colors)+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 16,face="bold"),
                               legend.text=element_text(size = 15),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm'),
                               strip.text = element_text(size = 16,face="bold"),
                               legend.position="bottom")+
                         guides(fill=guide_legend(title="Cell Type")))
  return(p)
}

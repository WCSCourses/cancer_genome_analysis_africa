## By Marco; Cristiane; Mariana Boroni on 28/02/2021
## Functions to plot CIBERSORTx results

plot.ciber.heat <- function(ciber.obj, ann_info, needHelp=F, is.Absolute=T, sample.column=1){
  
  if(needHelp==T){
    print("ciber.obj (required) is a df or tibble with cell types and sampleID in columns\n")
    print("ann_info is a data.frame with samples as rownames and infos to annotateas columns")
    print("is.Absolute assumes you ran CIBERSORTx in Absolute mode and normalize sample wise")
    print("sample.column define in which column is your sample ID")
  }
  require(data.table)
  require(RColorBrewer)
  require(dplyr)
  require(ggdendro)
  require(gridExtra)
  require(grid)
  require(ggplot2)
  require(cowplot)
  set.seed(123)
  
  # default values
  ciber.obj = as.data.frame(ciber.obj)
  sample.column=colnames(ciber.obj)[sample.column]
  value.columns <- colnames(ciber.obj)[!colnames(ciber.obj) %in% sample.column]
  
  ann.label1 <- colnames(ann_info)[1]
  if(ncol(ann_info)>1){ann.label2 <- colnames(ann_info)[2]} # to be tested
  ann_info[[sample.column]] <- rownames(ann_info)
  
  size.axis.Y = 12 # size do axis.text.y, mudar aqui para alterar no dend e barplot
  dend.Top = 50 # valor inversamente proporcional a altura do dendo, limite superior
  dend.Bot = 0 # valor do limite inferior da altura do dendo, alterar se for incluir multiplas labels
  height = F # T, F, plota valores de height pro dendo
  
  tema=list(theme_bw(),theme(axis.title=element_blank(), 
                             axis.text.y = element_text(size = size.axis.Y, vjust = 1), ### nao alterar
                             axis.text.x = element_blank(), #axis.text.y = element_blank(),
                             legend.position = "none",
                             plot.margin = unit(c(0,0,0,0),"mm")),
            scale_x_discrete(expand = c(0,0)))
  
  # normalize data
  if(is.Absolute==T){
    sum.ciber <- rowSums(ciber.obj %>% dplyr::select(-sample.column))
    ciber.obj[,value.columns] = ciber.obj[,value.columns]/sum.ciber*100
    scaleFUN <- function(x) sprintf("%.f", x*100) # , labels = scaleFUN
    scaleFUN.dend <- function(x) substr(round(x*10000, digits = 3),1,4) # , labels = scaleFUN
  }else{
    scaleFUN <- function(x) sprintf("%.f", x) # , labels = scaleFUN
    scaleFUN.dend <- function(x) round(x, digits = 1) # , labels = scaleFUN
  }
  
  # dendogram
  hc_complete <- hclust(as.dist(1-cor(t(ciber.obj[,value.columns]))), method="ward.D")
  ord <- hc_complete$order ## to change clustering order
  ciber.obj[[sample.column]] <- factor(ciber.obj[[sample.column]], levels = ciber.obj[[sample.column]][ord]) ## to change clustering order
  dendogram <- as.dendrogram(hc_complete)
  ddata <- ggdendro::dendro_data(dendogram, type = "rectangle")
  
  ## IMPROVE HERE - dend alignement left and right
  dend = ggplot(ggdendro::segment(ddata))+ 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    theme(axis.title = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = size.axis.Y, vjust = -0.2), axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(dend.Top,0.5,dend.Bot,0.5),"mm"))+
    scale_y_continuous(labels = scaleFUN.dend, expand = c(0,0), position = "left")+
    scale_x_continuous(expand = c(0,0))
  if(height == F){dend = dend+theme(axis.text.y = element_text(color = "white"), axis.ticks.y = element_blank())}
  
  ## barplot
  colors <- c(brewer.pal(12,"Paired"), brewer.pal(8,"Dark2"), brewer.pal(8,"Set2"), brewer.pal(8,"Set1"),  brewer.pal(8,"Set3"),  brewer.pal(8,"Accent"), brewer.pal(8,"Spectral"), "red") # suited for 20 cell types max
  p1=suppressWarnings(ciber.obj %>% melt() %>% 
                        ggplot(aes_string(x = sample.column, y = 'value', fill = 'variable'))+
                        geom_bar(width = 1, stat = "identity")+
                        labs(fill = "Subpopulation", y = "Relative percent")+
                        ylab("teste")+tema+
                        scale_y_continuous(expand = c(0,0), position = "left")+ ### nao alterar
                        guides(fill = guide_legend(ncol = 3))+
                        scale_fill_manual(values = colors))
  
  # barplot legend
  legend <- cowplot::get_legend(p1+theme(legend.text=element_text(size=8), legend.position ="bottom", legend.key.size = unit(0.5,"line")))
  
  
  ## ann_label
  plot.ann_label <- function(data, ann_label, palette){
    # check size of palette and group
    if(brewer.pal.info[palette,][1] < length(unique(data[[ann_label]]))){
      print("Provided more labels than colors in palette, Ramping palette instead")
      palette = colorRampPalette(brewer.pal(name=palette, n = 8))(length(unique(data[[ann_label]]))) # 8 cobre todas
      color=scale_fill_manual(values = palette)
    }else{color <- scale_fill_brewer(palette = palette)}
    
    ## to match clustering order
    data[[sample.column]] <- factor(data[[sample.column]], levels = ciber.obj[[sample.column]][ord])
    
    ggplot(data = data, aes_string(x=sample.column, y='value', fill=ann_label))+
      geom_col(width = 1,position = "fill")+tema+
      theme(axis.text.x = element_blank(), axis.ticks = element_blank(),axis.text.y = element_text(colour = "white"))+
      scale_y_continuous(expand=c(0,0), position = "left", labels = scaleFUN)+color
  }
  
  if(length(ann_info) >2){
    p.ann1=plot.ann_label(data = suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info), ann_label = ann.label1, palette = "Set3")
    p.ann2=plot.ann_label(data = suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info), ann_label = ann.label2, palette = "Paired")
    p.ann <- cowplot::plot_grid(p.ann2, p.ann1, ncol = 1)
    ann_legend1 = cowplot::get_legend(p.ann1+theme(legend.text=element_text(size=10), legend.position="right", legend.key.size = unit(0.5,"line"))+guides(fill=guide_legend(ncol=3)))
    ann_legend2 = cowplot::get_legend(p.ann2+theme(legend.text=element_text(size=10), legend.position="right", legend.key.size = unit(0.5,"line"))+guides(fill=guide_legend(ncol=3)))
    legends <- cowplot::plot_grid(ann_legend1, ann_legend2, legend, ncol=1, align = "v", rel_heights = c(2,2,2))
    
  }else{
    p.ann=plot.ann_label(data = suppressWarnings(ciber.obj %>% melt()) %>% inner_join(ann_info), ann_label = ann.label1, palette = "Set3")
    ann_legend <- cowplot::get_legend(p.ann+theme(legend.text=element_text(size=10), legend.position="right", legend.key.size = unit(0.8,"line"))+guides(fill=guide_legend(ncol=3)))
    legends <- cowplot::plot_grid(ann_legend, NULL, legend, ncol=1, align = "v", rel_heights = c(1,0,1))
  }
  
  ## Plot
  plots <- cowplot::plot_grid(dend, p.ann, p1, ncol=1, rel_heights = c(4,0.3,6))
  cowplot::plot_grid(plots, legends, rel_heights = c(2,1), ncol = 2, align = "hv")
}

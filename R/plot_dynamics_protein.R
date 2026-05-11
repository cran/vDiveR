#' Dynamics of Diversity Motifs (Protein) Plot
#'
#' This function compactly display the dynamics of diversity motifs (index and its variants: major, minor and unique)
#' in the form of dot plot(s) as well as violin plots for all the provided individual protein(s).
#'
#' @param df DiMA JSON converted csv file data
#' @param host number of host (1/2)
#' @param protein_order order of proteins displayed in plot
#' @param base_size base font size in plot
#' @param alpha any number from 0 (transparent) to 1 (opaque)
#' @param line_dot_size dot size in scatter plot
#' @param bw smoothing bandwidth of violin plot (default: ucv)
#' @param adjust adjust the width of violin plot (default: 1)
#' @return A plot
#' @examples plot_dynamics_protein(proteins_1host)
#' @importFrom cowplot plot_grid ggdraw draw_label draw_plot
#' @importFrom ggplot2 ggplot aes geom_point geom_violin geom_boxplot scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 labs theme_classic theme element_rect element_text element_blank guides guide_legend
#' @importFrom ggplot2 scale_colour_manual scale_fill_manual facet_grid ggtitle coord_cartesian vars xlab ylab
#' @export
plot_dynamics_protein<-function(df, 
                                host=1, 
                                protein_order=NULL, 
                                base_size=11, 
                                alpha = 1/3, 
                                line_dot_size = 3,
                                bw = "ucv",
                                adjust = 1){
    #single host
    if (host == 1){
        generate_protein_plots(data=df, protein_order=protein_order, base_size=base_size,alpha=alpha,  line_dot_size=line_dot_size, bw = bw, adjust = adjust)
    }else{ #multihost
        #split the data into multiple subsets (if multiple hosts detected)
        data_list<-split(df,df$host)
        multihost_plots <- lapply(data_list, function(df) {      
          generate_protein_plots(df, protein_order = protein_order, base_size = base_size, alpha = alpha,
                  line_dot_size = line_dot_size, bw = bw, adjust = adjust, host=host)
        })

        plot_grid(plotlist = multihost_plots,
              ncol = length(unique(df$host)))
    }
}

#' @importFrom ggpubr ggarrange
#' @importFrom cowplot plot_grid ggdraw draw_label draw_plot
generate_protein_plots<-function(data, protein_order=NULL,alpha=1/3, line_dot_size=3, base_size=11, host=1, bw = "ucv", adjust = 1){
    Total_Variants <- Incidence <- Group <- x <- proteinName <- entropy <- NULL

    scatter_plot_data<-data.frame()
    group_names<-c("Index",
                   "Major","Minor",
                   "Unique",
                   "Total variants","Distinct variants")

    for (i in 7:12){
        tmp<-data.frame(proteinName=data[1],position=data[2],incidence=data[i],total_variants=data[11],Group=group_names[i-6],Multiindex=data[13])

        names(tmp)[3]<-"Incidence"
        names(tmp)[4]<-"Total_Variants"
        scatter_plot_data<-rbind(scatter_plot_data,tmp)
    }

    scatter_plot_data$proteinName <- toupper(scatter_plot_data$proteinName)
    if (!is.null(protein_order) && protein_order != ""){
        #order the proteins based on user input
        protein_order <- toupper(trimws(protein_order))
        level<-strsplit(protein_order, ',')[[1]]
        level <- sapply(level, function(x) toupper(trimws(x)))
        #set protein order as factor
        
        scatter_plot_data$proteinName<-factor(scatter_plot_data$proteinName, levels=level)
        scatter_plot_data$size_f = factor(scatter_plot_data$proteinName,levels = level)
    }
    scatter_plot_data$Group<-factor(scatter_plot_data$Group, levels=c("Index","Total variants", "Major", "Minor", "Unique", "Distinct variants"))
    violin_plot_data<-scatter_plot_data

    #plot plot 4
    scatter_plots<-ggplot()+geom_point(scatter_plot_data,mapping=aes(x=Total_Variants,y=Incidence,color=Group),alpha=alpha,size=line_dot_size)+
        scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
        labs(y = NULL,x= NULL)+
        theme_classic(base_size = base_size)+
        theme(
            legend.background = element_rect(fill = "transparent"),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title   = element_text(size = base_size, face = "bold"),
            axis.title.x = element_text(size = base_size, face = "bold"),
            axis.title.y = element_text(size = base_size, face = "bold"),
            strip.text.x = element_text(size = base_size-1, face = "bold"),
            legend.position = "bottom"
        )+ guides(colour = guide_legend(override.aes = list(alpha = 1,size=2),keywidth = 1,keyheight = 0.1,nrow=1, byrow=TRUE))+
        scale_colour_manual('',values = c("Index"="black","Total variants"="#f7238a", "Major"="#37AFAF","Minor"="#42aaff","Unique"="#af10f1","Distinct variants"="#c2c7cb" ))
    scatter_plots<-scatter_plots+facet_grid(cols=vars(scatter_plot_data$proteinName))

    #host label
    if("host" %in% colnames(data)){
        scatter_plots<-scatter_plots+ggtitle(unique(data$host))+
            theme(plot.title = element_text(hjust = 0.5))
    }
    
    # Add y-label to scatter plot
    scatter_with_ylabel <- ggdraw() +
        draw_label("Incidence (%)",
                   x = 0.02, y = 0.5,
                   angle = 90,
                   fontface = "bold", size = base_size) +
        draw_plot(scatter_plots, x = 0.06, y = 0, width = 0.94, height = 1)
    
    if (length(unique(data$proteinName)) <=10){
      # Define violin plot theme
      theme_violin <- theme_classic(base_size = base_size) +
          theme(
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
              plot.title   = element_text(size = base_size, face = "bold", hjust = 0.5),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()
          )
      
      #prepare the data for each subplot of violin_plots
      index<-violin_plot_data[violin_plot_data$Group %in% c("Index"),]
      major<-violin_plot_data[violin_plot_data$Group %in% c("Major"),]
      minor<-violin_plot_data[violin_plot_data$Group %in% c("Minor"),]
      unique<-violin_plot_data[violin_plot_data$Group %in% c("Unique"),]
      nonatypes<-violin_plot_data[violin_plot_data$Group %in% c("Distinct variants"),]
      variants_max_yaxis<-ceiling((max(as.numeric(major$Incidence),as.numeric(minor$Incidence),as.numeric(unique$Incidence))/10))*10
  
      #plot 5
      violin_plot_index<-ggplot(index, aes(x=proteinName, y=Incidence))+
          geom_violin(fill="black",trim = TRUE, color="black",alpha=0.9, adjust=adjust, bw=bw)+
          geom_boxplot(outlier.shape = NA,width=0.08, color="white",alpha=0.20,fill="white")+
          coord_cartesian(ylim = c(0, 100))+
          ggtitle("Index k-mer (%)")+
          theme_violin
  
      violin_plot_tv<-ggplot(index, aes(x=proteinName, y=Total_Variants))+
          geom_violin(fill="#f7238a",trim = TRUE, color="#f7238a",alpha=0.9, adjust=adjust, bw=bw)+
          geom_boxplot(outlier.shape = NA,width=0.08, color="black",alpha=0.20,fill=NA)+
          coord_cartesian(ylim = c(0, 100))+
          ggtitle("Total variants (%)")+
          theme_violin
  
      violin_plot_major<-ggplot(major, aes(x=proteinName, y=Incidence)) +
          geom_violin(fill="#37AFAF",trim = TRUE, color="#37AFAF", alpha=0.9, adjust=adjust, bw=bw)+
          geom_boxplot(outlier.shape = NA,width=0.08, color="black", alpha=0.20,fill=NA)+
          coord_cartesian(ylim = c(0, variants_max_yaxis))+
          ggtitle("Major variants (%)")+
          theme_violin
  
      violin_plot_minor<-ggplot(minor, aes(x=proteinName, y=Incidence))+
          geom_violin(fill="#42aaff",trim = TRUE,color="#42aaff", alpha=0.9, adjust=adjust, bw=bw)+
          geom_boxplot(outlier.shape = NA,width=0.08, color="black", alpha=0.20,fill=NA)+
          coord_cartesian(ylim = c(0, variants_max_yaxis))+
          ggtitle("Minor variants (%)")+
          theme_violin
  
      violin_plot_unique<-ggplot(unique, aes(x=proteinName, y=Incidence)) +
          geom_violin(fill="#af10f1",trim = TRUE, color="#af10f1", alpha=0.9, adjust=adjust, bw=bw)+
          geom_boxplot(outlier.shape = NA,width=0.08, color="black", alpha=0.20,fill=NA)+
          coord_cartesian(ylim = c(0, variants_max_yaxis))+
          ggtitle("Unique variants (%)")+
          theme_violin
  
      violin_plot_nonatypes<-ggplot(nonatypes, aes(x=proteinName, y=Incidence)) +
          geom_violin(fill="#c2c7cb",trim = TRUE, color="#c2c7cb", alpha=0.9, adjust=adjust, bw=bw)+
          geom_boxplot(outlier.shape = NA,width=0.08, color="black", alpha=0.20,fill=NA) +
          coord_cartesian(ylim = c(0, 100))+
          ggtitle("Distinct variants (%)")+
          theme_violin
      
      violin_plots<-ggarrange(violin_plot_index,violin_plot_tv,violin_plot_nonatypes,violin_plot_major,violin_plot_minor,violin_plot_unique,ncol=3,nrow=2)
      
      # Add y-label to violin panel
      violin_with_ylabel <- ggdraw() +
          draw_label("Incidence (%)",
                     x = 0.02, y = 0.5,
                     angle = 90,
                     fontface = "bold", size = base_size) +
          draw_plot(violin_plots, x = 0.06, y = 0, width = 0.94, height = 1)
      
    } else {
      violin_plot_data$Group[violin_plot_data$Group == "Index"] <- "Index k-mer" 
      violin_plot_data$Group[violin_plot_data$Group == "Major"] <- "Major variant" 
      violin_plot_data$Group[violin_plot_data$Group == "Minor"] <- "Minor variants" 
      violin_plot_data$Group[violin_plot_data$Group == "Unique"] <- "Unique variants" 
      
      violin_plot_data$Group<-factor(violin_plot_data$Group, levels=c("Index k-mer","Total variants", "Distinct variants", "Major variant", "Minor variants", "Unique variants"))
      variants<-subset(violin_plot_data, Group=="Major variant" | Group=="Minor variants" | Group=="Unique variants")
      max_ylim<-ceiling((max(variants$Incidence)/10))*10
      
      breaks_fun <- function(x) {
        if (max(x)<= max_ylim){
          seq(0,max_ylim,10)
        }else{
          seq(0,100,20)
        }
      }
      
      limits_fun <- function(x) {
        if (max(x)<= max_ylim){
          c(0,max_ylim)
        }else{
          c(0,100)
        }
      }
      
      violin_plots<-ggplot()+
        geom_violin(data=violin_plot_data,aes(x=proteinName,y=Incidence, fill=Group, color=Group), trim=TRUE, adjust=adjust, bw=bw)+
        theme_classic(base_size = base_size)+xlab("Protein")+ylab("Incidence (%)\n")+
        theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
              legend.position="none")+
        scale_y_continuous(limits = limits_fun,breaks = breaks_fun)+
        facet_grid(rows = vars(Group),switch="y",scales = 'free')+
        scale_colour_manual('',values = c("Index k-mer"="black","Total variants"="#f7238a", "Major variant"="#37AFAF","Minor variants"="#42aaff","Unique variants"="#af10f1","Nonatypes"="#c2c7cb" ))+
        scale_fill_manual('',values = c("Index k-mer"="black","Total variants"="#f7238a", "Major variant"="#37AFAF","Minor variants"="#42aaff","Unique variants"="#af10f1","Nonatypes"="#c2c7cb" ))
      
      violin_with_ylabel <- violin_plots
    }
    
    plot_grid(scatter_with_ylabel, violin_with_ylabel, ncol = 1, rel_heights = c(1.6, 1))
}


#' Dynamics of Diversity Motifs (Protein) Plot
#'
#' This function compactly display the dynamics of diversity motifs (index and its variants: major, minor and unique)
#' in the form of dot plot(s) as well as violin plots for all the provided individual protein(s).
#'
#' @param df DiMA JSON converted csv file data
#' @param host number of host (1/2)
#' @param proteinOrder order of proteins displayed in plot
#' @param base_size base font size in plot
#' @param alpha any number from 0 (transparent) to 1 (opaque)
#' @param dot_size dot size in scatter plot
#' @return A plot
#' @examples plot_dynamics_protein(proteins_1host)
#' @importFrom gridExtra grid.arrange
#' @export
plot_dynamics_protein<-function(df, host=1, proteinOrder="", base_size=8, alpha = 1/3, dot_size = 3){
    #single host
    if (host == 1){
        plot4_5(df,proteinOrder, alpha, dot_size, base_size)
    }else{ #multihost
        #split the data into multiple subsets (if multiple hosts detected)
        plot4_list<-split(df,df$host)
        plot4_multihost<-lapply(plot4_list,plot4_5,proteinOrder, alpha, dot_size, base_size)

        #create spacing between multihost plots
        theme = theme(plot.margin = unit(c(0.5,1.0,0.1,0.5), "cm"))
        do.call("grid.arrange", c(grobs=lapply(plot4_multihost,"+",theme), ncol = length(unique(df$host))))
    }
}

#' @importFrom ggplot2 guides guide_legend scale_colour_manual ggtitle element_text
#' @importFrom ggplot2 geom_violin geom_boxplot ylim scale_color_grey margin element_line
#' @importFrom ggplot2 scale_fill_manual theme_bw facet_grid xlab ylab
#' @importFrom ggpubr annotate_figure ggarrange text_grob
plot4_5<-function(data, proteinOrder="",alpha=1/3, dot_size=3, base_size=8){
    Total_Variants <- Incidence <- Group <- x <- proteinName <- entropy <- NULL

    plot4_data<-data.frame()
    group_names<-c("Index","Major","Minor","Unique","Total variants","Distinct variants")

    for (i in 7:12){
        tmp<-data.frame(proteinName=data[1],position=data[2],incidence=data[i],total_variants=data[11],Group=group_names[i-6],Multiindex=data[13])
        names(tmp)[3]<-"Incidence"
        names(tmp)[4]<-"Total_Variants"
        plot4_data<-rbind(plot4_data,tmp)
    }

    if (proteinOrder !=""){
        #order the proteins based on user input
        level<-strsplit(proteinOrder, ',')[[1]]
        #set protein order as factor
        plot4_data$proteinName<-factor(plot4_data$proteinName, levels=level)
        plot4_data$size_f = factor(plot4_data$proteinName,levels = level)
    }

    plot5_data<-plot4_data

    #plot plot 4
    plot4<-ggplot()+geom_point(plot4_data,mapping=aes(x=Total_Variants,y=Incidence,color=Group),alpha=alpha,size=dot_size)+
        scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
        labs(y = "Incidence (%)",x= NULL)+
        theme_classic(base_size = base_size)+
        theme(
            legend.background = element_rect(fill = "transparent"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "bottom"
        )+ guides(colour = guide_legend(override.aes = list(alpha = 1,size=2),keywidth = 1,keyheight = 1,nrow=1,byrow=TRUE))+
        scale_colour_manual('',values = c("Index"="black","Total variants"="#f7238a", "Major"="#37AFAF","Minor"="#42aaff","Unique"="#af10f1","Distinct variants"="#c2c7cb" ))
    plot4<-plot4+facet_grid(cols=vars(plot4_data$proteinName))

    #host label
    if("host" %in% colnames(data)){
        plot4<-plot4+ggtitle(unique(data$host))+
            theme(plot.title = element_text(hjust = 0.5))
    }

    #prepare the data for each subplot of plot5
    index<-plot5_data[plot5_data$Group %in% c("Index"),]
    major<-plot5_data[plot5_data$Group %in% c("Major"),]
    minor<-plot5_data[plot5_data$Group %in% c("Minor"),]
    unique<-plot5_data[plot5_data$Group %in% c("Unique"),]
    nonatypes<-plot5_data[plot5_data$Group %in% c("Distinct variants"),]
    variants_max_yaxis<-ceiling((max(as.numeric(major$Incidence),as.numeric(minor$Incidence),as.numeric(unique$Incidence))/10))*10

    #plot 5
    plot5_index<-ggplot(index, aes(x=proteinName, y=Incidence))+
        geom_violin(fill="black",trim = FALSE, color="black",alpha=0.9)+ylim(0,100)+ylab("Index k-mer (%)")+xlab("") +theme_bw() +
        geom_boxplot(outlier.shape = NA,width=0.05, color="white",alpha=0.15,fill="white")+
        theme_classic(base_size = base_size)+
        theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks.x = element_blank())

    plot5_tv<-ggplot(index, aes(x=proteinName, y=Total_Variants))+
        geom_violin(fill="#f7238a",trim = FALSE, color="#f7238a",alpha=0.9)+ylim(0,100)+ylab("Total variant (%)")+xlab("") +theme_bw() +
        geom_boxplot(outlier.shape = NA,width=0.05, color="black",alpha=0.15,fill="white")+
        theme_classic(base_size = base_size)+
        theme(plot.margin = unit(c(0,0.1,0.1,0.1), "cm"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks.x = element_blank())

    plot5_major<-ggplot(major, aes(x=proteinName, y=Incidence)) +
        geom_violin(fill="#37AFAF",trim = FALSE, color="#37AFAF")+ylim(0,variants_max_yaxis)+ylab("Major variant (%)")+xlab("")+theme_bw() +
        geom_boxplot(outlier.shape = NA,width=0.04, color="black", alpha=0.15,fill="white")+
        theme_classic(base_size = base_size)+
        theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks.x = element_blank(),
              axis.text.y  = element_text(face="bold"))

    plot5_minor<-ggplot(minor, aes(x=proteinName, y=Incidence))+
        geom_violin(fill="#42aaff",trim = FALSE,color="#42aaff")+ylim(0,variants_max_yaxis)+ylab("Minor variants (%)")+xlab("") +theme_bw() +
        geom_boxplot(outlier.shape = NA,width=0.04, color="black", alpha=0.15,fill="white")+
        theme_classic(base_size = base_size)+
        theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks.x = element_blank(),
              axis.text.y  = element_text(face="bold"))

    plot5_unique<-ggplot(unique, aes(x=proteinName, y=Incidence)) +
        geom_violin(fill="#af10f1",trim = FALSE, color="#af10f1")+ylim(0,variants_max_yaxis)+ylab("Unique variants (%)")+xlab("")+theme_bw() +
        geom_boxplot(outlier.shape = NA,width=0.05, color="black", alpha=0.15,fill="white")+
        theme_classic(base_size = base_size)+
        theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks.x = element_blank(),
              axis.text.y  = element_text(face="bold"))

    plot5_nonatypes<-ggplot(nonatypes, aes(x=proteinName, y=Incidence)) +
        geom_violin(fill="#c2c7cb",trim = FALSE, color="#c2c7cb")+ylim(0,100)+ylab("Distinct variants (%)")+xlab("")+theme_bw()+
        geom_boxplot(outlier.shape = NA,width=0.05, color="black", alpha=0.15,fill="white") +
        theme_classic(base_size = base_size)+
        theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.ticks.x = element_blank())
    plot5<-ggarrange(plot5_index,plot5_tv,plot5_nonatypes,plot5_major,plot5_minor,plot5_unique,ncol=3,nrow=2)

    #plot4_5
    ggarrange(plot4,plot5,ncol=1,heights = c(1,0.5))
}





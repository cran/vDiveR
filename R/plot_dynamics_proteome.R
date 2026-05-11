#' Dynamics of Diversity Motifs (Proteome) Plot
#'
#' This function compactly display the dynamics of diversity motifs (index and its variants: major, minor and unique)
#' in the form of dot plot as well as violin plot for all the provided proteins at proteome level.
#'
#' @param df DiMA JSON converted csv file data
#' @param host number of host (1/2)
#' @param line_dot_size size of dot in plot
#' @param base_size word size in plot
#' @param alpha any number from 0 (transparent) to 1 (opaque)
#' @param bw smoothing bandwidth of violin plot (default: nrd0)
#' @param adjust adjust the width of violin plot (default: 1)
#' @return A plot
#' @examples plot_dynamics_proteome(proteins_1host)
#' @importFrom cowplot plot_grid
#' @export
plot_dynamics_proteome <- function(df,
                                   host=1,
                                   line_dot_size=2,
                                   base_size=10,
                                   alpha=1/3,
                                   bw = "nrd0",
                                   adjust = 1){
    #single host
    if (host == 1){
        generate_plots(data=df, base_size= base_size, alpha=alpha, line_dot_size=line_dot_size, bw = bw, adjust = adjust)
    }else{ #multihost
        #split the data into multiple subsets (if multiple hosts detected)
        data_list<-split(df,df$host)
        multihost_plots <- lapply(data_list, function(df) {
            generate_plots(df, line_dot_size = line_dot_size, base_size = base_size, 
                    alpha = alpha, bw = bw, adjust = adjust, host = host)
        })
        
        #create spacing between multihost plots
        theme = theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
        plot_grid(plotlist = lapply(multihost_plots, '+', theme),
                  ncol = length(unique(df$host)))
    }

}

#' @importFrom ggplot2 element_blank facet_wrap scale_colour_manual
#' @importFrom ggplot2 guides guide_legend scale_color_manual ggtitle element_text geom_violin geom_boxplot ylim scale_color_grey margin scale_fill_manual coord_cartesian
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#' @importFrom patchwork plot_layout wrap_plots
generate_plots<-function(data,
                        line_dot_size=2,
                        base_size=10,
                        host = 1,
                        alpha=1/3,
                        bw = "nrd0",
                        adjust = 1){

    Total_Variants <- Incidence <- Group <- x <- NULL
    df<-data.frame()
    group_names<-c("Index",
                   "Major", "Minor", "Unique",
                   "Total variants","Distinct variants")

    #transpose the data format
    for (i in 7:12){
        tmp<-data.frame(proteinName=data[1],
                        position=data[2],
                        incidence=data[i],
                        total_variants=data[11],
                        Group=group_names[i-6],
                        Multiindex=data[13])
        names(tmp)[3]<-"Incidence"
        names(tmp)[4]<-"Total_Variants"
        df<-rbind(df,tmp)
    }

    df_violin<-df
    minor<-rbind(df[df$Group == "Index",],df[df$Group == "Total variants",])
    uniq<-rbind(df[df$Group == "Index",],df[df$Group == "Total variants",])
    minor$motif<- "Minor"
    uniq$motif<-"Unique"

    df<-df%>%mutate(motif = case_when(
        df$Group == "Index" ~ "Major",
        df$Group == "Total variants" ~ "Major",
        df$Group == "Major" ~ "Major",
        df$Group == "Minor"  ~ "Minor",
        df$Group == "Unique" ~ "Unique",
        df$Group == "Distinct variants" ~ "Distinct variants"
    ))
    df<- rbind(df,minor,uniq)
    df$motif<-factor(df$motif,levels = c("Major","Minor","Unique","Distinct variants"))

    # Define common theme components
    common_left <- theme(
        axis.text.y  = element_text(size = base_size),
        axis.ticks.y = element_line(),
        plot.margin  = margin(t = 5, r = 5, b = 5, l = 5)
    )
    
    theme_clean <- theme_classic(base_size = base_size) +
        theme(
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title   = element_text(size = base_size, face = "bold", hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank()
        ) + common_left

    #plotting scatter plot
    proteins_point_plot<-ggplot()+
        geom_point(df,mapping=aes(x=Total_Variants,y=Incidence,color=Group),alpha=alpha,size= line_dot_size)+
        geom_point(df,mapping = aes(x =Total_Variants,y=Incidence),
                   col=ifelse(df$multiIndex== TRUE & df$Group== "Index", 'red', ifelse(df$multiIndex== FALSE, 'white', 'white')), 
                   alpha=ifelse(df$multiIndex ==TRUE & df$Group== "Index", 1, ifelse(df$multiIndex== TRUE, 0,0)),
                   pch=1,size=3,stroke=1.05)+ #multiIndex
        scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
        theme_classic(base_size = base_size) +
        theme(
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title   = element_text(size = base_size, face = "bold"),
            axis.title.x = element_text(size = base_size, face = "bold"),
            axis.title.y = element_text(size = base_size, face = "bold"),
            strip.text.x = element_blank(),
            legend.position="bottom"
        ) + common_left +
        labs(y= "Incidence (%)", x="\nTotal variants (%)")+
        facet_wrap(~ motif,ncol = 1)+
        guides(colour = guide_legend(override.aes = list(alpha = 1,size=2), nrow = host, byrow=T, keywidth = 1, keyheight = .1 ))+
        scale_colour_manual('',breaks=c("Index","Total variants","Major","Minor","Unique","Distinct variants"),
                            values = c("Index"="black", "Total variants"="#f7238a","MultiIndex"="red",
                                       "Major"="#37AFAF" , "Minor"="#42aaff","Unique"="#af10f1", "Distinct variants"="#c2c7cb"))
    #host label
    if("host" %in% colnames(data)){
        proteins_point_plot<-proteins_point_plot+ggtitle(unique(data$host))
    }

    # Create individual violin plots
    index<-df_violin[df_violin$Group %in% "Index",]
    total_var<-df_violin[df_violin$Group %in% "Total variants",]
    distinct_var<-df_violin[df_violin$Group %in% "Distinct variants",]
    major_var<-df_violin[df_violin$Group %in% "Major",]
    minor_var<-df_violin[df_violin$Group %in% "Minor",]
    unique_var<-df_violin[df_violin$Group %in% "Unique",]

    p_index <- ggplot(index, aes(x = "", y = Incidence)) +
        geom_violin(fill = "black", color = "black", bw = bw, adjust = adjust) +
        geom_boxplot(width=0.08, alpha=0.20, fill="white", outlier.shape=NA, color="white") +
        coord_cartesian(ylim = c(0, 100)) +
        ggtitle("Index") + theme_clean

    p_total <- ggplot(index, aes(x = "", y = Total_Variants)) +
        geom_violin(fill = "#f7238a", color = "#f7238a", bw = bw, adjust = adjust) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA) +
        coord_cartesian(ylim = c(0, 100)) +
        ggtitle("Total variants") + theme_clean

    p_distinct <- ggplot(distinct_var, aes(x = "", y = Incidence)) +
        geom_violin(fill = "#c2c7cb", color = "#c2c7cb", bw = bw, adjust = adjust) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA) +
        coord_cartesian(ylim = c(0, 100)) +
        ggtitle("Distinct variants") + theme_clean

    p_major <- ggplot(major_var, aes(x = "", y = Incidence)) +
        geom_violin(fill = "#37AFAF", color = "#37AFAF", bw = bw, adjust = adjust) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA) +
        coord_cartesian(ylim = c(0, 50)) +
        ggtitle("Major variants") + theme_clean

    p_minor <- ggplot(minor_var, aes(x = "", y = Incidence)) +
        geom_violin(fill = "#42aaff", color = "#42aaff", bw = bw, adjust = adjust) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA) +
        coord_cartesian(ylim = c(0, 50)) +
        ggtitle("Minor variants") + theme_clean

    p_unique <- ggplot(unique_var, aes(x = "", y = Incidence)) +
        geom_violin(fill = "#af10f1", color = "#af10f1", bw = bw, adjust = adjust) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fill = NA) +
        coord_cartesian(ylim = c(0, 50)) +
        ggtitle("Unique variants") + theme_clean

    # 2 × 3 layout: Index | Total | Distinct / Major | Minor | Unique
    violin_panel <- (p_index | p_total | p_distinct) / (p_major | p_minor | p_unique)

    # Combine scatter plot and violin panel
    combined_plot <- proteins_point_plot / violin_panel +
        plot_layout(heights = c(1.6, 1))

    return(combined_plot)
}






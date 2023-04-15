#' Conservation Levels Distribution Plot
#'
#' This function plots conservation levels distribution of k-mer positions, which consists of
#' completely conserved (black) (index incidence = 100\%), highly conserved (blue)
#' (90\% <= index incidence < 100\%), mixed variable (green) (20\% < index incidence <= 90\%),
#' highly diverse (purple) (10\% < index incidence <= 20\%) and
#' extremely diverse (pink) (index incidence <= 10\%).
#'
#' @param df DiMA JSON converted csv file data
#' @param protein_order order of proteins displayed in plot
#' @param conservation_label 0 (partial; show present conservation labels only) or 1 (full; show ALL conservation labels) in plot
#' @param host number of host (1/2)
#' @param base_size base font size in plot
#' @param label_size conservation labels font size
#' @param alpha any number from 0 (transparent) to 1 (opaque)
#' @examples plot_conservationLevel(proteins_1host, conservation_label = 1,alpha=0.8, base_size = 15)
#' @examples plot_conservationLevel(protein_2hosts, conservation_label = 0, host=2)
#' @return A plot
#' @importFrom dplyr case_when
#' @importFrom grid unit
#' @importFrom gridExtra grid.arrange
#' @export
plot_conservationLevel <- function(df, protein_order="",conservation_label=1,host=1, base_size = 11, label_size = 2.6, alpha=0.6){
    df<-df%>%mutate(ConservationLevel = case_when(
        df$index.incidence == 100 ~ "Completely conserved (CC)",
        df$index.incidence >= 90 ~ "Highly conserved (HC)",
        df$index.incidence >= 20 ~ "Mixed variable (MV)",
        df$index.incidence >= 10  ~ "Highly diverse (HD)",
        df$index.incidence < 10 ~ "Extremely diverse (ED)"
    ))

    #single host
    if (host == 1){
        plot_plot7(data = df,protein_order = protein_order,conservation_label = conservation_label, base_size = base_size, label_size = label_size)
    }else{ #multihost

        #split the data into multiple subsets (if multiple hosts detected)
        plot7_list<-split(df,df$host)
        plot7_multihost<-lapply(plot7_list,plot_plot7, protein_order,conservation_label, base_size, label_size)

        #create spacing between multihost plots
        theme = theme(plot.margin = unit(c(2.5,1.0,0.1,0.5), "cm"))
        do.call("grid.arrange", c(grobs=lapply(plot7_multihost,"+",theme), nrow = length(unique(df$host))))
    }
}

#' @importFrom plyr ddply .
#' @importFrom ggplot2 position_jitter scale_colour_manual
#' @importFrom ggplot2 position_dodge coord_cartesian
#' @importFrom gghalves geom_half_boxplot geom_half_point
#' @importFrom ggtext geom_richtext
#plotting function
plot_plot7<- function(data,protein_order="",conservation_label=1, base_size = 11, label_size = 2.6, alpha =0.6){
    proteinName <- Total <- index.incidence <- NULL
    Label <- ConservationLevel <- NULL
    #add word 'protein' in front of each protein name
    data$proteinName<-paste("Protein",data$proteinName)
    #create data for proteome bar "All" from existing data
    data1<-data
    data1$proteinName <- "All"
    data1$level <- "All"
    #set up the order of proteins in plot from left to right
    if (protein_order ==""){ #follow the default order in csv file
        level<-c("All",unique(data$proteinName))
    }else{ #order the proteins based on user input
        level<-c("All",strsplit(protein_order, ',')[[1]])
    }

    #determine the protein order
    data$level = factor(data$proteinName, levels=level)
    #combine proteome bar with protein bars
    data<-rbind(data1,data)
    #determine the protein order
    data$level = factor(data$proteinName, levels=level)
    #calculation for total and percentage of conservation levels for each protein
    #sum up the total positions for each conservation level of proteins
    plot7_data<-ddply(data,.(proteinName,ConservationLevel),nrow)
    names(plot7_data)[3]<-"Total"

    C_level<- c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)")
    #check the presence of conservation level: insert value 0 if it is absent
    if (conservation_label == 1){ #full label
        #check the presence of conservation level: insert value 0 if it is absent
        for ( conservation in C_level){ #conservation level
            for (name in level){ #proteinName
                if (!(conservation %in% plot7_data[plot7_data$proteinName==name,]$ConservationLevel)){
                    plot7_data<-rbind(plot7_data,c(name,conservation,0))
                }}}
    }

    #sort the dataframe
    plot7_data[order(plot7_data$proteinName),]
    plot7_data$Total<- as.integer(plot7_data$Total)
    #get the percentage of each conservation level for each protein
    plot7_data<-ddply(plot7_data,.(proteinName),transform, percent=Total/sum(Total)*100)

    #gather the protein label in multicolor
    plot7_data<-plot7_data%>%mutate(Label = case_when(
        plot7_data$ConservationLevel == "Completely conserved (CC)" ~ paste0(sprintf("<span style =
    'color:#000000;'>CC: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
        plot7_data$ConservationLevel == "Highly conserved (HC)" ~ paste0(sprintf("<span style =
    'color:#0057d1;'>HC: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
        plot7_data$ConservationLevel == "Mixed variable (MV)" ~ paste0(sprintf("<span style =
    'color:#02d57f;'>MV: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
        plot7_data$ConservationLevel == "Highly diverse (HD)" ~ paste0(sprintf("<span style =
    'color:#A022FF;'>HD: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
        plot7_data$ConservationLevel == "Extremely diverse (ED)" ~ paste0(sprintf("<span style =
    'color:#ff617d;'>ED: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
    ))

    #set conservation level in specific order (CC,HC,MV,HD,ED)
    plot7_data<-plot7_data[order(factor(plot7_data$ConservationLevel, levels=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"))),]

    #combine all conservation level labels into one for each protein
    Proteinlabel<- aggregate(Label~proteinName, plot7_data, paste, collapse="<br>")
    #get number of protein for labelling
    nProtein<-nrow(Proteinlabel)

    #plotting
    ggplot(data, aes(x=level,y=index.incidence)) +
        # gghalfves
        geom_half_boxplot(outlier.shape = NA) +
        geom_half_point(aes(col = ConservationLevel), side = "r",
                        position = position_jitter(width = 0, height=-0.7),alpha=alpha) +
        ylim(0,105) +
        labs(x=NULL, y="Index incidence (%)\n", fill="Conservation level")+
        theme_classic(base_size = base_size)+
        theme(
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.position = 'bottom',
            plot.margin = unit(c(5, 1, 1, 1), "lines"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)
        ) +
        scale_colour_manual('Conservation Level',
                            breaks = c("Completely conserved (CC)",
                                       "Highly conserved (HC)",
                                       "Mixed variable (MV)",
                                       "Highly diverse (HD)",
                                       "Extremely diverse (ED)"),
                            values = c("Completely conserved (CC)"="black",
                                       "Highly conserved (HC)"="#0057d1",
                                       "Mixed variable (MV)"="#02d57f",
                                       "Highly diverse (HD)"="#8722ff",
                                       "Extremely diverse (ED)"="#ff617d")) +
        geom_richtext(data = Proteinlabel,
                      aes(x=proteinName,label = Label, y=c(rep(105,nProtein)),
                          label.size=0, label.color="transparent"),
                      position = position_dodge(width=0.1),
                      size=label_size, color="black", hjust=0, angle=90) +
        guides(color = guide_legend(override.aes = list(size = 2), nrow=2))+
        coord_cartesian(clip = "off")+ #allow ggtext outside of the plot
        ggtitle(unique(data$host))

}









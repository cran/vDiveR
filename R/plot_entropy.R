#' Entropy plot
#'
#' This function plot entropy (black) of each
#' k-mer position across the studied proteins and highlight region with zero entropy in yellow.
#'
#' @param df DiMA JSON converted csv file data
#' @param host number of host (1/2)
#' @param proteinOrder order of proteins displayed in plot
#' @param kmer_size size of the k-mer window
#' @param ymax maximum y-axis
#' @param line_dot_size size of the line and dot in plot
#' @param wordsize size of the wordings in plot
#' @return A plot
#' @examples plot_entropy(proteins_1host)
#' @examples plot_entropy(protein_2hosts, host = 2)
#' @importFrom ggplot2 geom_rect geom_area geom_hline geom_line
#' @importFrom ggplot2 sec_axis element_line scale_colour_manual scale_linetype_manual
#' @export
plot_entropy <- function(df,host=1,proteinOrder="",kmer_size=9, ymax = 10,line_dot_size=2,wordsize=8){
    entropy <- end <- lowSupportPos <- totalVariants.incidence <- NULL
    #determine number of host
    #scale the amino acid position for each protein
    if (host ==1){ #single host
        if (proteinOrder ==""){
            a<-table(df$proteinName)
            proteinName<-as.vector(names(a))
            position<-as.vector(a)
            df$size_f = factor(df$proteinName,levels = proteinName)
            scales_x<-mapply(function(x,y){
                x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
            }, proteinName,position)
        }else{
            #order the proteins based on user input
            level<-strsplit(proteinOrder, ',')[[1]]
            position<-c()
            for (i in level){
                position<-append(position,table(df$proteinName)[names(table(df$proteinName)) == i])
            }
            df$size_f = factor(df$proteinName,levels = level)
            scales_x<-mapply(function(x,y){
                x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
            }, level,position)
        }
    }else{ # multihost
        #categorise data based on host
        df$host<- factor(df$host)
        #count the aa length for each proteins (each host is expected to have same number of proteins with same length)
        df_sub<-df[df$host==unique(df$host[1]),]
        if (proteinOrder ==""){
            a<-table(df_sub$proteinName)
            proteinName<-as.vector(names(a))
            position<-as.vector(a)
            df$size_f = factor(df$proteinName,levels = proteinName)
            scales_x<-mapply(function(x,y){
                x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
            }, proteinName,position)
        }else{
            #order the proteins based on user input
            level<-strsplit(proteinOrder, ',')[[1]]
            position<-c()
            for (i in level){
                position<-append(position,table(df_sub$proteinName)[names(table(df_sub$proteinName)) == i])
            }
            df$size_f = factor(df$proteinName,levels = level)
            scales_x<-mapply(function(x,y){
                x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
            }, level,position)
        }
    }

    #---------------set zero entropy region-----------------#
    #check if the kmer positions are zero entropy
    #Reference for getting all the positions of kmer based on kmer size and the starting position of kmer:
    #https://stackoverflow.com/questions/31120552/generate-a-sequence-of-numbers-with-repeated-intervals
    #https://stackoverflow.com/questions/41637518/adding-a-shaded-rectangle-to-an-existing-plot-using-geom-rect

    #identify the kmer position with zero entropy
    #position: startPosition with zero entropy
    #end: startPosition + kmersize - 1
    df_zeroEntropy<- df %>%
        dplyr::group_by(proteinName)%>%
        dplyr::summarize(
            end = position[which(entropy == min(entropy))+kmer_size-1], #end: startPosition + kmersize - 1
            position = position[which(entropy == min(entropy))] #extract those starting positions with zero entropy
        )%>%
        as.data.frame()

    #concatenate df with df_zeroEntropy
    df <-merge(df,df_zeroEntropy,id="position",all=T)
    #replace NAN in column "zero entropy" with FALSE
    df$end[is.na(df$end)]<- -1

    #---------------set maximum y limit-----------------#
    #if the max y value in data < 10 => ymax = 10
    #if max y value in data > 10 => ymax = ceiling(max y value)
    if (max(df$entropy) <= 10){
       ymax <-10.0
    }else if (max(df$entropy) > 10){
       ymax <- ceiling(max(df$entropy))
    }

    breaks_fun <- function(x) {
        seq(0,max(x),50)
    }

    limits_fun <- function(x) {
        c(0,max(x))
    }

    #----------------plotting--------------------#
    #detect if low support present
    if (TRUE %in% df$lowSupport){
        df$lowSupportPos <- -0.3
        df$lowSupportPos[df$lowSupport ==TRUE]<- -0.5
        plot1<-ggplot(df) +
            geom_rect(inherit.aes = FALSE,
                      aes(xmin=position, xmax=end, ymin=-Inf, ymax=+Inf),
                      fill='#FFECAF', alpha=ifelse(df$end == -1, 0, 0.5))+
            geom_area(mapping = aes(x = position, y = entropy,color= "k-mer Entropy", linetype="k-mer Entropy"), show.legend=F)+
            geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), linewidth = (line_dot_size/10))+
            geom_point(mapping = aes(x = position,y=lowSupportPos),col=ifelse(df$lowSupportPos==-0.5, 'black', ifelse(df$lowSupportPos==-0.3, 'white', 'white')), alpha=ifelse(df$lowSupportPos==-0.5, 1, ifelse(df$lowSupportPos==-0.3, 0,0)),pch=17, linewidth = (line_dot_size/10))+
            scale_x_continuous(limits = limits_fun,breaks = breaks_fun)+
            labs(y = "k-mer entropy (bits)\n",x= "\nk-mer position (aa)",color = "#f7238a")+
            ylim(0, ymax)+
            theme_classic() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.y.right = element_line(color = "#f7238a"),
                axis.ticks.y.right = element_line(color = "#f7238a"),
                axis.title.y.right =  element_text(color = "#f7238a"),
                legend.position="bottom"
            )+
            scale_colour_manual("",
                                values = c("Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"="black",
                                           "k-mer Entropy" = "black"),
                                guide = guide_legend(override.aes=aes(fill=NA)))+
            scale_linetype_manual("",values=c("Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                              "k-mer Entropy" = 1))

        #number of host
        if (host == 1){ #one host
            plot1 +facet_grid(cols=vars(df$size_f),scales = 'free',space = "free",switch = "x")
        }else{ # multi host
            plot1 +facet_grid(rows = vars(df$host),cols=vars(df$size_f),scales = 'free',space = "free",switch = "both")
        }

    }else{
        plot1<-ggplot(df) +
            geom_rect(inherit.aes = FALSE,
                      aes(xmin=position, xmax=end, ymin=-Inf, ymax=+Inf),
                      fill='#FFECAF', alpha=ifelse(df$end == -1, 0, 0.5))+
            geom_area(mapping = aes(x = position, y = entropy,color= "k-mer Entropy", linetype="k-mer Entropy"),show.legend = F)+
            geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), linewidth= (line_dot_size/10))+
            scale_x_continuous(limits = limits_fun,breaks = breaks_fun)+
            labs(y = "k-mer entropy (bits)\n",x= "\nk-mer position (aa)",color = "#f7238a")+
            ylim(0, ymax)+
            theme_classic(base_size = wordsize) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.y.right = element_line(color = "#f7238a"),
                axis.ticks.y.right = element_line(color = "#f7238a"),
                axis.title.y.right =  element_text(color = "#f7238a"),
                legend.position="bottom"
            )+
            scale_colour_manual("",
                                values = c("Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"="black",
                                           "k-mer Entropy" = "black"),
                                guide = guide_legend(override.aes=aes(fill=NA)))+
            scale_linetype_manual("",values=c("Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                              "k-mer Entropy" = 1))
        #number of host
        if (host == 1){ #one host
            plot1 +facet_grid(cols=vars(df$size_f),scales = 'free',space = "free",switch = "x")
        }else{ # multi host
            plot1 +facet_grid(rows = vars(df$host),cols=vars(df$size_f),scales = 'free',space = "free",switch = "both")
        }
    }


}

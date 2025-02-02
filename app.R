
#NOTES: Always save as app.R

##################################################################################################################################
#Load data and libraries
##################################################################################################################################

if(!require(BiocManager)) { install.packages("BiocManager"); library(BiocManager) }
if(!require(devtools)) { install.packages("devtools"); library(devtools) }

if(!require(shiny)) { install.packages("shiny"); library(shiny) }
if(!require(bslib)) { install.packages("bslib"); library(bslib) }
if(!require(phyloseq)) { BiocManager::install("phyloseq"); library(phyloseq) }
if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse) }
if(!require(ggplot2)) { install.packages("ggplot2"); library(ggplot2) }
if(!require(viridis)) { install.packages("viridis"); library(viridis) }
if(!require(kableExtra)) { install.packages("kableExtra"); library(kableExtra) }
if(!require(microViz)) {install.packages(BiocManager::install("microViz")); library(microViz) }
if(!require(microbiomeMarker)) { BiocManager::install("microbiomeMarker"); library(microbiomeMarker) }
if(!require(DT)) { install.packages("DT"); library(DT) }
if(!require(Maaslin2)) { install.packages("Maaslin2"); library(Maaslin2) }
if(!require(fantaxtic)) { devtools::install_github("gmteunisse/fantaxtic"); library(fantaxtic) }
if(!require(fantaxtic)) { devtools::install_github("gmteunisse/ggnested"); library(fantaxtic) }
if(!require(ggpubr)) { devtools::install_github("kassambara/ggpubr"); library(ggpubr) }


# Set a global theme for ggplot2 plots
theme_set(theme_bw() + theme(text = element_text(size = 18))) # Adjust the size as needed

#Load in data
# Load the saved Rds file - this comes from /CEDAR/Microbiome/oral_microbiome_exploratory/CEDAR_microbiome_exploratory on my local machine. This is the OSCC data from PRJNA822685, filtered for only matched samples with AN and Tumor reads.
phy_obj <- readRDS("00data.rds")
phy_obj_fs <-  prune_taxa(taxa_sums(phy_obj) > 1, phy_obj)
  
#Define global color
colors_all= c("Tumor" = "#FF495C", "AN"="#083D77",  "abnormal" = "#FF495C", "control"="#083D77")

##################################################################################################################################
#Define UI
##################################################################################################################################
# Define UI for application (this is just an HTML wrapper)
## This defines how the application looks - which tabs there are, how the sidebar is formatted, etc. 
ui <- fluidPage(

  titlePanel("Microbiome Browser"),
  sidebarLayout(
    
    #The input goes into the sidebarpanel
    sidebarPanel( width=2, #set the width of the sidebar
      
      selectInput("Taxa_input", "Taxanomic rank",
                  choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species")),
      
      numericInput("Num_top_taxa", "Taxa plots: Number of top taxa to display", value = 10, min = 0, max = 25),
      
      selectInput("PCoA_distance", "PCoA Distance",
                  choices = c("bray", "jaccard")),

      numericInput("ldacutoff", "Lefse: LDA Cutoff", value = 2, min = 0, max = 10),
      
      selectInput("Lefse_taxa", "Lefse: Taxa",
                  choices = c("all", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
      
      selectInput("Lefsenorm", "Lefse: Normalization (default: CPM)",
                  choices = c("CPM", "CLR", "RLE", "CSS", "none"))
    ),
  
    
    mainPanel(
      tabsetPanel(
        
          # Read Me Tab
    tabPanel("Read Me",
      fluidPage(
        titlePanel("Welcome to the Microbiome Browser"),
        
        # Add a formatted text section
        fluidRow(
          column(12, 
            p("This app is designed to demonstrate some basic microbiome analysis techniques on a real-world dataset. 
              The data used in this example was taken from PRJNA822685, a pubically avaliable dataset of oral sqaumous cell carcinoma tumor (abnormal) 
              and normal tumor-adjacent (control) tissue. This publication is avaliable at: https://pubmed.ncbi.nlm.nih.gov/38589501/ and the citation is below:" ),
            p("Cai L, Zhu H, Mou Q, Wong PY, Lan L, Ng CWK, Lei P, Cheung MK, Wang D, Wong EWY, Lau EHL, Yeung ZWC, Lai R, Meehan K, Fung S, Chan KCA, Lui VWY, Cheng ASL, Yu J, Chan PKS, Chan JYK, Chen Z. 
              Integrative analysis reveals associations between oral microbiota dysbiosis and host genetic and epigenetic aberrations in oral cavity squamous cell carcinoma.", 
              em("NPJ Biofilms Microbiomes."), "2024 Apr 8;10(1):39. doi: 10.1038/s41522-024-00511-x. PMID: 38589501; PMCID: PMC11001959."),
            
            p("The dataset was curated to include only tumor and normal tumor tissue from patients with matched samples."),
            p("To explore this dataset, click on the tabs. For analyses with adjustable parameters such as ", strong("PCoA"), "and", strong("LeFSe"), "adjust the sidebars accordingly."),
            hr(),
            p("Key Features:", style = "font-weight: bold; font-size: 16px;"),
            tags$ul(
              tags$li("Dynamically generated plots!"),
              tags$li("Interactive tables for results!"),
              tags$li("Customizable parameters!")
            ),
            hr(),
            p("For more details about sample processing, please visit Madeline Krieger's github page for this project ", 
              a("documentation", href = "hhttps://github.com/kriegerm/Microbiome_browser", target = "_blank"), ".")
          )
        )
      )
    ),

        
      #Tab outputs
      # You'll find the plotOutput() and tableOutput() functions in the "Render Output" section below.
      tabPanel("Average Taxa Plot", plotOutput("taxplot")),
      tabPanel("Taxa Plot by Sample", plotOutput("taxplot_sample")),
      tabPanel("Average Taxa Table", tableOutput("results")),
      tabPanel("Beta Diversity: PCA", plotOutput("PCA", height="1000px", width="800px"), plotOutput("PCoA", height="1000px", width="800px")),
      tabPanel("Beta Diversity: PCoA", plotOutput("PCoA", height="1000px", width="800px")),
      tabPanel("Alpha Diversity", plotOutput("Alphadiv")),
      tabPanel("LeFSe", plotOutput('lefse', height = "800px", width="800px"), dataTableOutput("lefse_results"))
               )
            )      
          )
        )


##################################################################################################################################
#Server
##################################################################################################################################

# Define server logic for creating output
# The Server function wraps all of the input/output, so the bottom bracket is waaaaay down at the end of the page above the "Run app" section.
server <- function(input, output) {
  
  
  ##################################################################################################################################
  #Define functions
  ##You'll use these functions for generating plots and stuff below in the "Render output" section
  ##################################################################################################################################
  
  # Plot PCA
  plot_pca <- function(phyloseq_obj, rank_transformation, variable, colors_list=NULL) {
    # Transform and calculate distance
    phylo_trans <- phyloseq_obj %>% tax_fix() %>% tax_transform(rank = rank_transformation, trans = "identity")
    dist_matrix <- dist_calc(phylo_trans, "euclidean")
    ord_res <- ord_calc(dist_matrix, "PCA")
    
    # Plot
    p <- ord_plot(ord_res, axes = c(1, 2), fill = variable, shape = variable, alpha = 0.8, size = 2) +
      ggtitle(label = paste0(
        rank_transformation)) +
      scale_shape_girafe_filled() +
      ggplot2::stat_ellipse(aes(color = !!sym(variable) ) )+
      theme(
        plot.title = element_text(face = "bold", size = 12, hjust = .5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 10, vjust = 0.5, hjust = 1)
      )
    
    # Add color scales if a custom colors_list is provided
    if (!is.null(colors_list)) {
      p <- p + scale_fill_manual(values = colors_list) +
        scale_color_manual(values = colors_list)
    }
    return(p)
  }
  
  
  #Plot PCoA
  plot_PCoA <- function(phyloseq_obj, rank_transformation, trans_type, dist_cal_type, ord_calc_method, variable, colors_list=NULL) {
    # Transform and calculate distance
    phylo_trans <- phyloseq_obj %>% tax_fix() %>%
      tax_transform(rank = rank_transformation, trans = trans_type)
    dist_matrix <- dist_calc(phylo_trans, dist_cal_type)
    ord_res <- ord_calc(dist_matrix, ord_calc_method)
    
    
    # Plot
    p <- ord_plot(ord_res, axes = c(1, 2), fill = variable, shape = variable, alpha = 0.8, size = 2, plot_taxa = 1:5, size = 2) +
      ggtitle(label = paste0(
        rank_transformation)) +
      scale_shape_girafe_filled() +
      ggplot2::stat_ellipse(aes(color =!!sym(variable) )) +
      theme(
        plot.title = element_text(face = "bold", size = 12, hjust = .5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 10, vjust = 0.5, hjust = 1)
      )
    
    # Add color scales if a custom colors_list is provided
    if (!is.null(colors_list)) {
      p <- p + scale_fill_manual(values = colors_list) +
        scale_color_manual(values = colors_list)
    }
    return(p)
  }
  
  
  ##################################################################################################################################
  #Define reactive variables
  ## These variables change as you modify the selection in the UI.
  ##################################################################################################################################

  #The averaged taxa table that will be displayed in the "Average Taxa Plot" tab.
  ave_filtered_for_taxplot <- reactive({
    phy_obj_fs %>%
      microbiomeMarker::normalize(., method="TSS") %>%
      tax_glom(., taxrank = input$Taxa_input, NArm = TRUE) %>%
      microViz::tax_fix(.) %>%
      fantaxtic::get_top_taxa(., input$Num_top_taxa, relative = TRUE, discard_other = TRUE) %>%
      psmelt() %>%
      dplyr::group_by(condition, !!sym(input$Taxa_input)) %>%
      dplyr::summarise(Average_Abundance = mean(Abundance))
  })
  
  #The averaged taxa table that will be displayed in the "Average Taxa Table" tab (the same as the one above, but with pivot_wider)
  ave_filtered_for_display <- reactive({
    taxa_input_sym <- sym(input$Taxa_input) # Convert to symbol
    
    phy_obj_fs %>%
      microbiomeMarker::normalize(., method="TSS") %>%
      tax_glom(., taxrank = input$Taxa_input, NArm = TRUE) %>%
      microViz::tax_fix(.) %>%
      fantaxtic::get_top_taxa(., input$Num_top_taxa, relative = TRUE, discard_other = TRUE) %>%
      psmelt() %>%
      dplyr::group_by(condition, !!sym(input$Taxa_input)) %>%
      dplyr::summarise(Average_Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = condition, values_from = Average_Abundance, names_prefix = "")
  })
  
  #This is the most basic taxa table generated - it is by sample, not averaged. Displayed in the "Taxa Plot by Sample" tab.
  filtered_for_taxplot <- reactive({
    
    phy_obj_fs %>%
      tax_glom(., taxrank = input$Taxa_input, NArm = TRUE) %>%
      microViz::tax_fix(.) %>%
      microbiomeMarker::normalize(., method="TSS") %>%
      get_top_taxa(., input$Num_top_taxa, relative = TRUE, discard_other = TRUE) %>%
      psmelt() %>%
      group_by(condition, !!sym(input$Taxa_input)) 
  })
  

  #Create a Lefse results object
  Lefse_results <- reactive({

        run_lefse(phy_obj_fs,
                               wilcoxon_cutoff = .05,
                               taxa_rank = input$Lefse_taxa,
                               norm = input$Lefsenorm, 
                               group = "condition",
                               kw_cutoff = .05,
                               multigrp_strat = TRUE,
                               lda_cutoff = input$ldacutoff)
    
          })
  
  ##################################################################################################################################
  #Render output
  ##################################################################################################################################
 
  # Render Average Taxa Table - have to use function() {} for kable styling
  output$results <- renderUI({
    filtered_data <- ave_filtered_for_display()  # Get the reactive filtered data
    HTML(
      filtered_data %>%
        knitr::kable("html") %>%
        kable_styling("striped", full_width = FALSE) %>%
        scroll_box(width = "1000px", height = "1000px")
    )})
    

  # Render Average Taxa Plot
  output$taxplot <- renderPlot({
    
    filtered_data <- ave_filtered_for_taxplot() # Get the reactive filtered data
    taxa_input_sym <- sym(input$Taxa_input) # Convert to symbol

    ggplot(filtered_data, aes(x = condition, y = Average_Abundance, fill = !!sym(input$Taxa_input))) +
      geom_bar(stat = "identity", position = "stack") + # Adjust position to "stack" for stacked bars
      scale_fill_viridis(option="magma", discrete=TRUE) +
      labs(fill = input$Taxa_input) # Label the legend with the input Taxa group
  })
  
  #Render Taxa Plot by Sample
  output$taxplot_sample <- renderPlot({
    
    filtered_data <- filtered_for_taxplot() # Get the reactive filtered data
    
    ggplot(filtered_data, aes(x = Sample, y = Abundance, fill = !!sym(input$Taxa_input))) +
      geom_bar(stat = "identity") + # Adjust position to "stack" for stacked bars
      scale_fill_viridis(option="magma", discrete=TRUE) +
      facet_wrap(vars(condition), scales = "free_x")+
      labs(fill = input$Taxa_input) + # Label the legend with the input Taxa group
      theme(axis.title.x=element_blank(), 
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  })

  
  #Render PCA
  output$PCA <- renderPlot({
      resolution <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
      plots <- lapply(resolution, function(rank) {
        plot_pca(phyloseq_obj = phy_obj, 
                 rank_transformation = rank, 
                 variable = "condition", 
                 colors_list = colors_all)})
      wrap_plots(plots, ncol = 2) & theme(legend.position = "bottom") })


  #Render PCoA 
  output$PCoA <- renderPlot ({
    resolution <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

    plots <- lapply(resolution, function(rank) {
      plot_PCoA(phyloseq_obj = phy_obj,
                rank_transformation = rank,
                trans_type = "identity",        
                dist_cal_type = input$PCoA_distance,   
                ord_calc_method = "NMDS",
                variable = "condition", 
                colors_list = colors_all)})
    
    wrap_plots(plots, ncol = 2) & theme(legend.position = "bottom") 
  })
  
  
  #Render Alpha Diversity
  output$Alphadiv <- renderPlot({

    comps <- list(c("AN", "Tumor"))

    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
    
    
    # Sample plots
    p_Shannon <- plot_richness(phy_obj, x= "condition", measures="Shannon", color = "condition") +
      geom_boxplot(alpha=0.6) + 
      theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
      stat_compare_means(method = "wilcox.test", comparisons = comps, label = "p.signif", symnum.args = symnum.args) +
      scale_color_manual(values = colors_all)
  
    p_Chao1 <- plot_richness(phy_obj, x= "condition", measures="Chao1", color = "condition") +
      geom_boxplot(alpha=0.6) + 
      theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
      stat_compare_means(method = "wilcox.test", comparisons = comps, label = "p.signif", symnum.args = symnum.args) +
      scale_color_manual(values = colors_all)
    
    p_Observed <- plot_richness(phy_obj, x= "condition", measures="Observed", color = "condition") +
      geom_boxplot(alpha=0.6) + 
      theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
      stat_compare_means(method = "wilcox.test", comparisons = comps, label = "p.signif", symnum.args = symnum.args) +
      scale_color_manual(values = colors_all)
    
    ggarrange(p_Shannon, p_Observed, p_Chao1, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
  })

  
  #Render lefse results plot
  output$lefse <- renderPlot({
    LEFSe_results <- Lefse_results()
    results <- as.data.frame(marker_table(LEFSe_results))
    number_of_results <- nrow(results)
    
    plot_ef_dot(LEFSe_results) +
      ggtitle("LefSE Abundance Graph") +
      theme(plot.title = element_text(face = "bold", vjust = 1, lineheight = 1)) +
      scale_color_manual(values = colors_all) 
    
  })
  
  
  #Render lefse result table
  output$lefse_results <- renderDataTable({
    LEFSe_results <- Lefse_results()
    results <- as.data.frame(marker_table(LEFSe_results))
    # Use DT::datatable for more interactive tables
    DT::datatable(results)
  })
  
  
}
##################################################################################################################################
#Run app
##################################################################################################################################
# Run the application - THIS IS THE LAST LINE OF CODE
shinyApp(ui, server)


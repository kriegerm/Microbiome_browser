
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
if(!require(fantaxtic)) { devtools::install_github("gmteunisse/fantaxtic"); library(fantaxtic) }
if(!require(fantaxtic)) { devtools::install_github("gmteunisse/ggnested"); library(fantaxtic) }
if(!require(ggpubr)) { devtools::install_github("kassambara/ggpubr"); library(ggpubr) }


# Set a global theme for ggplot2 plots
theme_set(theme_bw() + theme(text = element_text(size = 18))) # Adjust the size as needed

#Load in data
# Load the saved Rds file - this comes from /CEDAR/Microbiome/oral_microbiome_exploratory/CEDAR_microbiome_exploratory on my local machine. This is the OSCC data from PRJNA822685, filtered for only matched samples with AN and Tumor reads.
loaded_objects <- readRDS("00data.rds")
# Assign back the objects to the environment, if needed
list2env(loaded_objects, .GlobalEnv)

#Define global color
colors_all= c("Tumor" = "#FF495C", "AN"="#083D77" )

##################################################################################################################################
#Define UI
##################################################################################################################################
# Define UI for application (this is just an HTML wrapper)
ui <- fluidPage(

  titlePanel("Microbiome Browser"),
  sidebarLayout(
    
    #The input goes into the sidebarpanel
    sidebarPanel( width=2, #set the width of the sidebar
                  
      selectInput("ps_object", "Chose filtered singles (fs) or no",
          choices = c("phy_obj_oscc", "phy_obj_oscc_fs")),
      
      selectInput("Taxa_input", "Taxanomic rank",
                  choices = c("Phylum", "Genus", "Species")),
      
      
      selectInput("groupby", "Group by",
                  choices = c("Type", "condition")),
      
      numericInput("Num_top_taxa", "Taxa plots: Number of top taxa to display", value = 10, min = 0, max = 25),
      
      uiOutput("dynamicComparisons"), # Placeholder for the dynamic input
      
      numericInput("ldacutoff", "Lefse: LDA Cutoff", value = 2, min = 0, max = 10),
      
      selectInput("Lefse_taxa", "Lefse: Taxa",
                  choices = c("all", "Species", "Genus", "Phylum")),
      
      selectInput("Lefsenorm", "Lefse: Normalization (default: CPM)",
                  choices = c("CPM", "CLR", "RLE", "CSS", "none"))
    ),
    
    mainPanel(
      tabsetPanel(
        
      #Tab outputs
      tabPanel("Average Taxa Plot", plotOutput("taxplot")),
      tabPanel("Taxa Plot by Sample", plotOutput("taxplot_sample")),
      tabPanel("Average Taxa Table", tableOutput("results")),
      tabPanel("Beta Diversity", plotOutput("PCA")),
      tabPanel("Alpha Diversity", plotOutput("Alphadiv")),
      tabPanel("LeFSe", plotOutput('lefse'), dataTableOutput("lefse_results"))
               )
            )      
          )
        )


##################################################################################################################################
#Server
##################################################################################################################################

# Define server logic for creating output
server <- function(input, output) {
  
  #This code makes it possible to define comparisons between different subcategories or conditions
  
    output$dynamicComparisons <- renderUI({
      
    # Define choices based on the selection in 'groupby'
      # Can add more here, but for this example we just have An-Tumor
    choices <-list("AN vs Tumor" = "AN-Tumor")
    
    # Create the 'Comparisons' select input with dynamic choices
    selectInput("comparisons", "Lefse: Comparisons", choices = choices)
    })
    

##################################################################################################################################
#Define reactive variables
##################################################################################################################################
    # Get the name of the ps_obj you want to use
    ps_obj <- reactive({
      req(input$ps_object)  # Ensure input is not NULL
      get(input$ps_object, envir = .GlobalEnv)  # Get from the global environment
    })
    
    # Create a "filtered" ps object that is agglomerated at the taxa rank you want
    ps_filtered <- reactive ({
      ps_core <- ps_obj()
      ps_core %>%
        tax_transform(rank = input$Taxa_input, trans = "identity")
      })
  
  
  #The averaged taxa table that will be displayed in the "Average Taxa Plot" tab.
  ave_filtered_for_taxplot <- reactive({
    ps_core <- ps_obj()
    taxa_input_sym <- sym(input$Taxa_input) # Convert to symbol
    Group_by_input <- sym(input$groupby)
    
    ps_core %>%
      tax_glom(., taxrank = input$Taxa_input, NArm = TRUE) %>%
      transform_sample_counts(., function(x) (x)/sum(x)) %>%
      get_top_taxa(., input$Num_top_taxa, relative = TRUE, discard_other = TRUE) %>%
      psmelt() %>%
      group_by(!!Group_by_input, !!taxa_input_sym) %>%
      summarise(Average_Abundance = mean(Abundance))
  })
  
  #The averaged taxa table that will be displayed in the "Average Taxa Table" tab (the same as the one above, but with pivot_wider)
  ave_filtered_for_display <- reactive({
    ps_core <- ps_obj()
    taxa_input_sym <- sym(input$Taxa_input) # Convert to symbol
    Group_by_input <- sym(input$groupby)
    
    ps_core %>%
      tax_glom(., taxrank = input$Taxa_input, NArm = TRUE) %>%
      transform_sample_counts(., function(x) (x)/sum(x)) %>%
      psmelt() %>%
      group_by(!!Group_by_input, !!taxa_input_sym) %>%
      summarise(Average_Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = !!Group_by_input, values_from = Average_Abundance, names_prefix = "")
  })
  
  #This is the most basic taxa table generated - it is by sample, not averaged. Displayed in the "Taxa Plot by Sample" tab.
  filtered_for_taxplot <- reactive({
    ps_core <- ps_obj()
    taxa_input_sym <- sym(input$Taxa_input) # Convert to symbol
    Group_by_input <- sym(input$groupby)
    
    ps_core %>%
      tax_glom(., taxrank = input$Taxa_input, NArm = TRUE) %>%
      transform_sample_counts(., function(x) (x)/sum(x)) %>%
      get_top_taxa(., input$Num_top_taxa, relative = TRUE, discard_other = TRUE) %>%
      psmelt() %>%
      group_by(!!Group_by_input, !!taxa_input_sym) 
  })
  

  #Create a Lefse results object
  Lefse_results <- reactive({
    comparisons_vector <- strsplit(input$comparisons, "-")[[1]]
    Group_by_input <- sym(input$groupby)
    ps_core <- ps_obj()
    
    phylo_obj_fLefSe <- ps_core %>%
      ps_filter(., !!Group_by_input %in% comparisons_vector) 
    
    # Assuming run_lefse, marker_table, and plot_ef_dot are defined elsewhere and work as intended
    run_lefse(phylo_obj_fLefSe,
                               wilcoxon_cutoff = .05,
                               taxa_rank = input$Lefse_taxa,
                               norm = input$Lefsenorm, 
                               group = input$groupby,
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
    Group_by_input <- sym(input$groupby)
    
    ggplot(filtered_data, aes(x = !!sym(input$groupby), y = Average_Abundance, fill = !!sym(input$Taxa_input))) +
      geom_bar(stat = "identity", position = "stack") + # Adjust position to "stack" for stacked bars
      scale_fill_viridis(option="magma", discrete=TRUE) +
      labs(fill = input$Taxa_input) # Label the legend with the input Taxa group
  })
  
  #Render Taxa Plot by Sample
  output$taxplot_sample <- renderPlot({
    
    filtered_data <- filtered_for_taxplot() # Get the reactive filtered data
    Group_by_input <- sym(input$groupby)
    
    ggplot(filtered_data, aes(x = Sample, y = Abundance, fill = !!sym(input$Taxa_input))) +
      geom_bar(stat = "identity") + # Adjust position to "stack" for stacked bars
      scale_fill_viridis(option="magma", discrete=TRUE) +
      facet_wrap(vars(!!Group_by_input), scales = "free_x")+
      labs(fill = input$Taxa_input) # Label the legend with the input Taxa group
  })

  #Render PCA
  output$PCA <- renderPlot({
    
    Group_by_input <- sym(input$groupby)
    
    ps_filtered() %>% microViz::tax_fix() %>%
      ord_calc("PCA") %>% 
      ord_plot(axes = c(1, 2), fill = input$groupby, shape =input$groupby, alpha = 0.8, size = 2) + 
      ggtitle(label=paste0("PCA Plot")) +
      scale_shape_girafe_filled() +
      scale_fill_manual(values = colors_all) +
      scale_color_manual(values = colors_all) +
      ggplot2::stat_ellipse(aes(color =!!sym(input$groupby))) +
      theme(plot.title = element_text(face="bold", size=16, hjust=.5))+
      theme(plot.subtitle = element_text(face="italic", size=4, hjust=.5))+
      theme(legend.text=element_text(size=14))+ 
      theme(legend.title = element_text(size=14))+
      theme(axis.text=element_text(size=14, vjust = 0.5, hjust=1))
  })
  

  
  
  #Render Alpha Diversity
  output$Alphadiv <- renderPlot({

    comps <- list(c("AN", "Tumor"))

    
    Group_by_input <- sym(input$groupby)
    ps_core <- ps_obj()
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
    
    
    # Sample plots
    p_Shannon <- plot_richness(ps_core, x= input$groupby, measures="Shannon", color = input$groupby) +
      geom_boxplot(alpha=0.6) + 
      theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
      stat_compare_means(method = "wilcox.test", comparisons = comps, label = "p.signif", symnum.args = symnum.args) +
      scale_color_manual(values = colors_all)
  
    p_Chao1 <- plot_richness(ps_core, x= input$groupby, measures="Chao1", color = input$groupby) +
      geom_boxplot(alpha=0.6) + 
      theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
      stat_compare_means(method = "wilcox.test", comparisons = comps, label = "p.signif", symnum.args = symnum.args) +
      scale_color_manual(values = colors_all)
    
    p_Observed <- plot_richness(ps_core, x= input$groupby, measures="Observed", color = input$groupby) +
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


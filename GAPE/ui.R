#library(plyr)
library(dplyr)
library(DT)
library(shiny)
library(shinyjs)


appCSS <- "
#loading-content {
position: absolute;
background: #000000;
opacity: 0.8;
z-index: 100;
left: 0;
right: 0;
height: 100%;
text-align: center;
color: #FFFFFF;
}
"
##read in objects essential to Pseudomonas

PaLocusTag<-readRDS("Pa_LocusTagIntergenicSpotID_nocontrols_unzip.rds")
Pa_GSEnum<-readRDS("Pa_GSEnum_unzip.rds")

PaPathsGenes<-readRDS("Pa_KEGGPathGenes_unzip.rds")


####read in objects essential to Escherichia 
EcLocusTag<-readRDS("Ec_LocusTagIntergenicSpotID_nocontrols_unzip.rds")
Ec_GSEnum<-readRDS("Ec_GSEnum_unzip.rds")
EcPathsGenes<-readRDS("Ec_KEGGPathsGenes_unzip.rds")



fluidPage(
  useShinyjs(),
  inlineCSS(appCSS),
  # Loading message
  div(
    id = "loading-content",
    h1("Loading...")
  ),
  includeCSS('customcss.css'),
  titlePanel(h2("GAPE"),windowTitle = "PaPE"),
  HTML("<h3><strong>- <span style=\"text-decoration: underline;\">
       <strong>GA</strong></span>UGE-annotated <em><span style=\"text-decoration: underline;\"> <strong>P</strong></span>seudomonas <strong>a</strong>eruginosa &amp;&nbsp;<span style=\"text-decoration: underline;\">E</span>scherichia coli&nbsp;</em> transcriptomic ANOVA compendia for reanalysis</strong></h3>"),
  hr(),
hidden(
  div(
      id = "app-content",
  fixedRow(
    column(12, HTML("<h4 style=\"text-align: center;\"><span style=\"text-decoration: underline;\"><span style=\"color: #ff6600; text-decoration: underline;\"><strong>Step 1: select a species you want to analyze</strong></span></span></h4>"),
           column(4,
                  selectInput("species",NULL,multiple = FALSE, choices = c("","Pseudomonas aeruginosa","Escherichia coli"),selected = NULL,width = "100%")
           ),
           column(4,
                  actionButton("speciesButton","Submit",width = "100%")
                  
           ),
           column(4,
                  actionButton("reset_species","Reset Species",width = "100%")
                  
           )
    ),
    column(12,htmlOutput("chosen_compendium"))
  ),
  # when the species is chose and the submit button is hit, show the sidebarLayout
  conditionalPanel( condition =  "output.species=='Pseudomonas aeruginosa'|output.species=='Escherichia coli'",
      sidebarLayout(
            sidebarPanel(
                        width = 3,
                        tags$head(tags$style(type='text/css',".nav-tabs {font-size: 12px;font-weight: bold} ")),
                        wellPanel(
                          actionButton("instructionButton","Show User Guide",width = "100%"),
                          actionButton("browseButton","Browse Datasets",width = "100%")
                        ),
                        tabsetPanel(
                          tabPanel("Single Gene",
                                   wellPanel(
                                     conditionalPanel( condition =  "output.species=='Pseudomonas aeruginosa'",
                                                       
                                     selectInput("singleLocus","Search & select a PAO1 locus or an intergenic region on GPL84 (e.g. PA0001, Intergenic region...)",multiple = FALSE, choices = c("",PaLocusTag),selected = NULL)
                                   ),
                                   conditionalPanel( condition =  "output.species=='Escherichia coli'",
                                     selectInput("singleLocus_Ec","Search & select an E. coli locus or an intergenic region on GPL199 & GPL3154 (e.g. b2836, z0312, c4728, Ecs1315, IG_7960...)",multiple = FALSE, choices = c("",EcLocusTag),selected = NULL)
                                   ),               
                                     hr(),
                                     
                                     actionButton("singleLocusButton","Submit",width = "84%"),
                                     actionButton("reset_singleLocus","Reset",width = "84%"),
                                     
                                     hr(),
                                     
                                     uiOutput("download_singleLocus")
                                     
                                     
                                   ))
                          ,
                          tabPanel("Enrichment",
                                   tabsetPanel(
                                     tabPanel("KEGG",
                                              wellPanel(
                                                conditionalPanel( condition =  "output.species=='Pseudomonas aeruginosa'",# && input.speciesButton==true",
                                                    selectInput('KEGG', "Search & select a KEGG pathway", multiple = FALSE, choices = c("",names(PaPathsGenes)),selected = NULL)),
                                                conditionalPanel( condition =  "output.species=='Escherichia coli'",# && input.speciesButton==true",
                                                    selectInput('KEGG_Ec', "Search & select an E. coli KEGG pathway", multiple = FALSE, choices = c("",names(EcPathsGenes)),selected = NULL)),
                                                hr(),
                                                HTML("<p style=\"text-align: center;\"><span style=\"font-size: 15pt;\">or</span></p>"),
                                                conditionalPanel( condition =  "output.species=='Pseudomonas aeruginosa'",# && input.speciesButton==true",
                                                  selectInput('KEGG_selectedGSE', "Search & select a dataset", multiple = FALSE, choices = c("",Pa_GSEnum),selected = NULL)),
                                                conditionalPanel( condition =  "output.species=='Escherichia coli'",# && input.speciesButton==true",
                                                  selectInput('KEGG_selectedGSE_Ec', "Search & select a dataset", multiple = FALSE, choices = c("",Ec_GSEnum),selected = NULL)),
                                                hr(),
                                                actionButton("KEGGButton","Submit",width = "84%"),
                                                
                                                actionButton("reset_KEGG","Reset",width = "84%"),
                                                
                                                hr(),
                                                
                                                uiOutput("download_KEGG"),
                                                uiOutput("download_KEGG_selectedGSE")
                                                
                                                
                                              )),
                                      
                                     tabPanel("Custom List",        
                                            wellPanel(
                                              conditionalPanel( condition =  "output.species=='Pseudomonas aeruginosa'",# && input.speciesButton==true", 
                                                textAreaInput(inputId = "customlist",
                                                              label = "PAO1 locus IDs:",
                                                              value = NULL,
                                                              placeholder = "Paste a list of locus ID here.\nID separated by a new line.",
                                                              rows=20,
                                                              height = "100%",
                                                              width="100%"
                                                              )),
                                              conditionalPanel( condition =  "output.species=='Escherichia coli'",# && input.speciesButton==true", 
                                                                textAreaInput(inputId = "customlist_Ec",
                                                                              label = "E.coli locus IDs:",
                                                                              value = NULL,
                                                                              placeholder = "Paste a list of locus ID here.\nID separated by a new line.",
                                                                              rows=20,
                                                                              height = "100%",
                                                                              width="100%"
                                                                )),
                                                hr(),
                                                
                                                actionButton("customButton","Submit",width = "84%"),
                                                
                                                actionButton("reset_customlist","Reset",width = "84%"),
                                                
                                                hr(),
                                                
                                                uiOutput("download_custom")
                                                
                                              ))#close tabPanel custom list
                                     
                                     
                                   )# close tabsetpanel under enrichemnt 
                          )# close tabpanel of enrichment 
                        )# close tabpanel above enrichment
                        
                        
                      ),#close sidebarPanel
                      mainPanel(
                        width=9,
                        htmlOutput("instruction"),
                        htmlOutput("error"),
                        DT::dataTableOutput(outputId = "table")
                      )#close mainPanel
                    )#close sidebarLayout
  )#close conditional Panel
    )
  )
)#close fluidpage



function(input,output){
  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content", anim = TRUE, animType = "fade")    
  show("app-content")
  
  d <- reactiveValues(titleSummaries=NULL,ANOVAList=NULL,KEGGPathID=NULL,KEGGPathGenes=NULL)
  v <- reactiveValues(species=NULL,instruction=1,compendium=NULL,singleLocus=NULL
                      ,error=NULL,KEGG = NULL,KEGG_GSE=NULL,customlist=NULL)
  
##detect which species that user wants to analyze  
  observeEvent(input$speciesButton,{req(input$species!="")
    v$species<-input$species
    
    output$species<-renderText({v$species},quoted = TRUE)
    outputOptions(output, "species", suspendWhenHidden = FALSE)
    })
  
  #For Pseudomonas
  observeEvent(input$speciesButton,{
    req(v$species=="Pseudomonas aeruginosa")
      {
      d$titleSummaries<-readRDS("Pa_titleSummaries_unzip.rds")
      d$ANOVAList<-readRDS("Pa_GPL84_refine_ANOVA_List_unzip.rds")
      d$KEGGPathID<-readRDS("Pa_KEGGPathID_unzip.rds")
      d$KEGGPathGenes<-readRDS("Pa_KEGGPathGenes_unzip.rds")
    }
  })
    #For e. coli 
  observeEvent(input$speciesButton,{  
  req(v$species=="Escherichia coli")
    {
      d$titleSummaries<-readRDS("Ec_titleSummaries_unzip.rds")
      d$ANOVAList<-readRDS("Ecoli_twoPlatform_ANOVA_List_withGPL_unzip.rds")
      d$KEGGPathID<-readRDS("Ec_KEGGPathID_unzip.rds")
      d$KEGGPathGenes<-readRDS("Ec_KEGGPathsGenes_unzip.rds")
    }
    })  
##reset species when the species reset button is hit
  observeEvent(input$reset_species,{
    #shinyjs::js$refresh()
    reset("species")
    output$species<-NULL
    outputOptions(output, "species", suspendWhenHidden = FALSE)
    v$species<-NULL
    })

## for showing instruction  
  observeEvent(input$instructionButton,{
    v$instruction<-1
    v$compendium<-NULL
    v$singleLocus<-NULL
    v$error<-NULL
    v$KEGG<-NULL
    v$KEGG_GSE<-NULL
    v$customlist<-NULL
  })
## for browsing the compendium
  observeEvent(input$browseButton,{
    v$instruction<-NULL
    v$compendium<-{d$titleSummaries}
    v$singleLocus<-NULL
    v$error<-NULL
    v$KEGG<-NULL
    v$KEGG_GSE<-NULL
    v$customlist<-NULL
    })

  #for one single locus tag screening
  observeEvent(input$singleLocusButton,{
    v$instruction<-NULL
    if(input$singleLocus!=""){req(input$singleLocus)
                              input_singleLocus<-input$singleLocus}
    if(input$singleLocus_Ec!=""){req(input$singleLocus_Ec)
                              input_singleLocus<-input$singleLocus_Ec}
    v$singleLocus<-{
      nameofList<-names(d$ANOVAList)
      locustag<-input_singleLocus
      withProgress(message="Calculating differential expression",value = 0,{
        FishPv<-lapply(nameofList,function(n){
          print(n)
          num<-which(nameofList ==n)
          print(num)
          incProgress(1/length(nameofList), detail = paste("Doing dataset No.", num))
          y<-d$ANOVAList[[n]]
          #y<-apply(y,2,as.character)
          #y<-as.data.frame(y,stringsAsFactors = FALSE)
          finalMatrix<-y
          finalMatrix$LocusTagANDspotID<-paste(finalMatrix$LocusTag,finalMatrix$IntergenicSpotID,sep = "")
          
          if(sum(finalMatrix$LocusTagANDspotID==locustag)==1){
            output<-finalMatrix[finalMatrix$LocusTagANDspotID==locustag,]
            output<-c(n,output)}
          else{output<-NULL}
        })
        FishPv[sapply(FishPv, is.null)]<-NULL
        FishPv<-lapply(FishPv,function(x){
          y<-as.data.frame(x,stringsAsFactors = FALSE)
          y<-y[,c(1,2,3,4,10,6:9)]
          colnames(y)<-c("GSE","ProbeID","GeneTitle","Symbol","LocusTag\\SpotID"
                         ,"IntergenicSpotID","ANOVAPvalue","Log2FC","FDR" )
          y
        })
        bigdataframe<-ldply(FishPv)
        bigdataframe$FDR<-as.numeric(bigdataframe$FDR)
        bigdataframe$ANOVAPvalue<-as.numeric(bigdataframe$ANOVAPvalue)
        bigdataframe$Log2FC<-as.numeric(bigdataframe$Log2FC)
      })
      if(sum(grepl("GPL",bigdataframe$GSE)!=0)){bigdataframe$GSE<-gsub("_GPL3154","",bigdataframe$GSE)}
      bigdataframe<-merge(bigdataframe,d$titleSummaries,by="GSE")
      bigdataframe<-bigdataframe[order(bigdataframe$FDR,decreasing = FALSE),]
      
      
      bigdataframe}  
  })
  
  ### error message for selecting two inputs for KEGG enrichment analysis
  observeEvent(input$KEGGButton,{
    v$instruction<-NULL
     v$error<-{req(input$KEGG,input$KEGG_selectedGSE)
      1
     }
    })  
  observeEvent(input$KEGGButton,{
    v$instruction<-NULL
    v$error<-{req(input$KEGG_Ec,input$KEGG_selectedGSE_Ec)
      1
    }
  }) 
  
  ###for selecting one KEGG pathway for analysis 
  observeEvent(input$KEGGButton,{
    v$instruction<-NULL
    req(input$KEGG_selectedGSE=="",input$KEGG_selectedGSE_Ec=="")
    if(input$KEGG!="" & input$KEGG_selectedGSE==""){req(input$KEGG,input$KEGG_selectedGSE=="")
                                                        input_KEGG<-input$KEGG}
    if(input$KEGG_Ec!="" & input$KEGG_selectedGSE_Ec==""){req(input$KEGG_Ec,input$KEGG_selectedGSE_Ec=="")
                                                        input_KEGG<-input$KEGG_Ec}
    v$KEGG<-{
      
      nameofList<-names(d$ANOVAList)
      pathwayToTest<-as.character(unlist(d$KEGGPathGenes[[input_KEGG]]))
      withProgress(message="Calculating Gene set enrichment",value = 0,{
        FishPv<-lapply(nameofList,function(n){
          
          print(n)
          num<-which(nameofList ==n)
          print(num)
          incProgress(1/length(nameofList), detail = paste("Doing dataset No.", num))
          y<-d$ANOVAList[[n]]
          finalMatrix<-y
          str(finalMatrix)
          #make a new column that is a combination of LocusTag and IntergenicSpotID
          finalMatrix$LocusTagANDspotID<-paste(finalMatrix$LocusTag,finalMatrix$IntergenicSpotID,sep = "")
          
          #now every row is considered as a unique gene
          Genes<-(finalMatrix$LocusTagANDspotID)
          #get the significat differential expressed genes
          DE_matrix<-finalMatrix[finalMatrix$FDR<=0.05,]
          SigGenes<-(DE_matrix[,"LocusTagANDspotID"])
          Regulated <- SigGenes
          
          NotRegulated <- setdiff(Genes, SigGenes)
          
          OnPath <- intersect(Genes,pathwayToTest)
          
          NotOnPath <- setdiff(Genes, OnPath)
          
          
          vals <- c(length(intersect(OnPath,Regulated)), 
                    length(intersect(OnPath,NotRegulated)),  
                    length(intersect(NotOnPath,Regulated)),
                    length(intersect(NotOnPath,NotRegulated)))
          PathMatrix <- matrix(vals, nrow = 2, dimnames = list(
            Regulated=c("Yes", "No"),  
            OnPath = c("Yes", "No")
          ))
          PathMatrix
          p<-fisher.test(PathMatrix, alternative = "greater")$p.value
          #convert Onpath locus tag to symbol
          sym<-finalMatrix$Symbol[finalMatrix$LocusTag %in% intersect(OnPath,Regulated)]
          print(p)
          c(p,length(intersect(OnPath,Regulated)),paste(sym,collapse = ";"),length(Regulated))
        })
      })
      bigdataframe<-ldply(FishPv)
      colnames(bigdataframe)<-c("FDR","No. of DE genes on the pathway","DE genes on the pathway","No. of DE genes in the study")
      GSE<-grep("GSE|SRP",unlist(strsplit(nameofList,"_")),value=TRUE)
      bigdataframe<-cbind(GSE,bigdataframe)
      bigdataframe$FDR<-p.adjust(bigdataframe$FDR,method = "fdr")
      bigdataframe<-merge(bigdataframe,d$titleSummaries,by="GSE")
      bigdataframe<-bigdataframe[order(bigdataframe$FDR,decreasing = FALSE),]
      bigdataframe}
      })  
  
  ###for selecting one GSE to perform all KEGG pathways enrichment 
  observeEvent(input$KEGGButton,{
    v$instruction<-NULL
    req(input$KEGG=="",input$KEGG_Ec=="")
    if(input$KEGG=="" & input$KEGG_selectedGSE!=""){req(input$KEGG_selectedGSE,input$KEGG=="")
      input_KEGG_selectedGSE<-input$KEGG_selectedGSE}
    if(input$KEGG_Ec=="" & input$KEGG_selectedGSE_Ec!=""){req(input$KEGG_selectedGSE_Ec,input$KEGG_Ec=="")
      input_KEGG_selectedGSE<-input$KEGG_selectedGSE_Ec}
    v$KEGG_GSE<-{
     
      withProgress(message="Calculating Gene set enrichment",value = 0,{
        FishPv<-lapply(d$KEGGPathGenes,function(n){
          
          print(n)
          
          incProgress(1/length(names(d$KEGGPathGenes)))
          a<-which(names(d$ANOVAList)==input_KEGG_selectedGSE)
          y<-d$ANOVAList[[a]]
          
          #y<-apply(y,2,as.character)
          #y<-as.data.frame(y,stringsAsFactors = FALSE)
          finalMatrix<-y
          #make a new column that is a combination of LocusTag and IntergenicSpotID
          finalMatrix$LocusTagANDspotID<-paste(finalMatrix$LocusTag,finalMatrix$IntergenicSpotID,sep = "")
          #now every row is considered as a unique gene
          Genes<-(finalMatrix$LocusTagANDspotID)
          #get the significat differential expressed genes
          DE_matrix<-finalMatrix[finalMatrix$FDR<=0.05,]
          
          SigGenes<-(DE_matrix[,"LocusTagANDspotID"])
          Regulated <- SigGenes
          
          NotRegulated <- setdiff(Genes, SigGenes)
          
          OnPath <- intersect(Genes,n)
          
          NotOnPath <- setdiff(Genes, OnPath)
          
          
          vals <- c(length(intersect(OnPath,Regulated)), 
                    length(intersect(OnPath,NotRegulated)),  
                    length(intersect(NotOnPath,Regulated)),
                    length(intersect(NotOnPath,NotRegulated)))
          PathMatrix <- matrix(vals, nrow = 2, dimnames = list(
            Regulated=c("Yes", "No"),  
            OnPath = c("Yes", "No")
          ))
          PathMatrix
          p<-fisher.test(PathMatrix, alternative = "greater")$p.value
          #convert Onpath ENSEMBLID to symbol
          sym<-finalMatrix$Symbol[finalMatrix$LocusTag %in% intersect(OnPath,Regulated)]
          print(p)
          c(p,length(intersect(OnPath,Regulated)),paste(sym,collapse = ";"))
        })
      })
      FishPv
      bigdataframe<-ldply(FishPv)
      bigdataframe
      colnames(bigdataframe)<-c("KEGG pathway","FDR","No. of DE genes on the pathway","DE genes on the pathway")
      head(bigdataframe)
      bigdataframe$FDR<-p.adjust(bigdataframe$FDR,method = "fdr")
      #bigdataframe$FDR<-round(bigdataframe$FDR,digits = 5)
      bigdataframe$"KEGG"<-d$KEGGPathID
      bigdataframe<-bigdataframe[,c(5,1:4)]
      bigdataframe<-bigdataframe[order(bigdataframe$FDR,decreasing = FALSE),]
      bigdataframe
    }
  })
  
  observeEvent(input$customButton,{
    v$instruction<-NULL
    if(input$customlist!=""){req(input$customlist)
      input_customlist<-input$customlist}
    if(input$customlist_Ec!=""){req(input$customlist_Ec)
      input_customlist<-input$customlist_Ec}
    v$customlist<-{
      nameofList<-names(d$ANOVAList)
      pathwayToTest<-as.character(unlist(strsplit(input_customlist,"\n")))
      withProgress(message="Calculating Gene set enrichment",value = 0,{
        FishPv<-lapply(nameofList,function(n){
          
          print(n)
          num<-which(nameofList ==n)
          print(num)
          incProgress(1/length(nameofList), detail = paste("Doing dataset No.", num))
          y<-d$ANOVAList[[n]]
          
          #y<-apply(y,2,as.character)
          #y<-as.data.frame(y,stringsAsFactors = FALSE)
          finalMatrix<-y
          #make a new column that is a combination of LocusTag and IntergenicSpotID
          finalMatrix$LocusTagANDspotID<-paste(finalMatrix$LocusTag,finalMatrix$IntergenicSpotID,sep = "")
          #now every row is considered as a unique gene
          Genes<-(finalMatrix$LocusTagANDspotID)
          #get the significat differential expressed genes
          DE_matrix<-finalMatrix[finalMatrix$FDR<=0.05,]
          
          SigGenes<-(DE_matrix[,"LocusTagANDspotID"])
          Regulated <- SigGenes
          
          NotRegulated <- setdiff(Genes, SigGenes)
          
          OnPath <- intersect(Genes,pathwayToTest)
          
          NotOnPath <- setdiff(Genes, OnPath)
          
          
          vals <- c(length(intersect(OnPath,Regulated)), 
                    length(intersect(OnPath,NotRegulated)),  
                    length(intersect(NotOnPath,Regulated)),
                    length(intersect(NotOnPath,NotRegulated)))
          PathMatrix <- matrix(vals, nrow = 2, dimnames = list(
            Regulated=c("Yes", "No"),  
            OnPath = c("Yes", "No")
          ))
          PathMatrix
          p<-fisher.test(PathMatrix, alternative = "greater")$p.value
          #convert Onpath ENSEMBLID to symbol
          sym<-finalMatrix$Symbol[finalMatrix$LocusTag %in% intersect(OnPath,Regulated)]
          print(p)
          c(p,length(intersect(OnPath,Regulated)),paste(sym,collapse = ";"),length(Regulated))
        })
      })
      bigdataframe<-ldply(FishPv)
      colnames(bigdataframe)<-c("FDR","No. of DE genes on the List","DE genes on the List","No. of DE genes in the study")
      GSE<-grep("GSE|SRP",unlist(strsplit(nameofList,"_")),value=TRUE)
      bigdataframe<-cbind(GSE,bigdataframe)
      head(bigdataframe)
      bigdataframe$FDR<-p.adjust(bigdataframe$FDR,method = "fdr")
      #bigdataframe$FDR<-round(bigdataframe$FDR,digits = 5)
      bigdataframe<-merge(bigdataframe,d$titleSummaries,by="GSE")
      bigdataframe<-bigdataframe[order(bigdataframe$FDR,decreasing = FALSE),]
      bigdataframe}
  })
  
  
  

  observeEvent(input$reset_singleLocus, {
    v$instruction<-1
    v$compendium<-NULL
    v$singleLocus<-NULL
    v$error<-NULL
    v$KEGG<-NULL
    v$KEGG_GSE<-NULL
    v$customlist<-NULL
    reset("singleLocus")
    reset("singleLocus_Ec")
  })
  
  observeEvent(input$reset_KEGG, {
    v$instruction<-1
    v$compendium<-NULL
    v$singleLocus<-NULL
    v$error<-NULL
    v$KEGG<-NULL
    v$KEGG_GSE<-NULL
    v$customlist<-NULL
    reset("KEGG")
    reset("KEGG_selectedGSE")
    reset("KEGG_Ec")
    reset("KEGG_selectedGSE_Ec")
  })
 
  
   observeEvent(input$reset_customlist, {
     v$instruction<-1
     v$compendium<-NULL
     v$singleLocus<-NULL
     v$error<-NULL
     v$KEGG<-NULL
     v$KEGG_GSE<-NULL
     v$customlist<-NULL
    reset("customlist")
    reset("customlist_Ec")
  })

  
  
  output$table<-DT::renderDataTable({
    out<-NULL
    out_KEGG_GSE<-NULL
    if (!is.null(v$singleLocus)){out<-v$singleLocus[,c(1,5,4,3,8,9,10,11)]}
    else if (!is.null(v$compendium)){out<-v$compendium}
    else if (!is.null(v$KEGG)) {out<-v$KEGG[,c(-4,-5)]}
    else if(!is.null(v$KEGG_GSE)){out_KEGG_GSE<-v$KEGG_GSE}
    else if(!is.null(v$customlist)){out<-v$customlist[,c(-4,-5)]}
    else{out<-NULL
        out_KEGG_GSE<-NULL
        }
    if (!is.null(out)){   
    DT::datatable(data = out %>%
                  dplyr::mutate(URL = paste0("https://www.ncbi.nlm.nih.gov/search/all/?term=", GSE)) %>% ##paste two strings into an address 
                  dplyr::mutate(GSE = paste0("<a href='", URL,"'","target='_blank" ,"'>",GSE,"</a>"))%>% ##make this address into a hyperlink and show as A
                  dplyr::select(-c("URL")),
                  rownames = FALSE, escape = FALSE)
    }
    else if(!is.null(out_KEGG_GSE)){
      DT::datatable(data = out_KEGG_GSE %>%
                      dplyr::mutate(URL = paste0("https://www.genome.jp/dbget-bin/www_bget?pathway:", KEGG)) %>% ##paste two strings into an address 
                      dplyr::mutate(KEGG = paste0("<a href='", URL,"'","target='_blank" , "'>",KEGG,"</a>")) %>% ##make this address into a hyperlink and show as A
                      dplyr::select(-c("URL")),
                    rownames = FALSE, escape = FALSE)
      
    }
    
  })
  
  output$instruction<-renderUI({
    if(!is.null(v$instruction))
    {
      HTML("<h3 style=\"text-align: center;\"><strong>PaPE user guide&nbsp;</strong></h3>
<h4><strong>PaPE is a user-friendly Shiny App allows the user to investigate differential gene expression and gene set enrichment among auto-annotated <em>Pseudomonas aeruginosa</em> and <em>Escherichia coli</em> transcriptomic datasets.</strong></h4>
           <h4><strong>Step 1: Select and submit a species you want to analyze. Press \"Reset Species\" to restart the species selection&nbsp;</strong></h4>
           <ol>
           <li>Press \"Show User Guide\" when you need to revisit this guide.</strong></li>
           <li>Press \"Browse datasets\" to visit the information about the auto-annotated datasets. Use the top-right text input box to search globally to locate dataset of interest. Click the GSE number to visit the corresponding data website.&nbsp;&nbsp;</strong></li>
           </ol>
           <h4><strong>Step 2. Select the analysis function you need:</strong></h4>
           <ol>
           <li>For the differential expression analysis for the gene of interest across auto-annotated datasets, choose the tab &ldquo;Singel Gene&rdquo; on the left side panel. Go to step 3.1.</li>
           <li>For the gene set enrichment analysis, choose the tab &ldquo;Enrichment&rdquo; on the left side panel, another layer of tabs would appear.<br />
           <ul>
           <li>Choose the &ldquo;KEGG&rdquo; tab for performing KEGG pathways enrichment analysis. Go to step 3.2.</li>
           <li>Choose the &ldquo;Custom List&rdquo; tab for specifying your gene list. Go to step 3.3.</li>
           </ul>
           </li>
           </ol>
           <h4><strong>Step 3. &nbsp;Enter the gene/gene set/dataset of interest:</strong></h4>
           <p><strong>Note: All drop-down lists are searchable in this app</strong></p>
           <ol>
           <li>Select a single locus tag\\intergenic region of interest from the drop-down list and click \"Submit\". </li>
           <li>Under the &ldquo;KEGG&rdquo; tab, choose a KEGG pathway or a dataset of interest from the drop-down list and press &ldquo;Submit&rdquo;.</li>
           <li>Under the &ldquo;Custom List&rdquo; tab, enter a list of PA01 locus IDs separated by a new line and press &ldquo;Submit&rdquo;.</li>
           </ol>
           <h4><strong>Step 4.&nbsp; Inspect the result and press &ldquo;Download Report&rdquo; to download a CSV file containing detail result.</strong></h4>
           <h4><strong>Step 5. Press &ldquo;Reset&rdquo; or refresh your browser to initiate a new analysis.</strong>&nbsp;&nbsp;</h4>
           <h4>&nbsp;</h4>
           <h4><strong>*Contact Zhongyou Li (<a href=\"mailto:Zhongyou.li.gr@dartmouth.edu\">Zhongyou.li.gr@dartmouth.edu</a>) for any questions or suggestions</strong></h4>
           <h4><strong>*PaPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.</strong></h4>
           <h4><strong>*If this code is helpful to you, please cite our related publication:</strong></h4>"
        )
    }
      })
  output$error<-renderUI({
    if(!is.null(v$error))
    {
      HTML(
        "<h4 style=\"text-align: justify;\">&nbsp;</h4>
        <h4 style=\"text-align: justify;\">&nbsp;</h4>
        <h4 style=\"text-align: justify;\">&nbsp;</h4>
        <h4 style=\"text-align: justify;\"><strong>You have chosen both a pathway\\signature and a dataset for analysis, which is not allowed. 
        Please press the \"<span style=\"color: #ff0000;\">Reset</span>\" button and re-select <span style=\"text-decoration: underline;
        \">either a KEGG pathway\\signature or a dataset</span> to submit your analysis.</strong></h4>"
      )
    }
    })
  
  output$chosen_compendium<-renderUI({
    if(!is.null(v$species))
    {
      HTML(
        paste("<h3 style=\"text-align: center;\"><span style=\"color: #ff0000;\"><strong>You are analyzing <span style=\"text-decoration: underline;\"><em> ",v$species,"</em> </span> compendium</strong></span></h3>",sep="")
      )
    }
  })
  data<-reactive({
    if(!is.null(v$singleLocus)){v$singleLocus[,-6]}
    else if(!is.null(v$KEGG)){v$KEGG}
    else if(!is.null(v$KEGG_GSE)){v$KEGG_GSE}
    else if(!is.null(v$customlist)){v$customlist}
    })
  ##
  #Download file
  ## use downloadHandler to control the file name
  output$report <- downloadHandler(
    filename = "report.csv",
    content = function(file) { 
      write.csv(data(), file,row.names = FALSE) 
    },
    contentType = "csv")
  
  ## use renderUI to control whether show downloadButton or not
  output$download_singleLocus <- renderUI({
    if(!is.null(v$singleLocus)) {
      downloadButton('report', 'Download Report')
    }
  })
  output$download_custom <- renderUI({
    if(!is.null(v$customlist)) {
      downloadButton('report', 'Download Report')
    }
  })
  output$download_KEGG <- renderUI({
    if(!is.null(v$KEGG)) {
      downloadButton('report', 'Download Report')
    }
  })
  output$download_KEGG_selectedGSE <- renderUI({
    if(!is.null(v$KEGG_GSE)) {
      downloadButton('report', 'Download Report')
    }
  })
  }
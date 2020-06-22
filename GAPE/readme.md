# Source code for PaPE the Shiny app
PaPE is a user friendly Shiny app to perfome Parallel analysis of Pseudomonas aeruginosa & Escherichia coli auto-annotated 
transcriptomic compendia. The transcriptomic compendia was built using the automated sample group detection algorithm followed by ANOVA
analysis. Please refer the corresponding publication for further details.

The shiny app has deployed on the free shiny server, Shinyapps.io, at https://iamsoshiny.shinyapps.io/pape/ 
For users who prefer to run the app locally or modify the script, please follow the instruction below:

1. Clone this repository into your local directory
2. install required packages in your R environment:

```R
install.packages("plyr")
install.packages("dplyr")
install.packages("DT")
install.packages("shiny")
install.packages("shiyjs")
```
3. Launch PaPE the Shiny app
```R
library(shiny)
runApp(launch.browser = TRUE)
```
# Auto-annotated P. aeruginosa and E. coli ANOVA compendia
For user who wants to perform analysis other than PaPE provided, please download the auto-annotated P. aeruginosa and E. coli ANOVA compendia.
The Pa_GPL84_refine_ANOVA_List_unzip.rds and Ecoli_twoPlatform_ANOVA_List_withGPL_unzip.rds files store the R list objects of the compendia. Pa_titleSummaries_unzip.rds and Ec_titleSummaries_unzip.rds are the R objects for the study titles and summaries. 

Contact Zhongyou Li (Zhongyou.li.gr@dartmouth.edu) for any questions or suggestions

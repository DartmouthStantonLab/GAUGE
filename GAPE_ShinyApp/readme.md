# Source code for GAPE the Shiny app
GAPE is a user friendly Shiny app to perfome reanalysis of Pseudomonas aeruginosa & Escherichia coli GAUGE annotated 
transcriptomic compendia. The transcriptomic compendia was built using the automated sample group detection algorithm followed by ANOVA
analysis. Please refer the corresponding publication for further details.

The shiny app has deployed on the free shiny server, Shinyapps.io, at https://iamsoshiny.shinyapps.io/gape/ 
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
3. Launch GAPE the Shiny app
```R
library(shiny)
runApp(launch.browser = TRUE)
```

Contact Zhongyou Li (Zhongyou.li.gr@dartmouth.edu) for any questions or suggestions

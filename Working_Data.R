#Comments
install.packages("tidyverse")
library(tidyverse)
#in terminal copy data

#load data 
metadata <- read.table(file="Saanich.metadata.txt",header=TRUE,row.names=1,sep="\t", na.strings=c("NAN", "NA", "."))
OTU <- read.table(file="Saanich.OTU.txt",header=TRUE,row.names=1,sep="\t", na.strings=c("NAN", "NA", "."))

library(dplyr)

metadata %>%
  filter(Temperature_C <10) %>%
  select(matches("CH4|methane|nM"))

Convert all variables that are in nM to um output a talbe showing only the original nM and converted uM variables

```{r metadata} 
library(dplyr)

select(Saanich_metadata, starts_with("CH4"))
filter(Saanich_metadata, CH4_nM > 100 & Temperature_C < 10) %>% select(Depth_m)

```
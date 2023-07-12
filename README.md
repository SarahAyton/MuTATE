# MuTATE
Collection of functions to recursively partition datasets on binary splits across multiple targets.

The goal of MuTATE is to create automated, explainable, and comprehensive models across multiple dependent variables of interest.

It provides a collection of functions to recursively partition data on binary splits across multiple targets having different dependent variable types. This overcomes single-target limitations of traditional decision trees, without loosing model interpretability, and while handling continuous, categorical, count, and survival outcome variables. This suite of functions also includes a number of parameters for model customization, dependent variable weights, parameter tuning functions, and visualization tools.

# Installation

You can install the development version of MuTATE from [GitHub](https://github.com/) with:
      
``` {r}
# install.packages("devtools")
devtools::install_github("SarahAyton/MuTATE")
```

Once MuTATE is installed, it can be easily loaded:

```{r}
library(MuTATE)
```

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FSarahAyton%2FMuTATE%2F&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Page+Views&edge_flat=false)](https://hits.seeyoufarm.com)

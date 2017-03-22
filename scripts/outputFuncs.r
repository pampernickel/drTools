# PDobay for AG Bourquin, KISPI

#@...
# xtable output functions
#@...

library(xtable)

toTable <- function(cis){
  xtable(cis, label = 'tab:cis') -> cis.tab
  return(cis.tab)
}
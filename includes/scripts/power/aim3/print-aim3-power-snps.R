
```{r}
# print off sig snps from power-calcs-ind-assoc.Rmd
sig.snp = read.csv("~/dissertation/unc-dissertation-markdown/includes/scripts/power/aim3/sig-snp.csv",
         header=T)
kable(table.sig.snp)
```
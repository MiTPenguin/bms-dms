## Clonal Specificity in TYK2 Chunk 10

We have two clones, named `c1` and `c8`, containing the same TYK2 DMS library and assayed under identical conditions. We can handle this structure in two ways, either by analyzing the clones separately and comparing them post-hoc, or by combining them to increase our barcode coverage and overall power.

1. [Summary](#part1)
2. [Comparisons Between Clones](#part2)
3. [Impact of Clone Aggregation](#part3)

### Summary <a name="part1"></a>

-  The mutation effects quantified in both clones are strongly correlated, and there are no significant differences between clones when tested directly (lowest unadjusted p-value = 0.022, lowest adjusted p-value = 0.95).
-  More significant variants are identified in clone 1, likely due to the higher average number of barcodes per variant.
-  However, stop codons (and strong LoF in general) have a slightly larger effect size in clone 8 compared to clone 1.
-  When clones are combined, we observe a substantial reduction in standard errors and broadly consistent effect sizes or fold changes.
-  Combining clones increases the number of significant variants to 365, from 308 in clone 1 and 251 in clone 8 individually.


### Comparisons Between Clones <a name="part2"></a>

Since we already have summary statistics for each clone individually, the easiest starting point is to just compare the resulting fold changes and standard errors between the two clones at each variant. In the plot below, variants which are nominally significant or non-significant in both clones (unadjusted p-value <= 0.05) are colored red and purple, respectively. The other two colors indicate clones which are nominally significant in only one clone or the other but not both.


    
![png](OCNT-DMSLIB-1-Clonal-Var_files/OCNT-DMSLIB-1-Clonal-Var_7_0.png)
    


Clearly the effect sizes between clones correlate very well, although not completely linearly. Stop effects seem slightly stronger in clone 8 compared to clone 1, but the higher number of barcodes in clone 1 results in a larger number of clone1-only nominally significant variants compared to clone8-only.

We can also actually test whether the `clone1 - clone8` difference in variant effect is significantly different from zero for each variant. We apply a linear contrast analogous to previous comparisons and obtain summary statistics which we can compare with the originals. After correcting for multiple testing, none of the `clone1 - clone8` tests is close to significant (lowest unadjusted p-value = 0.026, lowest adjusted p-value = 0.95). In particular, the p-value histogram summarizes the overall result nicely and "looks", statistically speaking, like a complete absence of signal in the clone difference:


    
![png](OCNT-DMSLIB-1-Clonal-Var_files/OCNT-DMSLIB-1-Clonal-Var_11_0.png)
    


So, while we can gain something from combining clones and their barcodes, there is no statistically significant difference in variant effect profiles between them.

### Impact of Clone Aggregation <a name="part3"></a>

While we will typically have a single clone per chunk, the existence of two allows us to ask how our discovery set would change if the clones were combined. If the clones are highly consistent (i.e. capturing the same effects), we should observe a decrease in standard error leading to higher power and more significant variants identified. If they are not very consistent, we would expect to observe the opposite. We see a fairly large shift of all stop variants to become more negative, due to the expected decrease in standard errors:


    
![png](OCNT-DMSLIB-1-Clonal-Var_files/OCNT-DMSLIB-1-Clonal-Var_15_0.png)
    


Because of the relative scaling, the heatmap looks fairly similar to those we have already seen, except that the scale extends much further negative to about -50. However, comparing the stop variant distributions highlights the power increase:


    
![png](OCNT-DMSLIB-1-Clonal-Var_files/OCNT-DMSLIB-1-Clonal-Var_18_0.png)
    


## OCNT-DMSLIB-1 DMS Assay Run #2: IFN-beta and varying seeding cell densities

1. [Barcode Sequencing Quality Control](#part1)
2. [Summary Statistics](#part2)
    -  [Within Condition](#part2a)
    -  [Across Conditions](#part2b)
3. [Visualization and Interpretation](#part3) 

### Barcode Sequencing Quality Control <a name="part1"></a>

In this dataset, we assay a set of conditions in triplicate, and the set of conditions included form a tidy 3x3 grid. Thus, we have 3x3x3=27 total samples with the following design:

| Cell Density | Cytokine | Sample Group |
| :--- | :--- | :--- |
| 7 million | None | 1ABC |
| 7 million | IFN-beta | 2ABC |
| 7 million | IFN-alpha | 3ABC |
| 3 million | None | 4ABC |
| 3 million | IFN-beta | 5ABC |
| 3 million | IFN-alpha | 6ABC |
| 1 million | None | 7ABC |
| 1 million | IFN-beta | 8ABC |
| 1 million | IFN-alpha | 9ABC |

We consider IFN-beta analogously to forskolin in the MC4R DMS assays, and use the exact same model evaluated per-position across chunk 10. We can decompose the total depth from each sample into categories based on the barcodes each read originate from:


    
    
    |sample | total depth|  perfect| perfect, subthreshold| imperfect, previously observed| not previously observed|
    |:------|-----------:|--------:|---------------------:|------------------------------:|-----------------------:|
    |1A     |    30753237| 18858254|               4642656|                        6307039|                  945288|
    |1B     |    31470859| 19216148|               4755362|                        6535706|                  963643|
    |1C     |    30750445| 18938478|               4569751|                        6265922|                  976294|
    |2A     |    33377984| 20382899|               5089353|                        6919471|                  986261|
    |2B     |    34640783| 21097625|               5315695|                        7158793|                 1068670|
    |2C     |    34721558| 21182153|               5313254|                        7170219|                 1055932|
    |3A     |    37187425| 23737223|               5852696|                        6473273|                 1124233|
    |3B     |    30669239| 19509875|               4863295|                        5375824|                  920245|
    |3C     |    28951805| 18459528|               4582295|                        5033691|                  876291|
    |4A     |    32672796| 19873035|               5014150|                        6752741|                 1032870|
    |4B     |    29440263| 17989604|               4379096|                        6107643|                  963920|
    |4C     |    30787735| 18784978|               4684271|                        6326364|                  992122|
    |5A     |    33312859| 20330715|               5149844|                        6836485|                  995815|
    |5B     |    33381253| 20365997|               5151209|                        6884644|                  979403|
    |5C     |    33520294| 20457584|               5151164|                        6924904|                  986642|
    |6A     |    37280298| 23501353|               5946081|                        6609178|                 1223686|
    |6B     |    34802892| 22145070|               5485484|                        6142274|                 1030064|
    |6C     |    31159576| 19795009|               4899717|                        5531635|                  933215|
    |7A     |    20084094| 12272969|               3135613|                        4059598|                  615914|
    |7B     |    26485775| 16246334|               4050564|                        5372122|                  816755|
    |7C     |    29644058| 18266942|               4298582|                        6109910|                  968624|
    |8A     |    31133920| 18933139|               4906041|                        6354583|                  940157|
    |8B     |    27754816| 16996843|               4221511|                        5702727|                  833735|
    |8C     |    26750902| 16226844|               4174193|                        5538851|                  811014|
    |9A     |    31111400| 19725506|               4972947|                        5484341|                  928606|
    |9B     |    30694702| 19557947|               4957800|                        5267180|                  911775|
    |9C     |    29653466| 18829133|               4798483|                        5164372|                  861478|


Below, the same table expressed as percentages:


    
    
    |sample | total depth|  perfect| perfect, subthreshold| imperfect, previously observed| not previously observed|
    |:------|-----------:|--------:|---------------------:|------------------------------:|-----------------------:|
    |1A     |         100| 61.32120|              15.09648|                       20.50854|                3.073784|
    |1B     |         100| 61.06013|              15.11037|                       20.76749|                3.062017|
    |1C     |         100| 61.58766|              14.86076|                       20.37669|                3.174894|
    |2A     |         100| 61.06690|              15.24763|                       20.73064|                2.954825|
    |2B     |         100| 60.90401|              15.34519|                       20.66579|                3.085005|
    |2C     |         100| 61.00577|              15.30246|                       20.65063|                3.041142|
    |3A     |         100| 63.83132|              15.73837|                       17.40716|                3.023154|
    |3B     |         100| 63.61382|              15.85724|                       17.52839|                3.000547|
    |3C     |         100| 63.75951|              15.82732|                       17.38645|                3.026723|
    |4A     |         100| 60.82441|              15.34656|                       20.66778|                3.161254|
    |4B     |         100| 61.10545|              14.87451|                       20.74588|                3.274155|
    |4C     |         100| 61.01449|              15.21473|                       20.54833|                3.222459|
    |5A     |         100| 61.02963|              15.45903|                       20.52206|                2.989281|
    |5B     |         100| 61.01028|              15.43144|                       20.62428|                2.933991|
    |5C     |         100| 61.03044|              15.36730|                       20.65884|                2.943417|
    |6A     |         100| 63.03961|              15.94966|                       17.72834|                3.282393|
    |6B     |         100| 63.62997|              15.76158|                       17.64875|                2.959708|
    |6C     |         100| 63.52785|              15.72459|                       17.75260|                2.994954|
    |7A     |         100| 61.10790|              15.61242|                       20.21300|                3.066676|
    |7B     |         100| 61.33985|              15.29336|                       20.28305|                3.083750|
    |7C     |         100| 61.62092|              14.50065|                       20.61091|                3.267515|
    |8A     |         100| 60.81193|              15.75786|                       20.41048|                3.019719|
    |8B     |         100| 61.23926|              15.21001|                       20.54680|                3.003929|
    |8C     |         100| 60.65905|              15.60393|                       20.70529|                3.031726|
    |9A     |         100| 63.40282|              15.98432|                       17.62808|                2.984777|
    |9B     |         100| 63.71766|              16.15197|                       17.15990|                2.970464|
    |9C     |         100| 63.49724|              16.18186|                       17.41574|                2.905151|


This is generally consistent with the previously estimated "sequencing + synthesis" error rate of 35%, along with a small fraction (~3.3%) of barcodes which are observed for the first time in the new data. These are probably sequencing errors which occurred in the barcode sequencing data, but were not observed in any prior run.

Now, we consider only those barcodes which are in the `perfect` category in the above tables. Each of these barcodes is attached to a particular oligo, so we can view the spread of how many barcodes are attached to each codon or each amino-acid level sequence:


    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_13_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_13_1.png)
    


Each of the above points is a particular codon variant or residue variant, and the count is the number of unique barcodes for that variant. We can see the distributions a bit more clearly as a density:


    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_16_0.png)
    


Note that while these barcodes have only one read count in a given sample in the 1x26 dataset, the fact that it is in the barcode map means it was detected in at least two sequencing replicates with at least three reads each in the barcode mapping data. Nonetheless, we can verify that we still would have enough barcode coverage even if we removed "lowly expressed" barcodes - again, we will not do this for the actual analysis, but we can just check:


    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_19_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_19_1.png)
    


We can also assess coverage by plotting the numbers of barcodes per variant across each amino acid. Since the samples are fairly similar


    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_21_0.png)
    


If we collapse these counts into bars, we see the approximate effects of changing the number of assayed cells on total barcode diversity:


    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_23_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_24_0.png)
    


### Summary Statistics <a name="part2"></a>

Like for MC4R, we have a forskolin-like condition (IFNb), our main stimulatory condition (IFNa), and a none-like condition (none!). We can apply the same hierarchical mixed model to jointly infer mutant vs WT effects within each of these three conditions, and then compute the change in variant impact (the fold change of fold changes) across conditions. We start by considering the primary mutant vs WT effects within each treatment (IFNa, IFNb, or None) and seeing how many significant hits we detect in each at a 1% FDR:

#### Within Condition <a name="part2a"></a>


```R
sumstat_sig <- sumstats_primary %>%
    group_by(density, condition) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
           nom_sig = if_else(p.value <= 1e-2, "significant", "nonsignificant"),
           adj_sig = if_else(p.adj <= 1e-2, "significant", "nonsignificant"))
```


    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_30_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_30_1.png)
    


We can extract the summary statistics for the stop codons to get a sense of the effect of cell density on stop signal:

    Picking joint bandwidth of 0.477
    
    Picking joint bandwidth of 0.297
    
    Picking joint bandwidth of 0.301
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_32_1.png)
    


#### Across Conditions <a name="part2b"></a>

We can also generate the same plots for the differences between pairs of these conditions:

    Picking joint bandwidth of 0.401
    
    Picking joint bandwidth of 0.231
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_38_1.png)
    


Finally, we can compute the final difference: of (IFNa - IFNb) vs (None - IFNb). This comparison does not _quite_ simplify to (IFNa - None) only because different barcodes are shared or non-shared amongst different sample groups. This returns the following final comparison:

    Picking joint bandwidth of 0.278
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_40_1.png)
    


### Visualization and Interpretation <a name="part3"></a>

Finally, we can make the traditional heatmaps for each of the above comparisons. For IFNa, this will hopefully look pretty similar to the last assay, but we can also see what the IFNb-based signals look like:


    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_42_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_42_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_42_2.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_42_3.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_43_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run2_files/OCNT-DMSLIB-1-assay-run2_43_1.png)
    


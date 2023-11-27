TYK2 Full-Length FlowDMS Variant Interpretation
================

1.  [Allele Frequencies](#part1)
2.  [ClinVar](#part2)
3.  [EVE and ESM1b](#part3)
4.  [SIFT and Polyphen2](#part4)
5.  [AlphaMissense](#part5)

#### Allele Frequencies <a name="part1"></a>

Each variant has a population-level allele frequency, which we can plot
for each global population for which we have data from Gnomad. In the
plots below, red points are significant in the DMS data at an adjusted
p-value \< 0.01:

![](flowdms_variant_interpretation/unnamed-chunk-4-1.png)<!-- -->

#### ClinVar <a name="part2"></a>

There are 316 missense variants in TYK2 noted in ClinVar, which we can
compare to our data stratified by ClinVar classification from Benign to
Pathogenic. Notably, P1104A is classified as Benign/Likely Benign, so
there is not necessairly a strong expectation that our functional
effects should associate strongly with a ClinVar class:

![](flowdms_variant_interpretation/unnamed-chunk-6-1.png)<!-- -->

The vast majority of variants are either not in ClinVar, or are of
uncertain significance. So, if we extract all variants NOT in either of
those two categories, we can show them together in a table with the DMS
adjusted p-value:

| mutation   | ClinVar category                             | DMS adjusted p-value |
|:-----------|:---------------------------------------------|---------------------:|
| Leu1014Pro | Pathogenic                                   |            0.0000000 |
| Leu757Val  | Pathogenic                                   |            0.0448691 |
| Gly799Arg  | Pathogenic                                   |            0.0000000 |
| Arg701Thr  | Conflicting interpretations of pathogenicity |            0.3507835 |
| Gly512Arg  | Conflicting interpretations of pathogenicity |            0.2693573 |
| Arg118Gln  | Conflicting interpretations of pathogenicity |            0.7059977 |
| Pro871Ser  | Likely benign                                |            0.1832957 |
| Gly39Ser   | Likely benign                                |            0.8342422 |
| His993Tyr  | Likely benign                                |            0.3953974 |
| Gly761Val  | Likely benign                                |            0.0000000 |
| Gly634Glu  | Likely benign                                |            0.4104980 |
| Val237Ile  | Likely benign                                |            0.2777453 |
| Val15Ala   | Benign/Likely benign                         |            0.9734686 |
| Ala928Val  | Benign/Likely benign                         |            0.4628942 |
| Pro820His  | Benign/Likely benign                         |            0.2988644 |
| Pro1104Ala | Benign/Likely benign                         |            0.0000000 |
| Ala53Thr   | Benign                                       |            0.6765941 |
| Val362Phe  | Benign                                       |            0.0000223 |
| Gly363Ser  | Benign                                       |            0.4492083 |
| Glu1163Gly | Benign                                       |            0.0000000 |
| Arg197His  | Benign                                       |            0.2626720 |
| Ile684Ser  | Benign                                       |            0.0000000 |

#### EVE and ESM1b <a name="part3"></a>

EVE (both scores from the Marks lab) and ESB1b have extremely similar
patterns:

![](flowdms_variant_interpretation/unnamed-chunk-9-1.png)<!-- -->

![](flowdms_variant_interpretation/unnamed-chunk-11-1.png)<!-- -->

#### SIFT and PolyPhen2 <a name="part4"></a>

The older predictors also have L-shaped patterns, though with visibly
lower resolution to distinguish between similar variants:

![](flowdms_variant_interpretation/unnamed-chunk-13-1.png)<!-- -->

![](flowdms_variant_interpretation/unnamed-chunk-14-1.png)<!-- -->

#### AlphaMissense <a name="part5"></a>

Comparing AlphaMissense variant effect predictions with TYK2 VAMP-seq
for all 17 chunks.

![](flowdms_variant_interpretation/unnamed-chunk-18-1.png)<!-- -->

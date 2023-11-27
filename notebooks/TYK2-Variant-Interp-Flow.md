### TYK2 Full-Length FlowDMS Variant Interpretation

1. [Allele Frequencies](#part1)
2. [ClinVar](#part2)
3. [EVE and ESM1b](#part3)
4. [SIFT and Polyphen2](#part4)
5. [AlphaMissense](#part5)

#### Allele Frequencies <a name="part1"></a>

Each variant has a population-level allele frequency, which we can plot for each global population for which we have data from Gnomad. In the plots below, red points are significant in the DMS data at an adjusted p-value < 0.01:


    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_8_0.png)
    



    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_8_1.png)
    


#### ClinVar <a name="part2"></a>

There are 316 missense variants in TYK2 noted in ClinVar, which we can compare to our data stratified by ClinVar classification from Benign to Pathogenic. Notably, P1104A is classified as Benign/Likely Benign, so there is not necessairly a strong expectation that our functional effects should associate strongly with a ClinVar class:


    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_12_0.png)
    


The vast majority of variants are either not in ClinVar, or are of uncertain significance. So, if we extract all variants NOT in either of those two categories, we can show them together in a table with the DMS adjusted p-value:


<table class="dataframe">
<caption>A tibble: 22 × 4</caption>
<thead>
	<tr><th scope=col>mutation</th><th scope=col>ClinVar category</th><th scope=col>midpoint shift</th><th scope=col>DMS adjusted p-value</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Leu1014Pro</td><td>Pathogenic                                  </td><td>-0.398049762</td><td> 3.056952e-12</td></tr>
	<tr><td>Leu757Val </td><td>Pathogenic                                  </td><td>-0.308629769</td><td> 4.486910e-02</td></tr>
	<tr><td>Gly799Arg </td><td>Pathogenic                                  </td><td>-0.421945818</td><td> 1.730525e-16</td></tr>
	<tr><td>Arg701Thr </td><td>Conflicting interpretations of pathogenicity</td><td>-0.126109050</td><td> 3.507835e-01</td></tr>
	<tr><td>Gly512Arg </td><td>Conflicting interpretations of pathogenicity</td><td>-0.192161858</td><td> 2.693573e-01</td></tr>
	<tr><td>Arg118Gln </td><td>Conflicting interpretations of pathogenicity</td><td> 0.042268295</td><td> 7.059977e-01</td></tr>
	<tr><td>Pro871Ser </td><td>Likely benign                               </td><td>-0.193946073</td><td> 1.832957e-01</td></tr>
	<tr><td>Gly39Ser  </td><td>Likely benign                               </td><td>-0.034260258</td><td> 8.342422e-01</td></tr>
	<tr><td>His993Tyr </td><td>Likely benign                               </td><td>-0.110472522</td><td> 3.953974e-01</td></tr>
	<tr><td>Gly761Val </td><td>Likely benign                               </td><td>-0.464036766</td><td> 6.139221e-74</td></tr>
	<tr><td>Gly634Glu </td><td>Likely benign                               </td><td>-0.063984826</td><td> 4.104980e-01</td></tr>
	<tr><td>Val237Ile </td><td>Likely benign                               </td><td>-0.139147064</td><td> 2.777453e-01</td></tr>
	<tr><td>Val15Ala  </td><td>Benign/Likely benign                        </td><td> 0.004749789</td><td> 9.734686e-01</td></tr>
	<tr><td>Ala928Val </td><td>Benign/Likely benign                        </td><td>-0.061107045</td><td> 4.628942e-01</td></tr>
	<tr><td>Pro820His </td><td>Benign/Likely benign                        </td><td>-0.134179762</td><td> 2.988644e-01</td></tr>
	<tr><td>Pro1104Ala</td><td>Benign/Likely benign                        </td><td>-0.361154147</td><td> 3.462099e-81</td></tr>
	<tr><td>Ala53Thr  </td><td>Benign                                      </td><td>-0.070303125</td><td> 6.765941e-01</td></tr>
	<tr><td>Val362Phe </td><td>Benign                                      </td><td>-0.382455106</td><td> 2.232471e-05</td></tr>
	<tr><td>Gly363Ser </td><td>Benign                                      </td><td>-0.066182059</td><td> 4.492083e-01</td></tr>
	<tr><td>Glu1163Gly</td><td>Benign                                      </td><td>-0.319794600</td><td> 3.266925e-08</td></tr>
	<tr><td>Arg197His </td><td>Benign                                      </td><td>-0.192949230</td><td> 2.626720e-01</td></tr>
	<tr><td>Ile684Ser </td><td>Benign                                      </td><td>-0.453308901</td><td>4.103510e-111</td></tr>
</tbody>
</table>



#### EVE and ESM1b <a name="part3"></a>

EVE (both scores from the Marks lab) and ESB1b have extremely similar patterns:


    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_17_0.png)
    



    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_19_0.png)
    


#### SIFT and PolyPhen2 <a name="part4"></a>

The older predictors also have L-shaped patterns, though with visibly lower resolution to distinguish between similar variants:


    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_22_0.png)
    



    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_23_0.png)
    


#### AlphaMissense <a name="part5"></a>

Comparing AlphaMissense variant effect predictions with TYK2 VAMP-seq for chunk 10.


    
![png](TYK2-Variant-Interp-Flow_files/TYK2-Variant-Interp-Flow_27_0.png)
    


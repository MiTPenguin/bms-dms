## OCNT-DMSLIB-1 DMS Assay Run #6: TYK2 IL-23 Mini-DMS

This is the first DMS dataset derived from the IL-23 reporter system. There are 12 samples (4 replicates in each of 3 conditions), with the following structure:

| ID | Condition | Dosage | Replicates
| --- | ----------- | --- | ----------- |
| 1 | None | 0 | 4 |
| 2 | IL-23 | 100 | 4 |
| 3 | gDNA | NA | 4 |

Let's check some effects we expect to observe (mainly QC parameters and stop codons), then examine the global distribution of mutant vs WT effects across conditions.

1. [Barcode Sequencing Distributions](#part1)
2. [Inference and Stop Codon Effects](#part2)
3. [Visualizations](#part3)
4. [Raw Data Review](#part4)

### Barcode Sequencing Distributions <a name="part1"></a>


    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_5_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_5_1.png)
    


To get a sense of positional distribution, we can show the same data as lineplots across the length of TYK2. Below is an example using sample `1A` only; the remaining plots can be found [here](./coverage-plots):


    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_9_0.png)
    


### Inference and Stop Codon Effects <a name="part2"></a>


    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_13_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_14_0.png)
    



<table class="dataframe">
<caption>A tibble: 27 × 8</caption>
<thead>
	<tr><th scope=col>pos</th><th scope=col>condition</th><th scope=col>aa</th><th scope=col>estimate</th><th scope=col>std.error</th><th scope=col>statistic</th><th scope=col>p.value</th><th scope=col>p.adj</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 911</td><td>IL23100</td><td>*</td><td>-0.6933807</td><td>0.1586946</td><td>-4.369278</td><td>1.246582e-05</td><td>5.813229e-03</td></tr>
	<tr><td> 916</td><td>IL23100</td><td>C</td><td>-0.6552358</td><td>0.1522913</td><td>-4.302517</td><td>1.688685e-05</td><td>6.763842e-03</td></tr>
	<tr><td> 926</td><td>None0  </td><td>Y</td><td> 0.7069348</td><td>0.1652121</td><td> 4.278952</td><td>1.877749e-05</td><td>6.852969e-03</td></tr>
	<tr><td> 935</td><td>None0  </td><td>V</td><td> 0.9636704</td><td>0.1986476</td><td> 4.851156</td><td>1.227438e-06</td><td>1.470740e-03</td></tr>
	<tr><td> 936</td><td>None0  </td><td>Y</td><td> 1.0091754</td><td>0.1776301</td><td> 5.681332</td><td>1.336499e-08</td><td>3.739524e-05</td></tr>
	<tr><td> 945</td><td>IL23100</td><td>F</td><td>-0.7889376</td><td>0.1643836</td><td>-4.799369</td><td>1.591665e-06</td><td>1.484493e-03</td></tr>
	<tr><td> 954</td><td>IL23100</td><td>C</td><td>-0.8393561</td><td>0.1981266</td><td>-4.236463</td><td>2.270683e-05</td><td>7.941713e-03</td></tr>
	<tr><td> 954</td><td>IL23100</td><td>W</td><td>-0.8075816</td><td>0.1839288</td><td>-4.390728</td><td>1.129716e-05</td><td>5.578141e-03</td></tr>
	<tr><td> 955</td><td>None0  </td><td>F</td><td> 0.9665538</td><td>0.1979005</td><td> 4.884039</td><td>1.039343e-06</td><td>1.470740e-03</td></tr>
	<tr><td> 956</td><td>IL23100</td><td>Y</td><td> 0.9680600</td><td>0.1644993</td><td> 5.884887</td><td>3.983252e-09</td><td>1.671771e-05</td></tr>
	<tr><td> 960</td><td>IL23100</td><td>S</td><td> 0.8051737</td><td>0.1755737</td><td> 4.585957</td><td>4.519113e-06</td><td>2.709531e-03</td></tr>
	<tr><td> 973</td><td>None0  </td><td>I</td><td> 0.8069440</td><td>0.1837594</td><td> 4.391307</td><td>1.126711e-05</td><td>5.578141e-03</td></tr>
	<tr><td> 975</td><td>None0  </td><td>C</td><td> 0.8027506</td><td>0.1736246</td><td> 4.623485</td><td>3.773455e-06</td><td>2.535374e-03</td></tr>
	<tr><td>1051</td><td>None0  </td><td>C</td><td>-1.3912224</td><td>0.3048583</td><td>-4.563505</td><td>5.030651e-06</td><td>2.815152e-03</td></tr>
	<tr><td>1069</td><td>IL23100</td><td>I</td><td>-0.9280583</td><td>0.2166009</td><td>-4.284646</td><td>1.830302e-05</td><td>6.852969e-03</td></tr>
	<tr><td>1073</td><td>IL23100</td><td>K</td><td>-1.1283169</td><td>0.2322951</td><td>-4.857258</td><td>1.190227e-06</td><td>1.470740e-03</td></tr>
	<tr><td>1081</td><td>IL23100</td><td>E</td><td>-0.9138728</td><td>0.1925255</td><td>-4.746763</td><td>2.066981e-06</td><td>1.735024e-03</td></tr>
	<tr><td>1087</td><td>None0  </td><td>H</td><td>-1.0440507</td><td>0.2426862</td><td>-4.302061</td><td>1.692169e-05</td><td>6.763842e-03</td></tr>
	<tr><td>1103</td><td>IL23100</td><td>D</td><td>-0.9159106</td><td>0.2195300</td><td>-4.172144</td><td>3.017470e-05</td><td>9.723666e-03</td></tr>
	<tr><td>1104</td><td>None0  </td><td>C</td><td>-0.9240435</td><td>0.2144990</td><td>-4.307916</td><td>1.648000e-05</td><td>6.763842e-03</td></tr>
	<tr><td>1108</td><td>None0  </td><td>R</td><td>-1.5366967</td><td>0.3134930</td><td>-4.901854</td><td>9.493652e-07</td><td>1.470740e-03</td></tr>
	<tr><td>1109</td><td>IL23100</td><td>H</td><td>-1.1807548</td><td>0.2834792</td><td>-4.165226</td><td>3.110446e-05</td><td>9.723666e-03</td></tr>
	<tr><td>1111</td><td>IL23100</td><td>A</td><td>-0.9410562</td><td>0.2036239</td><td>-4.621542</td><td>3.808991e-06</td><td>2.535374e-03</td></tr>
	<tr><td>1112</td><td>IL23100</td><td>V</td><td>-1.2356740</td><td>0.2561109</td><td>-4.824762</td><td>1.401706e-06</td><td>1.470740e-03</td></tr>
	<tr><td>1114</td><td>IL23100</td><td>N</td><td>-1.0324866</td><td>0.2237129</td><td>-4.615230</td><td>3.926597e-06</td><td>2.535374e-03</td></tr>
	<tr><td>1114</td><td>None0  </td><td>Q</td><td>-1.0282692</td><td>0.2469448</td><td>-4.163964</td><td>3.127698e-05</td><td>9.723666e-03</td></tr>
	<tr><td>1115</td><td>None0  </td><td>*</td><td> 1.6849505</td><td>0.2711503</td><td> 6.214084</td><td>5.162493e-10</td><td>4.333397e-06</td></tr>
</tbody>
</table>



### Visualizations <a name="part3"></a>


    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_17_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_17_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_20_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_20_1.png)
    


### Raw Data Review <a name="part4"></a>

These data look very flat across conditions. To make sure we didn't miss anything computationally or statistically, we can examine the raw barcode count data of a few positions in very high resolution, and compare it to another dataset (here, run 3). Since it is a spike-in shared by both datasets, we choose position *930R* and compare it to WT for its segment (chunk 14). Let's just look at the raw counts:


    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_23_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_23_1.png)
    


To compare these distributions, we extract and plot the empirical cumulative density, which shows the varint effect size growing more negative with increasing IFNa concentration (but with no effect change in IL-23):


    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_25_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_25_1.png)
    


This is pretty clear, but we can get more precise. Since we have the barcode sequences associated with each count, we can subset each dataset to only include the barcodes which are detected in the _other_ dataset at least once.


    
    
    |group       | number of barcodes|
    |:-----------|------------------:|
    |One Assay   |              35713|
    |Both Assays |               7492|


Most barcodes are only detected in one assay, but 7492 barcodes are detected in both (though presumably not in all samples). We can subset to just these barcodes, and plot the same data views:


    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_30_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_30_1.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_31_0.png)
    



    
![png](OCNT-DMSLIB-1-assay-run6_files/OCNT-DMSLIB-1-assay-run6_31_1.png)
    


So, even for these specific ~7500 barcodes, there is no different in raw counts within each sample in a way analogous to that observed in assay 3 with an IFN-alpha dose response.

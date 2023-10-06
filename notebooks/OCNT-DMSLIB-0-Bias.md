## OCNT-DMSLIB-0 MC4R Gs/Gq Bias

We have three concentrations each of aMSH and THIQ in each of the Gs and Gq reporter systems, but not precisely the same concentrations. We also have a None condition in both systems, and Forskolin in Gs only. We generated Unnormalized summary statistics within each condition, and then can create the following concentration pairings between Gs and Gq:

| Comparison | Gs Concentration | Gq Concentration | 
| --- | ----------- | ----- |
| None | 0 | 0 |
| Low aMSH | 5e-10 | 2e-08 |
| Med aMSH | 5e-09 | 5e-08 |
| High aMSH | 2e-08 | 1e-06 |
| Low THIQ | 4e-10 | 3e-09 |
| Med THIQ | 4e-09 | 9e-09 |
| High THIQ | 1.2e-08 | 1e-07 | 

We compute each of these contrasts and examine the results below.

1. [Summary Statistics](#part1)
2. [Heatmaps](#part2)
3. [Biased Positions and Variants](#part3)
4. [Medium Comparisons](#part4)

### Summary Statistics <a name="part1"></a>

To get an initial intuition for the relationship between the summary statistics we have from Gs and Gq, we can plot the Z-statistics per variant from either dataset on the same plot, startified by the pairings defined above.


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_5_0.png)
    


We actually want to compute these differences, and then test for a significant effect at our typical thesholds of 1% and 5% FDR. Those results are now stored in the file called `MC4R-Bias.tsv` located [here](../sumstats/MC4R). We with our other contrast models, they seem fairly well calibrated:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_9_0.png)
    



    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_9_1.png)
    


The stop effects are informative, since we expect a stop to be equally damaging to both Gs and Gq. Thus, and differences in the stop profile are presumably due to different degrees of pathway activation by the agonist, not a difference in mutation effects on pathway activation.


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_11_0.png)
    



    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_11_1.png)
    


In each comparison, we can apply a 1% FDR threshold and count the number of significant variants in each direction, noting how many are stops in each case:


<table class="dataframe">
<caption>A tibble: 7 × 6</caption>
<thead>
	<tr><th scope=col>condition</th><th scope=col>Gs_Non-Stop</th><th scope=col>Gq_Non-Stop</th><th scope=col>Gs_Stop</th><th scope=col>Gq_Stop</th><th scope=col>total</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Gs_aMSH_Low - Gq_aMSH_Low  </td><td> 11</td><td>430</td><td> 0</td><td>46</td><td>487</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med  </td><td>199</td><td> 39</td><td> 2</td><td> 2</td><td>242</td></tr>
	<tr><td>Gs_aMSH_High - Gq_aMSH_High</td><td>830</td><td> 23</td><td>18</td><td> 0</td><td>871</td></tr>
	<tr><td>Gs_None_None - Gq_None_None</td><td>102</td><td>122</td><td> 0</td><td> 9</td><td>233</td></tr>
	<tr><td>Gs_THIQ_Low - Gq_THIQ_Low  </td><td> 15</td><td>309</td><td> 1</td><td>40</td><td>365</td></tr>
	<tr><td>Gs_THIQ_Med - Gq_THIQ_Med  </td><td> 66</td><td> 58</td><td> 0</td><td> 4</td><td>128</td></tr>
	<tr><td>Gs_THIQ_High - Gq_THIQ_High</td><td>140</td><td> 63</td><td> 8</td><td> 0</td><td>211</td></tr>
</tbody>
</table>



And then the same for the reversed comparisons (note that the medium comparison changes due to multiple testing correction being applied to a different total set of p-values, even though the raw p-values would be the same). In all conditions, there are fewer significant stops:


<table class="dataframe">
<caption>A tibble: 7 × 6</caption>
<thead>
	<tr><th scope=col>condition</th><th scope=col>Gs_Non-Stop</th><th scope=col>Gq_Non-Stop</th><th scope=col>Gs_Stop</th><th scope=col>Gq_Stop</th><th scope=col>total</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Gs_aMSH_Low - Gq_aMSH_High </td><td>182</td><td> 54</td><td>6</td><td>1</td><td>243</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med  </td><td>135</td><td> 29</td><td>2</td><td>1</td><td>167</td></tr>
	<tr><td>Gs_aMSH_High - Gq_aMSH_Low </td><td> 31</td><td> 38</td><td>0</td><td>2</td><td> 71</td></tr>
	<tr><td>Gs_None_None - Gq_None_None</td><td> 73</td><td> 92</td><td>0</td><td>6</td><td>171</td></tr>
	<tr><td>Gs_THIQ_Low - Gq_THIQ_High </td><td> 18</td><td>133</td><td>1</td><td>5</td><td>157</td></tr>
	<tr><td>Gs_THIQ_Med - Gq_THIQ_Med  </td><td> 37</td><td> 43</td><td>0</td><td>1</td><td> 81</td></tr>
	<tr><td>Gs_THIQ_High - Gq_THIQ_Low </td><td> 29</td><td> 30</td><td>1</td><td>0</td><td> 60</td></tr>
</tbody>
</table>



### Heatmaps <a name="part2"></a>

Next, we project these into heatmaps as usual. This Z-statistic is scaled in the range (-15 to +5) generally applied to other heatmaps:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_17_0.png)
    



    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_17_1.png)
    


When inverted, the Stop bands disappear almost entirely, which is an encouraging sign. We can explore these results in a bit more detail by looking at the top Gs and Gq biased positions across all comparisons, and then examining those positions.

### Biased Positions and Variants <a name="part3"></a>

From the set of significantly biased variants, we can extract the top 15 most significant for each comparison. Since None and Med seem the most well-matched, let's look at those first:


<table class="dataframe">
<caption>A tibble: 15 × 6</caption>
<thead>
	<tr><th scope=col>condition</th><th scope=col>pos</th><th scope=col>aa</th><th scope=col>estimate</th><th scope=col>std.error</th><th scope=col>p.adj</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>107</td><td>S</td><td>0.9310918</td><td>0.1428655</td><td>6.156698e-08</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>164</td><td>L</td><td>0.7992729</td><td>0.1332090</td><td>7.207639e-07</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>106</td><td>C</td><td>1.4110898</td><td>0.2417894</td><td>1.677164e-06</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>133</td><td>C</td><td>1.8388852</td><td>0.3232554</td><td>3.244243e-06</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>231</td><td>W</td><td>2.6472974</td><td>0.4700389</td><td>4.217819e-06</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>161</td><td>V</td><td>1.1088932</td><td>0.2061404</td><td>1.365518e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>216</td><td>H</td><td>1.7369203</td><td>0.3345990</td><td>3.015297e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>162</td><td>N</td><td>1.1323185</td><td>0.2192833</td><td>3.416490e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>208</td><td>S</td><td>1.3735765</td><td>0.2675047</td><td>3.823896e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>283</td><td>T</td><td>1.6009896</td><td>0.3151098</td><td>4.822350e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>132</td><td>G</td><td>1.2305014</td><td>0.2430200</td><td>5.196998e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>124</td><td>T</td><td>1.2668927</td><td>0.2524752</td><td>6.156492e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>205</td><td>I</td><td>1.5332576</td><td>0.3075878</td><td>6.863223e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>135</td><td>S</td><td>1.0223716</td><td>0.2067014</td><td>8.060000e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>190</td><td>G</td><td>1.5834051</td><td>0.3218140</td><td>8.958086e-05</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 15 × 6</caption>
<thead>
	<tr><th scope=col>condition</th><th scope=col>pos</th><th scope=col>aa</th><th scope=col>estimate</th><th scope=col>std.error</th><th scope=col>p.adj</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td> 79</td><td>R</td><td>-1.2806058</td><td>0.1907055</td><td>2.181565e-08</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>158</td><td>R</td><td>-1.2604460</td><td>0.1884001</td><td>2.462756e-08</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>254</td><td>P</td><td>-1.3015428</td><td>0.1952214</td><td>2.754899e-08</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>223</td><td>L</td><td>-0.8652945</td><td>0.1338754</td><td>7.907665e-08</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td> 79</td><td>S</td><td>-1.2543983</td><td>0.1981559</td><td>1.494266e-07</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>228</td><td>R</td><td>-0.9505143</td><td>0.1652512</td><td>2.497869e-06</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>152</td><td>R</td><td>-1.0361797</td><td>0.1812785</td><td>2.914687e-06</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td> 77</td><td>K</td><td>-2.0378593</td><td>0.3810963</td><td>1.571219e-05</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td> 77</td><td>R</td><td>-0.9997778</td><td>0.2062942</td><td>1.186649e-04</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>214</td><td>K</td><td>-1.6148591</td><td>0.3352080</td><td>1.313257e-04</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>230</td><td>S</td><td>-0.8206461</td><td>0.1721997</td><td>1.579628e-04</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>232</td><td>R</td><td>-0.7996768</td><td>0.1706802</td><td>2.100958e-04</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>239</td><td>C</td><td>-1.2797450</td><td>0.2841858</td><td>4.020729e-04</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td>214</td><td>A</td><td>-1.2279186</td><td>0.2819238</td><td>6.759191e-04</td></tr>
	<tr><td>Gs_aMSH_Med - Gq_aMSH_Med</td><td> 74</td><td>R</td><td>-0.9543596</td><td>0.2198311</td><td>7.105419e-04</td></tr>
</tbody>
</table>



Immediately, we can see some multi-variant positions, and we can extract and profile the position of the top hit in each direction:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_22_0.png)
    



    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_23_0.png)
    


The None condition is of particular interest, since it is the only condition where the concentration matching is absolutely precise. We can extract the top variants in both directions there as well:


<table class="dataframe">
<caption>A tibble: 15 × 6</caption>
<thead>
	<tr><th scope=col>pos</th><th scope=col>aa</th><th scope=col>estimate</th><th scope=col>std.error</th><th scope=col>condition</th><th scope=col>p.adj</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>140</td><td>V</td><td> 1.6394493</td><td>0.2421233</td><td>Gs_None_None - Gq_None_None</td><td>1.630696e-08</td></tr>
	<tr><td>140</td><td>N</td><td> 2.4913374</td><td>0.3733423</td><td>Gs_None_None - Gq_None_None</td><td>2.705063e-08</td></tr>
	<tr><td>223</td><td>S</td><td>-0.9563646</td><td>0.1530305</td><td>Gs_None_None - Gq_None_None</td><td>2.249252e-07</td></tr>
	<tr><td>140</td><td>Q</td><td> 2.0197471</td><td>0.3258109</td><td>Gs_None_None - Gq_None_None</td><td>2.865228e-07</td></tr>
	<tr><td>122</td><td>L</td><td> 1.3225628</td><td>0.2155991</td><td>Gs_None_None - Gq_None_None</td><td>3.854305e-07</td></tr>
	<tr><td>158</td><td>S</td><td> 1.2656483</td><td>0.2068701</td><td>Gs_None_None - Gq_None_None</td><td>4.109727e-07</td></tr>
	<tr><td>147</td><td>S</td><td>-1.2810549</td><td>0.2107587</td><td>Gs_None_None - Gq_None_None</td><td>4.947043e-07</td></tr>
	<tr><td>140</td><td>A</td><td> 1.4511983</td><td>0.2434871</td><td>Gs_None_None - Gq_None_None</td><td>8.736728e-07</td></tr>
	<tr><td>140</td><td>T</td><td> 1.4453240</td><td>0.2508869</td><td>Gs_None_None - Gq_None_None</td><td>2.398563e-06</td></tr>
	<tr><td>148</td><td>W</td><td> 2.4850969</td><td>0.4552593</td><td>Gs_None_None - Gq_None_None</td><td>9.479685e-06</td></tr>
	<tr><td>140</td><td>E</td><td> 1.8743433</td><td>0.3453761</td><td>Gs_None_None - Gq_None_None</td><td>1.099779e-05</td></tr>
	<tr><td>140</td><td>S</td><td> 1.8972448</td><td>0.3564399</td><td>Gs_None_None - Gq_None_None</td><td>1.743692e-05</td></tr>
	<tr><td> 79</td><td>R</td><td> 1.0647037</td><td>0.2007308</td><td>Gs_None_None - Gq_None_None</td><td>1.857352e-05</td></tr>
	<tr><td>128</td><td>W</td><td>-2.7489907</td><td>0.5228001</td><td>Gs_None_None - Gq_None_None</td><td>2.236487e-05</td></tr>
	<tr><td> 83</td><td>Y</td><td> 1.7611852</td><td>0.3351946</td><td>Gs_None_None - Gq_None_None</td><td>2.270466e-05</td></tr>
</tbody>
</table>



That is a lot of position 140:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_28_0.png)
    


And the most significant Gq-biased in the same comparison is at position 223. The stop effect here notably is positive in Gq, but not significantly so after adjustment (the most significant hit here is 223S):


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_30_0.png)
    


### Medium Comparisons <a name="part4"></a>

Since the medium concentration pair seems the best calibrated in terms of Stop effects, we can extract just that comparison and plot the Gs/Gq relationship in higher resolution:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_37_0.png)
    



    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_37_1.png)
    


There are three points that have FDR < 0.01 in both Gs and Gq based on the first plot, but in opposite directions. Those variants are:


<table class="dataframe">
<caption>A tibble: 3 × 5</caption>
<thead>
	<tr><th scope=col>pos</th><th scope=col>aa</th><th scope=col>statistic_Gs</th><th scope=col>statistic_Gq</th><th scope=col>Bias</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 48</td><td>Q</td><td> 3.130821</td><td>-3.391171</td><td>Biased</td></tr>
	<tr><td>214</td><td>K</td><td>-3.301927</td><td> 3.558693</td><td>Biased</td></tr>
	<tr><td>254</td><td>P</td><td>-4.411550</td><td> 5.000009</td><td>Biased</td></tr>
</tbody>
</table>



We can profile these positions accross both Gs and Gq, and including all variants at those positions for comparison:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_41_0.png)
    


Let's extract a more refined variant set with several requirements:

- Negative estimate (LoF) and Significant (FDR < 0.01) in Gq
- Non-Significant (FDR > 0.01) in Gs
- Significant (FDR < 0.01) in Gs - Gq

This returns a list of 166 variants, which we can rank by the `Gs - Gq` FDR (all of these are FDR < 0.01, but we can prioritize the most signifiant biases):


<table class="dataframe">
<caption>A tibble: 15 × 9</caption>
<thead>
	<tr><th scope=col>pos</th><th scope=col>aa</th><th scope=col>estimate_Gs</th><th scope=col>std.error_Gs</th><th scope=col>estimate_Gq</th><th scope=col>std.error_Gq</th><th scope=col>p.adj_Gs</th><th scope=col>p.adj_Gq</th><th scope=col>p.adj_bias</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>107</td><td>S</td><td> 0.05019250</td><td>0.08574044</td><td>-0.8808993</td><td>0.1142766</td><td>0.75112345</td><td>9.453057e-13</td><td>6.156698e-08</td></tr>
	<tr><td>164</td><td>L</td><td> 0.14259169</td><td>0.08740016</td><td>-0.6566812</td><td>0.1005279</td><td>0.25489846</td><td>2.197823e-09</td><td>7.207639e-07</td></tr>
	<tr><td>106</td><td>C</td><td> 0.32213328</td><td>0.14538607</td><td>-1.0889565</td><td>0.1931968</td><td>0.09137741</td><td>3.418358e-07</td><td>1.677164e-06</td></tr>
	<tr><td>133</td><td>C</td><td> 0.62766468</td><td>0.20178417</td><td>-1.2112206</td><td>0.2525414</td><td>0.01011171</td><td>2.022303e-05</td><td>3.244243e-06</td></tr>
	<tr><td>231</td><td>W</td><td> 0.72622662</td><td>0.24552341</td><td>-1.9210708</td><td>0.4008177</td><td>0.01556518</td><td>2.049254e-05</td><td>4.217819e-06</td></tr>
	<tr><td>161</td><td>V</td><td> 0.38688176</td><td>0.13985620</td><td>-0.7220114</td><td>0.1514401</td><td>0.02580560</td><td>2.293851e-05</td><td>1.365518e-05</td></tr>
	<tr><td>216</td><td>H</td><td> 0.01295013</td><td>0.18460890</td><td>-1.7239702</td><td>0.2790628</td><td>0.97512597</td><td>1.779972e-08</td><td>3.015297e-05</td></tr>
	<tr><td>162</td><td>N</td><td> 0.08642312</td><td>0.14297689</td><td>-1.0458954</td><td>0.1662612</td><td>0.74161960</td><td>9.285810e-09</td><td>3.416490e-05</td></tr>
	<tr><td>208</td><td>S</td><td> 0.19906606</td><td>0.14834827</td><td>-1.1745104</td><td>0.2226018</td><td>0.37479858</td><td>2.109963e-06</td><td>3.823896e-05</td></tr>
	<tr><td>283</td><td>T</td><td> 0.10881484</td><td>0.08537688</td><td>-1.4921748</td><td>0.3033232</td><td>0.40654953</td><td>1.155684e-05</td><td>4.822350e-05</td></tr>
	<tr><td>132</td><td>G</td><td> 0.31994454</td><td>0.15688135</td><td>-0.9105568</td><td>0.1855990</td><td>0.12886678</td><td>1.228792e-05</td><td>5.196998e-05</td></tr>
	<tr><td>124</td><td>T</td><td> 0.32667650</td><td>0.14969723</td><td>-0.9402162</td><td>0.2033088</td><td>0.09779677</td><td>4.281862e-05</td><td>6.156492e-05</td></tr>
	<tr><td>205</td><td>I</td><td>-0.16225090</td><td>0.15694485</td><td>-1.6955085</td><td>0.2645346</td><td>0.52316205</td><td>4.612841e-09</td><td>6.863223e-05</td></tr>
	<tr><td>135</td><td>S</td><td> 0.14756500</td><td>0.13412832</td><td>-0.8748066</td><td>0.1572738</td><td>0.49016815</td><td>5.024696e-07</td><td>8.060000e-05</td></tr>
	<tr><td>190</td><td>G</td><td> 0.29255983</td><td>0.16575652</td><td>-1.2908452</td><td>0.2758424</td><td>0.20738726</td><td>3.371521e-05</td><td>8.958086e-05</td></tr>
</tbody>
</table>



To verify the patterns we expect, we can extract some of the positions of the most significant variants and plot them as before:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_46_0.png)
    


These seem to consistently return cases where the stop effects are quite similar, but there are obvious individual variants which damage Gq but not Gs. We can count these by both position and by residue:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_48_0.png)
    



    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_49_0.png)
    


Upon examination, the following positions have 3 or 4 significant variants: 107, 124, 135, 137, 145, 148, 161, 162, 262, 274, 304. We can plot the Gs and Gq summary statistics for those positions in a harmonized format:


    
![png](OCNT-DMSLIB-0-Bias_files/OCNT-DMSLIB-0-Bias_53_0.png)
    


### TYK2 DMS/FlowDMS Quantitative Contrasts

  1.  [Error Propagation: FlowDMS](#part1)
  2.  [Error Propagation: DMS](#part2)
  3.  [Global Summary Statistics](#part3)
  4.  [Positional Distributions](#part4)

#### Error Propagation: FlowDMS <a name="part1"></a>

We start by finding the FlowDMS midpoints for chunks 3, 10, and 14 while propagating the per-bin errors. As before, we can extract the WT score estimate and standard error, with a separate estimate per chunk. The estimates for all three chunks are extremely similar, and the full midpoint summary statistics for the mini-FlowDMS are located [here](../sumstats/TYK2-VAMP/midpoint/).


    
    
    | chunk|  WT score| WT score standard error|
    |-----:|---------:|-----------------------:|
    |     3| 0.7037447|               0.0106709|
    |    10| 0.6756306|               0.0058669|
    |    14| 0.6874737|               0.0074624|


Next, we apply the transformation described in the [TYK2 calibration notebook](../TYK2-Calibration.md), where we use the median stop effect computed _within_ each chunk as with WT:

$$ Score = \frac{midpoint - median(stop)}{median(WT) - median(stop)} $$

We get the same median scores and standard errors for the set of stop effect per chunk:


    
    
    | chunk| Stop score| Stop score standard error|
    |-----:|----------:|-------------------------:|
    |     3|  0.2794127|                 0.0631790|
    |    10|  0.2995426|                 0.0477065|
    |    14|  0.3199645|                 0.0527895|



    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_9_0.png)
    


#### Error Propagation: DMS <a name="part2"></a>

Next, we transform the DMS summary statistics to the same scale and propagate the error with the following transformation:

$$ FC = 2^{Log2FoldChange} $$

$$ Score = \frac{FC - median(FC_{stop})}{1 - median(FC_{stop})} $$

As before, the subtraction and division steps are straightforward to propagate error through. However, the exponentiation involved in going from Log2FoldChange to FoldChange results in asymmetric confidence intervals and thus is not as simple but do-able by sampling. We can make an analogous distribution of stop effects as above on the new scale. Immediately, it is obvious that the range of each distribution is much wider than that of FlowDMS, which seems intuitive given the necessary FlowDMS binning:


    
    
    |chunk | median(FCstop)|
    |:-----|--------------:|
    |3     |      0.3266600|
    |10    |      0.3131436|
    |14    |      0.3187060|



    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_11_1.png)
    


#### Global Summary Statistics <a name="part3"></a>

Now, we want to characterize the difference between these two rescaled values across each assay. To evaluate this, we sample and transform the DMS summary statistics, including exponentiation, and compare it to an analogous sampling from the FlowDMS summary statistics. To do so, we sample "effects" from our summary statistics using the distributions:

$$ DMS:\ Normal(Log2FoldChange, Log2StandardError) $$

$$ FlowDMS:\ Normal(Midpoint, MidpointSE) $$

We sample from these distributions many times, perform the rescaling operations on each resampled point, and then obtain the resulting distributions after any transformations. The mean and standard deviation of this distribution (a _sampling distribution_) are our transformed estimates and standard errors, respectively. We can see what this looks like with an individual spike-in from chunk 10 (669P):


    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_13_0.png)
    


We can test the extent of overlap between the two distributions for each variant by approximating each as normal, or simply using an empirical FDR. Here, we use the pooled standard error to test for a non-zero mean difference with a Wald test. The first facet plot below shows the 0-1 scale estimates for each assay type, with Stops highlighted. The second facet plot shows the same, except with significant (FDR < 1%) cross-assay contrasts highlighted:


    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_17_0.png)
    



    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_18_0.png)
    


#### Positional Distributions <a name="part4"></a>
We can use these tests to identify the variants that are significantly LoF by signaling by with little protein stability effect. These are the red points in the lower facet plots above that are in the upper left quadrant:


    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_21_0.png)
    


Similarly, we can make heatmaps of specifically the contrast:


    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_24_0.png)
    



    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_24_1.png)
    



    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_24_2.png)
    



    
![png](TYK2-Signal-Stability-Contrast_files/TYK2-Signal-Stability-Contrast_24_3.png)
    


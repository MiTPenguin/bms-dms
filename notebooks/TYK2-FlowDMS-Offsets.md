### TYK2 FlowDMS Offsets

For the most recent full TYK2 FlowDMS, we realized that processing the count data per-chunk instead of all at once changed parts of the output, even though only counts within a given chunk are used for positions in that chunk. However, this change _does_ matter for computing the _offset_, which is taken as the `mean(log(count))` within each sample for FlowDMS.

In the plots below, the lines in blue indicate the WT level of each chunk, the red distributions are the stop midpoints of individual chunks, and the black distribution contains all non-stop, non-WT variants. Remember that each blue line is computed from the counts in each chunk - this only changes the offset term but not the input data to each regression:


    
![png](TYK2-FlowDMS-Offsets_files/TYK2-FlowDMS-Offsets_3_0.png)
    


Clearly, the per-sample offset is much closer to the peak in the black distribution than in the per-samplechunk offset. For each version, if you compute the Z-statistic of the difference between each variant (in black or red) and its correspoding blue WT midpoint, you get these distributions:


    
![png](TYK2-FlowDMS-Offsets_files/TYK2-FlowDMS-Offsets_5_0.png)
    


This looks decent in general (and the stops are in the Z range we expect of -10 to -5), but if we zoom in, we see a problem:


    
![png](TYK2-FlowDMS-Offsets_files/TYK2-FlowDMS-Offsets_7_0.png)
    


The Z-statsitics are systematically negative, when (in principle) the average variant should have no effect, and thus this distribution should peak at zero. However, it is clearly closer to zero per-sample instead of per-samplechunk.

To explore further, let's grab position 669, which contains a spike-in, and try different offsets to see what produces the best result.  We tried these seven models:

| Offset | Domain | 
| --- | --- |
| `log(sum(count))` | Per Sample |
| `log(sum(count))` | Per Sample-Chunk |
| `mean(log(count))` | Per Sample |
| `mean(log(count))` | Per Sample-Chunk |
| `log(sum(stop_count))` | Per Sample |
| `log(sum(stop_count))` | Per Sample-Chunk |
| None | NA |


    
    
    |offset                 |  midpoint|
    |:----------------------|---------:|
    |none                   | 0.6528855|
    |mean_log_count_sample  | 0.6636106|
    |log_stop_counts_sample | 0.6856211|
    |log_total_count_chunk  | 0.6982526|
    |log_stop_counts_chunk  | 0.7007422|
    |mean_log_count_chunk   | 0.7338009|
    |log_total_count_sample | 0.7940420|


This recapitulates what we initially observed, namely that `mean(log(count))` per-sample seems to produce the lower result in line with the observed non-WT distribution. Surprisingly, the no-offset model seems to produce an even slightly lower result. Taken together, this warrants further investigation into what offset structure is best for which experimental designs involving flow midpoints.

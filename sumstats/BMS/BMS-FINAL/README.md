### TYK2 DMS Datasets

Each enclosed dataset has the following columns:
  -  `chunk` indicates TYK2 segment
  -  `pos` indicates TYK2 amino acid position
  -  `condition` indicates the treatment condition as `Drug_Concentration`
  -  `aa` indicates the introduced amino acid
  -  `log2FoldChange` is the mutant vs wild-type effect size
  -  `log2StdError` is the mutant vs wild-type standard error
  -  `statistic` is the Wald test statistic defined as `log2FoldChange`/`log2StdError`
  -  `p.value` is the unadjusted p-value
  -  `version` is the analysis pipeline version

### Model Fitting Commands

#### IFN-alpha DMS
```
make RUN_ID=BMS/TYK2-run3 sumstats/TYK2-run3-cleaned.sumstats.tsv
make RUN_ID=BMS/TYK2-run4 sumstats/TYK2-run4-cleaned.sumstats.tsv
make RUN_ID=BMS/TYK2-run7 sumstats/TYK2-run7-cleaned.sumstats.tsv
```

#### IL-23 mini-DMS
```
make RUN_ID=BMS/TYK2-run10-DL41 sumstats/TYK2-run10-DL41-cleaned.sumstats.tsv
make RUN_ID=BMS/TYK2-run10-DL42 sumstats/TYK2-run10-DL42-cleaned.sumstats.tsv
make RUN_ID=BMS/TYK2-run10-DL7 sumstats/TYK2-run10-DL7-cleaned.sumstats.tsv
```

#### IFN-alpha FlowDMS
```
make RUN_ID=BMS/TYK2-FLOW sumstats/TYK2-FLOW-flow.sumstats.tsv
```
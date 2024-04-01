### TYK2 DMS Model Fitting Commands

#### IFN-alpha DMS
```
export RUN_ID=BMS/TYK2-run3
make sumstats/TYK2-run3-combined.sumstats.tsv

export RUN_ID=BMS/TYK2-run4
make sumstats/TYK2-run4-combined.sumstats.tsv

export RUN_ID=BMS/TYK2-run7
make sumstats/TYK2-run7-combined.sumstats.tsv
```

#### IL-23 mini-DMS
```
export RUN_ID=BMS/TYK2-run10-DL41
make sumstats/TYK2-run10-DL41-combined.sumstats.tsv

export RUN_ID=BMS/TYK2-run10-DL42
make sumstats/TYK2-run10-DL42-combined.sumstats.tsv

export RUN_ID=BMS/TYK2-run10-DL7
make sumstats/TYK2-run10-DL7-combined.sumstats.tsv
```

#### IFN-alpha FlowDMS
```
export RUN_ID=BMS/TYK2-FLOW
make sumstats/TYK2-FLOW-flow.sumstats.tsv
```
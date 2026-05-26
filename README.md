# Killifish Transcriptome Aging Clock

Python adaptation of the KillifishAtlas aging clock pipeline for [AAlab](https://www.age.mpg.de/antebi). 

Applies three transcriptomic aging clocks — **BayesAge 2.0**, **Elastic Net (EN)**, and **Principal Component Regression (PCR)** — to *N. furzeri* RNA-seq data and predicts transcriptomic age (tAge) in query samples, [Original work](https://github.com/emkcosta/KillifishAtlas).


---

## Setup
Download KilifishAtlas data from [GEO PRJNA1274512](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308970)

```bash
conda env create -f environment.yml
conda activate killifish-tx-clock
```

---

## Links

[Documents](https://iron-lion.github.io/Killifish-Tx-Clock/)



## Note
The refactoring has been done with the assistance of CLAUDE.

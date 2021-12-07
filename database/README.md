## Generate database for adaptive sequencing

Here we generate a (reduced) database of antimicrobial resistance genes for use in nanopore adaptive sequencing.

- Data source: [card.mcmaster.ca/download](https://card.mcmaster.ca/download)
- Data release metadata:

```
Data    July 2021 release - improved annotation of Neisseria gonorrhoeae macrolide and beta-lactam resistance mutations 3.1.3   JSON, TAB , FASTA   2021-07-05 16:10:04.04555
```


```bash
mmseqs easy-cluster --min-seq-id 0.95 -c 0.8 --cov-mode 1 nucleotide_fasta_protein_homolog_model.fasta cluster tmp

```

Reduces 2,979 sequences to 1,147.




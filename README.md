## README

Below we provide code for non-standard analyses. This repo also contains all data required to reproduce the analyses.

- Putative horizontal gene transfer
- Create database for adaptive sequencing
- Analysis of experiment that switched adaptive sequencing on and off
- Bayesian regression model of target count

First, download the `sequencing` data we [doposited](https://osf.io/wt7gc/) with the Open Science Framework (OSF) under project ID `wt7gc`, and unpack it.


### Samples

- `A2` .. _Raoultella ornithinolytica_
- `B1` .. _Citrobacter freundii_ (same as `A1`, not used here)
- `B2` .. _Citrobacter amalonaticus_


### Putative horizontal gene transfer

Pairwise align isolate genomes and MAGs and parse all genes longer than 1 kb and 99% sequence identity.

```bash
mkdir tmp
```


```python
from glob import glob
from itertools import combinations, islice
import subprocess
import os

import screed
from tqdm import tqdm


g1 = glob('assemblies/*.fasta')
d1 = {os.path.basename(i).split('_')[0]: i for i in g1}

g2 = glob('metagenome/MAGs_uncultured/*.fa')
d2 = {os.path.basename(i).replace('.fa', ''): i for i in g2}
d2 = {k.replace('.', '_'): v for k, v in d2.items()} 

d = d1 | d2 


# Pairwise genome alignment
for i, j in tqdm(combinations(d, 2)):
    name = f'{i}__{j}'  # dunder
    print(name)
    subprocess.run(['nucmer', '-p', f'tmp/{name}', d[i], d[j]])
```

Now retrieve all shared elements with a minimum:

- alignment nucleotide identity (0.999)
- alignment length (1 kb)


```bash
# stackoverflow.com/questions/965053
for i in $(find tmp -name '*.delta'); do
    filename=$(basename -- "$i")
    name="${filename%.*}"
    echo $name
    
    dnadiff -p tmp/${name} -d $i
    # Create a bed for all shared entries with a homology > 99.9%
    cut -f1,2,7,12 tmp/${name}.mcoords | \
    awk -v 'OFS=\t' '{print $4,$1,$2,$3}' | \
    awk '$4>=99.9' > tmp/${name}.bed
done
```


```python
from glob import glob
from itertools import product
import os

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt


# params
minlen = 1000


files = glob('tmp/*.bed')
m = {}  # matrix

for i in files:
    c1, c2 = os.path.basename(i).replace('.bed', '').split('__')
    try:
        df = pd.read_csv(i, sep='\t', header=None)
        cnt = 0
        for _, j in df.iterrows():
            if j[2] - j[1] > minlen:
                cnt += 1
        m[(c1, c2)] = cnt

    except pd.errors.EmptyDataError:
        m[(c1, c2)] = 0
        continue


ix = sorted(set([item for sublist in m.keys() for item in sublist]))
M = np.zeros([len(ix), len(ix)])

for n, i in enumerate(ix):
    for c1, c2 in product([i], ix):
        if not c1 == c2:
            try:
                M[n, ix.index(c2)] = m[(c1, c2)]
            except KeyError:
                M[n, ix.index(c2)] = m[(c2, c1)]


# Zero all entries below the diagonal, bc/ matrix is symmetric
M = np.triu(M)

colors = matplotlib.cm.get_cmap('gray_r')
labels = ['C. freundii', 'R. ornithinolytica', 'C. amalonaticus', 'E. faecium (MAG)', 'S. ureilytica (MAG)']


fig = plt.figure()
fig.set_size_inches(4, 2)
ax = fig.add_subplot(111)
cax = ax.matshow(M, interpolation='nearest', cmap=colors)
fig.colorbar(cax)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

ax.set_xticklabels([''])  # ix
ax.set_yticklabels([''] + labels)  # ix


# plt.show()
plt.savefig('matrix.png', dpi=300)
```


### Create database for adaptive sequencing

```bash
# For details see database/README.md
mmseqs easy-cluster --min-seq-id 0.95 -c 0.8 --cov-mode 1 database/nucleotide_fasta_protein_homolog_model.fasta cluster tmp
```


### Analysis of experiment that switched adaptive sequencing on and off

Aim: Extract all open reading frames marked as antimicrobial resistance gene. Then map reads from the experiment in both conditions ("adaptive", "standard") to them, compare.


```bash
mkdir aln
DB=database/protein_fasta_protein_homolog_model.fasta
for i in A2 B1 B2; do
    prodigal -i assemblies/${i}.fna -a assemblies/${i}.faa -d assemblies/${i}.coding.fna > /dev/null 
    mmseqs easy-search --threads 8 --max-accept 1 --min-seq-id 0.9 --cov-mode 1 -c 0.8 assemblies/${i}.faa $DB assemblies/${i}.m8 assemblies/tmp
    cut -f1 assemblies/${i}.m8 > select
    seqtk subseq assemblies/${i}.coding.fna select > assemblies/${i}.amr.fna
    rm select

    DATA=$(find sequencing/adaptive_on_off -name '*.fastq.gz' | grep $i)
    REF=assemblies/${i}.amr.fna

    mkdir aln/${i}
    for QRY in $DATA; do
        filename=`basename $QRY`
        minimap2 -ax map-ont --secondary=no -t 8 $REF $QRY | grep -v '^@' > aln/${i}/${filename}.sam
    done
done
```

How similar are the coding sequences we annotated as resistance genes to their respective matches in the reduced AMR database we used during adaptive sequencing?

```bash
DB=database/cluster_rep_seq.fasta
for i in A2 B1 B2; do
    QRY=assemblies/${i}.amr.fna
    # --search-type 3 means nucleotide search
    # Note: We search against reduced database, thus more lenient params
    mmseqs easy-search --search-type 3 --threads 8 --max-accept 1 --min-seq-id 0.5 --cov-mode 1 -c 0.5 $QRY $DB assemblies/${i}.amr.card_sim.m8 assemblies/tmp
done
# Sanity check:
# 44 grep '>' B2.amr.fna | wc -l
# 44 wc -l B2.amr.card_sim.m8
```

Integrate data.

```python
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd
import screed
from tqdm import tqdm


def parse_fp(fp):
    '''
    PosixPath('aln/B1/unrejected_B1_1_16h.fastq.gz.sam')

    parse_fp(fp)
    ('standard', 'B1', 1, 11)
    '''
    x = fp.name.strip('h.fastq.gz.sam')
    if 'unrejected' in x:
        group = 'adaptive'
        _, isolate, replicate, hours = x.split('_')
    else:
        group = 'standard'
        isolate, replicate, hours = x.split('_')
    return group, isolate, int(replicate), int(hours)


'''
Look at assembly graph (bandage, assemblies/*.gfa) and mark circular
chromosome, assuming all else is non-chromosomal.

Also, we note the coverage in the assembly, as a proxiy for relative copy number
of the contigs (for use later in a regression model, see misc/cn.json).
'''

# Which contigs belong to the chromosome?
chromosome = {
    'A2': 'contig_1',
    'B1': 'contig_1',
    'B2': 'contig_2',
}

# What is their copy number?
with open('cn.json', 'r') as file:
    cn = json.load(file)


# Collect information on which reads map to which ORF in which condition,
# how long the reads were, etc.
cnts, cnt_lens = {}, {}
read_lens = defaultdict(list)

for sample in ['A2', 'B1', 'B2']:
    cnt = defaultdict(dict)
    cnt_len = defaultdict(dict)
    df = pd.read_csv(f'assemblies/{sample}.m8', sep='\t', header=None)
    amr = set(df[0])  # AMR detected in isolate

    p = Path(f'aln/{sample}')
    for fp in tqdm(p.rglob('*.sam')):

        # How many reads map to AMR genes?
        with open(fp, 'r') as file:
            for line in file:
                l = line.strip().split('\t')
                contig = l[2]
                read_len = len(l[9])

                if (contig != '*') and (contig in amr):
                    
                    group, isolate, replicate, hours = parse_fp(fp)
                    # Read len distr chromosome vs plasmid
                    if '_'.join(contig.split('_')[:2]) == chromosome[sample]:
                        read_lens['chromosome'].append(read_len)
                    else:
                        read_lens['plasmid'].append(read_len)

                    try:
                        cnt[group][contig] += 1
                        cnt_len[group][contig].append(read_len)

                    except KeyError:
                        cnt[group][contig] = 1               # initialize
                        cnt_len[group][contig] = [read_len]  # initialize
    cnts[sample] = cnt
    cnt_lens[sample] = cnt_len


# Load similarity vals we calculated above for the ORFs.
sims = {}
for sample in ['A2', 'B1', 'B2']:
    sim = {}
    fp = f'assemblies/{sample}.amr.card_sim.m8'
    df = pd.read_csv(fp, sep='\t', header=None)
    for _, i in df.iterrows():
        sim[i[0]] = i[2]
    sims[sample] = sim


# Each line one gene (coding sequence) with properties of the same.
with open('misc/similarity.csv', 'w+') as out:
    out.write('sample,condition,contig,count,similarity,carrier,mu,coverage\n')
    
    # results = []
    for sample in ['A2', 'B1', 'B2']:
        for condition in ['standard', 'adaptive']:
            for k, v in cnts[sample][condition].items():
                # 'contig_2_3755': 43,
                coding = '_'.join(k.split('_')[:2])  
                # contig_2_34 to contig_2
                if coding == chromosome[sample]:
                    carrier = 'chromosome'
                else:
                    carrier = 'plasmid'

                cov = cn[sample][coding]  # coverage
                mu = np.median(cnt_lens[sample][condition][k])
                s = sims[sample][k]
                out.write(f'{sample},{condition},{k},{v},{s},{carrier},{mu},{cov}\n')
```

Visualize.

```r
library(ggplot2)
library(readr)
options(scipen=999)
library(dplyr)


df <- read_csv('misc/similarity.csv')
df$log_mu = log(df$mu)


p <- ggplot(df, aes(x=carrier, y=count, color=similarity)) + geom_jitter(size=.8) + theme_minimal() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + ylab('read count') + scale_colour_viridis_c() + facet_wrap(~condition)
ggsave('plot3.pdf', p, width=7, height=10, units='cm')
# Figure 2B

palette <- 'Accent'
p <- ggplot(df, aes(x=carrier, y=log_mu, color=condition)) + geom_jitter(size=.75, position = position_jitterdodge()) + geom_boxplot(alpha=0.75, outlier.shape=NA) + theme_minimal() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + ylab('log median read length') + scale_color_brewer(palette=palette)
# Figure 2D
```


### Bayesian regression model of target count

```bash
# Pull container with all dependencies, enter it and start R session
docker pull lcolling/brms
docker run --rm -it -v $PWD:/data lcolling/brms /bin/bash
R
```

Fun.

```r
library(brms)

library(ggplot2)
library(readr)
options(scipen=999)
library(dplyr)


df <- read_csv('/data/misc/similarity.csv')
df$adaptive <- ifelse(df$condition == 'adaptive', 1, 0)
df$plasmid <- ifelse(df$carrier == 'plasmid', 1, 0)
df$log_mu = log(df$mu)

# Standardize coverage predictor to make model fit easier
df2 <- df %>% mutate_at(c('coverage'), ~(scale(.) %>% as.vector))


m <- brm(
    data=df2,
    family=poisson(),
    count ~ 1 + similarity*adaptive + log_mu*plasmid + log_mu*adaptive + coverage
    )
# posterior_summary(m)


# Merge posterior samples with predictors for visualisation.
new <- fitted(m, newdata=df2) %>% as_tibble() %>% bind_cols(df2)


palette <- 'Accent'
p <- ggplot(new, aes(x=similarity, y=count, color=condition)) + geom_point(size=.8) +geom_smooth(data=new, aes(y=Estimate, ymin=Q2.5, ymax=Q97.5, fill=condition)) + theme_minimal() + scale_color_brewer(palette=palette) + scale_fill_brewer(palette=palette) + ylab('read count')
# Figure 2C

p <- ggplot(new, aes(x=log_mu, y=count, color=condition)) + geom_point(size=.8) +geom_smooth(data=new, aes(y=Estimate, ymin=Q2.5, ymax=Q97.5, fill=condition)) + theme_minimal() + xlab('log median read length') + scale_color_brewer(palette=palette) + scale_fill_brewer(palette=palette) + ylab('read count')
# Figure 2E
```

Model summary:

```
 Family: poisson
  Links: mu = log
Formula: count ~ 1 + similarity * adaptive + log_mu * plasmid + log_mu * adaptive + coverage
   Data: df2 (Number of observations: 238)
Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup samples = 4000

Population-Level Effects:
                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept               5.77      0.14     5.50     6.04 1.00     1999     2395
similarity             -0.09      0.12    -0.32     0.13 1.00     2407     2466
adaptive               -9.65      0.18   -10.01    -9.30 1.00     1540     1740
log_mu                 -0.03      0.01    -0.04    -0.02 1.00     1856     2444
plasmid                -0.41      0.07    -0.54    -0.28 1.00     1516     2096
coverage                0.05      0.00     0.04     0.05 1.00     4418     2890
similarity:adaptive    11.98      0.15    11.68    12.28 1.00     1665     1904
log_mu:plasmid          0.06      0.01     0.05     0.08 1.00     1456     1861
adaptive:log_mu        -0.22      0.01    -0.23    -0.20 1.00     2044     2000

Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

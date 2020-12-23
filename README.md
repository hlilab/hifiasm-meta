(WIP) A hifiasm fork (https://github.com/chhylp123/hifiasm) for metagenome assembly.

## Getting Started
```sh
# Install hifiasm-meta (g++ and zlib required)
git clone https://github.com/xfengnefx/hifiasm-meta.git
cd hifiasm-meta && make

# Run
hifiasm_meta -t32 -oasm reads.fq.gz 2>asm.log
hifiasm_meta -t32 -S -o asm reads.fq.gz 2>asm.log // if the dataset has high redundancy, or overlap & error correction takes way too long
```

## Output files

Raw unitig graph: asm.r\_utg\*.gfa

Cleaned unitig graph: asm.p\_utg\*.gfa 

Contig graph: asm.p\_ctg\*.gfa and asm.a\_ctg\*.gfa

Unitig/Contig naming: `^s[0-9]+\.[uc]tg[0-9]{6}[lc]` where the `s[0-9]+` is a subgraph label.

## Special Notes

Based on the limited available test data, real datasets are unlikely to require read selection; mock datasets, however, might need it.

Non-release commits may have extra debug outputs for dev/debug purposes, even without -V.

Bin file is one-way compatible with the stable hifiasm for now: stable hifiasm can use hifiasm\_meta's bin file, but not vice versa. Meta needs to store extra info from overlap & error correction step.

## Switches

See also README\_ha.md, the stable hifiasm doc.

```
#Interface
-B		Name of bin files. Allows to use bin files from other 
       		directories.

# Read selection
-S		Enable read selection.
--force-preovec Force kmer frequency-based read selection. 
                (otherwise if total number of read overlaps 
                 look realistic, won't do selection.)
--lowq-10Lower  10% quantile kmer frequency threshold, runtime. Lower value means less reads kept, if read selection is triggered. [150]

```

## Known issues

Read selection needs to speed up. Currently there's a blocking sequential step.

## Preliminary results

[Sheep fecal material dataset](https://www.ncbi.nlm.nih.gov/sra/SRX7628648[accn]): wall clock 17.2 h on 48 cpus, peak memory 183.3 GB.

A [Bandage](https://github.com/rrwick/Bandage) plot of primary contig graph:

<p align="center">
  <img src="https://user-images.githubusercontent.com/61363437/102846924-4f2c3f80-43df-11eb-8c00-e84b439a9398.png"/>
</p>

**<sub>Evaluation tool</sub>**|**<sub>Criteria</sub>**|**<sub>Count</sub>**|**<sub>Bases or relative ratio</sub>**
-----|-----|-----|-----
<sub>custom</sub>|<sub>contig >1Mb</sub>|<sub>320</sub>|<sub>638Mb</sub>
<sub></sub>|<sub>contig >100kb</sub>|<sub>3307</sub>|<sub>1.30Gb</sub>
<sub></sub>|<sub>circular contig >1Mb</sub>|<sub>147</sub>|<sub>340Mb</sub>
<sub>[Barrnap][ubarrnap]</sub>|<sub>contig with all three types of rRNA</sub>|<sub>1012</sub>|<sub>643Mb</sub>
<sub>[CheckM][ucheckm]</sub>|<sub>genome completeness >90%</sub>|<sub>177</sub>|<sub>423Mb</sub>
<sub></sub>|<sub>^~ contamintation  >5%</sub>|<sub>2</sub>|<sub>1.10%</sub>
<sub></sub>|<sub>genome completeness >25%</sub>|<sub>444</sub>|<sub>706MB</sub>
<sub></sub>|<sub>^~ contamintation  >5%</sub>|<sub>5</sub>|<sub>1.10%</sub>
<sub>[CheckV][ucheckv]</sub>|<sub>high quality virus genome\*</sub>|<sub>186</sub>|<sub>10Mb</sub>
<sub>[prodigal][uprodigal]</sub>|<sub>genes predicted</sub>|<sub>27039</sub>|<sub>-</sub>

[ubarrnap]: https://github.com/tseemann/barrnap
[ucheckm]: https://github.com/Ecogenomics/CheckM
[ucheckv]: https://bitbucket.org/berkeleylab/checkv/src
[uprodigal]: https://github.com/hyattpd/Prodigal

[Mock community ATCC MSA-1003](https://www.ncbi.nlm.nih.gov/sra/SRX8173258[accn]) (with -S --lowq-10 50)

**<sub>Strain</sub>**|**<sub>Abundance</sub>**|**<sub>Assembly status</sub>**| |**<sub>Strain</sub>**|**<sub>Abundance</sub>**|**<sub>Assembly status</sub>**
-----|-----|-----|-----|-----|-----|-----
<sub>Acinetobacter baumannii</sub>|<sub>0.18%</sub>|<sub>pass</sub>||<sub>Lactobacillus gasseri</sub>|<sub>0.18%</sub>|<sub>pass</sub>
<sub>Bacillus cereus</sub>|<sub>1.80%</sub>|<sub>pass</sub>||<sub>Neisseria meningitidis</sub>|<sub>0.18%</sub>|<sub>pass</sub>
<sub>Bacteroides vulgatus</sub>|<sub>0.02%</sub>|<sub>fragmented</sub>||<sub>Porphyromonas gingivalis</sub>|<sub>18.00%</sub>|<sub>pass\*</sub>
<sub>Bifidobacterium adolescentis</sub>|<sub>0.02%</sub>|<sub>lost</sub>||<sub>Pseudomonas aeruginosa</sub>|<sub>1.80%</sub>|<sub>pass</sub>
<sub>Clostridium beijerinckii</sub>|<sub>1.80%</sub>|<sub>pass</sub>||<sub>Rhodobacter sphaeroides</sub>|<sub>18.00%</sub>|<sub>pass</sub>
<sub>Cutibacterium acnes</sub>|<sub>0.18%</sub>|<sub>pass</sub>||<sub>Schaalia odontolytica</sub>|<sub>0.02%</sub>|<sub>lost</sub>
<sub>Deinococcus radiodurans</sub>|<sub>0.02%</sub>|<sub>fragmented</sub>||<sub>Staphylococcus aureus</sub>|<sub>1.80%</sub>|<sub>pass</sub>
<sub>Enterococcus faecalis</sub>|<sub>0.02%</sub>|<sub>fragmented</sub>||<sub>Staphylococcus epidermidis</sub>|<sub>18.00%</sub>|<sub>pass</sub>
<sub>Escherichia coli</sub>|<sub>18.00%</sub>|<sub>pass</sub>||<sub>Streptococcus agalactiae</sub>|<sub>1.80%</sub>|<sub>pass\*</sub>
<sub>Helicobacter pylori</sub>|<sub>0.18%</sub>|<sub>pass</sub>||<sub>Streptococcus mutans</sub>|<sub>18.00%</sub>|<sub>pass</sub>

\*: not one clean circle.

[Mock community Zymo D6331, standard input library](https://www.ncbi.nlm.nih.gov/sra/SRX9569057[accn]): wall clock 15.7 h on 32 cpus, peak memory 121.7 GB.

**<sub>Strains</sub>**|**<sub>Abundance</sub>**|**<sub>Assembly status</sub>**|**<sub></sub>**|**<sub>Strains</sub>**|**<sub>Abundance</sub>**|**<sub>Assembly status</sub>**
-----|-----|-----|-----|-----|-----|-----
<sub>Akkermansia muciniphila</sub>|<sub>1.36%</sub>|<sub>pass</sub>|<sub></sub>|<sub>Escherichia coli JM109</sub>|<sub>8.37%</sub>|<sub>unseparated\*</sub>
<sub>Bacteroides fragilis</sub>|<sub>13.13%</sub>|<sub>pass</sub>|<sub></sub>|<sub>Faecalibacterium prausnitzii</sub>|<sub>14.39%</sub>|<sub>pass</sub>
<sub>Bifidobacterium adolescentis</sub>|<sub>1.34%</sub>|<sub>pass</sub>|<sub></sub>|<sub>Fusobacterium nucleatum</sub>|<sub>3.78%</sub>|<sub>pass</sub>
<sub>Candida albican</sub>|<sub>1.61%</sub>|<sub>fragmented</sub>|<sub></sub>|<sub>Lactobacillus fermentum</sub>|<sub>0.86%</sub>|<sub>pass</sub>
<sub>Clostridioides difficile</sub>|<sub>1.83%</sub>|<sub>pass</sub>|<sub></sub>|<sub>Methanobrevibacter smithii</sub>|<sub>0.04%</sub>|<sub>3contigs</sub>
<sub>Clostridium perfringens</sub>|<sub>0.00%</sub>|<sub>lost</sub>|<sub></sub>|<sub>Prevotella corporis</sub>|<sub>5.37%</sub>|<sub>partial</sub>
<sub>Enterococcus faecalis</sub>|<sub>0.00%</sub>|<sub>lost</sub>|<sub></sub>|<sub>Roseburia hominis</sub>|<sub>3.88%</sub>|<sub>pass</sub>
<sub>Escherichia coli B1109</sub>|<sub>8.44%</sub>|<sub>unseparated\*</sub>|<sub></sub>|<sub>Saccharomyces cerevisiae</sub>|<sub>0.18%</sub>|<sub>fragmented</sub>
<sub>Escherichia coli b2207</sub>|<sub>8.32%</sub>|<sub>pass</sub>|<sub></sub>|<sub>Salmonella enterica</sub>|<sub>0.02%</sub>|<sub>unseparated\*</sub>
<sub>Escherichia coli B3008</sub>|<sub>8.25%</sub>|<sub>unseparated\*</sub>|<sub></sub>|<sub>Veillonella rogosae</sub>|<sub>11.02%</sub>|<sub>pass</sub>
<sub>Escherichia coli B766</sub>|<sub>7.83%</sub>|<sub>pass</sub>|<sub></sub>|<sub></sub>|<sub></sub>|<sub></sub>

\*: E.coli strains except B2207 and B766 presented in one subgraph.


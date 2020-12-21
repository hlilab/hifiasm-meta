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

## Current output files

Raw unitig graph: asm.r\_utg\*.gfa

Cleaned unitig graph: asm.p\_utg\*.gfa 

Contig graph: asm.p\_ctg\*.gfa and asm.a\_ctg\*.gfa

## Special Notes

Based on the limited available test data, real datasets are unlikely to require read selection; mock datasets, however, might need it.

Non-release commits may write excessive info to STDERR (and gfa files) for dev/debug purposes, even without -V, please pipe to gzip if this is a concern.

Bin file is one-way compatible with the stable hifiasm for now: stable hifiasm can use hifiasm\_meta's bin file, but not vice versa. Meta needs to store extra info from overlap & error correction step.

## Switches (see also README\_ha.md)

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


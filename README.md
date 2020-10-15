(WIP) A hifiasm fork (https://github.com/chhylp123/hifiasm) for metagenome assembly.

## Getting Started
```sh
# Install hifiasm-meta (g++ and zlib required)
git clone https://github.com/xfengnefx/hifiasm-meta.git
cd hifiasm-meta && make

# Run
hifiasm_meta -t32 --force-preovec -oasm reads.fq.gz 2>asm.log

# (debug) include previous commit ID to bin files
(after git commit, do)
./stamp.sh; make
```

## Special Notes

Bin file format has been slightly altered because of the need for auxiliary info etc, currently it's not compatible with the stable hifiasm release. This is temporary.

## Switches (see also README\_ha.md)

```
#Interface
-B		Name of bin files. Allows to use bin files from other 
       		directories (and write to the destination specified by -o). 
	        Use -o otherwise.

# Read selection
--preovec	Enable the current read selection strategy. Note that 
		if the total number of possible overlaps appears acceptable, 
		read selection will note be triggered. 
		To negate this check, use --force-preovec instead.
-X		Disable all read selection.
--force-preovec 	Force to do the read selection.
--preovec-coverage	Median kmer frequency threshold, runtime; semi deprecated. [150] 
--lowq-10		Lower 10% quantile kmer frequency threshold, runtime. [150]
--lowq-5		Lower 5% quantile kmer frequency threshold, runtime. [-1 (disabled)]
--lowq-3		Lower 3% quantile kmer frequency threshold, runtime. [-1 (disabled)]

# Graph cleaning
--exp-graph-cleaning	Enable experimental treatments (e.g. topo-aware coverage-based 
			arc drop). Basic routines are not affected by this switch.
```


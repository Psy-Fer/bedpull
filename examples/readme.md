

command for extracting subsequence from fasta and getting query coordinates


First, extract the paternal genome assembly from the diploid assembly
You can download the assembly from here

https://github.com/marbl/hg002


```
grep "_PATERNAL" hg002v1.1.fasta.fai | cut -f1 | xargs samtools faidx hg002v1.1.fasta > hg002.paternal_genome.fasta
```

```
bedpull --paf hg002pat_to_hs1.rfc1_only.paf --query_ref hg002.paternal_genome.fasta -r rfc1.bed -o pat_paf_extracted.fasta
```

You should get something like this

```
>chr4_PATERNAL|chr4:39318077-39318136|RFC1|chr4_PATERNAL:39438031-39438610|+
AAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAA
```

(the header format shown reflects the current 5-field structure; `--bed_out`, `--dedup`,
`--hap_split`, `--debug`, and CRAM mode were all added after this example was first written —
see the [mdBook](https://psy-fer.github.io/bedpull/) for those.)

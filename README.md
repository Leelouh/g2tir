# g2tir.pl

This script weas used in the publication (Helou _et al_, 2020) [work in progress].

This tool allow to search sequences flanked by TIRs from a genome/fasta file.

## Research of sequences flanked by TIRs

```
usage : g2tir.pl [--motif IUPAC] [--in file.fasta] [--out output_file]
```

**Required arguments**
```
  --motif IUPAC             Motif for the TIR in IUPAC notation (ex : TTAANNCMK)
  or
  --R1 motif1 --R2 motif2   Motif for TIR1 and TIR2 in iupac format if they are different
  
  --in IN                   Fasta file to search sequences flanked by TIR
  --out OUT                 Output file with the results in table format
```

**Optional arguments**
```
    --score INT       Score between two TIRs based on Levenshtein::XS (INT default : 10) (--score 0 for perfect TIR)
    --printScore      Print the score in the output file
    --maxL INT        Maximum distance between two TIRs (INT default : 6400 nt)
    --minL INT        Minimum distance between two TIRs (INT default : 43 nt)
    --pm              Print motifs of TIR founded in the output file
    --TSD             Print the TSD (4nt before each TIR) in the output file
    --degen           Degenerate the motif submitted (used the option --score to keep only similar TIRs)
    --fasta           Generate a fasta file with the sequences flanked by TIRs identified
    --config          Generate a config file with the option selected
```

**Exemple**
```
./g2tir.pl --score 0 --fasta --min 30  --printScore --config --pm --motif TTAANNN --in /data/hg38/chr/chr19.fa --out output_chr19.txt
```

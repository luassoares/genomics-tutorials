# 08 — Dfoil: direction of introgression

Test for introgression **and its direction** in five-taxon trees with the
D<sub>FOIL</sub> statistics, applied to all possible combinations of
individuals with the **ExDFOIL** scripts.

This tutorial follows ExDFOIL (Lambert et al. 2019), which was written for
RADseq data:
<https://github.com/slambert-git/ExDFOIL>
DFOIL: <https://github.com/jbpease/dfoil>

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/pick_individuals.py`](scripts/pick_individuals.py) | Chooses the N individuals with the least missing data in each population |
| [`scripts/dfoil_summary.R`](scripts/dfoil_summary.R) | Classifies the results by population and counts the tests in each class |

The ExDFOIL and DFOIL scripts are downloaded from their repositories (see
Install).

## Install

```bash
conda create -n dfoil -c conda-forge -c bioconda python=3.11 numpy scipy matplotlib \
    parallel perl vcftools perl-vcftools-vcf r-base r-ape r-phytools r-stringr
conda activate dfoil

cd ~/programs
git clone https://github.com/jbpease/dfoil.git
git clone https://github.com/slambert-git/ExDFOIL.git
```

- `perl-vcftools-vcf` provides `vcf-to-tab` (the `vcftools` package alone
  does not have it).
- The ExDFOIL scripts must be **in the folder where you run them** (the
  README of ExDFOIL says: use a clean folder). They are copied to
  `results/` in step 4.

## Prepare the folder

```
08_dfoil/
├── data/      my_data-out-miss0.9-maf0.01-th100.recode.vcf  (with outgroup, from 03_filtering)
│              popmap.txt
│              a rooted tree of the individuals with branch lengths (e.g. from 02_raxml)
├── scripts/   pick_individuals.py, dfoil_summary.R
└── results/
```

---

## 1. Background

The D statistic ([06_dsuite](../06_dsuite/)) uses four taxa and shows
**that** two lineages share an excess of alleles, but not the direction of
gene flow. D<sub>FOIL</sub> uses **five taxa** in a symmetric tree:

```
        ┌── P1
     ┌──┤
     │  └── P2        P1 and P2 split more recently than P3 and P4
  ┌──┤
  │  │  ┌── P3
  │  └──┤
  │     └── P4
  └──────── O  (outgroup)
```

It calculates four D statistics (D<sub>FO</sub>, D<sub>IL</sub>,
D<sub>FI</sub>, D<sub>OL</sub>). The combination of their signs (positive,
negative or zero) tells which taxa exchanged genes and in which direction:

| Result | Meaning |
|---|---|
| `13`, `14`, `23`, `24` | from P1 (or P2) into P3 (or P4), e.g. `13` = P1 → P3 |
| `31`, `41`, `32`, `42` | from P3 (or P4) into P1 (or P2), e.g. `31` = P3 → P1 |
| `123`, `124` | between the **ancestor** of P1 and P2 and P3 (or P4); direction not known |
| `none` | no significant introgression |
| `na` | the signs do not match any expected pattern (not interpretable) |

**ExDFOIL** runs DFOIL on all combinations of individuals that fit the
symmetric tree (taken from a rooted tree of the individuals) and then
summarizes how many tests show each type of introgression.

> **Rules of the method.** The P1–P2 split must be **younger** than the
> P3–P4 split. ExDFOIL chooses the combinations and the order of the taxa
> from the branch lengths of your tree, so the tree must be **rooted** and
> have **branch lengths**.

## 2. Choose the individuals

The number of tests grows very fast with the number of individuals: in the
original analysis, 132 individuals gave 3,873,154 combinations, and 64
individuals gave 181,919. Choose a few individuals per population (e.g. 2–5),
the ones with the least missing data:

```bash
cd results
vcftools --vcf ../data/my_data-out-miss0.9-maf0.01-th100.recode.vcf --missing-indv --out my_data
python3 ../scripts/pick_individuals.py ../data/popmap.txt my_data.imiss 4 > names.txt
```

The table printed on the screen shows the chosen individuals and their
missing data. `names.txt` has one individual per line, **without** the
outgroup.

The outgroup goes in a separate file with **one** outgroup individual (the
one with the least missing data):

```bash
echo "1201" > outgroup.txt      # the ID of one outgroup individual
```

## 3. FASTA alignment

The VCF is converted to a table and then to a FASTA alignment of the SNPs,
with the ExDFOIL version of the script `vcf_tab_to_fasta_alignment.pl` (it
reads missing data coded as `./.`):

```bash
cat ../data/my_data-out-miss0.9-maf0.01-th100.recode.vcf | vcf-to-tab > snps.tab
perl ~/programs/ExDFOIL/AppendixS2_Materials/SingleSite/vcf_tab_to_fasta_alignment_sml.pl \
    -i snps.tab > all_snps.fasta
grep -c ">" all_snps.fasta        # number of individuals
```

The heterozygous sites are written as IUPAC codes. The names in the FASTA
are the VCF sample IDs: the tree must use the **same names**.

## 4. Choose the tests

Copy the ExDFOIL and DFOIL scripts to `results/` and make them executable:

```bash
cp ~/programs/ExDFOIL/Scripts/* .
cp ~/programs/dfoil/*.py .
chmod +x *.sh *.pl *.py *.R
```

**The tree.** Use a tree of the individuals with branch lengths, e.g. the
RAxML-ng best tree ([02_raxml](../02_raxml/)), **rooted on the outgroup**.
To root it in R:

```r
library(ape)
tree <- read.tree("../../02_raxml/results/tree.raxml.bestTree")
tree <- root(tree, outgroup = c("1201", "1202", "1203"), resolve.root = TRUE)
write.tree(tree, "tree_rooted.nwk")
```

(The tree can have more individuals than `names.txt`: the others are
removed.)

Then:

```bash
Rscript DFOIL_Picker.R names.txt tree_rooted.nwk
wc -l myoutput.txt
head -3 myoutput.txt
```

`myoutput.txt` has one test per line, four individuals in the order
P1 P2 P3 P4 (the younger pair first):

```
D_1 E_1 A_1 C_1
D_2 E_1 A_1 C_1
```

> The ExDFOIL README and `run_ExDFOIL.sh` call this file `DFOIL_picked.txt`,
> but the current `DFOIL_Picker.R` writes **`myoutput.txt`**.

To test only some hypotheses, filter the lines, e.g. only the tests that
include an individual of `PopB`: `grep "PopB_" myoutput.txt > tests.txt`
(this works if the sample names include the population).

## 5. Run DFOIL

Run inside byobu ([00_linux_basics](../../00_linux_basics/)): with many
tests, this takes hours.

```bash
# 1. site-pattern counts for each test -> counts/
./fasta2foiler.sh myoutput.txt all_snps.fasta outgroup.txt

# 2. DFOIL on each count file (mode dfoilalt) -> dfoil/ and precheck/
./dfoiler_alt.sh

# 3. one summary line per test -> summary/
parallel './summarizer_alt.sh < <(echo {})' ::: counts/*.counts
cat summary/*.summary_alt > all.summary.alt.txt
wc -l all.summary.alt.txt
```

- `fasta2foiler.sh` uses all the processors. To limit them, edit the script
  and add `-j 16` after `parallel` (as in the original analysis).
- `dfoiler_alt.sh` runs DFOIL with `--mode dfoilalt`, which ignores the
  singleton patterns. The authors recommend it because uneven error or
  rates between samples can create false signals in the singletons.
- `SyntaxWarning: "is" with a literal` messages from `fasta2dfoil.py` can be
  ignored.

Each line of `all.summary.alt.txt` has: the four individuals, four counts,
the result (`none`, `123`, ...), and each statistic with its p-value:

```
A_1 A_2 D_1 E_1 145 126 332 308 none -0.038 0.69 0.077 0.43 -0.077 0.58 0.15 0.27
C_1 C_2 A_1 B_2 68 67 216 120 124 ...
```

## 6. Summarize by population

Add the population of each individual with `associate.sh`. It needs a
sample information file with three columns, separated by **spaces**:
`<sample ID> <population> <species>` (if you have no species, repeat the
population):

```bash
awk '{print $2, $1, $1}' ../data/popmap.txt > sampleinfo.txt
./associate.sh all.summary.alt.txt sampleinfo.txt
ls all.summary.alt_appended.txt
```

Then classify and count:

```bash
Rscript ../scripts/dfoil_summary.R all.summary.alt_appended.txt
```

→ `freq_introgression.csv` and a table on the screen. The classes read
as:

| Class | Meaning |
|---|---|
| `PopA -> PopC` | gene flow from PopA into PopC (the arrow gives the direction) |
| `PopA / PopB <-> PopC` | between the ancestor of PopA and PopB (P1 and P2) and PopC |
| `none` | no introgression (this includes `na`) |

Example (original analysis, 64 individuals, 181,919 tests):

| Class | Number of tests |
|---|---|
| none | 178,921 |
| Psch2 / Pgua <-> Psaxi | 577 |
| Pgua2 / Pgua <-> Psaxi | 546 |
| Pgua2 / Psch2 <-> Psaxi | 366 |
| Pgua2 / Pgua <-> Preit | 299 |
| Psch2 / Pgua <-> Preit | 241 |
| ... | |

Only 2,998 tests (1.6%) showed introgression; the most frequent classes
involve the ancestor of the Pgua / Pgua2 / Psch2 clade and Psaxi or Preit.

## 7. Interpretation

- **Most tests are `none`.** With RADseq SNPs, each test has few
  informative sites (look at the counts in columns 5–8), so the power of a
  single test is low. Look at the **patterns across many tests**: which
  pairs of populations come out again and again, and in which direction.
- **Ancestral signals** (`123`, `124`) are common when gene flow happened
  before the P1–P2 split.
- **Check the counts.** Results from tests with very few informative sites
  are not reliable.
- Lambert et al. (2019) show the power and the pitfalls of DFOIL with RAD
  data (missing data, batch effects, singletons): read it before
  interpreting the results.
- Compare with the other introgression tests: HyDe ([07](../07_hyde/)),
  Dsuite ([06](../06_dsuite/)), ADMIXTOOLS ([09](../09_admixtools/)).

## Next step

[09 — ADMIXTOOLS 2](../09_admixtools/): admixture graphs.

## References

- Pease JB, Hahn MW (2015) Detection and polarization of introgression in a
  five-taxon phylogeny. *Systematic Biology* 64: 651–662.
- Lambert SM, Streicher JW, Fisher-Reid MC, et al. (2019) Inferring
  introgression using RADseq and D<sub>FOIL</sub>: power and pitfalls
  revealed in a case study of spiny lizards (*Sceloporus*). *Molecular
  Ecology Resources* 19: 818–837.

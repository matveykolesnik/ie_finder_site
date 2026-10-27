# ie_finder_site

A reduced, self-contained copy of the MGE_finder discovery workflow. The integrase HMMs, the search scripts, and the conda environment are in this repository. `UPSTREAM` records which commit of [MGE_finder](https://github.com/rljech13/MGE_finder) the copy was taken from. A checkout of that repository is not required.

The pipeline reads a genome assembly in FASTA and writes attL/attR, the integrase, and the element span in the coordinates of that assembly.

## What an element is

A confident element is a phage-type tyrosine integrase and a tRNA, with a direct repeat of that tRNA 3′ end. The integrase sits at one end of the element: next to the tRNA on the opposite strand, their 3′ ends at most `trna_max_distance_bp` (700 bp; the published finder used 500) apart, or next to the repeat (attL) on the tRNA's strand, its 3′ end within `trna_max_distance_bp` of attL — the layout of the tailless IN93/P23-77-like Thermus phages, which the published rule never makes candidates. Setting `integrase_at: trna` keeps only the first, as the published finder did. An integrase at the attL end is paired with the nearest tRNA on its 3′ side whose 3′ end it finds copied next to itself (`scripts/attl_side.py`); after that, the candidate goes through the same filters as any other. The integrase protein is longer than 300 aa, attL does not overlap a CDS, and the repeat is either an exact run of at least 14 bp or a BLAST hit of at least 17 bp. Alignment gaps in the repeat are allowed. Setting `attl_reject_gapped: true` rejects them, as the published finder did. Among the BLAST hits anchored at the tRNA 3′ end, attL is the one with the highest bitscore, so a long but poor match far from the tRNA does not win over a near-perfect repeat. Setting `attl_select_by: length` takes the longest hit instead, as the published finder did. An ORF is an integrase when it reaches the Pfam gathering threshold of PF00589 or PF22022. Setting `integrase_hmm_cutoff: none` keeps every hit HMMER reports instead, as the published finder did. Integrase length is the predicted protein without its stop codon; setting `integrase_len_from: span` divides the ORF span by 3 instead, which counts the stop codon, as the published finder did. Those thresholds are in `search_params.yaml`. Cohort deduplication from the main repository is not part of this copy.

## Tandems

Elements can integrate one after another at the same tRNA. Each integration leaves a copy of the tRNA 3′ end at the element's far end, and the next element takes the same site, so a tandem reads, outwards from the tRNA:

    tRNA - element 1 - copy - element 2 - copy ...

**How the cascade alone gets tandems wrong.** attL is the best-scoring copy, so a tandem comes out one of two ways:

- **The outer copy wins.** Elements 1 and 2 are called as one.
- **The inner copy wins.** Element 2 is left out. Its integrase sits next to a copy, not a tRNA gene, so it is never paired. At these inner junctions bakta often annotates a pseudo-tRNA, which Aragorn does not call.

Tandems whose elements were each found from an integrase of their own were already resolved (`attR_after`).

**How the tandem search fixes it** (`tandem_search: true`, the default; `scripts/tandems.py`):

1. **Find the copies.** For each tRNA with confident elements, the copies of its 3′ end are found again, from the tRNA to `tandem_window_bp` (60 kb) past the outermost attL.
   - A copy must pass the attL gates of the cascade: 3′-anchored, exact run or BLAST length, intergenic, gap rule.
   - It must also be at least `tandem_min_att_bp` (30) long and `tandem_min_att_pident` (80%) identical to the tRNA.
   - A copy within 1 kb of another junction is a duplicated att, not a junction.
2. **Find the members.** A stretch between two junctions is an element when it holds an element's integrase, or an integrase that faces one of its ends as the cascade requires:
   - the near copy, from the opposite strand, 3′ ends within `trna_max_distance_bp`;
   - with `integrase_at: either`, the far copy, from the tRNA's strand.
3. **Split and add.**
   - An element is split only where both parts are elements. A stretch without an integrase stays with the element before it.
   - Past the outermost attL, stretches are added while each is an element.

**How members are bounded.** Each member runs from the copy before it (attR; for the first, the part of the tRNA aligned to its attL) to its own copy (attL), so neighbours share a copy.

`tandem_search: false` gives the output of the finder without the search, byte for byte.

## Requirements

Linux. The first `./run.sh` builds the conda environment named in `envs/IE_finder.yaml` (`IE_finder_site`): Python 3.10, Biopython, pandas, PyYAML, pyrodigal (Prodigal's gene finder as a Python library), HMMER, Aragorn, and BLAST+. Later runs reuse that environment. `run.sh` only sets up that environment; the pipeline itself is `scripts/ie_finder.py`. The first build needs a network connection.

If `conda` is already installed, the script uses it. Otherwise it installs Miniforge into `~/miniforge3` (`MINIFORGE_ROOT` changes that directory) and creates the environment there.

## Install

```bash
git clone https://github.com/rljech13/ie_finder_site.git
cd ie_finder_site
```

## Run

One assembly, or a directory of assemblies. The sample name is the file name without its last suffix, so name the file after the strain.

```bash
./run.sh TTHB27c.fasta outdir
./run.sh genomes_dir outdir
```

Accepted suffixes in a directory are `.fna`, `.fa`, `.fasta`, and `.fsa`. Two files that share a stem are an error.

`outdir/` receives two files per assembly:

| file | role |
|---|---|
| `TTHB27c.ie.gff3` | features on the input contigs, for a browser track or a download |
| `TTHB27c.ie.report.txt` | one report: counts, why a candidate was rejected, the published coordinates, and any step log that was not empty |

With `tandem_search: true` every element is also written as a sequence of its own, sorted by whether it is in a tandem, and once more, all together, in `all/`:

```
outdir/
├── TTHB27c.ie.gff3, TTHB27c.ie.report.txt
├── single_elements/              elements that are not in a tandem
│   └── TTHB27c_IE1.fna, TTHB27c_IE1.gbk, ...
├── tandem_elements/
│   ├── full/                     each tandem whole, from the tRNA to its last member's attL
│   │   └── TTHB27c_TA1.fna, TTHB27c_TA1.gbk, ...
│   └── split/                    each member of a tandem alone
│       └── TTHB27c_IE2.fna, TTHB27c_IE2.gbk, ...
├── all/                          every element above: single ones and tandem members
│   └── TTHB27c_IE1.gbk, TTHB27c_IE2.gbk, ...
└── truncated_elements/           with TRUNCATED_ELEMENTS=1: elements whose attL search window ran into a contig end
    ├── single_elements/
    │   └── TTHB27c_IE3.fna, TTHB27c_IE3.gbk, ...
    └── tandem_elements/
        ├── full/
        └── split/
```

- **One FASTA and one GenBank file per sequence**, named after the sample and the element or tandem. Genomes share these directories, so the sample name keeps the files apart.
- **`all/` holds GenBank files only**: the files of `single_elements/` and `tandem_elements/split/` again, one per element, for tools that take a directory of elements. Whole tandems are not in it, nor anything from `truncated_elements/`.
- **Elements whose attL search window ran into a contig end** get no sequence files by default. The window is `attl_window_bp` on the tRNA's 3′ side; the element itself is bounded by its attR and attL, but a longer one, or a further member of its tandem, could lie past that end, on another contig. These elements stay in the GFF3, with `window_truncated=true`, and in the report. With `TRUNCATED_ELEMENTS=1` their sequences are written to `truncated_elements/`, in the same layout as above, a tandem whole when its tRNA's window was cut. On a complete circular chromosome the first and last base count as ends too: the search does not wrap around them.
- **The GenBank files carry features.** Each marks its element (`mobile_element`; in a whole tandem, one per member), attR and attL (`misc_feature`) and the integrase (`CDS`).
- **The FASTA headers** give the contig span and the tRNA. A whole tandem's header also gives each member's span on it; a member's header gives its tandem and place (k/n).
- **Orientation.** Every sequence is read outwards from the tRNA, so one at a minus-strand tRNA is reverse-complemented, with its features turned with it.
- **Absent directories.** A directory is missing when no genome has anything for it, and none of them is written with `tandem_search: false`.

Coordinates are 1-based and inclusive. They refer to the contigs of the FASTA you passed in. The first column of the GFF3 is that FASTA record id.

Each element has four features. All four are on the strand of the integrase: opposite the tRNA when the integrase sits next to it, on the tRNA's strand when it sits at attL (`integrase_at=attL`):

| GFF3 type | GenBank key, in the element files | meaning |
|---|---|---|
| `mobile_genetic_element` | `mobile_element` | span from attR to attL |
| `attachment_site` with `Note=attR` | `misc_feature` / `attR` | the part of the tRNA, the integration site, that aligns with attL; `product` names the tRNA |
| `attachment_site` with `Note=attL` | `misc_feature` / `attL` | the direct repeat of the tRNA 3′ end |
| `CDS` | `CDS` | the tyrosine integrase |

`confidence=confident` passed the filter. A run with no confident element still writes a valid GFF3 header, plus a report that says so.

Attributes that appear on an element only where they apply:

- `integrase_at=attL` when its integrase sits at the attL end.
- `attR_after=IE<n>` when it integrated into the copy of the tRNA 3′ end that element IE<n> left. Elements integrated one after another into one tRNA are delimited by consecutive copies, so its attR is that copy, not the tRNA.
- `tandem=TA<n>` and `tandem_member=k/m` on each member of a tandem, k counted outwards from the tRNA.
- `found_by=tandem` on a member that only the tandem search found.
- `window_truncated=true` when its attL search window ran into the end of the contig; its sequences are written only with `TRUNCATED_ELEMENTS=1`, to `truncated_elements/`.

An element found from both ends (an integrase at each) is written once.

Intermediate tables and BLAST output are written to a temporary directory and removed when the run succeeds.

## Options

```bash
ANNOTATE_ALL=1 ./run.sh TTHB27c.fasta outdir
KEEP_WORK=1 ./run.sh TTHB27c.fasta outdir
TRUNCATED_ELEMENTS=1 ./run.sh TTHB27c.fasta outdir
SEARCH_PARAMS=strict.yaml ./run.sh TTHB27c.fasta outdir
./run.sh TTHB27c.fasta outdir --cores 8
```

No activation step. `./run.sh` puts `IE_finder_site` on `PATH` for that run. After a change to `envs/IE_finder.yaml`, the next run updates the environment to match. The environment keeps a checksum of the file it was built from.

`ANNOTATE_ALL=1` also writes candidates that have an attL coordinate but failed a later filter. In the GFF3 and the report their `confidence` is `candidate`.

`KEEP_WORK=1` leaves the temporary directory. The path is printed only when a run fails or when this flag is set. Use it to see the raw tables behind the report.

`TRUNCATED_ELEMENTS=1` also writes the sequences of elements whose attL search window ran into a contig end, into `truncated_elements/` (see above). Without it they are only in the GFF3 and the report.

`SEARCH_PARAMS=strict.yaml` reads search thresholds from that file on top of `search_params.yaml`. Keys set in it win, and keys left out keep their values from `search_params.yaml`. An unknown key or a bad value stops the run before any step starts. For example, a file with these four lines brings back the published rules:

```yaml
integrase_hmm_cutoff: none
integrase_len_from: span
attl_select_by: length
attl_reject_gapped: true
```

`--cores N` sets how many genomes are processed at the same time; the default is one per CPU. Each genome runs its steps in order. Arguments after the input are passed to `scripts/ie_finder.py`, which also takes the options above as `--all-candidates`, `--keep-work`, `--truncated-elements` and `--search-params FILE` (see `scripts/ie_finder.py --help`).

If a step fails for one genome, the other genomes still finish and get their files. The failed genome gets a report that shows how far it got, the failing step's log is printed, the temporary directory is kept, and the run exits with status 1.

## The report

`strain.ie.report.txt` is the log for that genome. It has five blocks:

- `search_params` — the thresholds this run used, after any `SEARCH_PARAMS` override
- `counts` — integrase HMM hits, integrase–tRNA pairs, attL BLAST hits, how many passed, and how many confident and rejected candidates had their attL search window cut by a contig end. A rejected candidate with a cut window may be an element that continues on another contig, not a pair without one
- `candidates` — one line each: `confident` or `rejected` plus the first failing reason (`exact_run_too_short`, a CDS overlap, an integrase shorter than 300 aa, a gapped alignment); a candidate whose window was cut says so, with how many bases were searched
- `published_coordinates` and `audit` — the full tables, so a rejected locus can be reconstructed without the temporary files; with the tandem search, `published_coordinates` also gives each element's name and tandem
- `tandems` — with the tandem search, one line per tandem: contig, span, tRNA strand and type, and its members outwards from the tRNA
- `logs` — stdout of the steps that printed something. Empty logs are omitted

If the run crashes, the report is still written from whichever steps finished.

# Evidence and dataset decisions

Audit date: 2026-10-29. Dates record provenance only. Entries summarise public records; no summary-statistics file or expression matrix has been downloaded or opened yet.

## Why this audit

- Human genetic support is associated with higher clinical success for drug targets (Minikel et al., Nature 2024), which motivates checking genetic target claims carefully.
- Several 2024–2025 MR studies propose endometriosis targets, often from FinnGen alone, for example:
  - EPHB4 and others from a druggable-genome MR using FinnGen R10 (Hereditas 2025): https://pmc.ncbi.nlm.nih.gov/articles/PMC12126875/
  - RSPO3 and others from plasma/CSF protein MR (Front. Endocrinol. 2025): https://pmc.ncbi.nlm.nih.gov/articles/PMC11794050/
  - RSPO3 and FLT1 from proteome and metabolite MR (Front. Genet. 2025): https://www.frontiersin.org/articles/10.3389/fgene.2025.1631446/full
  - FSHB, RSPO3, SEZ6L2, EPHB4 from a plasma-proteome screen (2024): https://pmc.ncbi.nlm.nih.gov/articles/PMC11257892/
- Clinical context: hormonal suppression dominates current treatment and the pipeline. Non-hormonal candidates are early, for example vipoglanstat (mPGES-1 inhibitor) entered Phase 2 in 2026 (company and press reports, no efficacy data found).

## Genetic outcome data

- **FinnGen R9, R10, R12 N14_ENDOMETRIOSIS**: public downloads verified reachable 2026-10-29 with MD5 from storage headers. GRCh38. Exact case/control counts per release still to verify; citing papers disagree.
- **Rahmioglu et al. 2023 (GCST90205183)**: meta-analysis excluding 23andMe, GRCh37, public. The paper reports 24 datasets with 58,961 cases and 700,345 controls including 23andMe; fixed-effect inverse-variance; study-level genomic control. Supplementary Table 1 (read 2026-10-31 from a copy Artem downloaded; PubMed Central and nature.com block scripted access) lists the 24 datasets, and their case counts add up exactly to the paper's 58,961 and to the public file's 21,779 European cases. **FinnGen is not one of them**: a Supplementary Table 20 footnote says FinnGen (6,502 cases) contributed only sub-phenotype and adenomyosis results. The public file does contain UK Biobank (6,611 cases), Estonian Biobank (EGCUT, 1,716), deCODE (1,857), BioBank Japan (1,423) and the Northern Finland Birth Cohort (167 cases, 5,235 controls). NFBC1966 and NFBC1986 are FinnGen legacy cohorts ([FinnGen documentation](https://docs.finngen.fi/finngen-data-specifics/red-library-data-individual-level-data/genotype-data/affymetrix-chip-and-its-design/legacy-cohorts-and-chips)), so the file shares about 0.8% of its effective sample with FinnGen. Cohort list: `data/reference/rahmioglu2023_cohorts.json`. https://doi.org/10.1038/s41588-023-01323-z
- **Pan-UKBB phecode 615 (females)**: 4,354 European cases and 217,500 controls; GRCh37; manifest lists no population passing Pan-UKBB QC (`num_pops_pass_qc = 0`), so positive controls must be checked. Independent of FinnGen. Power will be limited for small effects.
- **Multi-ancestry GWAS (Guare et al., Penn, Global Biobank Meta-analysis Initiative)**: 14 biobanks, about 44,000 cases (broad definition) in the preprint, 58 loci in the latest version. No public summary-statistics download or GWAS Catalog accession found (search 2026-10-29). It pools biobanks that include FinnGen and UK Biobank, so it would not be independent of claims from either. https://www.medrxiv.org/content/10.1101/2024.11.26.24316723.full.pdf

- **Pujol Gualdo et al. 2025 (GCST90454213)**: Estonian Biobank plus FinnGen, 233,257 women, GRCh37, public. Found because one included paper shipped it as a supplement. Not independent of FinnGen claims; a possible replication source for the 31 claims discovered in UK Biobank alone.

## Exposure data used to re-derive instruments (downloaded 2026-10-31)

Listed in `data/datasets.json` (role `exposure_for_rederivation`); provenance with SHA-256 in `data/upstream/exposure/*.provenance.json`.

- **Zhao et al. 2023 Olink inflammation proteins**: 13 per-protein GWAS Catalog files, GRCh37, source MD5 matched. The F statistics Gao 2025 prints for 11 instruments are reproduced exactly from these files, which confirms the data and allele coding. The 98 variants Guo 2025 lists for 5 proteins are all present; their exposure p-values all lie below 1e-5, so Guo apparently used P < 1e-5 for proteins, which the paper never states.
- **Sun et al. 2018 INTERVAL cathepsin H**: two aptamers (CTSH.8465.52.3, CTSH.3737.6.3); Aru 2025 does not say which. Neither reproduces Aru's 10 instruments (11 and 17).
- **Ahola-Olli et al. 2017 TRAIL** (GCST004424): no published checksum; whole gzip stream read without error.
- **GLGC 2021 log triglycerides, European**: no published checksum; whole gzip stream read without error. Li 2026's instrument counts cannot be reproduced under any reading of its window (as printed ±100 Mb, ±50 Mb, ±100 kb). Its counts are also internally inconsistent: SCN3A (36) and ITGAV (3) lie 21 Mb apart, so a 100 Mb window would give them nearly the same variants. The methods text cites Willer 2013 (188,577 people), while Table 1 names GLGC 2021. The paper may have used the older, smaller GWAS; this was not tested.
- **GTEx v8 significant cis-eQTL pairs** (whole blood, uterus): source MD5 matched. 5 of the 6 SMR top SNPs Zeng 2024 lists are GTEx's own top variant for the gene; SURF6 (uterus) is not (GTEx top rs1141748, listed rs886090, exposure p 7e-10).
- **1000 Genomes phase 3 European LD reference** (MRC IEU `1kg.v3`): plain http only, no published checksum (approved by Artem). Size and gzip CRC matched, the download resumed once after the server closed the connection, and its structure is correct (503 people, 8,550,156 variants). All 199 of 200 seeded random variants that Ensembl knows agree with Ensembl's 1000 Genomes European frequencies within 0.01 (`results/ld_reference_check.json`).
- **eQTLGen**: still blocked. Rechecked 2026-10-31: www.eqtlgen.org has a valid certificate again, but every download link still points to download.gcc.rug.nl, whose certificate has expired (as has molgenis26.gcc.rug.nl). This is not bypassed. 61 claims depend on it, plus the eQTLGen claims added in batch G.

## Single-cell data

- **GSE266265** (public 2026-03-20): eutopic endometrium, 35 endometriosis and 25 control donors without exogenous hormones, 466,371 cells reported. One 39 GB Seurat RDS plus a 191 MB epithelial object. Raw counts, donor IDs and cycle phase inside the object are unverified. Paper: https://doi.org/10.1038/s41467-026-73020-4. This machine has 63 GB RAM, likely too little to load the 39 GB compressed object whole; look for an alternative format (for example a cellxgene copy) or a chunked reader before Gate D. R 4.3.2 is installed at `C:\Program Files\R\R-4.3.2` but not on PATH; Seurat packages were not installed as of the celiac audit.
- **GSE179640** (Tan et al. 2022): eutopic and lesion tissue; 759 MB raw archive; includes organoid samples.
- **GSE213216** (Fonseca et al. 2023): multiple lesion types; 16 GB raw archive.
- **Human Endometrial Cell Atlas** (Mareckova et al., Nat Genet 2024): annotation reference; integrates GWAS and implicates decidualised stromal cells and macrophages.

## Novelty checkpoint

One web search (2026-10-29) for systematic replication audits of endometriosis MR drug-target claims found none. This is not enough to claim novelty; repeat the check at Gate E, including recent preprints.

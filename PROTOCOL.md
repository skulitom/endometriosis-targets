# Protocol, version 0.1 (DRAFT, not frozen)

Status: written 2026-10-29, before any replication lookup and before any claim extraction was reviewed. This is a working protocol, not a registered preregistration. It becomes binding only when frozen: the SHA-256 of this file and of the claim table are written to `results/freeze.json` before the replication script touches any genomic region near a claimed gene. Every later change goes in the change log with its timing relative to outcome inspection.

## Question

Among genes or proteins that published studies propose as genetically supported drug targets for endometriosis, which claims (1) replicate in an outcome GWAS that shares no participants with the discovery GWAS, (2) show evidence that the molecular signal and the endometriosis signal share a causal variant, and (3) concern genes expressed in endometriosis-relevant cell types in independent patient single-cell data?

Estimand per claim: the cis-MR estimate of the effect of a genetically predicted increase in the gene's expression or protein level (or in the drug-target proxy) on endometriosis risk, using the claim's own instruments, in an independent outcome GWAS. This is a lifetime genetic effect on disease risk. It is not a drug effect, not an effect on pain or progression, and not an effect in women who already have the disease.

## 1. Claim inventory (Gate A)

### Eligibility

Include peer-reviewed articles and preprints that:
- use genetic instruments for a gene- or protein-level exposure: expression (eQTL, SMR), splicing (sQTL), circulating or tissue protein level (pQTL, including inflammatory-protein panels), or a drug-target proxy (variants in or near a drug target gene), and
- use endometriosis (any case definition) as the outcome, and
- report an effect estimate for that gene or protein.

Instrument scope is recorded per claim (`cis`, `cis_and_trans`, `genome_wide`, `unclear`). Claims whose instruments are not all near the gene are kept, because the audit asks whether published claims hold up, but they are summarised separately: trans instruments are more prone to pleiotropy, so they are weaker target evidence.

Exclude: MR of exposures that are not gene products (BMI, steroid hormone levels, lipid traits with genome-wide instruments, immune-cell counts or phenotypes, microbiota, metabolites, behaviours, other diseases) unless the paper also reports a gene- or protein-specific instrument; TWAS/FUSION without MR; colocalisation-only studies (logged separately as context); endometriosis as exposure; reviews; abstracts without estimates; retracted papers.

Change made 2026-10-29 before screening began: broadened from cis-only instruments to gene- or protein-level exposures with recorded instrument scope.

**Clinical biomarkers and peptide hormones** (rule added 2026-10-29 after title/abstract screening, before full-text review; screeners flagged TSH, IGF-1, leptin, insulin, PTH, SHBG, CRP, apolipoproteins as ambiguous). Such a protein is eligible only if at least one holds: (a) it was measured on a proteomic platform (Olink, SomaScan or similar) as part of a protein panel; (b) it was instrumented with variants in or near its own gene; (c) the paper explicitly proposes it, or the drug acting on it, as a therapeutic target. Otherwise, for example SHBG or IGF-1 from UK Biobank clinical-chemistry GWAS with genome-wide instruments, it is treated as a physiological trait and excluded (code E2). Rationale: genome-wide instruments for such traits mostly capture whole hormonal or metabolic axes rather than evidence about one gene as a target.

**Adenomyosis and conference abstracts** (clarifications made 2026-10-31 while screening the supplementary search, before any replication lookup; provisional until Artem confirms them):
- A study whose only endometriosis-related outcome is adenomyosis (FinnGen `N14_ENDOMETRIOSIS_UTERUS`, ICD-10 N80.0) is excluded with code E3. Adenomyosis is diagnosed and analysed as a separate condition, for example in Rahmioglu 2023. No claim in the table so far rests on an adenomyosis-only outcome.
- Conference abstract books are excluded with code E4, because eligibility covers peer-reviewed articles and preprints.
- FinnGen "female infertility associated with endometriosis" is an infertility endpoint (infertile women with endometriosis against women without infertility), not an endometriosis case definition, so it gives no claim (E3 if it is a paper's only endometriosis-related outcome). This was added before the batch G extraction.

Preprints whose published version is also in the records are merged into one paper at full-text review; the published version is the source of record, and disagreements between versions are noted.

### Unit

One claim = paper × gene × exposure dataset (tissue/platform) × outcome GWAS. A claim is "headline" if the abstract or conclusion names the gene as a candidate target or causal gene. Claims from the same paper are linked by paper ID; the same gene claimed by several papers stays as separate claims and is also summarised at gene level.

### Search

PubMed and Europe PMC (Europe PMC adds preprints). Exact query strings, run date and hit counts are stored by `src/literature.py` in `data/literature/`. Search cut-off is the run date. Reference lists of included papers are scanned for missed studies; additions are logged with their source.

**Reference-list scan** (run 2026-10-31, before any replication lookup; `src/refscan.py`, `data/literature/refscan/`). Reference lists come from Europe PMC, Crossref and the PMC full text (JATS), combined. A cited study is screened if it is not in the search records and either (a) its title or abstract mentions endometriosis together with a genetic-causal term (Mendelian, SMR, summary-data-based, QTL, druggable, drug target, genetically predicted/proxied/determined, instrumental variable, causal, colocalisation), or (b) its title has such a term and Europe PMC's full-text index finds endometriosis in it, or (c) its title or abstract mentions endometriosis and its open full text mentions MR, SMR or a drug target. Rules (a) and (b) were fixed before any reference list was fetched. Rule (c) was added after the first match run and before screening, because Rahmioglu 2023 reports gene-level SMR only in its full text and supplement. Candidates are screened with the same codes as the database records; full-text passages used as evidence are kept in `fulltext_snippets.json`. One included paper (an Authorea preprint) has no reference list available by any open route; its publisher blocks scripted access.

**Supplementary full-text search** (approved by Artem 2026-10-31, after the reference-list scan, before any replication lookup). The scan showed that the title/abstract search misses papers reporting SMR only in their full text. Europe PMC is searched once more, matching full texts as well as titles and abstracts: `"endometriosis" AND ("summary-data-based" OR "summary data-based" OR HEIDI OR "mendelian randomization" OR "mendelian randomisation")`. The query, run date and hits are stored by `src/literature.py supplementary` as `search_europepmc_fulltext_<date>.json`, and the search cut-off is that run date. An exploratory count with a similar query (`refscan/fulltext_search_probe.json`) was run first to size the gap; its titles were seen, but nothing was screened from it. Only open full texts are indexed, so papers whose SMR appears only in a paywalled full text can still be missed; this is reported as a limitation. Screening:
- Records not already in the search records are screened with the same codes. The evidence is the title, the abstract and the full-text passages that name MR, SMR or HEIDI (`refscan/supplementary_passages.json`).
- Screening is AI-assisted in batches. The lead re-checks every include and unsure. The lead also re-checks every exclusion whose passages describe the paper's own MR or SMR analysis ("we performed", "we conducted" or similar) and mention endometriosis.
- A random 20% of the remaining exclusions (minimum 10) is offered for human review, separate from the database-search sample.

Result (run 2026-10-31):
- The search returned 968 records. 717 were in neither the search records nor the reference-list candidates (`records_supplementary.json`).
- Eight screeners each handled one batch of about 90 records (`screening/decisions_supp_1..8.json`). The lead re-checked 56 records and changed 22 decisions (`screening/lead_checks_supp.json`).
- The re-check covered every include and unsure, the 13 exclusions flagged by the own-analysis rule, and 17 exclusions decided on incomplete excerpts.
- The incomplete excerpts came from 303 records fetched by an earlier run without the endometriosis sentences. They were refetched before the re-check, and the old excerpts are kept as `refscan/supplementary_passages_oldformat.json`.
- Europe PMC serves no open XML for 118 records, so their excerpts come from NCBI PMC. For 97 records no full text could be read (publisher-restricted, or not in PMC), and these were screened on title and abstract.
- 13 records passed to full-text review, alongside the 6 from reference lists. 132 unchecked exclusions form the human-review sample (`exclusion_review_sample_supplementary.json`).

### Screening and extraction

- Title/abstract screening, then full-text screening; each exclusion gets a recorded reason.
- Extraction fields are defined in `src/claims.py`. Every number carries its source location (table, figure, supplement, page). Where only the abstract is accessible, the claim is marked `full_text_accessed=false`.
- Screening and extraction are AI-assisted. Every headline claim gets an independent second extraction; disagreements are resolved against the source and logged. Human (Artem) verification status is recorded per claim. A random sample of exclusions (at least 20%, minimum 10) is offered for human review.

### Independence map

For each claim, record the cohorts in the discovery outcome GWAS. A replication source is eligible for that claim only if it shares no participants with the discovery outcome, with one accepted exception: the Northern Finland Birth Cohort, which is both a FinnGen legacy cohort and one of the Rahmioglu 2023 datasets (167 cases, 5,235 controls; 0.8% of Rahmioglu's effective sample; estimated correlation between estimates about 0.009). Results that rest on this overlap are labelled "near-independent (NFBC overlap)" (`src/claims.py`, `independence_note`). Accepted by Artem on 2026-10-31, before any replication lookup. Later FinnGen releases contain earlier releases' participants, so FinnGen R12 is never an independent replication of a FinnGen R5–R11 discovery; it may be reported as a labelled, non-independent update.

## 2. Data and sanity checks (Gate B)

### Candidate replication outcome sources

| Source | Cases / controls (European) | Independent of FinnGen? | Notes |
|---|---|---|---|
| Pan-UKBB phecode 615, females | 4,354 / 217,500 | Yes | Pan-UKBB manifest lists `num_pops_pass_qc = 0`; must pass positive controls before use |
| Rahmioglu 2023 meta-analysis excluding 23andMe (GCST90205183) | 21,779 / 449,087 European plus 1,713 / 1,581 Japanese (GWAS Catalog) | Near-independent: shares only NFBC (accepted) | Supplementary Table 1 read 2026-10-31 (file supplied by Artem; `data/reference/rahmioglu2023_cohorts.json`): 24 datasets, whose case counts add up exactly to the paper's and the GWAS Catalog's totals. FinnGen is not one of them (Supplementary Table 20 footnote: FinnGen contributed only sub-phenotype and adenomyosis results). The public file contains UK Biobank (6,611 cases), Estonian Biobank (1,716), deCODE (1,857), BioBank Japan (1,423) and the Northern Finland Birth Cohort (167 cases, 5,235 controls). NFBC1966/1986 are FinnGen legacy cohorts (FinnGen documentation), so the file shares those participants with FinnGen: 0.8% of its effective sample size. Artem accepted this overlap on 2026-10-31 (section 1, independence map), so it is eligible for FinnGen claims, labelled near-independent |
| FinnGen R12 | To verify | n/a | Independent only of claims whose discovery outcome had no FinnGen participants |

If Rahmioglu 2023 includes FinnGen, removing FinnGen's contribution from the fixed-effect meta-analysis is considered only if the exact release used is public, per-study genomic-control handling can be matched, and positive controls pass. This decision is made and recorded before the freeze.

### Exposure-side overlap

If a claim's exposure QTL data come from UK Biobank participants (for example UKB-PPP) and the replication outcome is UK Biobank, report the overlap and prefer a non-UK Biobank exposure source for the same gene when one exists. Record instrument strength (F statistic).

### Checks before freezing

1. **Checksums** for every input (`src/fetch.py`).
2. **Harmonisation**: genome build per file (FinnGen R9–R12: GRCh38; Pan-UKBB and GCST90205183: GRCh37); match on rsID with position cross-check; align effect alleles; drop palindromic SNPs with minor allele frequency above 0.42 unless strand can be resolved.
3. **Positive controls**: lead variants of established endometriosis GWAS loci (Rahmioglu 2023) that lie more than 1 Mb from every claimed gene. Each replication source should show mostly concordant directions; report the concordance rate and its binomial p-value. A source that fails is not used.
4. **Reproduction**: re-compute one published MR estimate from its stated discovery outcome data and instruments; report agreement or the discrepancy.

### Positive-control and reproduction rules (fixed 2026-10-30, before any GWAS file was read)

- **Lead variants**: from the Rahmioglu 2023 file read through `src/blinding.py` (so regions within 1 Mb of any claimed gene are already removed), take variants with p < 5e-8, then keep the lowest-p variant and drop all others within 1 Mb, repeating until none remain.
- **Matching**: by chromosome, position and alleles for GRCh37 sources (Pan-UKBB, Estonia + FinnGen); by rsID for FinnGen (GRCh38), with rsIDs taken from Ensembl GRCh37 at the lead positions. Effects are aligned to the Rahmioglu effect allele; palindromic variants with allele frequency 0.42 to 0.58 are dropped.
- **Columns used**: Pan-UKBB European columns (`beta_EUR`, `se_EUR`, effect = `alt`), excluding rows with `low_confidence_EUR = true`; FinnGen `beta`/`sebeta` (effect = `alt`); Estonia + FinnGen and Rahmioglu `beta`/`standard_error` (effect = `effect_allele`). rsIDs for FinnGen matching come first from the Estonia + FinnGen file's `rs_id` at the same GRCh37 position, otherwise from Ensembl GRCh37.
- **Pass criterion per source**: at least 70% of matched lead variants have the same direction as in Rahmioglu 2023, and the one-sided binomial p-value for concordance is below 0.01. A source that fails is not used for replication. Concordance is expected to be inflated where a source overlaps Rahmioglu 2023 (for example UK Biobank), so this checks file handling and allele coding, not independent replication.
- **Reproduction target**: among claims with extracted instruments, a headline claim whose discovery outcome is FinnGen R10 (the file held), restricted to genes for which every claim in the table includes FinnGen in its outcome (so FinnGen data near that gene can never be an independent replication source). Choose the one with the most instruments; ties by claim_id. Its gene is unmasked in the FinnGen R10 file only, through a logged exemption. Agreement: recomputed log-odds within 0.02 or 5% of the reported value.
  - Clarification (2026-10-30, after applying the rule to the claim table, before reading FinnGen R10): "stated discovery outcome data" means the paper names the endpoint held (N14_ENDOMETRIOSIS). The rule's first pick, Xie2026FI-1 (LRP1, 3 genome-wide instruments), does not name its endpoint, so the primary reproduction is the next pick, Tao2024-1 (EPHB4). Xie2026FI-1 is also run and reported as secondary, with that caveat.
  - Clarification (2026-10-30, after a first run returned no rows because the instruments were masked by neighbouring claimed genes): other claimed genes within 1 Mb of a target's instrument are also exempted in FinnGen R10, but only if every claim of each such gene includes FinnGen in its outcome. Otherwise the instrument stays masked and the target is reported as not reproducible. For Tao2024-1 this adds FIS1 and GIGYF1. No outcome value was seen in the first run (the log shows 2 rows returned, both for LRP1, and none were printed).

### Blinding rule

Before the freeze, no replication-source lookup within 1 Mb of any claimed gene, except positive-control loci selected by the rule above. All endometriosis GWAS reads go through `src/blinding.py`, which drops masked rows before returning anything and appends each read to `results/blinding_log.jsonl`. Exposure files (gene expression, protein and lipid GWAS, the LD reference) contain no endometriosis data and are read directly (section 3, re-derivation). Gene positions come from Ensembl for both builds (`data/reference/claimed_genes.json`); the reader refuses to run if any claimed gene lacks a position.

## 3. Primary replication analysis (Gate C)

- **Instruments**: the variants and exposure effect sizes the original paper reports for the claim. If not reported, re-derive them as below when the exposure data are openly available; otherwise the claim is "not testable".

### Re-derivation of instruments not published per variant (plan added 2026-10-31, before any replication lookup)

For the 177 primary claims without published per-variant effect sizes:

- **Eligible**: the paper names its exposure dataset precisely enough to find the original public file (94 claims across eQTLGen, GTEx, Olink 91 proteins from Zhao 2023, INTERVAL proteins, Ahola-Olli 2017 cytokines and GLGC lipids). Each paper's stated selection rule is recorded in `data/literature/instruments/rederive/rules_*.json` with its source location. Not eligible: claims whose exposure data need an account or accepted terms (UKB-PPP, deCODE: 39 claims) unless Artem decides otherwise, and claims with no stated exposure dataset (about 44); these stay "not testable".
- **Route 1, known variants** (preferred): when the paper lists the variant IDs but not their effects (tables of rsIDs and F statistics, top SNPs of SMR, figure labels), take those variants' exposure effects from the stated dataset. Check: every listed variant is found; reported F statistics, where given, match (beta/se)² within 10%.
- **Route 2, stated rule**: otherwise apply the paper's own rule to the stated dataset in this order: restrict to the stated window if any; p-value threshold; minor allele frequency filter if stated; LD clumping with the stated r² and distance against the stated reference (1000 Genomes European if unstated); F statistic filter if stated. A stated fallback (for example a relaxed threshold when too few variants pass) is applied as written. Check: the number of variants equals the paper's reported number.
- **Classification**: re-derived claims are analysed and reported separately from claims with published instruments. A route-2 claim whose re-derived count differs from the reported count, or a route-1 claim with missing variants or an F-statistic mismatch, is analysed as a sensitivity result only.
- Exposure files contain no endometriosis outcome data, so reading them does not break the blinding rule; outcome data are read only after the freeze.
- If a dataset cannot be obtained over a verified connection, its claims stay "not testable" unless Artem decides otherwise.
- **Fixed details** (2026-10-31, written before any exposure file was opened):
  - *Routes*: Route 1 for Guo2025Med-1..5 (rsIDs in Stable8), Gao2025Med-1..3 (rsIDs and F in Table 1) and Zeng2024Gene-1,2,3,9,10,11 (SMR top SNP and probe in Table S1; GTEx v8 whole blood or uterus as the table names). Route 2 for WeiY2024ITT-1..6, Li2026tg-1..7, Aru2025cath-1 and Xu2025JCMM-1. (The plan of the same day said 15 and 14; Xu2025JCMM-1 lists no variants, so it belongs to Route 2: 14 and 15.)
  - *Unstated parameters*: LD reference 1000 Genomes phase 3 European (MRC IEU `1kg.v3`, 503 people), the panel TwoSampleMR uses; clumping distance 10,000 kb (TwoSampleMR default) where only r² is stated; variants absent from the reference are dropped, as TwoSampleMR clumping does. Clumping is greedy by ascending p (ties by position): each index variant removes the variants within the distance whose r² with it exceeds the threshold. r² is the squared Pearson correlation of allele dosages in the 503 people. This approximates PLINK and the OpenGWAS server and can differ for variants near the threshold.
  - *Palindromic and outcome-stage steps*: the count check is made at the exposure stage. Steps that need the outcome file (dropping variants missing from the outcome, palindromic variants whose strand cannot be inferred) happen only after the freeze, and are applied with the replication outcome. Aru2025cath-1's exposure-side rule (drop palindromic variants with 0.42 < effect allele frequency < 0.58) is applied at the exposure stage. Unitemised removals (PhenoScanner, MR-PRESSO) cannot be reproduced; a count that differs because of them still counts as a mismatch.
  - *Ambiguous wording*: where the paper admits more than one reading, each reading is rebuilt. If exactly one reproduces the reported count, it is used; otherwise the reading listed first below is used, and the claim is sensitivity-only. Readings fixed now: Li2026tg window "within a 100 Mb window surrounding the genomic region" = (1) gene region ± 100 Mb as printed, (2) ± 50 Mb (100 Mb total), (3) ± 100 kb (a kb/Mb misprint). Aru2025cath-1 aptamer, not stated: (1) CTSH.8465.52.3 (GCST90240624), (2) CTSH.3737.6.3 (GCST90240625).
  - *Gene regions*: as printed in the paper's own files where given (Li2026tg S4, read as GRCh37); otherwise Ensembl GRCh37 gene start and end (`data/reference/claimed_genes.json`).
  - *p-value column*: the file's main p-value; for GLGC 2021 the uncorrected `pvalue_neg_log10`, not the genomic-control column (the paper does not say which it used). For the reported count (Li2026tg-4: Fig 5 says 235, S1 Table 250), the figure or table the claim's estimate came from is used (Fig 5).
  - *Route 1 matching*: variants matched by rsID. Where the exposure file has positions only, rsID ↔ position comes from the 1000 Genomes reference (GRCh37) and alleles must agree. GTEx v8 (GRCh38) positions for the Zeng top SNPs come from Ensembl. A listed variant absent from the exposure file gets one LD proxy (r² ≥ 0.8) and is marked as a proxy. SMR top SNPs must have exposure p < 5e-8 (the SMR default); otherwise the variant is flagged. Gao2025Med-3 (SLAMF1) is sensitivity-only from the start: Table 1 lists 7 variants but the paper's IVW tables report 3 (published list does not match the analysed set; amendment of 2026-10-30).
  - *1000 Genomes file check* (it comes over plain http with no published checksum): its structure must match (`.bed` size = 3 + variants × ⌈people/4⌉, 503 European people). A seeded random sample of 200 European variants must agree with Ensembl's 1000 Genomes phase 3 European allele frequencies (over https) within 0.01 for at least 98% of the variants found in Ensembl, with matching alleles. If this fails, the file is not used and Route 2 claims stay not testable.
- **Missing variants**: one LD proxy (r² ≥ 0.8, 1000 Genomes European) per missing instrument; record proxy use.
- **Estimator**: Wald ratio for one variant; inverse-variance weighted (fixed effect) for two or three; random effects if four or more and Cochran's Q p < 0.05. Report odds ratio per unit of exposure as in the original paper, with 95% CI.
- **Power**: power at the strict threshold to detect (a) the discovery estimate and (b) half of it, to allow for winner's curse.
- **Choice of replication source** (added 2026-10-31, before any replication lookup): eligible sources are those that passed the positive controls and share no participants with the claim's discovery outcome (`independent`; the NFBC overlap is accepted, section 1). Rahmioglu 2023 supplied the positive-control variants, so it counts as passing: the agreement of the four other sources with it checks its allele coding too. Case counts used: Rahmioglu 2023 public file 23,492 (21,779 European plus 1,713 Japanese); FinnGen R12 20,190 (R12 manifest); Estonia + FinnGen 19,588; Pan-UKBB 4,354. Applied in `src/tier_dryrun.py`. If several are eligible and they share participants with each other (for example FinnGen R12 and Estonia + FinnGen), the one with the most endometriosis cases, by the source's own documentation, gives the primary estimate and the others are reported as supportive. If eligible sources share no participants with each other, their estimates are combined (fixed-effect inverse-variance) as the primary estimate. A FinnGen release later than the discovery release is never eligible (it is reported as a labelled, non-independent update only).
- **Classification** (strict threshold α = 0.05 / number of distinct testable genes):
  - *Replicated*: same direction as discovery and p < α.
  - *Directionally consistent*: same direction and p < 0.05, but not p < α.
  - *Not replicated*: 95% CI excludes half the discovery log-odds ratio (the true effect is likely much smaller than claimed, or opposite).
  - *Inconclusive*: none of the above.
  - *Not testable*: no usable instruments, no independent source, or failed harmonisation.
- Report heterogeneity between discovery and replication estimates (Cochran's Q).

### Colocalisation (secondary)

coloc.abf with priors p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, window ±500 kb around the lead cis variant. Run when full regional exposure summary statistics are openly available (eQTL Catalogue, GTEx, eQTLGen). Gated protein datasets are used only if Artem decides to accept their terms. PP.H4 ≥ 0.8 supports a shared variant; PP.H3 ≥ 0.8 suggests distinct variants (possible confounding by linkage); otherwise inconclusive. Run in each eligible replication source and, labelled as non-independent, in the largest available outcome GWAS. SuSiE-based colocalisation is a sensitivity analysis when an LD reference is available.

### Sensitivity

Where the primary estimate is near-independent (NFBC overlap), the estimate from the largest strictly independent eligible source is also reported (for FinnGen-discovered claims, Pan-UKBB). Steiger directionality test; flag claims inside the MHC region; alternative endometriosis definitions where available (for example FinnGen sub-endpoints); leave-one-variant-out for multi-variant instruments.

## 4. Cell-type layer (Gate D; detailed plan frozen separately before looking at disease contrasts)

- Candidate cohorts: GSE266265 (eutopic endometrium, 35 endometriosis and 25 control donors without hormone treatment), GSE179640 (Tan et al. 2022, eutopic and lesion tissue), GSE213216 (Fonseca et al. 2023, lesion types). Human Endometrial Cell Atlas as an annotation reference.
- Aggregate raw counts per donor and cell type (pseudobulk). Never treat cells as independent patients.
- Primary output: where each claimed gene is expressed (detection fraction and normalised level by cell type), consistent across cohorts. Secondary: endometriosis-versus-control differences in eutopic tissue with cycle phase as covariate where recorded; lesion-versus-eutopic within donor where paired.
- Expression in a relevant cell type supports plausibility only. No differential expression does not refute a genetic target.

## 5. Interpretation (Gate E)

- Drug cross-reference from Open Targets Platform and ChEMBL: existing drugs, mechanism, development stage. Check direction: if higher genetically predicted level increases risk, an inhibitor is the compatible drug class, and vice versa.
- Flag fertility, pregnancy and safety considerations relevant to the patient population.
- Tiers combine replication class and colocalisation; cell-type evidence and drugs are annotations (rules below). The tier rules are frozen with this protocol.

### Tier rules (drafted and approved by Artem 2026-10-31, before any replication lookup; frozen with this protocol)

Tiers are assigned per gene from that gene's primary claims (mQTL claims are reported separately and never set a tier). Claim-level classes are those of section 3, from each claim's primary replication source.

A claim is a **qualifying replication** when all hold:
1. it is classified *Replicated*;
2. its instruments are published and verified, or re-derived with the check passed (`analysis_role = primary_rederived`); claims marked sensitivity-only do not qualify;
3. its instrument scope is `cis`;
4. no instrument lies in the MHC region (GRCh37 chr6:28,477,797–33,448,354), or, for a claim without positions, the gene does not;
5. the Steiger test does not point the other way (outcome variance explained exceeding exposure variance, p < 0.05);
6. colocalisation in the replication source, where run, does not give PP.H3 ≥ 0.8;
7. the gene does not lie within 1 Mb of a Gate B positive-control lead variant (`results/gateb_positive_controls.json`, 14 leads, GRCh37 gene coordinates).

**Claims at Gate B positive-control loci** (decided by Artem 2026-10-31, before any replication lookup). Gate B read every replication source at 14 lead variants chosen more than 1 Mb from the claimed genes known then. Genes added later from reference lists and the supplementary search can lie near those variants: 29 claims on 13 genes at 6 leads, listed in `results/blinding_impact_additions.json`. Their replication is not blind, and at a known GWAS locus it is expected and cannot tell neighbouring genes apart. These claims are therefore:
- still tested in the frozen replication;
- labelled "not blind: Gate B positive-control locus";
- reported in a separate table;
- capped at Tier 2 (condition 7).

The rule is applied in code (`src/tier_dryrun.py`, `gate_b_lead_near`) to every claim, not to a fixed list. Claims found earlier cannot meet it, because the leads were chosen away from their genes.

**Conflicting evidence** for a gene: another primary claim for the same gene and the same exposure type (eQTL, pQTL, sQTL or drug-target proxy) that is *Not replicated*, or *Replicated* in the opposite direction. Claims of different exposure types are shown side by side and do not count as conflict, because expression in one tissue and circulating protein can legitimately differ.

**Direction scale** (added on approval, 2026-10-31): directions are compared on one scale. For expression, splicing and protein claims, that is per higher genetically predicted level of the gene product. For drug-target proxies, it is per higher level of the proxied biomarker (for example LDL-C, triglycerides) through variants at the target gene. The field `effect_scale` records each claim's convention (`src/claims.py`, `common_direction`):
- *per_lower*: the paper reports the effect per lower level or per inhibition. Flipped before conflict checks and drug matching: Aru2025pcsk9-1 (LDL-lowering), Xue2025ugcg-1 (inverted expression OR).
- *unclear*: the scale is not stated. Left out of conflict checks and drug-direction matching, and flagged: Zhou2024lipid-1..3, Xiang2024bmd-1.
- absent: per higher level (all other claims).

Replication estimates are computed on each paper's own scale, so the discovery and replication estimates of a claim are always compared like for like; the shared scale is used only when claims are compared with each other or with a drug.

Gene tiers, applied in this order:
- **Tier 1, holds up**: at least one qualifying replication and no conflicting evidence.
  - *1A*: colocalisation PP.H4 ≥ 0.8 in an independent source as well.
  - *1B*: colocalisation not available or inconclusive.
- **Tier 2, partial support**: not Tier 1, and at least one claim *Replicated* or *Directionally consistent*. This includes replications that fail a qualifying condition (trans or genome-wide instruments, sensitivity-only instruments, MHC, Gate B positive-control locus, Steiger, PP.H3, or conflict). The reason is listed for every Tier 2 gene.
- **Not replicated**: none of the above, and at least one primary claim *Not replicated* (its 95% CI excludes half the discovery effect).
- **Inconclusive**: at least one testable claim, but none of the above. Power at half the discovery effect is shown, so underpowered tests are not read as negative.
- **Not testable**: no claim had usable instruments and an eligible source.

Annotations, which never change a tier:
- **Cell type**: expressed in an endometriosis-relevant cell type in at least two independent single-cell cohorts (yes / no / not assessed), with the cell types; definitions follow the Gate D plan, frozen before disease contrasts. Not tier-setting because section 4 treats expression as plausibility only, and circulating proteins can be made outside the endometrium.
- **Drugs**: highest development phase in Open Targets / ChEMBL and the mechanism; whether the mechanism matches the genetic direction on the shared scale (inhibitor, or a drug that lowers the proxied biomarker, when higher level raises risk; activator when it lowers risk).
- **Context**: number of papers and claims for the gene, re-derived instruments, fertility and pregnancy considerations, and the MR assumptions that remain untested.

Wording: Tier 1 means "genetic association replicated in independent data", not "causal" or "a treatment". Every tier is a hypothesis for laboratory testing.
- Novelty check: search for prior replication audits of endometriosis MR targets before claiming the audit is new.
- MR assumptions (relevance, independence, exclusion restriction) are stated per tier; a replicated MR estimate is a hypothesis for laboratory testing.

## Stop or pivot conditions

Fewer than five testable claims; every replication source fails positive controls; instruments unavailable for most claims; unresolvable sample overlap. In each case write a closeout report instead of producing an optimistic target list.

## Change log

- 0.1 (2026-10-29): first draft, before claim extraction review and before any replication lookup.
- 0.1 amendments (2026-10-29, all before any replication lookup):
  - Eligibility broadened to gene- or protein-level exposures with recorded instrument scope (before screening).
  - Clinical-biomarker and peptide-hormone rule added (after title/abstract screening, before full text).
  - DNA-methylation QTL claims labelled `mQTL` and made secondary-only (after first extraction): CpG-to-gene mapping is indirect and methylation is not itself a drug target.
  - Second-extraction agreement rules fixed in `src/second.py` before comparing (rounding: absolute 0.006 or 1% for estimates and CIs; factor 1.15 for p-values).
  - `instrument_scope = cis` may rest on observed variant positions in a supplement (all within about 1 Mb of the gene) when the methods do not state a window; the evidence is kept in the adjudication record.
  - Instrument extraction (2026-10-30, before any replication lookup): variants are copied from supplement tables by `src/instruments.py` using per-paper specs that select rows by labels only. Extraction is verified by recomputing the reported estimate where outcome effects are listed, by an independent spreadsheet reader, and by a hand check of the rest. Claims whose published variant list does not match the analysed set (`instruments_check.flags`, `source_issues`) are analysed with their published list and reported as a sensitivity result, not a primary one.
  - Re-derivation plan (2026-10-31, before any replication lookup): section 3 "Re-derivation of instruments", approved by Artem the same day. Fixed details (routes, unstated parameters, ambiguous readings, 1000 Genomes file check) written before any exposure file was opened. The blinding rule now says explicitly that it covers endometriosis files; exposure files never fell under it (section 3).
  - Rahmioglu 2023 overlap (2026-10-31, before any replication lookup): Estonian Biobank, deCODE and "other" were added to the overlap map. Supplementary Table 1, read later the same day from a copy Artem downloaded (the site's bot check blocks scripts), confirms those cohorts and shows FinnGen is not included. The only FinnGen link is NFBC, a FinnGen legacy cohort (167 cases, 5,235 controls; estimated correlation between estimates about 0.009). Decision (Artem, 2026-10-31, before any replication lookup): the NFBC overlap is accepted. Every result that rests on it is labelled "near-independent (NFBC overlap)", and the strictly independent estimate is reported alongside where one exists (section 3, sensitivity). Effect on the dry run: 98 FinnGen-discovered claims move from Pan-UKBB to Rahmioglu 2023 as primary source, and 26 Rahmioglu-discovered claims become testable in FinnGen R12 (`results/tier_ceiling_dryrun.json`).
  - Replication-source choice and tier rules (2026-10-31, before any replication lookup): section 3 "Choice of replication source"; section 5 "Tier rules". Cell-type evidence changed from tier-setting to an annotation, for the reason given there.
  - Tier rules approved by Artem (2026-10-31, before any replication lookup), with one addition: the shared direction scale (section 5). Six claims' conventions recorded in `corrections_lead.json` (2 per lower level, 4 unclear).
  - Reference-list scan (2026-10-31, before any replication lookup): candidate rules written into section 1 (Search). Rule (c), the full-text MR/SMR check, was added after the first match run, before screening. Two matching bugs found in the first run (italic tags in titles, DOIs inside free-text citations) were fixed before screening; the first run's files are kept with the suffix `_run1_buggy`. Result: 29 candidates screened, 6 included for full text and claim extraction, all reporting gene- or methylation-level SMR on endometriosis only in their full text or supplement (`data/literature/refscan/additions.json`).
  - Supplementary full-text search (2026-10-31, approved by Artem, before any replication lookup): rule and result in section 1 (Search). Two screening clarifications were made while screening it, both provisional until Artem confirms them: adenomyosis-only outcomes are E3, and conference abstract books are E4 (section 1, Eligibility). Records from other sources are kept in `records_supplementary.json` with their own decisions and review sample, so the database-search records, decisions and review sample are unchanged (checked byte for byte).
  - Batch G extraction (2026-10-31, before any replication lookup).
    - **Inputs:** supplements for the 19 additions were downloaded through the Europe PMC archive route, as approved by Artem (95 MB, 17 papers). The Rahmioglu 2023 supplement is Artem's manual copy. The 2026 Nature Genetics GWAS supplement is not on Europe PMC, so its claims come from the main text. Open full texts are saved as text by `src/fulltext.py`.
    - **Extraction:** four extractors followed `extraction/instructions_G.md` and produced 133 claims from 15 papers; 4 papers were excluded at full text. A blind second extraction of the 15 headline claims was compared under the fixed rules in `src/second.py`: 10 agreed and 6 disagreed on one field each. The 6 were adjudicated (`adjudication/adjudicated_4.json`). The 174 earlier comparisons and corrections were regenerated unchanged.
    - **Blinding consequence** (`results/blinding_impact_additions.json`): 13 of the 63 new genes lie within 1 Mb of 6 of the 14 Gate B positive-control leads, which carry 29 of the new claims. The replication sources were therefore read at those variants before the genes became claims. Decided by Artem the same day, before any replication lookup: they are tested, labelled "not blind: Gate B positive-control locus", reported separately and capped at Tier 2 (section 5, qualifying condition 7). The blinding reader now masks all 332 resolved genes, and refuses to run until the last 2 unresolved symbols (HCG2P7, ENSG00000289405) have positions.

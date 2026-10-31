# Claim extraction guide

Used for every included paper. One JSON object per claim, validated by `python src/claims.py validate <file>`.

## What counts as one claim

One paper × gene (or protein) × exposure dataset (tissue or platform) × outcome GWAS. Extract a claim for every gene/protein the paper reports as significant for endometriosis in its main or primary MR analysis (after the paper's own multiple-testing correction if it applies one, otherwise its own stated significance threshold). Proteome-wide or transcriptome-wide screens: extract the significant hits only, not every protein tested, and record the number tested in `notes`.

Set `headline = true` when the abstract or conclusion names the gene as a candidate drug target, therapeutic target, or causal gene.

## Fields

| Field | Rule |
|---|---|
| `claim_id` | `<paper short id>-<n>`, for example `Li2025-1` |
| `paper_id` | `doi:<doi>` if a DOI exists, else `pmid:<pmid>` |
| `first_author`, `year`, `journal`, `preprint` | From the article |
| `gene_symbol` | HGNC symbol. Protein names are mapped to their gene (e.g. "IL-6 receptor" → `IL6R`). Complexes: one claim per gene only if the paper reports genes separately, else the main subunit with a note |
| `exposure_type` | `pQTL`, `eQTL`, `sQTL`, `SMR_eQTL`, `SMR_pQTL`, `drug_target_proxy`, or `other` |
| `exposure_source` | Dataset and first author/year, e.g. `deCODE (Ferkingstad 2021)`, `UKB-PPP (Sun 2023)`, `eQTLGen (Vosa 2021)`, `GTEx v8`, `Zhao 2023 Olink inflammatory proteins` |
| `exposure_tissue` | Tissue for eQTL/sQTL (e.g. `Uterus`, `Whole blood`), `plasma` for most pQTL |
| `instrument_scope` | `cis` (all instruments near the gene, usually ±500 kb or ±1 Mb), `cis_and_trans`, `genome_wide` (selected genome-wide without regard to gene location), or `unclear` |
| `outcome_source` | As described, e.g. `FinnGen R9 N14_ENDOMETRIOSIS`, `IEU ukb-b-N80`, `Rahmioglu 2023 meta-analysis` |
| `outcome_release` | FinnGen release or GWAS ID if given |
| `outcome_cases`, `outcome_controls` | Integers as reported, else `null` |
| `outcome_cohorts` | List from: `FinnGen`, `UK Biobank`, `23andMe`, `BioBank Japan`, `Rahmioglu2023_meta`, `deCODE`, `Estonian Biobank`, `Million Veteran Program`, `All of Us`, `other`, `unknown`. Use `Rahmioglu2023_meta` for the Rahmioglu/IEC meta-analysis, not its parts |
| `n_instruments` | Number of variants used for this gene's estimate |
| `instruments` | List of `{rsid, effect_allele, other_allele, beta_exposure, se_exposure}` only when the paper reports them per gene; else `null`. Add `beta_outcome`, `se_outcome` if reported |
| `instruments_location` | Where the variant list is (e.g. `Supplementary Table S2`), even if not extracted |
| `method` | The primary method for the reported estimate, e.g. `Wald ratio`, `IVW`, `SMR` |
| `estimate_scale` | `OR` (most), `beta`, or `HR` |
| `estimate`, `ci_low`, `ci_high`, `p_value` | Primary estimate as reported (per unit or per SD of exposure). If only beta and SE are given with OR scale implied, convert OR = exp(beta), CI = exp(beta ± 1.96·SE) and say so in `notes` |
| `direction` | `risk_increasing` if higher exposure → higher endometriosis risk; `protective` if lower; `null_or_unclear` otherwise |
| `coloc_reported`, `coloc_pph4` | Whether colocalisation was run for this gene; PP.H4 if given |
| `replication_reported`, `replication_details` | Whether the paper tested the gene in a second outcome GWAS, and what it found (source, estimate, p) in ≤ 30 words |
| `experimental_validation` | Short description of any wet-lab or patient-sample check (e.g. `ELISA serum, 30 cases vs 30 controls`), else `null` |
| `source_location` | Where the primary estimate appears, e.g. `Table 2`, `Figure 3B`, `Supplementary Table S5`, `Abstract` |
| `full_text_accessed` | `true` only if the main text (not just the abstract) was read |
| `extraction_by` | Who extracted, e.g. `claude-opus-5-5 agent, 2026-10-29` |
| `second_extraction` | Filled during the independent second extraction: `agree` or a note of disagreements |
| `human_verified` | `false` until Artem checks it |
| `notes` | Anything that affects interpretation: number of proteins tested, thresholds, sensitivity results, conversions, ambiguities |

## Rules

- Copy numbers exactly; never compute or infer a value that is not in the paper except the documented OR conversion.
- When the abstract and a table disagree, use the table and note the disagreement.
- If full text is not accessible, extract what the abstract reports, set `full_text_accessed = false`, and leave missing fields `null`.
- Do not judge whether a claim is plausible. Extraction records what the paper says.
- Papers with no eligible claim after full-text reading are listed with an exclusion reason instead.

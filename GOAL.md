# Endometriosis drug-target check

## Goal

Complete a reproducible, pre-specified check of published genetic drug-target claims for endometriosis. For each claim, establish whether it replicates in genetic data independent of the data it was found in, whether the gene signal and the disease signal share a causal variant, whether the gene is expressed in disease-relevant cells at the patient level, and whether an existing drug acts on it. Deliver a tiered target table that includes claims that fail or are inconclusive.

Why: a cure probably needs a non-hormonal, disease-modifying drug. Targets with human genetic support are more likely to succeed in trials, but the endometriosis literature now contains many small Mendelian randomization (MR) studies that name different genes, mostly from one cohort (FinnGen) and mostly without independent replication. Labs cannot test them all. A careful triage helps them choose.

Completion is defined by evidence and deliverables, not a date. A robust negative result (most claims do not hold up) or a documented feasibility limit is a valid outcome. Finding a cure, a positive signal, or journal acceptance is not a completion condition we can promise.

## Completion criteria

1. **Claim inventory complete.** Documented, dated literature search; screening with recorded reasons; extraction of every eligible claim with source locations; independent second extraction of headline claims; per-claim map of which cohorts the discovery result used.
2. **Data feasibility and sanity checks passed.** Replication summary statistics retrieved with verified checksums; genome builds and alleles harmonised; positive-control loci behave as expected in each replication source; one published MR estimate reproduced from its original discovery data, or the failure documented.
3. **Frozen replication analysis completed.** Protocol and claim table frozen (hashes recorded) before any replication lookup near a claimed gene. Per-claim estimates, uncertainty, power and a pre-defined classification: replicated, directionally consistent, not replicated, inconclusive, or not testable. Colocalisation where open exposure data allow.
4. **Cell-type evidence completed, or infeasibility established.** Donor-level (not cell-level) expression of claimed genes across cell types in at least two independent single-cell cohorts, with a protocol frozen before looking at disease contrasts.
5. **Interpretation and reviewable package completed.** Tiered target table with drug cross-reference, direction-of-effect check and limitations; associations kept separate from causal or therapeutic claims; contradictory and null results reported. Code, input manifest with checksums, full results, figures, and a manuscript draft if findings justify one, otherwise a closeout report naming the next useful experiment.

## Decision gates

| Gate | Evidence needed to proceed | Current state |
|---|---|---|
| A. Claim inventory | Search log, screening decisions, validated claim table, independence map | In progress: 488 claims; 189 headline claims double-extracted and adjudicated; instruments extracted and verified for 128 claims, re-derived from public exposure files for 29 (17 match the paper), unavailable for 148. Reference-list scan and supplementary full-text search done: 19 papers added, 15 included, 133 new claims (488 total), headline claims double-extracted. Claims near Gate B positive-control loci are labelled not blind and capped at Tier 2 (decided by Artem). Outstanding: instrument specs for the new claims, 2 unresolved gene positions, provisional screening rules, human review, eQTLGen access (61 claims), 9 missing full texts |
| B. Data and sanity checks | Checksummed replication data, harmonisation checks, positive controls, one reproduction | Mostly done: 5 outcome files and 19 exposure files verified; blinded reader in place; all 4 replication sources pass positive controls; EPHB4 estimate reproduced (OR 1.398 vs 1.40); LD reference checked against Ensembl. Rahmioglu 2023 cohort list verified: FinnGen not included; the shared NFBC legacy cohort (about 1%) accepted as near-independent |
| C. Frozen replication | Freeze record, per-claim results, power, classifications | Not started; tier rules and source choice approved by Artem for the freeze |
| D. Cell-type layer | Frozen single-cell plan, donor-level results in two cohorts | Not started |
| E. Interpretation | Tiered table, drug cross-reference, novelty check, proportional conclusions | Not started |
| F. Reproducible release | Clean rerun, reviewed report or manuscript, release package | Not started |

If a gate fails, document why and choose between a narrower question, other data, or a completed negative/feasibility report. Do not keep re-analysing until something becomes significant.

## Working boundaries

- Work is local and unpublished. No repository, preprint, submission or researcher outreach without Artem's explicit go-ahead.
- Public, de-identified summary data and public single-cell deposits only. Datasets that need an account, a signed data-use agreement or accepted terms (for example deCODE or UKB-PPP protein data, 23andMe, UK Biobank individual data) are not assumed available; Artem decides whether to accept any terms.
- Never attempt to identify individuals. Give no treatment recommendations. A replicated genetic association is a hypothesis for laboratory testing, not evidence that a drug helps patients.
- Record AI assistance (literature screening, extraction, code) and which items a human has verified, for eventual disclosure.
- Original authors own their data and findings; link sources rather than redistributing them.

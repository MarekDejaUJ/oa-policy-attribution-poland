# Open access because of policy: data and code

Aggregate data, analysis script and result tables for the article *Open access because of policy: Probabilities of causation for Polish transformative agreements and the open-access mandate of the National Science Centre* by Marek Deja (Institute of Information Studies, Jagiellonian University).

The study estimates the share of articles by authors from Polish universities that were licensed at the publisher because of the national transformative agreements and because of the open-access mandate of the National Science Centre (NCN). Probabilities of causation (PNS and PN) are bounded under counterfactual comparisons, negative-control margins, a publisher-choice region and a sampling projection.

## Source data

The counts derive from OpenAlex records retrieved on 9 October 2026 for 419,157 articles and reviews published in 2014–2025 with at most 50 authors and at least one author at one of the 20 universities eligible for the first competition of the Excellence Initiative – Research University programme. OpenAlex data are released under CC0. Only aggregate counts are included; article identifiers, author identifiers and provider responses are not part of the dataset.

Access routes follow the OpenAlex classification: `closed` (subscription only), `green` (free only in a repository), `bronze` (free at the publisher without a licence) and `licensed` (open at the publisher under a licence). The mandate design splits the licensed route into `ccby` (CC BY at the publisher) and `other_licence`. Shares are proportions between 0 and 1.

## Reproduction

Requirements: R (4.5 or later) with the packages `pnsbib` (version 0.2.0, https://github.com/MarekDejaUJ/pnsbib), `digest` and `jsonlite`. From the top folder of the dataset:

```
Rscript analysis.R data/aggregates /path/to/new-output-folder
```

The script refuses an existing output folder, computes every interval, writes the result tables, figures, `checks.tsv`, `manifest.json` and `sessionInfo.txt`, and stops with an error unless every check passes. Its result tables are identical to those in `data/results`.

## Inputs (`data/aggregates/`)

| File | Content |
|---|---|
| `agreement_counts.tsv` | Articles of 2016–2023 in journals outside the Directory of Open Access Journals published by the agreement publishers, by `publisher`, `agreement_group` (`2019` or `2020`), `corresponding_author` (`PL` or `foreign`), `year` and `agreement_year` (1 when the publisher had a national agreement in the publication year), with counts in the four routes and the total `n`; 50,798 articles |
| `venue_counts.tsv` | Journal articles of 2016–2023 by corresponding-author group and year: counts at the agreement publishers of each group, in DOAJ journals and elsewhere |
| `venue_outside_open.tsv` | Articles with Polish corresponding authors outside the agreement publishers, and the number open at the publisher, by year |
| `mandate_counts.tsv` | NCN-funded articles of 2022–2025 by `year`, grant `cohort` (`M` call rounds 38–41, under the mandate; `C` rounds 34–37, comparison; `P` rounds 30–33, earlier grants), NCN `panel` and coverage by an agreement, with counts in the five routes; 7,267 articles |
| `mandate_opus.tsv` | Articles and CC BY counts for articles whose only NCN grants come from OPUS 18, OPUS 19 or OPUS 20 |
| `mandate_licences.tsv` | Licences of the articles open at the publisher, by cohort |
| `measurement_check.tsv` | Agreement of the licensed route with a Creative Commons licence registered in Crossref, and of any open access with the Scopus open-access flag, on a random sample of 300 exposed articles |
| `agreement_ncn_funding.tsv` | Licensed shares of exposed and pre-agreement articles with and without NCN funding |

## Results (`data/results/`)

| File | Content |
|---|---|
| `agreement_inputs.tsv` | For each agreement group and year: publishers, exposed articles `n`, licensed count, exposed route shares `p_*`, counterfactual route shares under method (a) `qa_*` (level before the agreement) and method (b) `qb_*` (with the trend among foreign corresponding authors), negative-control margins `margin_a` and `margin_b`, the control year, the foreign change in the licensed share, the attracted share `attracted` and `qb_clipped` (TRUE when a route share of method (b) is truncated at zero) |
| `mandate_inputs.tsv` | For each publication year: exposed, comparison and earlier-grant articles, CC BY count, articles standardized across panels (`fallback_articles`), exposed shares `p_*`, counterfactual shares `qa_*` and `qb_*`, the earlier-grant CC BY share and the margin |
| `attribution_bounds.tsv` | PNS and PN intervals at point inputs from the exact linear program (`exact_*`), the Li–Pearl formulas (`lipearl_*`, with the theorem number and its sharpness label) and the Shu–Wang–Li closed form (`shu_*`) |
| `attribution_regions.tsv` | PNS and PN intervals for each policy, year and method, by `region`: `point` (identification interval at the counterfactual estimate), `margin` (outer interval over the negative-control margin), `margin_choice` (margin and publisher-choice region), `sampling` and `sampling_*` (simultaneous 95% sampling projection added to the region) |
| `pooled.tsv` | The same regions for the pooled exposed articles of each policy, with years as disjoint strata; `opened_lower` and `opened_upper` give the endpoints as numbers of articles |
| `breakdown.tsv` | Breakdown values: `breakdown_margin` (largest margin with a positive lower PNS endpoint), `ratio` (to the calibrated margin), `breakdown_attracted` (attracted share at which the lower endpoint reaches zero), Clopper–Pearson limits `p_lo` and `p_hi` of the exposed event share and the Hoeffding radius `q_halfwidth` of the counterfactual share |
| `agreement_monotone.tsv`, `mandate_monotone.tsv` | Values under monotonicity: compatibility of the binary and categorical models with the inputs, PNS and PN, the shares opened from each route (`from_*`), the observed change in each route (`change_*`) and the share with the event also without the policy (`retained`) |
| `flow_bounds.tsv` | Exact and Li–Pearl intervals for the share of exposed articles opened from each route, without monotonicity |
| `agreement_checks.tsv`, `mandate_checks.tsv`, `agreement_publisher_choice.tsv` | Checks fixed before the outcome shares were computed, and the publisher-choice measurement |
| `checks.tsv` | Verification checks of the calculation |
| `fig-design.data.tsv`, `fig-attribution.data.tsv`, `fig-flows.data.tsv` | Values drawn in the figures |
| `manifest.json`, `sessionInfo.txt` | Checksums of inputs and outputs, software versions |

## Licence

Data, script and documentation are released under the Creative Commons Attribution 4.0 International licence (`LICENSE`).

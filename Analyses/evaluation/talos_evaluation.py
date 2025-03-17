#! /usr/bin/env python3

"""
Script to evaluate the performance of Talos and Exomiser against a "Gold Standard" truth set.
"""

import argparse
import csv
import json
from typing import Optional

from cloudpathlib import CloudPath
from pydantic import BaseModel, Field, field_validator
from talos_models import ParticipantResults, ReportVariant, StructuralVariant

# Constants
VARIANT_TYPES_BASIC = {"SNV_INDEL", "CNV_SV"}
VARIANT_TYPES_ALL = {"SNV_INDEL", "CNV_SV", "STR", "MITO", "MITO_SV"}


class CausativeVariant(BaseModel):
    variant_id: str
    gene: str
    variant_type: str
    variant_description: Optional[str] = None
    talos_hit: Optional[ReportVariant] = None
    exomiser_hit: Optional[dict] = None


class Family(BaseModel):
    family_id: str
    individual_id: str
    sample_id: str
    external_id: str
    subcohort: str = ""
    trio: bool
    solved: bool
    causative_gene: Optional[str] = ""
    mosaic: Optional[bool] = False
    incomplete_penetrance: Optional[bool] = False
    intergenic: bool = False
    genotyping_error: bool = False
    variant1: Optional[CausativeVariant] = None
    variant2: Optional[CausativeVariant] = None

    analysed_by_talos: bool = False
    talos_results: list[ReportVariant] = Field(default_factory=list)
    _talos_phenotype_match_found: Optional[bool] = None

    exomiser_results: dict = {}
    exomiser_rank: Optional[int] = None

    # will be overwritten if results are found
    exomiser_results_exists: bool = False
    solved_by_exomiser: bool = False

    in_scope_variant_types: set = VARIANT_TYPES_BASIC
    in_scope_mosaics: bool = False
    in_scope_incomplete_penetrance: bool = False
    in_scope_intergenic: bool = False

    @field_validator("mosaic", "incomplete_penetrance", mode="before")
    def empty_str_to_false(cls, v: str) -> str | bool:
        if not v:
            return False
        return v

    @property
    def solved_by_talos(self):
        """Has this family been solved by talos?"""
        if not self.solved:
            return None
        if self.variant1 and not self.variant1.talos_hit:
            return False
        if self.variant2 and not self.variant2.talos_hit:
            return False
        return True

    @property
    def solved_in_scope(self):
        """Is this family considered solved AND should it be considered in scope?"""
        if not self.solved:
            return False
        if self.mosaic and not self.in_scope_mosaics:
            return False
        if self.incomplete_penetrance and not self.in_scope_incomplete_penetrance:
            return False
        if self.intergenic and not self.in_scope_intergenic:
            return False
        if self.has_out_of_scope_variant_type:
            return False
        return True

    @property
    def has_out_of_scope_variant_type(self):
        """Does this family have variants tyoes that are out of scope?"""
        if (
            self.variant1
            and self.variant1.variant_type not in self.in_scope_variant_types
        ):
            return True
        if (
            self.variant2
            and self.variant2.variant_type not in self.in_scope_variant_types
        ):
            return True
        return False

    @property
    def solved_by_talos_and_in_scope(self):
        """Is this family solved by talos and in scope"""
        return self.solved_by_talos and self.solved_in_scope

    @property
    def talos_phenotype_match_found(self):
        """Return a list of talos results that are in scope"""
        if self._talos_phenotype_match_found is not None:
            return self._talos_phenotype_match_found

        if self.variant1 and (talos_hit := self.variant1.talos_hit):
            if (
                talos_hit.phenotype_labels
                or talos_hit.panels.forced
                or talos_hit.panels.matched
            ):
                self._talos_phenotype_match_found = True
            else:
                self._talos_phenotype_match_found = False

        return self._talos_phenotype_match_found

    @property
    def talos_candidate_count(self):
        """Return the number of talos candidates"""
        if not self.talos_results:
            return None
        return len(self.talos_results.variants)

    @property
    def talos_candidate_count_w_phe_match(self):
        """Return the number of talos candidates with a phenotype match"""

        if not self.talos_results:
            return None

        return len(
            [
                r
                for r in self.talos_results.variants
                if r.phenotype_labels or r.panels.forced or r.panels.matched
            ]
        )

    @property
    def talos_candidate_gene_count(self):
        """Return the number of talos candidate GENES"""
        if not self.talos_results:
            return None

        return len(set([r.gene for r in self.talos_results.variants]))

    @property
    def talos_candidate_gene_count_w_phe_match(self):
        """Return the number of talos candidate GENES with a phenotype match"""
        if not self.talos_results:
            return None
        return len(
            set(
                [
                    r.gene
                    for r in self.talos_results.variants
                    if r.phenotype_labels or r.panels.forced or r.panels.matched
                ]
            )
        )

    @property
    def variant_types(self):
        """Return a set of variant types that are causative in this family"""
        return {
            x.variant_type
            for x in [self.variant1, self.variant2]
            if x and x.variant_type
        }

    def find_causative_variants_in_talos_resuts(self):
        """For each causative variant, look for a matching variant in the talos candidate list

        sets talos_hit attribute on the causative variant/s
        """
        if self.talos_results and self.talos_results.variants:
            if self.variant1:
                self.variant1.talos_hit = self.find_causative_variant_in_talos(
                    self.variant1
                )
            if self.variant2:
                self.variant2.talos_hit = self.find_causative_variant_in_talos(
                    self.variant2
                )

    def find_causative_variant_in_talos(
        self, causative_variant: CausativeVariant
    ) -> Optional[ReportVariant]:
        """
        For a given causative variant, find a matching variant in the talos candidate list
        """
        for r in self.talos_results.variants:
            if r.var_data.coordinates.string_format == causative_variant.variant_id:
                return r

            if causative_variant.variant_type == "CNV_SV" and isinstance(
                r.var_data, StructuralVariant
            ):
                # CNV_SV match based on gene symbol in the predicted_lof field
                if causative_variant.gene in r.var_data.info["predicted_lof"].split(
                    ","
                ):
                    return r
        return None

    def find_exomiser_results(self, exomiser_results: dict) -> None:
        """
        Look for an exomiser result for this family. If found, find ranks of causative variants

        If exomiser result file is found self.exomiser_results_exists is set to True

        If the causative variant/s are found in the exomiser results, the self.exomiser_rank
        is set to the highest rank.

        If only one causative variant is found for an AR condition, the family is considered
        unsolved by exomiser and the rank is set to None
        """

        self.exomiser_results = {}

        # result blocks in this file look like
        # "chr1:1234567:G:C": [
        #     {
        #         "rank": 1,
        #         "moi": "AD"
        #     },
        #     {
        #         "rank": 2,
        #         "moi": "AR"
        #     },
        # ],
        for key, values in exomiser_results.get(self.sample_id, {}).items():
            new_key = key.replace("chr", "").replace(":", "-")
            # Variants can appear multiple times in the exomiser results (eg AD/AR)
            # For this analysis we are only interested in presence/absense so we
            # will only keep the highest ranking variant
            self.exomiser_results[new_key] = min([entry["rank"] for entry in values])

        if self.exomiser_results:
            self.exomiser_results_exists = True

        # Find the rank of the causative variants in exomiser results
        if self.variant1:
            rank_1 = self.exomiser_results.get(self.variant1.variant_id)
        else:
            rank_1 = None

        if self.variant2:
            rank_2 = self.exomiser_results.get(self.variant2.variant_id)
        else:
            rank_2 = None

        # Set exomiser_rank if ALL causative variants are found
        # if two variants are found, use the highest rank
        if self.variant1 and self.variant2:
            if rank_1 is not None and rank_2 is not None:
                self.exomiser_rank = max(rank_1, rank_2)
                self.solved_by_exomiser = True
        elif self.variant1 and rank_1 is not None:
            self.exomiser_rank = rank_1
            self.solved_by_exomiser = True


def parse_truth_data(
    truth_tsv_path: str,
    process_only_trios: bool,
    subcohort_label: str,
    in_scope_variant_types: set,
    in_scope_mosaics: bool,
    in_scope_incomplete_penetrance: bool,
    in_scope_intergenic: bool,
):
    """
    Parse the truth data from tsv file

    Returns a list of Family objects
    """
    families = []
    for fam in csv.DictReader(CloudPath(truth_tsv_path).open(), delimiter="\t"):
        if fam["variant_1_type"]:
            v1 = CausativeVariant(
                variant_id=fam["variant_1"]
                .strip()
                .removeprefix("chr")
                .replace(":", "-"),
                gene=fam["causative_gene"],
                variant_type=fam["variant_1_type"],
                variant_description=fam["variant_1_description"],
            )
        else:
            v1 = None

        if fam["variant_2_type"]:
            v2 = CausativeVariant(
                variant_id=fam["variant_2"]
                .strip()
                .removeprefix("chr")
                .replace(":", "-"),
                gene=fam["causative_gene"],
                variant_type=fam["variant_2_type"],
                variant_description=fam["variant_2_description"],
            )
        else:
            v2 = None
        family = Family(**fam, variant1=v1, variant2=v2, subcohort=subcohort_label)

        family.in_scope_variant_types = in_scope_variant_types
        family.in_scope_mosaics = in_scope_mosaics
        family.in_scope_incomplete_penetrance = in_scope_incomplete_penetrance
        family.in_scope_intergenic = in_scope_intergenic

        if process_only_trios and not family.trio:
            continue

        # Add scope
        families.append(family)

    return families


def generate_summary_stats(families):
    total_families = len(families)
    families_analysed_by_talos = len([f for f in families if f.analysed_by_talos])
    solved_families = len([f for f in families if f.solved])
    solved_families_analysed_by_talos = len([f for f in families if f.solved and f.analysed_by_talos])
    solved_in_scope_families = len([f for f in families if f.solved_in_scope])
    solved_in_scope_families_analysed_by_talos = len([f for f in families if f.solved_in_scope and f.analysed_by_talos])
    solved_by_talos = len([f for f in families if f.solved_by_talos_and_in_scope])
    solved_by_talos_w_phen_match = len(
        [
            f
            for f in families
            if f.solved_by_talos_and_in_scope and f.talos_phenotype_match_found
        ]
    )

    pct_talos_solved_all = solved_by_talos / solved_families_analysed_by_talos * 100
    pct_talos_solved_in_scope = solved_by_talos / solved_in_scope_families * 100
    pct_talos_solved_w_phen_match = solved_by_talos_w_phen_match / solved_families_analysed_by_talos * 100
    pct_talos_solved_w_phen_match_in_scope = (
        solved_by_talos_w_phen_match / solved_in_scope_families * 100
    )

    summary_stats = f"""
        Total families: {total_families}
        Families considered solved: {solved_families} ({solved_families/total_families*100:.1f}%)
        Families considered solved and in scope: {solved_in_scope_families} ({solved_in_scope_families/solved_families*100:.1f}%

        Solved by talos: {solved_by_talos} ({pct_talos_solved_in_scope:.1f}% of in scope, {pct_talos_solved_all:.1f}% of all solved)
        Solved by talos with phenotype match: {solved_by_talos_w_phen_match} ({pct_talos_solved_w_phen_match_in_scope:.1f}% of in scope, {pct_talos_solved_w_phen_match:.1f}% of all solved)

        Average number of talos candidates per family: {sum([f.talos_candidate_count for f in families if f.talos_results]) / families_analysed_by_talos:.1f}
        Average number of talos candidates per family with phenotype match: {sum([f.talos_candidate_count_w_phe_match for f in families if f.talos_results]) / families_analysed_by_talos:.1f}
        Average number of talos candidate genes per family: {sum([f.talos_candidate_gene_count for f in families if f.talos_results]) / families_analysed_by_talos:.1f}
        Average number of talos candidate genes per family with phenotype match: {sum([f.talos_candidate_gene_count_w_phe_match for f in families if f.talos_results]) / families_analysed_by_talos:.1f}
    """

    # Exomiser stats
    total_solved_and_run_exomiser = len([f for f in families if (f.solved and f.exomiser_results_exists)])

    # bail if exomiser results not provided
    if not total_solved_and_run_exomiser:
        return summary_stats, ""

    total_solved_and_run_exomiser_snv_only = len(
        [
            f
            for f in families
            if f.solved
            and f.exomiser_results_exists
            and f.variant_types == {"SNV_INDEL"}
        ]
    )
    solved_by_exomiser = len([f for f in families if f.solved_by_exomiser])
    solved_by_exomiser_snv = len([f for f in families if f.solved_by_exomiser and f.variant_types == {"SNV_INDEL"}])
    solved_by_exomiser_top1 = len(
        [f for f in families if f.solved_by_exomiser and f.exomiser_rank == 1]
    )
    solved_by_exomiser_top5 = len(
        [f for f in families if f.solved_by_exomiser and f.exomiser_rank <= 5]
    )
    solved_by_exomiser_top10 = len(
        [f for f in families if f.solved_by_exomiser and f.exomiser_rank <= 10]
    )

    solved_by_exomiser_set = {f.family_id for f in families if f.solved_by_exomiser}
    solved_by_talos_set = {
        f.family_id
        for f in families
        if f.solved_by_talos and f.exomiser_results_exists
    }

    exomiser_summary = f"""
        Solved by exomiser: {solved_by_exomiser} ({solved_by_exomiser/total_solved_and_run_exomiser*100:.1f}%, {solved_by_exomiser_snv/total_solved_and_run_exomiser_snv_only*100:.1f}% of SNV only)
        Talos solved {len([f for f in families if (f.solved_by_exomiser and f.solved_by_talos)])} of these ({len([f for f in families if f.solved_by_exomiser is not None and f.solved_by_talos and f.variant_types == set(["SNV_INDEL"]) ])} SNV only)
        Solved by exomiser (top 1): {solved_by_exomiser_top1} ({solved_by_exomiser_top1/total_solved_and_run_exomiser*100:.1f}%, {solved_by_exomiser_top1/total_solved_and_run_exomiser_snv_only*100:.1f}% of SNV only)
        Solved by exomiser (top 5): {solved_by_exomiser_top5} ({solved_by_exomiser_top5/total_solved_and_run_exomiser*100:.1f}%), {solved_by_exomiser_top5/total_solved_and_run_exomiser_snv_only*100:.1f}% of SNV only)
        Solved by exomiser (top 10): {solved_by_exomiser_top10} ({solved_by_exomiser_top10/total_solved_and_run_exomiser*100:.1f}%, {solved_by_exomiser_top10/total_solved_and_run_exomiser_snv_only*100:.1f}% of SNV only)
        Exomiser solves not found by Talos: {len(solved_by_exomiser_set - solved_by_talos_set)}: {", ".join(solved_by_exomiser_set - solved_by_talos_set)}
        Solved by talos and not exomiser (SNV only): {len(solved_by_talos_set - solved_by_exomiser_set)}: {", ".join(solved_by_talos_set - solved_by_exomiser_set)}'
    """
    return summary_stats, exomiser_summary


def write_per_family_results(families, out_tsv):
    """
    Write a sumamry tsv with one line per family with all features. Allows easy
    filtering and sorting in a excel etc
    """

    fieldnames = [
        "family_id",
        "individual_id",
        "subcohort",
        "trio",
        "solved",
        "in_scope",
        "mosaic",
        "incomplete_penetrance",
        "intergenic",
        "causative_gene",
        "causative_variants",
        "causative_types",
        "solved_by_talos",
        "talos_phenotype_match_found",
        "notes",
        "exomiser_rank",
        "solved_by_exomiser",
        "talos_candidate_count",
        "talos_candidate_count_w_phe_match",
        "talos_candidate_gene_count",
        "talos_candidate_gene_count_w_phe_match",
    ]

    # with open(out_tsv, 'w') as output:
    writer = csv.DictWriter(out_tsv, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for family in families:
        writer.writerow(
            {
                "family_id": family.family_id,
                "individual_id": family.individual_id,
                "subcohort": family.subcohort,
                "trio": family.trio,
                "solved": family.solved,
                "in_scope": family.solved_in_scope,
                "mosaic": family.mosaic,
                "incomplete_penetrance": family.incomplete_penetrance,
                "intergenic": family.intergenic,
                "causative_gene": family.causative_gene,
                "causative_variants": ",".join(
                    [
                        x
                        for x in [
                            family.variant1.variant_id if family.variant1 else None,
                            family.variant2.variant_id if family.variant2 else None,
                        ]
                        if x
                    ]
                ),
                "causative_types": ",".join(
                    [
                        x
                        for x in [
                            family.variant1.variant_type if family.variant1 else None,
                            family.variant2.variant_type if family.variant2 else None,
                        ]
                        if x
                    ]
                ),
                "solved_by_talos": family.solved_by_talos,
                "talos_phenotype_match_found": family.talos_phenotype_match_found,
                "notes": "| ".join(
                    [
                        x
                        for x in [
                            (
                                family.variant1.variant_description
                                if family.variant1
                                else None
                            ),
                            (
                                family.variant2.variant_description
                                if family.variant2
                                else None
                            ),
                        ]
                        if x
                    ]
                ),
                "exomiser_rank": family.exomiser_rank,
                "solved_by_exomiser": family.solved_by_exomiser,
                "talos_candidate_count": family.talos_candidate_count,
                "talos_candidate_count_w_phe_match": family.talos_candidate_count_w_phe_match,
                "talos_candidate_gene_count": family.talos_candidate_gene_count,
                "talos_candidate_gene_count_w_phe_match": family.talos_candidate_gene_count_w_phe_match,
            }
        )


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_scope_variants",
        help="Variant types to include in the evaluation",
        choices=["all", "core"],
        default="core",
    )
    parser.add_argument(
        "--incomplete_penetrance",
        help="Consider incomplete penetrance variants as in scope",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--intergenic_in_scope",
        help="Consider intergenic variants as in scope",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--mosaic_in_scope",
        help="Consider mosaic variants as in scope",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--process_full_trios_only",
        help="process only trios",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--summary_tsv", help="Output path for summary tsv", type=argparse.FileType("w")
    )

    args, unknown = parser.parse_known_args()

    if args.in_scope_variants == "all":
        in_scope_variants = VARIANT_TYPES_ALL
    else:
        in_scope_variants = VARIANT_TYPES_BASIC

    if unknown:
        raise ValueError(unknown)
    main(
        in_scope_variants=in_scope_variants,
        process_only_trios=args.process_full_trios_only,
        in_scope_incomplete_penetrance=args.incomplete_penetrance,
        in_scope_intergenic=args.intergenic_in_scope,
        in_scope_mosaics=args.mosaic_in_scope,
        summary_tsv=args.summary_tsv,
    )


def main(
    in_scope_variants: set,
    process_only_trios: bool,
    in_scope_incomplete_penetrance: bool,
    in_scope_intergenic: bool,
    in_scope_mosaics: bool,
    summary_tsv,
):
    """ """

    # Should pass this as a config file but...
    sub_cohorts = {
        "genome": {
            "talos_results": "gs://cpg-acute-care-main/reanalysis/2025-01-22/pheno_annotated_report.json",
            "truth_tsv_path": "gs://cpg-acute-care-main-upload/talos_truth_data/240829_acute_care-genome-gold_std.tsv",
            "exomiser_results": "gs://cpg-acute-care-main-analysis/39c12fb9076d3e45a7b6a9c09aed7512dc2491_2405/exomiser_variant_results.json",
        },
        "exome": {
            "talos_results": "gs://cpg-acute-care-main/exome/reanalysis/2025-01-22/pheno_annotated_report.json",
            "truth_tsv_path": "gs://cpg-acute-care-main-upload/talos_truth_data/240822_acute_care-exome-gold_std.tsv",
            "exomiser_results": "gs://cpg-acute-care-main-analysis/exome/d671000b77331661dd38ec58250671b6807863_342/exomiser_variant_results.json",
        },
        # "rgp": {
        #     "talos_results": "gs://cpg-broad-rgp-main-upload/talos_truth_data/pheno_annotated_report_2024_10_29.json",
        #     "truth_tsv_path": "gs://cpg-broad-rgp-main-upload/talos_truth_data/241213_RGP_Data_For_Talos_Paper.tsv",
        # },
    }

    # Parse families from truth data
    families = []
    for subcohort_label, subcohort_dict in sub_cohorts.items():
        talos_results_json = json.load(
            CloudPath(subcohort_dict["talos_results"]).open()
        )

        sub_cohort_families = parse_truth_data(
            truth_tsv_path=subcohort_dict["truth_tsv_path"],
            process_only_trios=process_only_trios,
            subcohort_label=subcohort_label,
            in_scope_variant_types=in_scope_variants,
            in_scope_mosaics=in_scope_mosaics,
            in_scope_incomplete_penetrance=in_scope_incomplete_penetrance,
            in_scope_intergenic=in_scope_intergenic,
        )

        # look for exomiser results
        if exomiser_results_path := subcohort_dict.get("exomiser_results"):
            # single aggregate of all exomiser results
            exomiser_results = json.load(CloudPath(exomiser_results_path).open())
        else:
            exomiser_results = None

        # Annotate families with talos results
        for family in sub_cohort_families:
            # Allow for sample_id or individual_id to be used as the key in the talos results
            if family.sample_id in talos_results_json["results"]:
                use_id = family.sample_id
            elif family.individual_id in talos_results_json["results"]:
                use_id = family.individual_id
            else:
                family.talos_results = None
                continue

            family.analysed_by_talos = True

            # Add talos results to family
            family.talos_results = ParticipantResults.model_validate(
                talos_results_json["results"][use_id]
            )
            family.find_causative_variants_in_talos_resuts()

            # look for exomiser results
            if exomiser_results is not None:
                family.find_exomiser_results(exomiser_results)

        families.extend(sub_cohort_families)

    # Write summary stats to stdout
    summary_stats, exomiser_summary = generate_summary_stats(families)
    print(summary_stats)
    if exomiser_summary:
        print("\n\n", "#Exomiser summary:\n", exomiser_summary)

    # Write per family results to tsv
    if summary_tsv:
        write_per_family_results(families, summary_tsv)


if __name__ == "__main__":
    cli_main()

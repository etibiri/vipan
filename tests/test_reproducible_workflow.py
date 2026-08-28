from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_script(name: str, *arguments: object) -> subprocess.CompletedProcess[str]:
    command = [sys.executable, str(ROOT / "scripts" / name), *map(str, arguments)]
    return subprocess.run(command, capture_output=True, text=True, check=False)


class ReproducibleWorkflowTests(unittest.TestCase):
    def test_prepare_sequences_audits_duplicates_and_short_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            raw = directory / "raw.fasta"
            raw.write_text(
                ">A description\nTAATATTACAAAA\n"
                ">B\nTAATATTACAAAA\n"
                ">C\nTAATATTACAAAT\n"
                ">D\nTAATATTACAATA\n"
                ">E\nTAAT\n",
                encoding="utf-8",
            )
            result = run_script(
                "prepare_sequences.py",
                "--input", raw,
                "--output-fasta", directory / "kept.fasta",
                "--excluded-fasta", directory / "excluded.fasta",
                "--manifest", directory / "manifest.tsv",
                "--min-length", 10,
                "--orient-by-motif",
                "--deduplicate-exact",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (directory / "manifest.tsv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            reasons = {row["accession"]: row["reason"] for row in rows}
            self.assertEqual(reasons["B"], "exact_duplicate")
            self.assertEqual(reasons["E"], "sequence_too_short")

    def test_alignment_filter_removes_a_distant_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            alignment = directory / "alignment.fasta"
            alignment.write_text(
                ">A\nAAAAAAAAAACCCCCCCCCC\n"
                ">B\nAAAAAAAAATCCCCCCCCCC\n"
                ">C\nAAAAAAAATTCCCCCCCCCC\n"
                ">D\nTTTTTTTTTTGGGGGGGGGG\n",
                encoding="utf-8",
            )
            result = run_script(
                "filter_alignment_outliers.py",
                "--alignment", alignment,
                "--kept-fasta", directory / "kept.fasta",
                "--kept-aligned-fasta", directory / "kept_aligned.fasta",
                "--excluded-fasta", directory / "excluded.fasta",
                "--metrics", directory / "metrics.tsv",
                "--pairwise-distances", directory / "distances.tsv",
                "--summary-json", directory / "summary.json",
                "--plot-prefix", directory / "filter",
                "--mode", "medoid",
                "--max-p-distance", 0.20,
                "--min-comparable-sites", 10,
                "--min-overlap-fraction", 0.5,
                "--min-valid-bases", 10,
                "--max-gap-fraction", 0.5,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (directory / "metrics.tsv").open(newline="", encoding="utf-8") as handle:
                rows = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["D"]["status"], "excluded")
            self.assertEqual(rows["D"]["reason"], "too_distant_from_anchor")
            self.assertEqual(sum(row["status"] == "kept" for row in rows.values()), 3)

    def test_rdp_parser_expands_continuation_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            alignment = directory / "alignment.fasta"
            accessions = ["NC_000001.1", "NC_000002.1", "NC_000003.1", "NC_000004.1"]
            alignment.write_text(
                "".join(f">{accession}\n{'A' * 20}\n" for accession in accessions),
                encoding="utf-8",
            )
            report = directory / "rdp.csv"
            rows = [
                ["", "", "", "Breakpoint Positions"],
                ["", "", "In Alignment"],
                [
                    "Recombination Event Number", "Number In .RDP File", "Begin", "End",
                    "Begin", "End", "Begin", "End", "Recombinant Sequence(s)",
                    "Minor Parental Sequence(s)", "Major Parental Sequence(s)",
                    "RDP", "GENECONV", "Bootscan", "Maxchi", "Chimaera", "SiSscan",
                    "PhylPro", "LARD", "3Seq",
                ],
                [
                    "1", "1", "3", "15", "3", "15", "3", "15",
                    "NC_000001.1", "NC_000002.1", "NC_000003.1",
                    "1e-5", "1e-4", "1e-3", "NS", "NS", "NS", "NS", "NS", "1e-6", "",
                ],
                ["1", "", "", "", "", "", "", "", "NC_000004.1", "", ""],
                [],
                ["Table Key:"],
            ]
            with report.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle).writerows(rows)
            result = run_script(
                "parse_rdp_csv.py",
                "--input", report,
                "--alignment", alignment,
                "--output-all", directory / "all.tsv",
                "--output-robust", directory / "robust.tsv",
                "--min-methods", 3,
                "--min-tested-methods", 7,
                "--min-support-fraction", 0.3,
                "--min-fragment-length", 5,
                "--require-known-parents",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (directory / "robust.tsv").open(newline="", encoding="utf-8") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(parsed), 1)
            self.assertEqual(parsed[0]["number_recombinants"], "2")
            self.assertIn("NC_000004.1", parsed[0]["recombinant_accessions"])

    def test_strict_classifier_uses_bipartition_size(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            calls = directory / "calls.tsv"
            fields = [
                "Evidence_tier", "Recombinant_accession", "3SEQ_Event_ID", "RDP_Event_ID",
                "Breakpoint_scheme", "Different_nearest_parents", "Expected_role_switch",
                "Role_switch_evidence",
            ]
            with calls.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "Evidence_tier": "Tier_B", "Recombinant_accession": "NC_000001.1",
                        "3SEQ_Event_ID": "3S1", "RDP_Event_ID": "R1", "Breakpoint_scheme": "RDP",
                        "Different_nearest_parents": "True", "Expected_role_switch": "True",
                        "Role_switch_evidence": "RDP",
                    }
                )
            trees = directory / "trees.tsv"
            tree_fields = [
                "Recombinant_accession", "3SEQ_Event_ID", "RDP_Event_ID", "Breakpoint_scheme",
                "Segment", "Nearest_parent_SH_aLRT", "Nearest_parent_UFBoot",
                "Nearest_parent_supported_bipartition_size",
            ]
            with trees.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=tree_fields, delimiter="\t")
                writer.writeheader()
                for segment, size in (("A_between_breakpoints", 3), ("B_complement", 2)):
                    writer.writerow(
                        {
                            "Recombinant_accession": "NC_000001.1", "3SEQ_Event_ID": "3S1",
                            "RDP_Event_ID": "R1", "Breakpoint_scheme": "RDP", "Segment": segment,
                            "Nearest_parent_SH_aLRT": 99, "Nearest_parent_UFBoot": 100,
                            "Nearest_parent_supported_bipartition_size": size,
                        }
                    )
            result = run_script(
                "classify_validation_calls.py",
                "--calls", calls,
                "--tree-summary", trees,
                "--output-calls", directory / "strict.tsv",
                "--output-events", directory / "events.tsv",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (directory / "strict.tsv").open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["Strict_validation_class"], "A_local_supported_parent_switch")

    def test_filtering_strategy_comparison_flags_discordant_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            fields = [
                "accession", "status", "reason", "anchor_accession",
                "distance_to_anchor", "comparable_sites_to_anchor",
            ]
            medoid = directory / "medoid.tsv"
            reference = directory / "reference.tsv"
            with medoid.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "accession": "A", "status": "kept", "reason": "passed_all_filters",
                            "anchor_accession": "A", "distance_to_anchor": "0.0",
                            "comparable_sites_to_anchor": "2800",
                        },
                        {
                            "accession": "B", "status": "kept", "reason": "passed_all_filters",
                            "anchor_accession": "A", "distance_to_anchor": "0.10",
                            "comparable_sites_to_anchor": "2750",
                        },
                        {
                            "accession": "C", "status": "excluded", "reason": "too_distant_from_anchor",
                            "anchor_accession": "A", "distance_to_anchor": "0.30",
                            "comparable_sites_to_anchor": "2700",
                        },
                    ]
                )
            with reference.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "accession": "A", "status": "kept", "reason": "passed_all_filters",
                            "anchor_accession": "NC_004650.1", "distance_to_anchor": "0.10",
                            "comparable_sites_to_anchor": "2800",
                        },
                        {
                            "accession": "B", "status": "excluded", "reason": "too_distant_from_anchor",
                            "anchor_accession": "NC_004650.1", "distance_to_anchor": "0.30",
                            "comparable_sites_to_anchor": "2750",
                        },
                        {
                            "accession": "C", "status": "excluded", "reason": "too_distant_from_anchor",
                            "anchor_accession": "NC_004650.1", "distance_to_anchor": "0.35",
                            "comparable_sites_to_anchor": "2700",
                        },
                    ]
                )
            output = directory / "comparison"
            result = run_script(
                "compare_filtering_strategies.py",
                "--medoid", medoid,
                "--reference", reference,
                "--output-dir", output,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (output / "discordant_accessions.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["Accession"], "B")
            self.assertEqual(rows[0]["Decision_concordance"], "medoid_only")

    def test_validation_strategy_comparison_detects_class_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            fields = [
                "Recombinant_accession",
                "Best_breakpoint_scheme",
                "Best_strict_validation_class",
            ]
            medoid = directory / "medoid_strict.tsv"
            reference = directory / "reference_strict.tsv"
            with medoid.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "Recombinant_accession": "A",
                        "Best_breakpoint_scheme": "RDP",
                        "Best_strict_validation_class": "A_local_supported_parent_switch",
                    }
                )
            with reference.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "Recombinant_accession": "A",
                        "Best_breakpoint_scheme": "3SEQ",
                        "Best_strict_validation_class": "D_expected_switch_insufficient_branch_support",
                    }
                )
            output = directory / "validation_comparison"
            result = run_script(
                "compare_validation_strategies.py",
                "--medoid", medoid,
                "--reference", reference,
                "--output-dir", output,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (output / "strict_validation_comparison.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["Validation_concordance"], "different_class")


if __name__ == "__main__":
    unittest.main()

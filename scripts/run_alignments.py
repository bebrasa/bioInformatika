from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO, pairwise2


ROOT = Path(__file__).resolve().parent.parent
FASTA_DIR = ROOT / "data" / "fasta"
OUT_DIR = ROOT / "results" / "alignments"


@dataclass
class PairJob:
    gene: str
    human_fasta: Path
    mouse_fasta: Path


def read_single_fasta(path: Path) -> tuple[str, str]:
    record = next(SeqIO.parse(path, "fasta"))
    return record.id, str(record.seq).upper()


def calc_metrics(aln_a: str, aln_b: str) -> dict[str, float]:
    matches = sum(1 for a, b in zip(aln_a, aln_b) if a == b and a != "-" and b != "-")
    aligned_non_gap = sum(1 for a, b in zip(aln_a, aln_b) if a != "-" and b != "-")
    gaps = sum(1 for a, b in zip(aln_a, aln_b) if a == "-" or b == "-")
    identity = (matches / aligned_non_gap * 100) if aligned_non_gap else 0.0
    return {
        "alignment_length": float(len(aln_a)),
        "aligned_non_gap": float(aligned_non_gap),
        "matches": float(matches),
        "gaps": float(gaps),
        "identity_percent": identity,
    }


def write_result(
    out_path: Path,
    method_name: str,
    score: float,
    human_id: str,
    mouse_id: str,
    aln_human: str,
    aln_mouse: str,
) -> None:
    metrics = calc_metrics(aln_human, aln_mouse)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    snippet_len = 700
    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"Method: {method_name}\n")
        f.write(f"Score: {score}\n")
        f.write(f"Human: {human_id}\n")
        f.write(f"Mouse: {mouse_id}\n\n")
        f.write(f"Alignment length: {int(metrics['alignment_length'])}\n")
        f.write(f"Aligned non-gap positions: {int(metrics['aligned_non_gap'])}\n")
        f.write(f"Matches: {int(metrics['matches'])}\n")
        f.write(f"Gaps: {int(metrics['gaps'])}\n")
        f.write(f"Identity (%): {metrics['identity_percent']:.2f}\n\n")
        f.write("Alignment snippet (first 700 columns):\n")
        f.write(aln_human[:snippet_len] + "\n")
        f.write(aln_mouse[:snippet_len] + "\n")


def main() -> None:
    jobs = [
        PairJob(
            gene="BRCA1",
            human_fasta=FASTA_DIR / "BRCA1_Homo_sapiens_NM_007294.4.fasta",
            mouse_fasta=FASTA_DIR / "Brca1_Mus_musculus_NM_009764.fasta",
        ),
        PairJob(
            gene="BRCA2",
            human_fasta=FASTA_DIR / "BRCA2_Homo_sapiens_NM_000059.4.fasta",
            mouse_fasta=FASTA_DIR / "Brca2_Mus_musculus_NM_009765.fasta",
        ),
    ]

    for job in jobs:
        human_id, human_seq = read_single_fasta(job.human_fasta)
        mouse_id, mouse_seq = read_single_fasta(job.mouse_fasta)

        global_aln = pairwise2.align.globalms(
            human_seq,
            mouse_seq,
            2,
            -1,
            -5,
            -0.5,
            one_alignment_only=True,
        )[0]
        write_result(
            OUT_DIR / f"{job.gene}_global_needleman_wunsch.txt",
            "Needleman-Wunsch (global)",
            global_aln.score,
            human_id,
            mouse_id,
            global_aln.seqA,
            global_aln.seqB,
        )

        local_aln = pairwise2.align.localms(
            human_seq,
            mouse_seq,
            2,
            -1,
            -5,
            -0.5,
            one_alignment_only=True,
        )[0]
        write_result(
            OUT_DIR / f"{job.gene}_local_smith_waterman.txt",
            "Smith-Waterman (local)",
            local_aln.score,
            human_id,
            mouse_id,
            local_aln.seqA,
            local_aln.seqB,
        )


if __name__ == "__main__":
    main()

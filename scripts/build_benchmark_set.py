"""Write examples/clinvar_missense_benchmark.vcf: a fixed, dated benchmark set.

    python scripts/build_benchmark_set.py

For each gene, the best-reviewed ClinVar missense variants on each side of the
classification (Pathogenic/Likely pathogenic vs Benign/Likely benign), up to
PER_CLASS each. `var2fold.clinvar` ranks by review stars and spreads ties
along the gene. Variants without a clear aggregate call are dropped. The VCF
keeps ClinVar's classification in INFO/CLNSIG, so the benchmark notebook reads
its labels from the file and the set does not drift as ClinVar updates.
"""

from __future__ import annotations

import logging
from pathlib import Path

from var2fold.clinvar import ClinVarClient, write_vcf
from var2fold.validation import clinvar_label

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "examples" / "clinvar_missense_benchmark.vcf"
GENES = ("BRCA1", "BRCA2", "TP53", "MSH2", "MLH1", "PTEN")
PER_CLASS = 25
SIDES = {1: ("pathogenic", "likely_pathogenic"), 0: ("benign", "likely_benign")}


def main() -> None:
    logging.basicConfig(level=logging.WARNING)
    client = ClinVarClient()
    selected = []
    for gene in GENES:
        counts = {}
        for label, significance in SIDES.items():
            found = client.variants_for_gene(
                gene, consequence="missense", significance=significance, max_variants=None
            )
            # Keep only variants whose aggregate call is on the queried side.
            usable = [cv for cv in found if clinvar_label(cv.classification) == label]
            pick = _take(usable, PER_CLASS)
            counts[label] = len(pick)
            selected.extend(pick)
        print(f"{gene}: {counts[1]} pathogenic, {counts[0]} benign")
    write_vcf(selected, OUT, source="ClinVar missense benchmark")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(selected)} variants)")


def _take(variants, n):
    """The first n in ClinVar ranking order (stars first), spread along the gene."""
    if len(variants) <= n:
        return variants
    top_stars = sorted(variants, key=lambda cv: -cv.stars)
    cutoff = top_stars[n - 1].stars
    sure = [cv for cv in top_stars if cv.stars > cutoff]
    ties = sorted((cv for cv in variants if cv.stars == cutoff), key=lambda cv: cv.variant.position)
    room = n - len(sure)
    step = len(ties) / room
    return sure + [ties[int(i * step + step / 2)] for i in range(room)]


if __name__ == "__main__":
    main()

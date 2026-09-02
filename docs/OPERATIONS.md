# Operations

## ClinVar mirror v2 rebuild

The workers read the mirror path from `CLINVAR_MIRROR_DB`.

A schema v1 file makes the screen use the bundled panel.

The provider status then states that schema v2 requires a rebuild.

Run `python -m geneask.annotators.clinvar_mirror` on alien02.

Copy the completed database to the NAS hub and alien03.

Follow the copy procedure in `memory/dna_report_gnomad_mirror_nearly_lost.md`.

Do not rebuild the mirror over the WAN.

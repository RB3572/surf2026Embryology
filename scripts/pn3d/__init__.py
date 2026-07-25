"""
3D Pronuclear Pseudotime Calibration.

A structure-aware pipeline that estimates a calibrated developmental pseudotime
tau in [0,1] for a static 3D zygote stack:

    image  ->  audited 3D segmentation (cell body, PN1, PN2, polar body)
           ->  dimensionless geometry
           ->  monotonic probabilistic clock (calibrated on live-imaging tau)
           ->  tau posterior + interval + QC / OOD.

The mounted MEFISH dataset is FIXED MERFISH (target imaging domain) with 3D
Slicer segmentations; it carries no independent time labels. True developmental
time comes from the external Scheffler 2021 live-imaging trajectories. Segment
biological meaning is established GEOMETRICALLY, never from the (generic) file /
segment names.
"""

__version__ = "pn3d-0.1.0"
SCHEMA_VERSION = 1

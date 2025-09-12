# tiger-paw
# Tandem and Intricate GEnomic Repeat - Puzzle piece Assembly and annotation Workflow

Workflow for very accurately closing gaps caused by large tandem repeats (e.g. 400 copies of 10,000 bp repeat) in eukaryotic genome assemblies.

Used to assemble the two Nucleolus Organizer Regions (NORs) in [The Telomere to Telomere Gapless Lettuce (_lactuca sativa c. salinas_) Genome Assembly](https://kittishgames.com/pounce/).

**Background**: 
* In large tandem repeats, real biological variation between adjacent repeats is often vastly outweighed by sequencing errors, rendering the regions unresolveable by classical genome assembly algorithms (DeBrujin Graphs, String Graphs).
* In many modern genome assemblies, these regions are left unresolved and gaps are filled with "N"s.
* This workflow employs custom methods converting this problem from "impossible" to "very difficult".

**Notes**:
This workflow is intended to produce a near-perfect sequence from an input of ONT reads in combination with HiFi or Illumina reads and a large time investment in manual alignment of "pseudoreads".


**Steps**: 

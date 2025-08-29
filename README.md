# traw - tandem repeat assembly workflow
Workflow for closing large gaps caused by large tandem repeats in eukaryotic genome assemblies.

(...link to paper...)
Used to assemble the two Nucleolus Organizer Regions in Telomere to Telomere Gapless Lettuce (_lactuca sativa c. salinas_) Genome.

**Background**: 
* In large tandem repeats, real biological variation between adjacent repeats is often vastly outweighed by sequencing errors, rendering the regions unresolveable by classical genome assembly algorithms (DeBrujin Graphs, String Graphs).
* In many modern genome assemblies, these regions are left unresolved and gaps are filled with "N"s.
* This workflow employs custom methods converting this problem from "impossible" to "very difficult".

**Notes**:
This workflow is intended to produce a near-perfect sequence from an input of ONT reads in combination with HiFi or Illumina reads and a large time investment in manual alignment of "pseudoreads".


**Steps**: 

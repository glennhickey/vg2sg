# vg2sg
Prototype code for converting [VG](https://github.com/ekg/vg) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)

(c) 2015 Glenn Hickey. See [LICENSE](https://github.com/glennhickey/hal2sg/blob/development/LICENSE) for details.

See also:
* [sg2vg](https://github.com/glennhickey/sg2vg): Convert [Global Alliance (Side Graph) Server](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format) to [VG](https://github.com/ekg/vg)
* [hal2sg](https://github.com/glennhickey/hal2sg): Convert [HAL](https://github.com/glennhickey/hal) (output by [Cactus](https://github.com/glennhickey/progressiveCactus) and[CAMEL](https://github.com/adamnovak/sequence-graphs)) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)

## Algorithm

Iteratatively add VG paths to side graph.  Consecutive VG nodes will be merged greedily when possible.  Option to generate paths covering all edges in VG input not already in a path to ensure all nodes and edges get converted. 

**Constraints**
1. Edits in paths are not supported.  Any non-trivial (snp/indel) edits will result in an error. 

## Important

This is organized as a stand-alone exectuable (ie reads VG protbuf directly) mostly so I can start developing on my mac (where VG won't build).  Chances are, some logic from VG will be needed eventually. At that point will have to look at either integrating into VG or linking against it (and its millions of deps)...

## Instructions

**Cloning:** Don't forget to clone submodules:

     git clone https://github.com/glennhickey/vg2sg.git --recursive

**Note** Start by verifying that the unit tests all pass:

	  make test

To run the converter:

	  vg2sg input.vg output.fa output.sql

`input.vg` Input variant graph to convert

`output.fa` Output fasta file of all Side Graph sequences

`output.sql` Output text file listing INSERT commands for Sequences, Joins and Paths (for each input sequence) in the graph.

To see all the options, run with no args or use `--help`.

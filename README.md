# vg2sg
Prototype code for converting [VG](https://github.com/ekg/vg) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)

(c) 2015 Glenn Hickey. See [LICENSE](https://github.com/glennhickey/hal2sg/blob/development/LICENSE) for details.

## Algorithm

Iteratatively add VG paths to side graph.  Consecutive VG nodes will be merged greedily when possible. 

## Instructions

**Dependencies:**    [VG](https://github.com/ekg/vg).  Expected to be in sister directory to vg2sg but can be changed in include.mk. 

**Cloning:** Don't forget to clone submodules:

     git clone git@github.com:glennhickey/vg2sg.git --recursive
or
     git clone https://github.com/glennhickey/vg2sg.git --recursive


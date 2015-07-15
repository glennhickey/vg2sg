# vg2sg
Prototype code for converting [VG](https://github.com/ekg/vg) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)

(c) 2015 Glenn Hickey. See [LICENSE](https://github.com/glennhickey/hal2sg/blob/development/LICENSE) for details.

## Algorithm

Iteratatively add VG paths to side graph.  Consecutive VG nodes will be merged greedily when possible. 

## Important

This is organized as a stand-alone exectuable (ie reads VG protbuf directly) mostly so I can start developing on my mac (where VG won't build).  Chances are, some logic from VG will be needed eventually (from what I understand from vg.proto, graphs > 64MB need to be merged, for example). At that point will have to look at either integrating into VG or linking against it (and its millions of deps)...

## Instructions

**Cloning:** Don't forget to clone submodules:

     git clone git@github.com:glennhickey/vg2sg.git --recursive
or
     git clone https://github.com/glennhickey/vg2sg.git --recursive


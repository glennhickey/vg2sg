/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _PATHMAPPER_H
#define _PATHMAPPER_H

#include <string>
#include <map>

#include "vglight.h"
#include "sglookup.h"
#include "sidegraph.h"

/** iteratively map vg paths onto a sidegraph
 */
class PathMapper
{
public:
   PathMapper();
   ~PathMapper();

   /** load the vg */
   void init(const VGLight* vg);

   /** get the side graph we made */
   const SideGraph* getSideGraph() const;

   /** get a chunk of DNA sequence from the side graph */
   std::string getSideGraphDNA(sg_int_t seqID, sg_int_t offset = 0,
                               sg_int_t length = -1, bool reversed = false);

   /** add a path by name (leave control of order of addition to 
    * calling code) */
   void addPath(const std::string& name);

   /** copied from halCommon.h -- dont want hal dep just for this*/
   static char reverseComplement(char c);
   static void reverseComplement(std::string& s);

protected:

   void addSegment(sg_int_t pathID, sg_int_t pathPos,
                   const vg::Position& pos, bool reversed,
                   sg_int_t segLength);

   std::string makeSeqName(sg_int_t pathID, sg_int_t pathPos);

   const std::string& getPathName(sg_int_t id) const;
   sg_int_t getPathID(const std::string& name) const;

   SideGraph* _sg;
   SGLookup* _lookup;
   const VGLight* _vg;
   std::vector<std::string> _pathNames;
   std::map<std::string, sg_int_t> _pathIDs;
   std::vector<std::string> _seqStrings;
   // make sure node ids are in range [0, numNodes)
   // note to self- some of the maps should be hash tables
   std::map<int64_t, sg_int_t> _nodeIDMap;
   SGSequence* _curSeq;
};

inline const SideGraph* PathMapper::getSideGraph() const
{
  return _sg;
}

inline const std::string& PathMapper::getPathName(sg_int_t id) const
{
  assert(id >=0 && id < _pathNames.size());
  return _pathNames[id];
}

inline sg_int_t PathMapper::getPathID(const std::string& name) const
{
  std::map<std::string, sg_int_t>::const_iterator i = _pathIDs.find(name);
  assert(i != _pathIDs.end());
  return i->second;
}

#endif

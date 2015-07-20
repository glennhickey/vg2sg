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
                               sg_int_t length = -1, bool reversed = false)
     const;

   /** get the DNA sequence of a path in the *side graph* */
   std::string getSideGraphPathDNA(const std::string& pathName) const;

   /** get the path in the Side Graph that corresponds to an added VG path
    */
   const std::vector<SGSegment>& getSideGraphPath(
     const std::string& vgPathName) const;

   /** get the name of the VG path from which a Side Graph sequence was
    * derived */
   const std::string& getVGPathName(const SGSequence* seq) const;

   /** add a path by name (leave control of order of addition to 
    * calling code) */
   void addPath(const std::string& name);

   /** get name of path (index is the order in which the VG path
    * was added) */
   const std::string& getPathName(sg_int_t id) const;

   /** number of paths that were added using addPath */
   size_t getNumPaths() const;

   /** throw an exception if side graph path's dna doesn't jive with
    * vg path's dna*/
   void verifyPaths() const;

protected:

   /** add a segment corresponding to an input node */
   void addSegment(sg_int_t pathID, sg_int_t pathPos,
                   const vg::Position& pos, bool reversed,
                   sg_int_t segLength);

   /** second pass of input path to compute joins and side graph
    * segments. */
   void addPathJoins(const std::string& name,
                     const VGLight::MappingList& mappings);

   /** append path onto the end of prevPath, merging the last segment
    * of prevPath with first segment of nextPath if possible */
   void mergePaths(std::vector<SGSegment>& prevPath,
                   const std::vector<SGSegment>& path) const;


   std::string makeSeqName(sg_int_t pathID, sg_int_t pathPos);


   sg_int_t getPathID(const std::string& name) const;

   SideGraph* _sg;
   SGLookup* _lookup;
   const VGLight* _vg;
   std::vector<std::string> _pathNames;
   std::map<std::string, sg_int_t> _pathIDs;
   std::vector<std::string> _seqStrings;
   std::vector<std::vector<SGSegment> > _sgPaths;
   // make sure node ids are in range [0, numNodes)
   // note to self- some of the maps should be hash tables
   std::map<int64_t, sg_int_t> _nodeIDMap;
   SGSequence* _curSeq;
   std::vector<sg_int_t> _sgSeqToVGPathID;
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

inline size_t PathMapper::getNumPaths() const
{
  return _pathNames.size();
}

inline sg_int_t PathMapper::getPathID(const std::string& name) const
{
  std::map<std::string, sg_int_t>::const_iterator i = _pathIDs.find(name);
  assert(i != _pathIDs.end());
  return i->second;
}

inline std::ostream& operator<<(std::ostream& os,
                                const std::vector<SGSegment>& segs)
{
  os << "(";
  for (int i = 0; i < segs.size(); ++i)
  {
    os << segs[i] << (i < segs.size() - 1 ? "," : ")");
  }
  return os;
}

#endif

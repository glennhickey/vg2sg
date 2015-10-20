/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _VGLIGHT_H
#define _VGLIGHT_H

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <stdexcept>
#include "vg.pb.h"

/*
 * Lightweight wrapper to get a VG graph out of a protobuf stream.  Written
 * for prototyping only, and with the intention of replacing with actual
 * VG class from vg.hpp in near term
 */
class VGLight
{
public:
   VGLight();
   virtual ~VGLight();
   
   /** Read a graph in from a protobuf stream
    */
   void loadGraph(std::istream& inStream);

   /** Copy in a graph
    */
   void loadGraph(const vg::Graph& graph);

   struct NodePtrLess {
      bool operator()(const vg::Node* node1, const vg::Node* node2) const;
   };
   typedef std::list<vg::Mapping> MappingList;
   typedef std::set<const vg::Node*, NodePtrLess> NodeSet;
   typedef std::multimap<int64_t, const vg::Edge*> EdgeMap;
   typedef std::map<std::string, MappingList> PathMap;

   const NodeSet& getNodeSet() const;
   const PathMap& getPathMap() const;
   const EdgeMap& getFromEdgeMap() const;
   const EdgeMap& getToEdgeMap() const;
   size_t getNumEdges() const;

   const vg::Node* getNode(int64_t id) const;
   const vg::Edge* getEdge(int64_t from_id, int64_t to_id,
                           bool from_start, bool to_end) const;
   const MappingList& getPath(const std::string& name) const;
   void removePath(const std::string& name);

   /* all edges that touch node in either direction */
   void getInEdges(const vg::Node* node,
                         std::vector<const vg::Edge*>& ins) const;
   void getOutEdges(const vg::Node* node,
                         std::vector<const vg::Edge*>& outs) const;

   /** get string for a VG path */
   void getPathDNA(const std::string& name, std::string& outDNA) const;
   void getPathDNA(const MappingList& mappingList, std::string& outDNA) const;

   /** get length of a path segment */
   int64_t getSegmentLength(const vg::Mapping& mapping) const;

   /** copied from halCommon.h -- dont want hal dep just for this*/
   static char reverseComplement(char c);
   static void reverseComplement(std::string& s);
   
protected:

   /** Convert _graphs into _nodes/_edges/_paths */
   void mergeGraphs();
       
   /** Graphs read direct from protobuf.  They are unmerged and we keep
    * them around just for storage */
   std::vector<vg::Graph> _graphs;


   /** Split out graph compoments into sets to (hopefully) handle merging. 
    * these are what we access during conversion */
   NodeSet _nodes;
   EdgeMap _fromEdges;
   EdgeMap _toEdges;
   PathMap _paths;
   size_t _numEdges;
};

inline bool VGLight::NodePtrLess::operator()(const vg::Node* node1,
                                             const vg::Node* node2) const
{
  return node1->id() < node2->id();
}

inline const VGLight::NodeSet& VGLight::getNodeSet() const
{
  return _nodes;
}

inline const VGLight::EdgeMap& VGLight::getFromEdgeMap() const
{
  return _fromEdges;
}

inline const VGLight::EdgeMap& VGLight::getToEdgeMap() const
{
  return _toEdges;
}

inline size_t VGLight::getNumEdges() const
{
  return _numEdges;
}

inline const VGLight::PathMap& VGLight::getPathMap() const
{
  return _paths;
}

inline const vg::Node* VGLight::getNode(int64_t id) const
{
  vg::Node node;
  node.set_id(id);
  NodeSet::const_iterator i = _nodes.find(&node);
  return i != _nodes.end() ? *i : NULL;
}

inline const vg::Edge* VGLight::getEdge(int64_t from_id, int64_t to_id,
                                        bool from_start, bool to_end) const 
{
  std::pair<EdgeMap::const_iterator, EdgeMap::const_iterator> ret =
     _fromEdges.equal_range(from_id);
  for (EdgeMap::const_iterator i = ret.first; i != ret.second; ++i)
  {
    if (i->second->to() == to_id &&
        i->second->from_start() == from_start &&
        i->second->to_end() == to_end)
    {
      return i->second;
    }
  }
  return NULL;
}

inline const VGLight::MappingList& VGLight::getPath(const std::string& name)
  const
{
  return _paths.find(name)->second;
}

inline void VGLight::removePath(const std::string& name)
{
  assert(_paths.find(name) != _paths.end());
  _paths.erase(_paths.find(name));
}


#endif

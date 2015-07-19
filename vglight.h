/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _VGLIGHT_H
#define _VGLIGHT_H

#include <string>
#include <map>
#include <list>
#include <iostream>
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
   struct EdgePtrLess {
      bool operator()(const vg::Edge* edge, const vg::Edge* edge2) const;
   };
   typedef std::list<const vg::Mapping> MappingList;
   typedef std::set<const vg::Node*, NodePtrLess> NodeSet;
   typedef std::set<const vg::Edge*, EdgePtrLess> EdgeSet;
   typedef std::map<std::string, MappingList> PathMap;

   const NodeSet& getNodeSet() const;
   const PathMap& getPathMap() const;
   const EdgeSet& getEdgeSet() const;

   const vg::Node* getNode(int64_t id) const;
   const vg::Edge* getEdge(int64_t from_id, int64_t to_id,
                           bool from_start, bool to_end);

   /** get string for a VG path */
   void getPathDNA(const std::string& name, std::string& outDNA) const;

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
   EdgeSet _edges;
   PathMap _paths;
};

inline bool VGLight::NodePtrLess::operator()(const vg::Node* node1,
                                             const vg::Node* node2) const
{
  return node1->id() < node2->id();
}

inline bool VGLight::EdgePtrLess::operator()(const vg::Edge* edge1,
                                             const vg::Edge* edge2) const
{
  if (edge1->from() < edge2->from())
  {
    return true;
  }
  else if (edge1->from() == edge2->from())
  {
    if (edge1->to() < edge2->to())
    {
      return true;
    }
    else if (edge1->to() == edge2->to())
    {
      if (edge1->from_start() < edge2->from_start())
      {
        return true;
      }
      else
      {
        return edge1->to_end() < edge2->to_end();
      }
    }
  }
  return false;
}

inline const VGLight::NodeSet& VGLight::getNodeSet() const
{
  return _nodes;
}

inline const VGLight::EdgeSet& VGLight::getEdgeSet() const
{
  return _edges;
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
                                        bool from_start, bool to_end)
{
  vg::Edge edge;
  edge.set_from(from_id);
  edge.set_to(to_id);
  edge.set_from_start(from_start);
  edge.set_to_end(to_end);
  EdgeSet::const_iterator i = _edges.find(&edge);
  return i != _edges.end() ? *i : NULL;
}
   


#endif

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

protected:

   /** Convert _graphs into _nodes/_edges/_paths */
   void mergeGraphs();
       

   /** Graphs read direct from protobuf.  They are unmerged and we keep
    * them around just for storage */
   std::vector<vg::Graph> _graphs;

   struct NodePtrLess {
      bool operator()(const vg::Node* node1, const vg::Node* node2) const;
   };
   struct EdgePtrLess {
      bool operator()(const vg::Edge* edge, const vg::Edge* edge2) const;
   };
   typedef std::set<const vg::Node*, NodePtrLess> NodeSet;
   typedef std::set<const vg::Edge*, EdgePtrLess> EdgeSet;
   typedef std::list<const vg::Mapping> MappingList;
   typedef std::map<std::string, MappingList> PathMap;

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


#endif

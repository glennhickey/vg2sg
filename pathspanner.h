/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _PATHSPANNER_H
#define _PATHSPANNER_H

#include <string>
#include <map>
#include <vector>
#include <set>

#include "vglight.h"

/** simple method to generate a bunch of paths that cover
 * ever edge in a VG (not already covered in path)
 */
class PathSpanner
{
public:
   PathSpanner();
   ~PathSpanner();

   /** load vg paths.     */
   void init(const VGLight* vg);

   /** can we get another path with getnextpath() ? */
   bool hasNextPath() const;

   /** get a directed path of uncovered vg edges */
   void getNextPath(VGLight::MappingList& mappings);
   
protected:

   // replace ptr comp to make determinstic and help debug
   struct EdgePtrLess {
      bool operator()(const vg::Edge* e1, const vg::Edge* e2) const;
   };
   typedef std::set<const vg::Edge*, EdgePtrLess> EdgeSet;
   
   const VGLight* _vg;
   EdgeSet _uncovered;
};


#endif

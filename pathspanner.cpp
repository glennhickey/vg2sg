/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include <vector>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "pathspanner.h"

using namespace std;
using namespace vg;

PathSpanner::PathSpanner() : _vg(0)
{
}

PathSpanner::~PathSpanner()
{
  _covered.clear();
  _uncovered.clear();
}

void PathSpanner::init(const VGLight* vg)
{
  _vg = vg;
  _covered.clear();
  _uncovered.clear();
  
  // get edges from existing paths and mark them covered
  const VGLight::PathMap& pathMap = _vg->getPathMap();
  for (VGLight::PathMap::const_iterator i = pathMap.begin();
       i != pathMap.end(); ++i)
  {
    const VGLight::MappingList& mappings = i->second;
    if (mappings.size() > 0)
    {
      VGLight::MappingList::const_iterator prev = mappings.begin();
      VGLight::MappingList::const_iterator cur = prev;;
      for (++cur; cur != mappings.end(); ++cur)
      {
        int64_t from = prev->position().node_id();
        int64_t to = cur->position().node_id();
        bool from_start = prev->is_reverse();
        bool to_end = cur->is_reverse();
        const Edge* edge = _vg->getEdge(from, to, from_start, to_end);
        if (edge == NULL)
        {
          stringstream ss;
          ss << "Can't find edge (" << from << "," << to <<") from_start="
             << from_start << ", to_end=" << to_end << " implied by path "
             << i->first << ".  This means the path is invalid or, likely "
             << "I've made a wrong assumption abot the reversal flags";
          throw runtime_error(ss.str());
        }
        _covered.insert(edge);
      }
    }
  }

  // mark other edges as uncovered
  const VGLight::NodeSet& nodeSet = _vg->getNodeSet();
  vector<const Edge*> edges;
  for (VGLight::NodeSet::const_iterator i = nodeSet.begin();
       i != nodeSet.end(); ++i)
  {
    _vg->getOutEdges(*i, edges);
    for (vector<const Edge*>::iterator j = edges.begin(); j != edges.end(); ++j)
    {
      if (_covered.find(*j) != _covered.end())
      {
        _uncovered.insert(*j);
      }
    }
  }
}


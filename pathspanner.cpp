/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include <vector>
#include <deque>
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
  _uncovered.clear();
}

void PathSpanner::init(const VGLight* vg)
{
  _vg = vg;
  _uncovered.clear();
  EdgeSet covered;
  
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
        bool from_start = prev->position().is_reverse();
        bool to_end = cur->position().is_reverse();
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
        covered.insert(edge);
        prev = cur;
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
      if (covered.find(*j) == covered.end())
      {
        _uncovered.insert(*j);
      }
    }
  }
  
}

bool PathSpanner::hasNextPath() const
{
  return _uncovered.size() > 0;
}

void PathSpanner::getNextPath(VGLight::MappingList& mappings)
{
  mappings.clear();
  assert(hasNextPath() == true);
  const Edge* edge = *_uncovered.begin();
  _uncovered.erase(_uncovered.begin());
  deque<const Edge*> pathEdges;
  vector<const Edge*> nextEdges;
  pathEdges.push_back(edge);

  // extend right
  for (int prevSize = 0; prevSize < nextEdges.size(); ++prevSize)
  {
    _vg->getOutEdges(_vg->getNode(pathEdges.back()->to()), nextEdges);
    // want to find an uncovered edge thats not on the to_end.
    for (int j = 0; j < nextEdges.size(); ++j)
    {
      if (nextEdges[j]->from_start() == pathEdges.back()->to_end())
      {
        EdgeSet::iterator setIt = _uncovered.find(nextEdges[j]);
        if (setIt != _uncovered.end())
        {
          pathEdges.push_back(edge);
          _uncovered.erase(setIt);
        }
      }
    }
  }
  // extend left
  for (int prevSize = nextEdges.size() - 1; prevSize < nextEdges.size();
       ++prevSize)
  {
    _vg->getInEdges(_vg->getNode(pathEdges.front()->from()), nextEdges);
    // want to find an uncovered edge thats not on the to_end.
    for (int j = 0; j < nextEdges.size(); ++j)
    {
      if (nextEdges[j]->to_end() == pathEdges.back()->from_start())
      {
        EdgeSet::iterator setIt = _uncovered.find(nextEdges[j]);
        if (setIt != _uncovered.end())
        {
          pathEdges.push_front(edge);
          _uncovered.erase(setIt);
        }
      }
    }
  }

  // convert path edges to mapping list based on from node
  for (int i = 0; i < pathEdges.size(); ++i)
  {
    edge = pathEdges[i];
    Mapping mapping;
    Position* position = mapping.mutable_position();
    position->set_is_reverse(edge->from_start());
    position->set_node_id(edge->from());
    int64_t nodeLen = _vg->getNode(position->node_id())->sequence().length();
    position->set_offset(!mapping.position().is_reverse() ? 0 : nodeLen - 1);
    mappings.push_back(mapping);
  }
  // pop on last to node
  if (pathEdges.size() > 0)
  {
    edge = pathEdges.back();
    Mapping mapping;
    Position* position = mapping.mutable_position();
    position->set_is_reverse(edge->to_end());
    position->set_node_id(edge->to());
    int64_t nodeLen = _vg->getNode(position->node_id())->sequence().length();
    position->set_offset(!mapping.position().is_reverse() ? 0 : nodeLen - 1);
    mappings.push_back(mapping);
  }
}

bool PathSpanner::EdgePtrLess::operator()(
  const vg::Edge* e1, const vg::Edge* e2) const
{
  if (e1 == NULL && e2 != NULL)
  {
    return true;
  }
  if (e1 != NULL && e2 == NULL)
  {
    return false;
  }

  if (e1->from() < e2->from())
  {
    return true;
  }
  else if (e1->from() == e2->from())
  {
    if (e1->to() < e2->to())
    {
      return true;
    }
    else if (e1->to() == e2->to())
    {
      if (e1->from_start() < e2->from_start())
      {
        return true;
      }
      else if (e1->from_start() == e2->from_start())
      {
        if (e1->to_end() < e2->to_end())
        {
          return true;
        }
      }
    }
  }
  
  return false;
}

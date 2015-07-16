/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include <vector>
#include <sstream>
#include <cassert>
#include "pathmapper.h"

using namespace std;
using namespace vg;

PathMapper::PathMapper() : _sg(0), _lookup(0), _vg(0)
{
}

PathMapper::~PathMapper()
{
  delete _sg;
  delete _lookup;
}

void PathMapper::init(const VGLight* vg)
{
  assert(vg != NULL);
  _vg = vg;
  
  delete _sg;
  _sg = new SideGraph();

  delete _lookup;
  _lookup = new SGLookup();
  _nodeIDMap.clear();
  // lookup structure uses string names (relic from hal2sg sequences)
  // here we are mapping node coordinates, so just use nodeId
  // (note these strings aren't really used for much)
  vector<string> nodeNames;
  const VGLight::NodeSet& nodeSet = _vg->getNodeSet();
  for (VGLight::NodeSet::const_iterator i = nodeSet.begin();
       i != nodeSet.end(); ++i)
  {
    stringstream ss;
    ss << (*i)->id();
    nodeNames.push_back(ss.str());
    // make sure we can index our nodes with some number <= numNodes
    _nodeIDMap.insert(pair<int64_t, sg_int_t>((*i)->id(), _nodeIDMap.size()));
  }
  _lookup->init(nodeNames);

  // keep all vg paths indexed by name and id
  _pathNames.clear();
  _pathIDs.clear();
}

void PathMapper::addPath(const std::string& pathName)
{
  assert(_pathIDs.find(pathName) == _pathIDs.end());
  const VGLight::PathMap& pathMap = _vg->getPathMap();
  assert(pathMap.find(pathName) != pathMap.end());
  const VGLight::MappingList& mappings = pathMap.find(pathName)->second;

  _pathNames.push_back(pathName);
  _pathIDs.insert(pair<string, sg_int_t>(pathName, _pathIDs.size()));

  _curSeq = NULL;
  sg_int_t pathID = getPathID(pathName);
  sg_int_t pathPos = 0;
  size_t mappingCount = 0;
  for (VGLight::MappingList::const_iterator i = mappings.begin();
       i != mappings.end(); ++i, ++mappingCount)
  {
    const Position& pos = i->position();
    bool reversed = i->is_reverse();
    sg_int_t segmentLength = 0;
    int numEdits = i->edit_size();
    assert(numEdits > 0);
    for (int j = 0; j < numEdits; ++i)
    {
      const Edit& edit = i->edit(j);
      assert(edit.from_length() > 0);
      if (edit.from_length() != edit.to_length())
      {
        stringstream msg;
        msg << "Nontrivial edit found: " << j <<"th Edit of " << mappingCount
            << "th Mapping of Path " << pathName
            << ": to_length != from_length";
        throw runtime_error(msg.str());
      }
      else
      {
        const string& cigar = edit.sequence();
        stringstream expected;
        expected << edit.from_length();
        expected << "M";
        if (cigar.length() > 0 && cigar != expected.str())
        {
          stringstream msg;
          msg << "Nontrivial edit found: " << j <<"th Edit of " << mappingCount
              << "th Mapping of Path " << pathName
              << ": sequence=" << cigar << " != expected=" << expected.str();
          throw runtime_error(msg.str());          
        }
      }
      segmentLength += edit.from_length();
      pathPos += edit.from_length();
    }
    addSegment(pathID, pathPos, pos, reversed, segmentLength);    
  }
}

void PathMapper::addSegment(sg_int_t pathID, sg_int_t pathPos,
                            const Position& pos, bool reversed,
                            sg_int_t segLength)
{
  // when mapping, our "from" coordinate is node-relative
  sg_int_t sgNodeID = _nodeIDMap.find(pos.node_id())->second;
  SGPosition sgPos(sgNodeID, pos.offset());

  SGSide mapResult = _lookup->mapPosition(sgPos);
  bool found = mapResult.getBase() != SideGraph::NullPos;

  if (!found)
  {
    // CASE 1) : Create new SG Sequence
    if (_curSeq == NULL)
    {
      _curSeq = new SGSequence(_sg->getNumSequences(), segLength,
                               makeSeqName(pathID, pathPos));
      assert(_curSeq->getID() == _seqStrings.size());
      _seqStrings.push_back("");
    }
    // CASE 2) : Extend existing SG Sequence
    assert(_curSeq->getID() < _seqStrings.size());
    const Node* node = _vg->getNode(pos.node_id());
    string dna = node->sequence();
    if (reversed == true)
    {
      reverseComplement(dna);
    }
    assert(pos.offset() + segLength <= dna.length());
    // add dna string to the sequence
    _seqStrings[_curSeq->getID()].append(dna.substr(pos.offset(), segLength));
    _curSeq->setLength(_curSeq->getLength() + segLength);
    assert(_curSeq->getLength() == _seqStrings[_curSeq->getID()].length());

    // update lookup
  }
  else
  {
    // CASE 3) : Segment maps to existing SG Sequence
    _curSeq = NULL;
  }

  // add joins based on lookup
}

string PathMapper::makeSeqName(sg_int_t pathID, sg_int_t pathPos)
{
  stringstream ss;
  ss << getPathName(pathID) << "_" << pathPos;
  return ss.str();
}
  
char PathMapper::reverseComplement(char c)
{
  switch (c)
  {
  case 'A' : return 'T'; 
  case 'a' : return 't'; 
  case 'C' : return 'G'; 
  case 'c' : return 'g';
  case 'G' : return 'C';
  case 'g' : return 'c';
  case 'T' : return 'A';
  case 't' : return 'a';
  default : break;
  }
  return c;
}

void PathMapper::reverseComplement(std::string& s)
{
  if (!s.empty())
  {
    size_t j = s.length() - 1;
    size_t i = 0;
    char buf;
    do
    {
      while (j > 0 && s[j] == '-')
      {
        --j;
      }
      while (i < s.length() - 1 && s[i] == '-')
      {
        ++i;
      }
      
      if (i >= j || s[i] == '-' || s[j] == '-')
      {
        if (i == j && s[i] != '-')
        {
          s[i] = reverseComplement(s[i]);
        }
        break;
      }

      buf = reverseComplement(s[i]);
      s[i] = reverseComplement(s[j]);
      s[j] = buf;

      ++i;
      --j;
    } while (true);
  }
}

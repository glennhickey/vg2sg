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

  _seqStrings.clear();
  _sgPaths.clear();
  _sgSeqToVGPathID.clear();
  
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

string PathMapper::getSideGraphDNA(sg_int_t seqID, sg_int_t offset,
                                   sg_int_t length, bool reversed) const
{
  assert(seqID >= 0 && seqID < _seqStrings.size());
  if (length == -1)
  {
    length = _seqStrings[seqID].length();
  }
  string dna = _seqStrings[seqID].substr(offset, length);
  if (reversed)
  {
    VGLight::reverseComplement(dna);
  }
  return dna;
}

string PathMapper::getSideGraphPathDNA(const string& pathName) const
{
  string outString;
  const vector<SGSegment>& path = getSideGraphPath(pathName);
  size_t pathLen = 0;
  (void)pathLen;
  for (size_t i = 0; i < path.size(); ++i)
  {
    outString += getSideGraphDNA(path[i].getSide().getBase().getSeqID(),
                                 path[i].getMinPos().getPos(),
                                 path[i].getLength(),
                                 !path[i].getSide().getForward());
    pathLen += path[i].getLength();
  }
  assert(pathLen == outString.length());
  return outString;
}

const vector<SGSegment>& PathMapper::getSideGraphPath(const string& pathName)
  const
{
  return _sgPaths[getPathID(pathName)];
}

const string& PathMapper::getVGPathName(const SGSequence* seq) const
{
  return _pathNames[_sgSeqToVGPathID[seq->getID()]];
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
    sg_int_t segmentLength = _vg->getSegmentLength(*i);
    addSegment(pathID, pathPos, pos, reversed, segmentLength);
    pathPos += segmentLength;
  }
  if (_curSeq != NULL)
  {
    _sg->addSequence(_curSeq);
  }
  _curSeq = NULL;
  addPathJoins(pathName, mappings);
}

void PathMapper::verifyPaths() const
{
  for (int i = 0; i < _pathNames.size(); ++i)
  {
    string sgDNA = getSideGraphPathDNA(_pathNames[i]);
    string vgDNA;
    _vg->getPathDNA(_pathNames[i], vgDNA);
    if (vgDNA != sgDNA)
    {
      stringstream ss;
      ss << "Verification failed for VG path " << _pathNames[i] << ": "
         << "output DNA differs from input.  Please report this bug";
      throw runtime_error(ss.str());
    }
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
      _curSeq = new SGSequence(_sg->getNumSequences(), 0,
                               makeSeqName(pathID, pathPos));
      assert(_curSeq->getID() == _seqStrings.size());
      _seqStrings.push_back("");
      _sgSeqToVGPathID.push_back(_pathNames.size()-1);
    }
    // CASE 2) : Extend existing SG Sequence
    sg_int_t curSeqLen = _curSeq->getLength();
    assert(_curSeq->getID() < _seqStrings.size());
    const Node* node = _vg->getNode(pos.node_id());
    string dna = node->sequence();
    if (reversed == true)
    {
      VGLight::reverseComplement(dna);
    }
    assert(pos.offset() + segLength <= dna.length());
    // add dna string to the sequence
    _seqStrings[_curSeq->getID()].append(dna);
    _curSeq->setLength(_curSeq->getLength() + segLength);
    assert(_curSeq->getLength() == _seqStrings[_curSeq->getID()].length());

    // update lookup
    SGPosition toPos(_curSeq->getID(), curSeqLen);
    // todo: reverse mapping shift
    // todo: offset shift
    
    _lookup->addInterval(sgPos, toPos, segLength, reversed);
  }
  else
  {
    // CASE 3) : Segment maps to existing SG Sequence
    if (_curSeq != NULL)
    {
      _sg->addSequence(_curSeq);
    }
    _curSeq = NULL;
  }
}

// doesn't really need a second pass but whatever
void PathMapper::addPathJoins(const string& name,
                              const VGLight::MappingList& mappings)
{
  _sgPaths.resize(_sgPaths.size() + 1);
  vector<SGSegment>& sgPath = _sgPaths.back();
  int mappingCount = 0;
  for (VGLight::MappingList::const_iterator i = mappings.begin();
       i != mappings.end(); ++i, ++mappingCount)
  {
    const Position& pos = i->position();
    bool reversed = i->is_reverse();
    const Node* node = _vg->getNode(pos.node_id());
    int64_t segmentLength = _vg->getSegmentLength(*i);
    int64_t offset = pos.offset();
    assert(offset >= 0);
    if (offset > 0 && mappingCount != 0)
    {
      stringstream ss;
      ss << "Offset > 0 found on " << mappingCount << "th mappping of "
         << "path " << name << ".  Offsets only permitted on 0th mapping";
      throw runtime_error(ss.str());
    }
    sg_int_t nodeID = _nodeIDMap.find(node->id())->second;
    SGPosition start(nodeID, offset);
    SGPosition end(nodeID, offset + segmentLength - 1);
    if (reversed)
    {
      swap(start, end);
    }
    vector<SGSegment> nextPath;
    _lookup->getPath(start, end, nextPath);
    if (reversed && start == end)
    {
      // weird. doesn't seem to be issue in hal2sg. should figure out why,
      // but just hack here for now (case where 1-base path not reversed)
      assert(nextPath.size() == 1);
      nextPath[0] = SGSegment(SGSide(nextPath[0].getSide().getBase(),
                                     !nextPath[0].getSide().getForward()),
                              nextPath[0].getLength());
    }
    mergePaths(sgPath, nextPath);
  }

  for (size_t i = 1; i < sgPath.size(); ++i)
  {
    SGSide srcSide = sgPath[i-1].getOutSide();
    SGSide tgtSide = sgPath[i].getInSide();
    SGJoin* join = new SGJoin(srcSide, tgtSide);
          
    if (!join->isTrivial())
    {
      _sg->addJoin(join);
    }
    else
    {
      delete join;
    }
  }
}

void PathMapper::mergePaths(vector<SGSegment>& prevPath,
                            const vector<SGSegment>& path) const
{
  assert(path.size() > 0);
  const SGSegment& nextSeg = path[0];
  int i = 0;
  if (!prevPath.empty() &&
      prevPath.back().getSide().getForward() ==
      nextSeg.getSide().getForward() &&
      prevPath.back().getOutSide().lengthTo(nextSeg.getInSide()) == 0)
  {
    prevPath.back().setLength(prevPath.back().getLength() +
                              nextSeg.getLength());
    ++i;
  }
  for (; i < path.size(); ++i)
  {
    prevPath.push_back(path[i]);
  }
}

string PathMapper::makeSeqName(sg_int_t pathID, sg_int_t pathPos)
{
  stringstream ss;
  ss << getPathName(pathID) << "_" << pathPos;
  return ss.str();
}
  

/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include <fstream>
#include <sstream>
#include "google/protobuf/stubs/common.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/io/gzip_stream.h"
#include "google/protobuf/io/coded_stream.h"

#include "vglight.h"

using namespace std;
using namespace vg;
using namespace google::protobuf::io;

VGLight::VGLight() : _numEdges(0)
{
}

VGLight::~VGLight()
{
}

/** Most code copy-pasted from VG stream constructor (vg.cpp)  and 
 * the functions it calls
 */
void VGLight::loadGraph(istream& in)
{
  ZeroCopyInputStream *raw_in = new IstreamInputStream(&in);
  GzipInputStream *gzip_in = new GzipInputStream(raw_in);
  CodedInputStream *coded_in = new CodedInputStream(gzip_in);
 
  uint64_t count;
  coded_in->ReadVarint64(&count);
  if (!count)
  {
    throw runtime_error("Empty stream");
  }

  
  _graphs.clear();
  
  do
  {
    std::string s;
    for (uint64_t i = 0; i < count; ++i)
    {
      uint32_t msgSize = 0;
      delete coded_in;
      coded_in = new CodedInputStream(gzip_in);
      // the messages are prefixed by their size
      coded_in->ReadVarint32(&msgSize);
      if ((msgSize > 0) &&
          (coded_in->ReadString(&s, msgSize)))
      {
        Graph object; 
        object.ParseFromString(s);
        _graphs.push_back(object);
      }
    }   
  } while (coded_in->ReadVarint64(&count));

  mergeGraphs();

  delete coded_in;
  delete gzip_in;
  delete raw_in;
}

void VGLight::loadGraph(const Graph& graph)
{
  _graphs.resize(1);
  _graphs[0] = graph;
  mergeGraphs();
}

void VGLight::mergeGraphs()
{
  _paths.clear();
  _nodes.clear();
  _fromEdges.clear();
  _toEdges.clear();
  _numEdges = 0;

  for (size_t i = 0; i < _graphs.size(); ++i)
  {
    for (size_t j = 0; j < _graphs[i].node_size(); ++j)
    {
      _nodes.insert(_graphs[i].mutable_node(j));
    }
    for (size_t j = 0; j < _graphs[i].edge_size(); ++j)
    {
      const Edge* edge = _graphs[i].mutable_edge(j);
      _fromEdges.insert(pair<int64_t, const Edge*>(edge->from(), edge));
      _toEdges.insert(pair<int64_t, const Edge*>(edge->to(), edge));
      ++_numEdges;
    }
    for (size_t j = 0; j < _graphs[i].path_size(); ++j)
    {
      Path* path = _graphs[i].mutable_path(j);
      pair<PathMap::iterator, bool> ret = _paths.insert(
        pair<string, MappingList>(path->name(), MappingList()));
      for (size_t k = 0; k < path->mapping_size(); ++k)
      {
        ret.first->second.push_back(path->mapping(k));
      }
      // sort by rank and remove dupes
      set<Mapping, MappingRankLess> mappingSet(ret.first->second.begin(),
                                               ret.first->second.end());
      ret.first->second = MappingList(mappingSet.begin(), mappingSet.end());
    }
  }
}

void VGLight::getInEdges(const Node* node, 
                         vector<const Edge*>& ins) const
{
  ins.clear();
  pair<EdgeMap::const_iterator, EdgeMap::const_iterator> ret;

  ret = _toEdges.equal_range(node->id());
  for (EdgeMap::const_iterator i = ret.first; i != ret.second; ++i)
  {
    ins.push_back(i->second);
  }
}

void VGLight::getOutEdges(const Node* node, 
                         vector<const Edge*>& outs) const
{
  outs.clear();
  pair<EdgeMap::const_iterator, EdgeMap::const_iterator> ret;

  ret = _fromEdges.equal_range(node->id());
  for (EdgeMap::const_iterator i = ret.first; i != ret.second; ++i)
  {
    outs.push_back(i->second);
  }
}

void VGLight::getPathDNA(const string& pathName, string& outDNA) const
{
  assert(_paths.find(pathName) != _paths.end());
  const MappingList& mappingList = _paths.find(pathName)->second;
  getPathDNA(mappingList, outDNA);
}

void VGLight::getPathDNA(const MappingList& mappingList, string& outDNA) const
{
  outDNA.erase();
  for (MappingList::const_iterator i = mappingList.begin();
       i != mappingList.end(); ++i)
  {
    const Position& pos = i->position();
    bool reversed = i->is_reverse();
    const Node* node = getNode(pos.node_id());
    int64_t segmentLength = getSegmentLength(*i);
    int64_t offset = pos.offset();
    if (reversed)
    {
      offset -= segmentLength - 1;
    }
    string dna = node->sequence().substr(offset, segmentLength);
    if (reversed)
    {
      reverseComplement(dna);
    }
    outDNA += dna;
  }
}

int64_t VGLight::getSegmentLength(const Mapping& mapping) const
{
  int64_t segmentLength = 0;
  // no edits: take distance from offset to end of node
  if (mapping.edit_size() == 0)
  {
    const Position& pos = mapping.position();
    const Node* node = getNode(pos.node_id());
    if (!mapping.is_reverse())
    {
      segmentLength = node->sequence().length() - pos.offset();
    }
    else
    {
      segmentLength = pos.offset() + 1;
    }
  }
  // has edits: take total length of edits (???)
  else
  {
    int numEdits = mapping.edit_size();
    for (int i = 0; i < numEdits; ++i)
    {
      const Edit& edit = mapping.edit(i);
      if (edit.from_length() != edit.to_length())
      {
        stringstream msg;
        msg << "Nontrivial edit found: to_length != from_length";
        throw runtime_error(msg.str());
      }
      else if (edit.sequence().length() > 0)
      {
        stringstream msg;
        msg << "Nontrivial edit found: sequence=" << edit.sequence()
            << " (only empty sequence supported)";
        throw runtime_error(msg.str());          
      }
      segmentLength += mapping.edit(i).from_length();
    }
  }
  return segmentLength;
}

char VGLight::reverseComplement(char c)
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

void VGLight::reverseComplement(std::string& s)
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

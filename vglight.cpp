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

VGLight::VGLight()
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
  _edges.clear();

  for (size_t i = 0; i < _graphs.size(); ++i)
  {
    for (size_t j = 0; j < _graphs[i].node_size(); ++j)
    {
      _nodes.insert(_graphs[i].mutable_node(j));
    }
    for (size_t j = 0; j < _graphs[i].edge_size(); ++j)
    {
      _edges.insert(_graphs[i].mutable_edge(j));
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
    }
  }
}

void VGLight::getPathDNA(const string& pathName, string& outDNA) const
{
  outDNA.erase();
  assert(_paths.find(pathName) != _paths.end());
  const MappingList& mappingList = _paths.find(pathName)->second;
  for (MappingList::const_iterator i = mappingList.begin();
       i != mappingList.end(); ++i)
  {
    const Position& pos = i->position();
    bool reversed = i->is_reverse();
    const Node* node = getNode(pos.node_id());
    int64_t segmentLength = getSegmentLength(*i);
    int64_t offset = pos.offset();
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
    segmentLength = node->sequence().length() - pos.offset();
  }
  // has edits: take total length of edits (???)
  else
  {
    int numEdits = mapping.edit_size();
    for (int i = 0; i < numEdits; ++i)
    {
      const Edit& edit = mapping.edit(i);
      assert(edit.from_length() > 0);
      if (edit.from_length() != edit.to_length())
      {
        stringstream msg;
        msg << "Nontrivial edit found: : to_length != from_length";
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
          msg << "Nontrivial edit found: sequence=" << cigar
              << " != expected=" << expected.str();
          throw runtime_error(msg.str());          
        }
        segmentLength += mapping.edit(i).from_length();
      }
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

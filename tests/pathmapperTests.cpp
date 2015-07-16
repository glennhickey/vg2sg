/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <sstream>
#include "unitTests.h"
#include "pathmapper.h"

using namespace std;
using namespace vg;

///////////////////////////////////////////////////////////
//  Helper stuff
///////////////////////////////////////////////////////////

static string randDNA(size_t len)
{
  string out;
  for (size_t i = 0; i < len; ++i)
  {
    long x = rand() % 4;
    char c;
    switch(x)
    {
    case 0 : c = 'a'; break;
    case 1 : c = 'c'; break;
    case 2 : c = 'g'; break;
    default:
      c = 't';
    }
    out += c;
  }
  return out;
}

static Node* makeNode(Graph& graph, int64_t id, const string& dna)
{
  Node* node = graph.add_node();
  node->set_id(id);
  node->set_sequence(dna);
  return node;
}

static Edge* makeEdge(Graph& graph, const int64_t from, const int64_t to,
                      bool from_start = false, bool to_start = false)
{
  Edge* edge = graph.add_edge();
  edge->set_from(from);
  edge->set_to(to);
  edge->set_from_start(from_start);
  edge->set_to_end(to_start);
  return edge;
}

// simple path where a node is cmpletely spanned by 1 edit.
static Path* makePath(Graph& graph, const string& name,
                      const vector<const Node*>& nodes,
                      const vector<bool>& flips)
{
  assert(flips.size() == nodes.size());
  Path* path = graph.add_path();
  path->set_name(name);
  bool forward = true;
  for (int i = 0; i < nodes.size(); ++i)
  {
    Mapping* mapping = path->add_mapping();
    Position* position = mapping->mutable_position();
    position->set_node_id(nodes[i]->id());
    position->set_offset(0);
    if (flips[i] == true)
    {
      forward = !forward;
    }
    mapping->set_is_reverse(!forward);
    if (i < nodes.size() - 1)
    {
      bool from_start = !forward;
      bool to_start = forward && flips[i+1];
      makeEdge(graph, nodes[i]->id(), nodes[i+1]->id(), from_start, to_start);
    }
    Edit* edit = mapping->add_edit();
    edit->set_from_length(nodes[i]->sequence().length());
    edit->set_to_length(nodes[i]->sequence().length());
    stringstream ss;
    ss << edit->to_length() << "M";
    edit->set_sequence(ss.str());
  }
  return path;
}

///////////////////////////////////////////////////////////
//  Simple Test
//    - one forward path through a few nodes
///////////////////////////////////////////////////////////
void simpleTest(CuTest *testCase)
{
  CuAssertTrue(testCase, true);
  string dna = "ACAAACACACAGGGTACACGTACAGACCGACTTAGCAGAGAT";
  
  Graph graph;
  vector<const Node*> path;
  vector<bool> flips(3, false);
  path.push_back(makeNode(graph, 0, dna.substr(0, 6)));
  path.push_back(makeNode(graph, 2, dna.substr(6, 20)));
  path.push_back(makeNode(graph, 1, dna.substr(26, 16)));

  makePath(graph, "path", path, flips);
  VGLight vg;
  vg.loadGraph(graph);
  PathMapper pm;
  pm.init(&vg);
  pm.addPath("path");
  const SideGraph* sg = pm.getSideGraph();
  CuAssertTrue(testCase, sg->getNumSequences() == 1);
  const SGSequence* seq = sg->getSequence(0);
  CuAssertTrue(testCase, seq->getLength() == dna.length());
  string seqDNA = pm.getSideGraphDNA(seq->getID());
  CuAssertTrue(testCase, seqDNA == dna);
}

CuSuite* pathMapperTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, simpleTest);
  return suite;
}

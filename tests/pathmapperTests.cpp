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
//    - also one back path 
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
  const vector<SGSegment>& sgPath = pm.getSideGraphPath("path");
  CuAssertTrue(testCase, sgPath.size() == 1);
  CuAssertTrue(testCase, sgPath[0].getMinPos() == SGPosition(0, 0));
  CuAssertTrue(testCase,
               sgPath[0].getMaxPos() == SGPosition(0, dna.length() -1));
  CuAssertTrue(testCase, sgPath[0].getSide().getForward() == true);
  std::string vgPath;
  vg.getPathDNA("path", vgPath);
  CuAssertTrue(testCase, vgPath == dna);
  CuAssertTrue(testCase, vgPath == pm.getSideGraphPathDNA("path"));

  vector<const Node*> htap;
  vector<bool> spilf(3, false);
  spilf[0] = true;
  htap.push_back(path[2]);
  htap.push_back(path[1]);
  htap.push_back(path[0]);
  makePath(graph, "htap", htap, spilf);
  vg.loadGraph(graph);
  pm.init(&vg);
  pm.addPath("htap");
  sg = pm.getSideGraph();
  CuAssertTrue(testCase, sg->getNumSequences() == 1);
  seq = sg->getSequence(0);
  CuAssertTrue(testCase, seq->getLength() == dna.length());
  seqDNA = pm.getSideGraphDNA(seq->getID());
  VGLight::reverseComplement(seqDNA);
  CuAssertTrue(testCase, seqDNA == dna);
  const vector<SGSegment>& sgHtap = pm.getSideGraphPath("htap");
  CuAssertTrue(testCase, sgHtap.size() == 1);
  CuAssertTrue(testCase, sgHtap[0].getMinPos() == SGPosition(0, 0));
  CuAssertTrue(testCase,
               sgHtap[0].getMaxPos() == SGPosition(0, dna.length() -1));
  CuAssertTrue(testCase, sgHtap[0].getSide().getForward() == true);
  std::string vgHtap;
  vg.getPathDNA("htap", vgHtap);
  seqDNA = pm.getSideGraphDNA(seq->getID());
  CuAssertTrue(testCase, vgHtap == seqDNA);
  CuAssertTrue(testCase, vgHtap == pm.getSideGraphPathDNA("htap"));

  try {
    pm.verifyPaths();
  }
  catch(...)
  {
    CuAssertTrue(testCase, false);
  }    
}

///////////////////////////////////////////////////////////
//  Simple Inversion Test
//    - one forward path through a few nodes with inversion
//      in middle
///////////////////////////////////////////////////////////
void inversionTest(CuTest *testCase)
{
  CuAssertTrue(testCase, true);
  string dna = "ACAAACACACAGGGTACACGTACAGACCGACTTAGCAGAGAT";
  
  Graph graph;
  vector<const Node*> path;
  vector<bool> flips(3, false);
  flips[1] = true; // flips strand coming into 2nd node
  flips[2] = true; // flips strand coming into 3rd node
  path.push_back(makeNode(graph, 0, dna.substr(0, 6)));
  string node2dna = dna.substr(6, 20);
  VGLight::reverseComplement(node2dna);
  path.push_back(makeNode(graph, 2, node2dna));

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
  string pathDNA = pm.getSideGraphPathDNA("path");
  CuAssertTrue(testCase, pathDNA == dna);
  vector<SGSegment> sgPath = pm.getSideGraphPath("path");
  CuAssertTrue(testCase, sgPath.size() == 1);

  try {
    pm.verifyPaths();
  }
  catch(...)
  {
    CuAssertTrue(testCase, false);
  }    
}

///////////////////////////////////////////////////////////
//  Simple Overlap Test
//    - paths that share some nodes
///////////////////////////////////////////////////////////
void overlapTest(CuTest *testCase)
{
  CuAssertTrue(testCase, true);
  Graph graph;

  // "path1" simple forward path;
  vector<const Node*> path1;
  vector<bool> flips1(4, false);
  int nc = 0;
  path1.push_back(makeNode(graph, nc++, randDNA(5)));
  path1.push_back(makeNode(graph, nc++, randDNA(7)));
  path1.push_back(makeNode(graph, nc++, "A"));
  path1.push_back(makeNode(graph, nc++, randDNA(8)));
  makePath(graph, "path1", path1, flips1);
  
  string trueSeq0 = path1[0]->sequence() + path1[1]->sequence() +
     path1[2]->sequence() + path1[3]->sequence();

  // "path2" single snp at 3rd node
  vector<const Node*> path2;
  vector<bool> flips2(4, false);
  path2.push_back(path1[0]);
  path2.push_back(path1[1]);
  path2.push_back(makeNode(graph, nc++, "G"));
  path2.push_back(path1[3]);
  makePath(graph, "path2", path2, flips2);

  string trueSeq1 = path2[2]->sequence();
  SGJoin trueJoin1(SGSide(SGPosition(0, 5 + 7 - 1), false),
                   SGSide(SGPosition(1, 0), true));
  SGJoin trueJoin2(SGSide(SGPosition(1, 0), false),
                   SGSide(SGPosition(0, 5 + 7 + 1), true));

  // "path3" start backwards at last two nodes of path2,
  // and continue backwards over two new nodes
  vector<const Node*> path3;
  vector<bool> flips3(4, false);
  flips3[0] = true;
  path3.push_back(path2[3]);
  path3.push_back(path2[2]);
  path3.push_back(makeNode(graph, nc++, randDNA(10)));
  path3.push_back(makeNode(graph, nc++, randDNA(2)));
  makePath(graph, "path3", path3, flips3);

  string trueSeq2 = path3[3]->sequence() + path3[2]->sequence();
  VGLight::reverseComplement(trueSeq2);
  SGJoin trueJoin3(SGSide(SGPosition(1, 0), true),
                   SGSide(SGPosition(2, 0), true));
  
  // "path4" path1 but delete second node and duplicated 3rd
  vector<const Node*> path4;
  vector<bool> flips4(4, false);
  path4.push_back(path1[0]);
  path4.push_back(path1[2]);
  path4.push_back(path1[2]);
  path4.push_back(path1[3]);
  makePath(graph, "path4", path4, flips4);

  SGJoin trueJoin4(SGSide(SGPosition(0, 4), false),
                   SGSide(SGPosition(0, 12), true));
  SGJoin trueJoin5(SGSide(SGPosition(0, 12), false),
                   SGSide(SGPosition(0, 12), true));
  
  VGLight vg;
  vg.loadGraph(graph);
  PathMapper pm;
  pm.init(&vg);
  pm.addPath("path1");
  pm.addPath("path2");
  pm.addPath("path3");
  pm.addPath("path4");
  const SideGraph* sg = pm.getSideGraph();
  CuAssertTrue(testCase, sg->getNumSequences() == 3);
  const SGSequence* seq = sg->getSequence(0);
  string seqDNA = pm.getSideGraphDNA(seq->getID());
  CuAssertTrue(testCase, seq->getLength() == seqDNA.length());
  CuAssertTrue(testCase, seqDNA == trueSeq0);

  seq = sg->getSequence(1);
  seqDNA = pm.getSideGraphDNA(seq->getID());
  CuAssertTrue(testCase, seq->getLength() == seqDNA.length());
  CuAssertTrue(testCase, seqDNA == trueSeq1);

  seq = sg->getSequence(2);
  seqDNA = pm.getSideGraphDNA(seq->getID());
  CuAssertTrue(testCase, seq->getLength() == seqDNA.length());
  CuAssertTrue(testCase, seqDNA == trueSeq2);

  string sgPath1 = pm.getSideGraphPathDNA("path1");
  string vgPath1;
  vg.getPathDNA("path1", vgPath1);
  CuAssertTrue(testCase, sgPath1 == vgPath1);

  string sgPath2 = pm.getSideGraphPathDNA("path2");
  string vgPath2;
  vg.getPathDNA("path2", vgPath2);
  CuAssertTrue(testCase, sgPath2 == vgPath2);

  string sgPath3 = pm.getSideGraphPathDNA("path3");
  string vgPath3;
  vg.getPathDNA("path3", vgPath3);

  string temp = path2[3]->sequence();
  VGLight::reverseComplement(temp);

  CuAssertTrue(testCase, sgPath3 == vgPath3);

  string sgPath4 = pm.getSideGraphPathDNA("path4");
  string vgPath4;
  vg.getPathDNA("path4", vgPath4);
  CuAssertTrue(testCase, sgPath4 == vgPath4);

  CuAssertTrue(testCase, sg->getJoin(&trueJoin1) != NULL);
  CuAssertTrue(testCase, sg->getJoin(&trueJoin2) != NULL);
  CuAssertTrue(testCase, sg->getJoin(&trueJoin3) != NULL);
  CuAssertTrue(testCase, sg->getJoin(&trueJoin4) != NULL);
  CuAssertTrue(testCase, sg->getJoin(&trueJoin5) != NULL);
  try {
    pm.verifyPaths();
  }
  catch(...)
  {
    CuAssertTrue(testCase, false);
  }    
}


CuSuite* pathMapperTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, simpleTest);
  SUITE_ADD_TEST(suite, inversionTest);
  SUITE_ADD_TEST(suite, overlapTest);
  return suite;
}

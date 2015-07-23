/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <string>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <getopt.h>

#include "pathmapper.h"
#include "vgsgsql.h"

using namespace std;

void help(char** argv)
{
  cerr << "usage: " << argv[0] << " <graph.vg> <out.fa> <out.sql> [options]\n"
       << "args:\n"
       << "    graph.vg:  Input VG graph to convert\n"
       << "    out.fa  :  Output Side Graph sequences file in FASTA format\n"
       << "    out.sql :  Output Side Graph SQL inserts file\n"
       << "options:\n"
       << "    -h, --help         \n"
       << "    -p, --primaryPath  Primary path name\n"
       << "    -s, --span         Create a dummy path that spans all edges\n"
       << "                       to make sure entire graph gets converted.\n"
       << endl;
}

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    help(argv);
    return 1;
  }

  string primaryPathName;
  bool span = false;
  optind = 4;
  while (true)
  {
    static struct option long_options[] =
       {
         {"help", no_argument, 0, 'h'},
         {"primaryPath", required_argument, 0, 'p'},
         {"span", no_argument, 0, 's'}
       };
    int option_index = 0;
    int c = getopt_long(argc, argv, "hp:s", long_options, &option_index);

    if (c == -1)
    {
      break;
    }
    
    switch(c)
    {
    case 'h':
    case '?':
      help(argv);
      exit(1);
    case 'p':
      primaryPathName = optarg;
      break;
    case 's':
      span = true;
      break;
    default:
      abort();
    }
  }
  
  string vgPath = argv[1];
  string outFaPath = argv[2];
  string outSQLPath = argv[3];

  ifstream vgStream(vgPath.c_str());
  if (!vgStream)
  {
    throw runtime_error(string("Error opening " + vgPath));
  }
  
  VGLight vglight;
  cout << "Reading input graph from disk" << endl;
  vglight.loadGraph(vgStream);
  cout << "Graph has " << vglight.getNodeSet().size() << " nodes, "
       << vglight.getNumEdges() << " edges and "
       << vglight.getPathMap().size() << " paths";
  size_t numMappings = 0;
  for (VGLight::PathMap::const_iterator i = vglight.getPathMap().begin();
       i != vglight.getPathMap().end(); ++i)
  {
    numMappings += i->second.size();
  }
  cout << " with a total of " << numMappings << " mappings." << endl;

  const VGLight::PathMap& paths = vglight.getPathMap();
  
  if (primaryPathName.length() > 0)
  {
    if (paths.find(primaryPathName) == paths.end())
    {
      throw runtime_error(string("Primary path ") + primaryPathName +
                          string(" not found in vg"));
    }
  }
  else
  {
    if (paths.empty())
    {
      throw runtime_error("No paths in input vg");
    }
    primaryPathName = paths.begin()->first;
  }

  PathMapper pm;
  pm.init(&vglight);
  cout << "Adding (primary) VG path: " << primaryPathName << endl;
  pm.addPath(primaryPathName, vglight.getPath(primaryPathName));
  for (VGLight::PathMap::const_iterator i = paths.begin(); i != paths.end();
       ++i)
  {
    if (i->first != primaryPathName)
    {
      cout << "Adding VG path: " << i->first << endl;
      pm.addPath(i->first, i->second);
    }
  }
  if (span == true)
  {
    cout << "Adding dummy path that spans all VG edges" << endl;
    pm.addSpanningPaths();
  }
  pm.verifyPaths();


  VGSGSQL sqlWriter;
  sqlWriter.exportGraph(&pm, outSQLPath, outFaPath, vgPath);

  cout << "side graph = " << *pm.getSideGraph() << endl;
  
  
}

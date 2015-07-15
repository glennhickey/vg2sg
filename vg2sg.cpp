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

#include "vglight.h"

using namespace std;

void help(char** argv)
{
  cerr << "usage: " << argv[0] << " <graph.vg> <out.fa> <out.sql> [options]\n"
       << "   graph.vg:  Input VG graph to convert\n"
       << "   out.fa  :  Output Side Graph sequences file in FASTA format\n"
       << "   out.sql :  Output Side Graph SQL inserts file\n"
       << endl;
}

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    help(argv);
    return 1;
  }

  optind = 4;
  while (true)
  {
    static struct option long_options[] =
       {
         {"help", no_argument, 0, 'h'}
       };
    int option_index = 0;
    int c = getopt_long(argc, argv, "h", long_options, &option_index);

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
    default:
      abort();
    }
  }
  
  string vgPath = argv[1];
  string outFaPath = argv[2];
  string outSQLPath = argv[3];

  ifstream vgStream(vgPath);
  if (!vgStream)
  {
    throw runtime_error(string("Error opening " + vgPath));
  }
  
  VGLight vglight;
  vglight.loadGraph(vgStream);

  // debugf
  cout << "numNodes " << vglight.getNodeSet().size() << endl
       << "numEdges " << vglight.getEdgeSet().size() << endl
       << "numPaths " << vglight.getPathMap().size() << endl;
  
  
}

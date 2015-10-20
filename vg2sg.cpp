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
#include <sstream>
#include <stdexcept>

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
       << "    -s, --span         Create a path set that spans all edges\n"
       << "                       to make sure entire graph gets converted.\n"
       << endl;
}

/** Check if a path has edits.  Spit warning to stderr and return false
 *  if it does */
static bool checkPath(const VGLight& vglight,
                      const std::string& name,
                      const VGLight::MappingList& mappings,
                      bool span);

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    help(argv);
    return 1;
  }

  string primaryPathName;
  bool span = false;
  optind = 1;
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
  
  string vgPath = argv[optind++];
  string outFaPath = argv[optind++];
  string outSQLPath = argv[optind];

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
  else if (!paths.empty())
  {
    primaryPathName = paths.begin()->first;
  }
  
  PathMapper pm;
  pm.init(&vglight);

  if (!primaryPathName.empty())
  {
    cout << "Adding (primary) VG path: " << primaryPathName << endl;
    if (checkPath(vglight, primaryPathName,
                  vglight.getPath(primaryPathName), span))
    {
      pm.addPath(primaryPathName, vglight.getPath(primaryPathName));
    }
  }
  else if (!span)
  {
    assert(paths.empty());
    throw runtime_error("No paths to convert using default logic.  Use "
                        "--span option to convert entire graph with inferred"
                        " spanning paths.");
  }
  
  for (VGLight::PathMap::const_iterator i = paths.begin(); i != paths.end();
       ++i)
  {
    if (i->first != primaryPathName)
    {
      cout << "Adding VG path: " << i->first << endl;
      if (checkPath(vglight, i->first, i->second, span))
      {
        pm.addPath(i->first, i->second);
      }
    }
  }
  if (span == true)
  {
    cout << "Adding set of paths that span all remaining VG edges" << endl;
    pm.addSpanningPaths();
  }
  pm.verifyPaths();


  VGSGSQL sqlWriter;
  sqlWriter.exportGraph(&pm, outSQLPath, outFaPath, vgPath);

  //cout << "side graph = " << *pm.getSideGraph() << endl;
  
  
}

/** Check if a path has edits.  Spit warning to stderr and return false
 *  if it does */
bool checkPath(const VGLight& vglight,
               const string& name,
               const VGLight::MappingList& mappings,
               bool span)
{
  try
  {
    // this is rather wasteful, but does the job.
    // can make smarter check function if performance comes up
    // as issue...
    string buffer;
    vglight.getPathDNA(mappings, buffer);
  }
  catch(runtime_error e)
  {
    if (span == true)
    {
      cerr << "Warning: Skipping path " << name << " because: "
           << e.what() << endl;
    }
    else
    {
      stringstream msg;
      msg << e.what() << " -- NOTE: This error can be turned into a "
          << "warning by using the --span option.";
      throw runtime_error(msg.str());
    }
    return false;
  }
  return true;
}

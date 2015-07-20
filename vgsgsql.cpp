/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "md5.h"
#include "vgsgsql.h"

using namespace std;
using namespace vg;

VGSGSQL::VGSGSQL() : SGSQL(), _pm(0)
{
}

VGSGSQL::~VGSGSQL()
{
}

void VGSGSQL::exportGraph(const PathMapper* pm,
                          const string& sqlInsertPath,
                          const string& fastaPath, const string& vgPath)
{
  _pm = pm;
  _halPath = vgPath;

  writeDb(pm->getSideGraph(), sqlInsertPath, fastaPath);
}

void VGSGSQL::getSequenceString(const SGSequence* seq,
                                 string& outString) const
{
  outString = _pm->getSideGraphDNA(seq->getID());
}

string VGSGSQL::getOriginName(const SGSequence* seq) const
{
  return _pm->getVGPathName(seq);
}

string VGSGSQL::getPrimaryOriginName() const
{
  return _pm->getPathName(0);
}

string VGSGSQL::getDescription() const
{
  return string("vg2sg ") + _halPath;
}

/*
CREATE TABLE VariantSet (ID INTEGER PRIMARY KEY,
	referenceSetID INTEGER NOT NULL REFERENCES ReferenceSet(ID),
	name TEXT);

CREATE TABLE Allele (ID INTEGER PRIMARY KEY, 
	variantSetID INTEGER REFERENCES VariantSet(ID), 
	name TEXT); -- Naming the allele is optional
--
CREATE TABLE AllelePathItem (alleleID INTEGER REFERENCES allele(ID), 
	pathItemIndex INTEGER NOT NULL, -- zero-based index of this pathItem within the entire path
	sequenceID INTEGER NOT NULL REFERENCES Sequence(ID), 
	start INTEGER NOT NULL,
	length INTEGER NOT NULL, 
	strandIsForward BOOLEAN NOT NULL,
	PRIMARY KEY(alleleID, pathItemIndex));
*/
void VGSGSQL::writePathInserts()
{
  // create single Variant set
  _outStream << "INSERT INTO VariantSet VALUES ("
             << 0 << ", "
             << 0 << ", "
             << "'" << "vg2sg" << "');\n";

  _outStream << endl;

  // For every VG path, we create one
  // Allele, and a list of Allele Path Items.
  for (size_t i = 0; i < _pm->getNumPaths(); ++i)
  {
    _outStream << "INSERT INTO Allele VALUES ("
               << i << ", "
               << 0 << ", "
               << "'" << _pm->getPathName(i) << "'" 
               << ");\n";
  }
  _outStream << endl;

  // create a path (AellePathItem) for every sequence
  for (size_t i = 0; i < _pm->getNumPaths(); ++i)
  {
    _outStream << "-- PATH for VG input sequence "
               << _pm->getPathName(i) << "\n";
    const vector<SGSegment>& path = _pm->getSideGraphPath(_pm->getPathName(i));
    for (size_t j = 0; j < path.size(); ++j)
    {
      _outStream << "INSERT INTO AllelePathItem VALUES ("
                 << i << ", "
                 << j << ", "
                 << path[j].getSide().getBase().getSeqID() << ", "
                 << path[j].getSide().getBase().getPos() << ", "
                 << path[j].getLength() << ", "
                 << (path[j].getSide().getForward() ? "\'TRUE\'" : "\'FALSE\'")
                 << ");\n";
      
    }
    _outStream <<endl;
  }

}


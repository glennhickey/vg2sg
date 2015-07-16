/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include <sstream>
#include "unitTests.h"
#include "sidegraph.h"

using namespace std;

void dummyTest(CuTest *testCase)
{
  CuAssertTrue(testCase, true);
}

CuSuite* pathMapperTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, dummyTest);
  return suite;
}

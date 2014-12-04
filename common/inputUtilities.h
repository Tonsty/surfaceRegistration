#ifndef __SurfaceRegistration__inputUtilities__
#define __SurfaceRegistration__inputUtilities__

#include <stdio.h>
#include <iostream>
#include <string>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>

#include <vtkExtractEdges.h>
#include <vtkLine.h>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> BGLUndirectedGraph;
typedef std::pair < int, int > Edge;

BGLUndirectedGraph vtkGetInput (std::string path);
BGLUndirectedGraph vtpGetInput (std::string path);

#endif

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

struct vertex_info {
    double x;
    double y;
    double z;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, vertex_info> BGLUndirectedGraph;

BGLUndirectedGraph vtkGetInput (std::string path);
BGLUndirectedGraph vtpGetInput (std::string path);

#endif

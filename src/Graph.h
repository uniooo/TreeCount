#pragma once
#define _CRT_SECURE_NO_DEPRECATE

#include <assert.h>

#include <map>
#include <utility>
#include <vector>

using std::vector;
using std::map;

class Graph {
 public:
  int n_, m_;  // number of nodes & edges
  vector<map<int, int> >
      E_;               // record edges and edge weight: u--->(v, dis)
  vector<int> D_;  // degree

 public:
  Graph();
  Graph(const char* file);
  void ReadGraph(const char* file);
  void ReadWeightedGraph(const char* file);
  bool isEdgeExist(int u, int v);
  void insertEdge(int u, int v, int w);
  void deleteEdge(int u, int v);

  inline int n() const { return n_; }
  inline int m() const { return m_; }
};

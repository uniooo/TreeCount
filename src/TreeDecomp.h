#pragma once
#define _CRT_SECURE_NO_DEPRECATE

#include <xmmintrin.h>

#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <chrono>

#include "Graph.h"
#include "Util.h"

using std::map;
using std::set;
using std::sort;
using std::vector;
using std::pair;

struct Node;

class TreeDecomp {
 protected:
  Graph G_;

  vector<vector<int> > neighbor_, length_;
  vector<vector<count_t> >
      edge_count_;  // record the number of shortest paths over short-cut edges.
                    // aligned with neighbor_, length_.

  vector<map<int, count_t> > C_;  // count the number of shortest paths.

  vector<int> delete_vertex_order_;  // the order vertex is deleted.

  int tree_width_ = 0;  // Tree width is derived from deleting node,i.e. max-deg
                        // vertex when deleting.
  int tree_height_ = 0;

  int* rank_ = NULL;  // the rank is the actual position in
                      // the tree_ vector. i.e. vertex id ---> tree_pos

  vector<Node> tree_;  // the decomposed tree structure
  int root_ = 0;       // root of the decompose tree, the position in tree_

  // RMQ Index for finding LCA
  vector<int> euler_tour_;      // The Euler Tour of the decomposed tree:
                                // euler_tour_pos ---> tree_pos
  int* euler_tour_pos_ = NULL;  // tree_pos ---> position in the euler_tour_
  vector<vector<int> > rmq_index_;

  // IO
  FILE* fout_;

  const int INF = 999999999;

 public:
  TreeDecomp();
  TreeDecomp(const char* filein, const char* fileout);
  ~TreeDecomp() {
    if (fout_ != NULL) fclose(fout_);
    if (rank_ != NULL) free(rank_);
  }

  //---
  void Reduce();
  void MakeTree();
  //void MakeIndexLocal();
  void MakeIndex();

  void EulerTour(int p, int height);  // generate euler_tour list

  static const int SIZEOFINT = sizeof(int);
  void PrintIntArray(int* a, int n);
  void PrintIntVector(vector<int>& a);
  void PrintLLVector(vector<count_t>& a);
  void PrintIndex();
  int GetParentNode(int vid, vector<int>& neighbors);

  void MakeRMQIndex();
  int DistanceQueryAncestorToPosterity(int p, int q);
  void makeIndexDFS(int p, vector<int>& list, int* toList);
  void MakeCountIndex();
  int DistanceCountQuery(int p, int q, count_t& count);

 private:
  void InsertCount(int u, int v, count_t count);
  void UpdateCount(int u, int v, count_t count);
  int LCAQuery(int _p, int _q);
};

struct Node {
  int vertex;
  vector<int> vertices, vertices_length;
  vector<int> pos, dis;
  vector<count_t> cnt;  // count the number to each ancestor.
  vector<int> child;
  int parent;
  int height;

  Node() {
    vertices.clear();
    vertices_length.clear();
    pos.clear();
    dis.clear();
    cnt.clear();
    vertex = -1;
    parent = -1;
    height = 0;
  }
};

class TreeQuery {
 private:
  FILE* findex_ = NULL;
  static const int SIZEOFINT = sizeof(int);
  int n_ = 0;                     // number of vertices.
  int tree_size_ = 0;             // length of the tree
  int* tree_node_height_ = NULL;  //
  int* rank_ = NULL;              // vertex id --> position in tree
  int* euler_tour_pos_ = NULL;    // tree_pos --> postion in the euler tour
  int rmq_index_size_ = 0;
  int euler_tour_size_ = 0;

  int tree_height_ = 0;
  int tree_width_ = 0;

  int** rmq_index_ = NULL;

  int root_ = 0;  // root of the tree_, typically 0
  int* tree_child_size_ =
      NULL;  // tree_child_size_[i]: child size of each tree node i
  int** tree_child_ = NULL;  // tree_child[i]: children of each tree node i

  int* pos_size_ = NULL;
  int** pos_ = NULL;

  int* self_pos_ = NULL;
  int** dis_ = NULL;
  count_t** cnt_ = NULL;

  const int INF = 999999999;
  //const int INF = MAXDIS_;


  int *LOG2 = NULL, *LOGD = NULL;

  // for testing
  long long query_count_ = 0;

 public:
  TreeQuery();
  ~TreeQuery();
  TreeQuery(const char* index_file);
  void ReadIndex(const char* index_file);
  inline int LCAQuery(int _p, int _q);
  int DistanceQuery(int p,
                    int q);  // return the distance between vid:p and vid:q

  int DistanceCountQuery(int p, int q, count_t& count);

  void HopSizeCount(const char* queryfile);

  void QueryTest(const char* queryfile, const char* resultfile);
  void PrintIndex(const char* index_file);

 private:
  void ScanIntArray(int* a, int n);
  int ScanLLVector(count_t*& a);
  int ScanIntVector(int*& a);  // return the length of vector.
  void InitLOG();
};

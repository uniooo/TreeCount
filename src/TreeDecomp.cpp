#include "TreeDecomp.h"

TreeDecomp::TreeDecomp() {
  delete_vertex_order_.clear();
  length_.clear();
  neighbor_.clear();
  edge_count_.clear();
  fout_ = NULL;
}

TreeDecomp::TreeDecomp(const char* filein, const char* fileout) : TreeDecomp() {
  G_ = Graph();

#ifndef WEIGHTED_GRAPH_
  G_.ReadGraph(filein);
#else
  G_.ReadWeightedGraph(filein);
#endif

  fout_ = fopen(fileout, "wb");

  vector<int> vectmp;
  vectmp.clear();
  vector<count_t> vecLLtmp;
  for (int i = 0; i <= G_.n(); ++i) {
    neighbor_.push_back(vectmp);  // record neighbors of each vertex, including
                                  // short-cut neighbors
    length_.push_back(vectmp);    // record the edge length between each vertex
                                  // and its neighbors
    edge_count_.push_back(vecLLtmp);  // record the count of shortest paths who
                                      // have the edge length
  }

  for (int i = 0; i <= G_.n(); ++i) {
    map<int, count_t> tmp;
    tmp.clear();
    C_.push_back(tmp);
  }
  for (int i = 0; i <= G_.n(); ++i) {
    for (auto v : G_.E_[i]) {
      InsertCount(i, v.first, 1);
    }
  }
}

void TreeDecomp::Reduce() {
  // To implement node selection in smallest-degree-first order.
  static int* DD;
  int* tmp_DD;

  // An extra order based on the sum of original degree and increased degree,
  // that is, do not count the removed edges
  static int* DD2;
  int* tmp_DD2;

  struct SelectElement {
    int x;
    SelectElement(int x) { this->x = x; }
    bool operator<(const SelectElement se) const {
      if (DD[x] != DD[se.x]) return DD[x] < DD[se.x];
      if (DD2[x] != DD2[se.x]) return DD2[x] < DD2[se.x];
      return x < se.x;
    }
  };

  set<SelectElement> deg;
  DD = (int*)malloc(sizeof(int) * (G_.n() + 1));
  tmp_DD = (int*)malloc(sizeof(int) * (G_.n() + 1));
  DD2 = (int*)malloc(sizeof(int) * (G_.n() + 1));
  tmp_DD2 = (int*)malloc(sizeof(int) * (G_.n() + 1));
  bool* changed = (bool*)malloc(sizeof(int) * (G_.n() + 1));

  deg.clear();

  for (int i = 1; i <= G_.n(); ++i) {
    DD[i] = G_.D_[i];
    tmp_DD[i] = G_.D_[i];
    DD2[i] = G_.D_[i];
    tmp_DD2[i] = G_.D_[i];
    changed[i] = false;
    deg.insert(SelectElement(i));
  }

  bool* exist = (bool*)malloc(sizeof(bool) * (G_.n() + 1));
  for (int i = 1; i <= G_.n(); ++i) exist[i] = true;
  delete_vertex_order_.clear();
  int delete_cnt = 0;

  while (!deg.empty()) {
    int x = (*deg.begin()).x;  // select x with the smallest degree
    while (changed[x]) {       // reinsert the element whose degree has changed
                               // to the right place in the set.
      deg.erase(SelectElement(x));
      DD[x] = tmp_DD[x];
      DD2[x] = tmp_DD2[x];
      deg.insert(SelectElement(x));
      changed[x] = false;
      x = (*deg.begin()).x;
    }

#ifdef PRINT_REDUCE_
    printf("deleted %d vertices. node degree = %d.\n", ++delete_cnt, G_.D_[x]);
#endif

#ifdef PRINT_REDUCE_LARGE_DEGREE_
    if (unlikely(G_.D_[x] > LARGE_DEGREE_))
      printf("deleted %d vertices. node degree = %d.\n", ++delete_cnt,
             G_.D_[x]);
    else
      ++delete_cnt;
#endif
    // correct x is selected
    delete_vertex_order_.push_back(x);
    deg.erase(deg.begin());
    exist[x] = false;  // mark vertex x as deleted;

    vector<int> tmp_neighbor, tmp_length;
    tmp_neighbor.clear();
    tmp_length.clear();

    vector<count_t> tmp_edge_count;
    tmp_edge_count.clear();

    // fetch x's neighbors y
    for (auto y : G_.E_[x]) {
      if (exist[y.first]) {
        tmp_neighbor.push_back(y.first);
        tmp_length.push_back(y.second);
        tmp_edge_count.push_back(C_[x][y.first]);
      }
    }

    // delete edges between x,y
    for (int y : tmp_neighbor) {
      G_.deleteEdge(x, y);
      tmp_DD[y] = G_.D_[y];  // here do not reduce DD2
      changed[y] = true;
    }
    // each pair u,v insert or update edge
    for (int pu = 0; pu < tmp_neighbor.size(); ++pu) {
      for (int pv = pu + 1; pv < tmp_neighbor.size(); ++pv) {
        if (pu != pv) {
          int u = tmp_neighbor[pu], v = tmp_neighbor[pv];
          if (G_.isEdgeExist(u, v)) {
            int newlength = tmp_length[pu] + tmp_length[pv];
            if (G_.E_[u][v] == newlength) {
              UpdateCount(u, v,
                          C_[u][v] + tmp_edge_count[pu] * tmp_edge_count[pv]);
            }
            if (G_.E_[u][v] > newlength) {
              G_.E_[u][v] = newlength;
              G_.E_[v][u] = newlength;
              UpdateCount(u, v, tmp_edge_count[pu] * tmp_edge_count[pv]);
            }
          } else {
            G_.insertEdge(u, v, tmp_length[pu] + tmp_length[pv]);
            InsertCount(u, v, tmp_edge_count[pu] * tmp_edge_count[pv]);
            tmp_DD[u] = G_.D_[u];
            tmp_DD[v] = G_.D_[v];
            ++tmp_DD2[u];
            ++tmp_DD2[v];
            changed[u] = changed[v] = true;
          }
        }
      }
    }

    if (tmp_neighbor.size() > tree_width_) tree_width_ = tmp_neighbor.size();
    neighbor_[x] = tmp_neighbor;
    length_[x] = tmp_length;
    edge_count_[x] = tmp_edge_count;
  }
  free(DD);
  free(tmp_DD);
  free(DD2);
  free(tmp_DD2);
  free(changed);
  free(exist);
}

void TreeDecomp::MakeTree() {
  rank_ = (int*)malloc(sizeof(int) * (G_.n() + 1));
  if (rank_ == NULL) return;
  tree_height_ = 0;
  for (int len = delete_vertex_order_.size() - 1; len >= 0; --len) {
    int vid = delete_vertex_order_[len];
    Node nod;
    nod.vertex = vid;
    nod.vertices = neighbor_[vid];
    nod.vertices_length = length_[vid];
    rank_[vid] = tree_.size();
    if (len == G_.n() - 1) {
      nod.parent = -1;
      nod.height = 1;
      tree_height_ = 1;
    } else {
      int parent_position = GetParentNode(vid, neighbor_[vid]);
      tree_[parent_position].child.push_back(tree_.size());
      nod.parent = parent_position;
      nod.height = tree_[parent_position].height + 1;
      if (nod.height > tree_height_) tree_height_ = nod.height;
    }
    tree_.push_back(nod);
  }
  // free(rank_);
}

void TreeDecomp::MakeIndex() {
  MakeRMQIndex();
  PrintIndex();

  MakeCountIndex();
}

void TreeDecomp::EulerTour(int p, int height) {  // p: position in tree_
  euler_tour_pos_[p] = euler_tour_.size();
  euler_tour_.push_back(p);
  for (int i = 0; i < tree_[p].child.size(); ++i) {
    EulerTour(tree_[p].child[i], height + 1);
    euler_tour_.push_back(
        p);  // add the parent vertex everytime when backtracking
  }
}

void TreeDecomp::PrintIntArray(int* a, int n) {
  fwrite(a, SIZEOFINT, n, fout_);
}

void TreeDecomp::PrintIntVector(vector<int>& a) {
  if (a.size() == 0) {
    int x = 0;
    fwrite(&x, SIZEOFINT, 1, fout_);
    return;
  }
  int x = a.size();
  fwrite(&x, SIZEOFINT, 1, fout_);
  for (int i = 0; i < a.size(); i++) {
    fwrite(&a[i], SIZEOFINT, 1, fout_);
  }
}

void TreeDecomp::PrintLLVector(vector<count_t>& a) {
  if (a.size() == 0) {
    int x = 0;
    fwrite(&x, SIZEOFINT, 1, fout_);
    return;
  }
  int x = a.size();
  fwrite(&x, SIZEOFINT, 1, fout_);
  for (int i = 0; i < a.size(); i++) {
    fwrite(&a[i], sizeof(count_t), 1, fout_);
  }
}

void TreeDecomp::PrintIndex() {
  // 1. G_.n()
  int vertex_size = G_.n();
  fwrite(&vertex_size, SIZEOFINT, 1, fout_);
  // 2. tree_.size() tree_.height
  int x = tree_.size();
  fwrite(&x, SIZEOFINT, 1, fout_);
  for (int i = 0; i < tree_.size(); i++) {
    fwrite(&tree_[i].height, SIZEOFINT, 1, fout_);
  }
  // rank_
  PrintIntArray(rank_, vertex_size + 1);
  // LCA - toRMQ - RMQIndex
  PrintIntArray(euler_tour_pos_, vertex_size + 1);
  x = rmq_index_.size();
  fwrite(&x, SIZEOFINT, 1, fout_);
  x = euler_tour_.size();
  fwrite(&x, SIZEOFINT, 1, fout_);
  for (int i = 0; i < rmq_index_.size(); i++) PrintIntVector(rmq_index_[i]);
  // rootDistance
  fwrite(&root_, SIZEOFINT, 1, fout_);

  for (int i = 0; i < tree_.size(); i++) {
    int t = tree_[i].child.size();
    fwrite(&t, SIZEOFINT, 1, fout_);
    for (int j = 0; j < t; j++) fwrite(&tree_[i].child[j], SIZEOFINT, 1, fout_);
  }
}

void TreeDecomp::MakeRMQIndex() {
  euler_tour_.clear();
  euler_tour_pos_ = (int*)malloc(sizeof(int) * (G_.n() + 1));
  EulerTour(root_, 1);
  rmq_index_.clear();
  rmq_index_.push_back(euler_tour_);
  int m = euler_tour_.size();
  for (int i = 2, k = 1; i < m; i *= 2, ++k) {
    vector<int> tmp;
    tmp.clear();
    tmp.resize(euler_tour_.size());
    for (int j = 0; j < m - i; ++j) {
      int x = rmq_index_[k - 1][j], y = rmq_index_[k - 1][j + i / 2];
      if (tree_[x].height < tree_[y].height)
        tmp[j] = x;
      else
        tmp[j] = y;
    }
    rmq_index_.push_back(tmp);
  }
}

int TreeDecomp::DistanceQueryAncestorToPosterity(int p, int q) {
  if (p == q) return 0;
  int x = rank_[p], y = rank_[q];
  return tree_[y].dis[tree_[x].pos[tree_[x].pos.size() - 1]];
}

void TreeDecomp::makeIndexDFS(int p, vector<int>& list, int* toList) {
  tree_[p].pos.resize(tree_[p].vertices.size() + 1);
  tree_[p].dis.resize(list.size());
  tree_[p].cnt.resize(list.size());

  for (int i = 0; i < tree_[p].vertices.size(); i++) {
    int j;
    for (j = 0; j < list.size(); j++)
      if (list[j] == tree_[p].vertices[i]) break;
    tree_[p].pos[i] = j;
  }
  tree_[p].pos[tree_[p].vertices.size()] = list.size();

  for (int i = 0; i < list.size(); i++) {
    tree_[p].dis[i] = INF;
    tree_[p].cnt[i] = 0;
  }
  int x = tree_[p].vertex;

  for (int i = 0; i < tree_[p].vertices.size(); i++) {
    // copy the known distance from vertices_length[i] to dis
    if (tree_[p].dis[toList[tree_[p].vertices[i]]] >
        tree_[p].vertices_length[i]) {
      tree_[p].dis[toList[tree_[p].vertices[i]]] = tree_[p].vertices_length[i];
      tree_[p].cnt[toList[tree_[p].vertices[i]]] =
          edge_count_[tree_[p].vertex][i];
    }
  }

  for (int i = 0; i < tree_[p].vertices.size(); ++i) {
    int x = tree_[p].vertices[i];  // x: vid of i-th neighbor of p
    int k;
    for (k = 0; k < list.size(); k++)
      if (list[k] == x) break;
    for (int j = 0; j < list.size(); j++) {
      if (k == j) continue;
      int y = list[j];
      int z, dis_new;
      count_t route_new;
#ifndef GLOBAL_DIS_

      if (k < j) {
        dis_new = INF;
        route_new = 0;
      } else if (k > j) {
        z = DistanceQueryAncestorToPosterity(y, x);

        // Here we use the origional count of the neighboring vertices.
        dis_new = tree_[p].vertices_length[i] + z;
        route_new = edge_count_[tree_[p].vertex][i] * tree_[rank_[x]].cnt[j];
      }

#else
      if (k < j) {
        z = DistanceQueryAncestorToPosterity(x, y);
        dis_new = tree_[p].dis[k] + z;
        route_new = 0;
      } else if (k > j) {
        z = DistanceQueryAncestorToPosterity(y, x);

        // Here we use the origional count of the neighboring vertices.
        dis_new = tree_[p].vertices_length[i] + z;
        route_new = edge_count_[tree_[p].vertex][i] * tree_[rank_[x]].cnt[j];
      }

#endif  // GLOBAL_DIS_

      if (tree_[p].dis[j] > dis_new) {
        tree_[p].dis[j] = dis_new;
        tree_[p].cnt[j] = route_new;
      } else if (tree_[p].dis[j] == dis_new) {
        tree_[p].cnt[j] = tree_[p].cnt[j] + route_new;
      }

#ifdef DEBUG_
      assert(j == toList[y]);
      assert(k == toList[tree_[p].vertices[i]]);
#endif  // DEBUG_
    }
  }
  toList[tree_[p].vertex] = list.size();
  list.push_back(tree_[p].vertex);
  for (int i = 0; i < tree_[p].child.size(); i++) {
    makeIndexDFS(tree_[p].child[i], list, toList);
  }
  list.pop_back();

  // sort(tree_[p].pos.begin(), tree_[p].pos.end());

  fwrite(&p, SIZEOFINT, 1, fout_);
#ifndef OMIT_POS_LIST_
  PrintIntVector(tree_[p].pos);
#else
  fwrite(&tree_[p].pos[tree_[p].pos.size()-1], SIZEOFINT, 1, fout_);
#endif
  PrintIntVector(tree_[p].dis);
  PrintLLVector(tree_[p].cnt);

#ifdef RELEASE_MEM_IMIDIATELY_AFTER_WRITTEN_
  vector<int> v1;
  v1.clear();
  tree_[p].dis.swap(v1);
  vector<count_t> v2;
  v2.clear();
  tree_[p].cnt.swap(v2);
#endif
}

void TreeDecomp::MakeCountIndex() {
  vector<int> list;
  list.clear();
  int* toList;
  toList = (int*)malloc(sizeof(int) * (G_.n() + 1));
  if (toList == NULL) {
    printf("toList is NULL");
    return;
  }

  tree_[root_].pos.clear();

  toList[tree_[root_].vertex] = 0;
  list.push_back(tree_[root_].vertex);
  tree_[root_].pos.push_back(0);

  fwrite(&root_, SIZEOFINT, 1, fout_);

  // fwrite the node's own position instead of it's pos list.
#ifndef OMIT_POS_LIST_
  PrintIntVector(tree_[root_].pos);
#else
  fwrite(&tree_[root_].pos[tree_[root_].pos.size()-1], SIZEOFINT, 1, fout_);
#endif
  
  PrintIntVector(tree_[root_].dis);
  PrintLLVector(tree_[root_].cnt);
  for (int i = 0; i < tree_[root_].child.size(); i++)
    makeIndexDFS(tree_[root_].child[i], list, toList);

  free(toList);
}

void TreeDecomp::InsertCount(int u, int v, count_t count) {
  if (C_[u].find(v) != C_[u].end()) return;
  C_[u].insert(std::make_pair(v, count));
  C_[v].insert(std::make_pair(u, count));
}

void TreeDecomp::UpdateCount(int u, int v, count_t count) {
  if (C_[u].find(v) == C_[u].end()) return;
  C_[u][v] = count;
  C_[v][u] = count;
}

inline int TreeDecomp::LCAQuery(int _p, int _q) {
  int p = euler_tour_pos_[_p], q = euler_tour_pos_[_q];

  if (p > q) {
    int x = p;
    p = q;
    q = x;
  }
  int len = q - p + 1;
  int* LOG2 = (int*)malloc(sizeof(int) * (G_.n_ * 2 + 10));
  int* LOGD = (int*)malloc(sizeof(int) * (G_.n_ * 2 + 10));

  int k = 0, j = 1;
  for (int i = 0; i < G_.n_ * 2 + 10; i++) {
    if (i > j * 2) {
      j *= 2;
      k++;
    }
    LOG2[i] = k;
    LOGD[i] = j;
  }

  int i = LOGD[len];
  k = LOG2[len];

  q = q - i + 1;
  if (tree_[rmq_index_[k][p]].height < tree_[rmq_index_[k][q]].height)
    return rmq_index_[k][p];
  else
    return rmq_index_[k][q];
}

int TreeDecomp::DistanceCountQuery(int p, int q, count_t& count) {
  if (p == q) return 0;
  int x = rank_[p], y = rank_[q];
  int lca = LCAQuery(x, y);
  int dis = INF;
  count_t cnt = 0;
  int ps = tree_[lca].pos.size();

  if (lca == x || lca == y) {
    if (lca == y) {
      int v = y;
      y = x;
      x = v;
      v = p;
      p = q;
      q = v;
    }
    dis = tree_[y].dis[tree_[x].pos[tree_[x].pos.size() - 1]];
    count = tree_[y].cnt[tree_[x].pos[tree_[x].pos.size() - 1]];
    ps--;
  }
  int tmp_dis = INF;
  vector<int>&dx = tree_[x].dis, &dy = (tree_[y].dis), &p2 = tree_[lca].pos;
  vector<count_t>&cx = tree_[x].cnt, &cy = tree_[y].cnt;
  for (int i = 0; i < ps; i++) {
    int tmp = dx[p2[i]] + dy[p2[i]];
    if (tmp_dis > tmp) {
      tmp_dis = tmp;
      cnt = cx[p2[i]] * cy[p2[i]];
    } else if (tmp_dis == tmp) {
      cnt += cx[p2[i]] * cy[p2[i]];
    }
  }

  if (dis > tmp_dis) {
    dis = tmp_dis;
    count = cnt;
  } else if (dis == tmp_dis) {
    count += cnt;
  }
  return dis;
}

int TreeDecomp::GetParentNode(int vid, vector<int>& neighbors) {
  int largest_rank_neighbor = neighbors[0];

  for (auto neighbor : neighbors) {
    if (rank_[neighbor] > rank_[largest_rank_neighbor]) {
      largest_rank_neighbor = neighbor;
    }
  }

#ifdef CHECK_NEIGHBORS 
  int p = rank_[largest_rank_neighbor];
  vector<int> a = tree_[p].vertices;
  assert(tree_[p].vertex >= 0);
  a.push_back(tree_[p].vertex);
  sort(a.begin(), a.end());
  sort(neighbors.begin(), neighbors.end());

  int i = 0, j = 0;
  for (; i < neighbors.size() && j < a.size();) {
    if (neighbors[i] == a[j]) {
      ++i;
      ++j;
    } else if (neighbors[i] < a[j])
      break;
    else
      ++j;
  }
  if (i >= neighbors.size()) return p;
  printf("wrong!!!!");
#else
  return rank_[largest_rank_neighbor];  
#endif
}

TreeQuery::TreeQuery() {}

TreeQuery::TreeQuery(const char* index_file) {
  ReadIndex(index_file);
  InitLOG();
}

TreeQuery::~TreeQuery() {
  if (findex_ != NULL) fclose(findex_);
  if (tree_node_height_ != NULL) {
    free(tree_node_height_);
    tree_node_height_ = NULL;
  }

  if (rank_ != NULL) {
    free(rank_);
    rank_ = NULL;
  }

  if (euler_tour_pos_ != NULL) {
    free(euler_tour_pos_);
    euler_tour_pos_ = NULL;
  }
  if (rmq_index_ != NULL) {
    for (int i = 0; i < rmq_index_size_; ++i) {
      free(rmq_index_[i]);
    }
    free(rmq_index_);
  }

  if (tree_child_size_ != NULL) {
    free(tree_child_size_);
    tree_child_size_ = NULL;
  }

  if (tree_child_ != NULL) {
    for (int i = 0; i < tree_size_; ++i) free(tree_child_[i]);
    free(tree_child_);
    tree_child_ = NULL;
  }

  if (pos_size_ != NULL) {
    free(pos_size_);
    pos_size_ = NULL;
  }

  if (pos_ != NULL) {
    for (int i = 0; i < tree_size_; ++i) free(pos_[i]);
    free(pos_);
    pos_ = NULL;
  }

  if (dis_ != NULL) {
    for (int i = 0; i < tree_size_; ++i) free(dis_[i]);
    free(dis_);
    dis_ = NULL;
  }

  if (cnt_ != NULL) {
    for (int i = 0; i < tree_size_; ++i) free(cnt_[i]);
    free(cnt_);
    cnt_ = NULL;
  }

  if (LOG2 != NULL) {
    free(LOG2);
    LOG2 = NULL;
  }

  if (LOGD != NULL) {
    free(LOGD);
    LOGD = NULL;
  }
}

void TreeQuery::ReadIndex(const char* index_file) {
  findex_ = fopen(index_file, "rb");

  const auto beg = std::chrono::steady_clock::now();

  fread(&n_, SIZEOFINT, 1, findex_);
  fread(&tree_size_, SIZEOFINT, 1, findex_);

  if (tree_size_ < 1) {
    printf("tree_size_ < 1");
    return;
  }

  // read node_height, array of size (tree_size_ + 1)
  if (tree_node_height_ == NULL)
    tree_node_height_ = (int*)malloc(sizeof(int) * (tree_size_ + 1));
  if (tree_node_height_ == NULL) {
    printf("tree_node_height == NULL");
    return;
  }

  ScanIntArray(tree_node_height_, tree_size_);

  // read rank_, array of size (n_ + 1)
  if (rank_ == NULL) rank_ = (int*)malloc(SIZEOFINT * (n_ + 1));
  if (rank_ == NULL) {
    printf("rank_ == NULL");
    return;
  }
  ScanIntArray(rank_, n_ + 1);

  // read euler_tour_pos_
  if (euler_tour_pos_ == NULL)
    euler_tour_pos_ = (int*)malloc(SIZEOFINT * (n_ + 1));
  if (euler_tour_pos_ == NULL) {
    printf("euler_tour_pos_ == NULL");
    return;
  }

  ScanIntArray(euler_tour_pos_, n_ + 1);

  // read rmq_index_size_, euler_tour_size_
  fread(&rmq_index_size_, SIZEOFINT, 1, findex_);
  fread(&euler_tour_size_, SIZEOFINT, 1, findex_);

  // read rmq_index_
  if (rmq_index_size_ < 1) {
    printf("rmq_index_size_ < 1");
    return;
  }
  if (rmq_index_ == NULL)
    rmq_index_ = (int**)malloc(sizeof(int*) * rmq_index_size_);
  if (rmq_index_ == NULL) {
    printf("rmq_index_ == NULL");
    return;
  }
  for (int i = 0; i < rmq_index_size_; ++i) {
    ScanIntVector(rmq_index_[i]);
  }

  // read root
  fread(&root_, SIZEOFINT, 1, findex_);

  // read tree_child_size, tree_child
  if (tree_child_size_ == NULL) {
    tree_child_size_ = (int*)malloc(SIZEOFINT * tree_size_);
  }
  if (tree_child_size_ == NULL) {
    printf("tree_child_size_ == NULL ");
    return;
  }

  if (tree_child_ == NULL) {
    tree_child_ = (int**)malloc(sizeof(int*) * tree_size_);
  }
  if (tree_child_ == NULL) {
    printf("tree_child_ == NULL");
    return;
  }
  for (int i = 0; i < tree_size_; ++i) {
    int child_size = ScanIntVector(tree_child_[i]);
    tree_child_size_[i] = child_size;
  }

  // read pos_size_, pos_, dis_;
  #ifndef OMIT_POS_LIST_
  if (pos_size_ == NULL) pos_size_ = (int*)malloc(SIZEOFINT * tree_size_);
  if (pos_size_ == NULL) {
    printf("pos_size_ == NULL");
    return;
  }

  if (pos_ == NULL) pos_ = (int**)malloc(sizeof(int*) * tree_size_);
  if (pos_ == NULL) {
    printf("pos_ == NULL");
    return;
  }
  #else
  if (self_pos_ == NULL) self_pos_ = (int*)malloc(SIZEOFINT * tree_size_);
  if (self_pos_ == NULL) {
    printf("self_pos_ == NULL");
    return;
  }
  #endif

  if (dis_ == NULL) dis_ = (int**)malloc(sizeof(int*) * tree_size_);
  if (dis_ == NULL) {
    printf("dis_ == NULL");
    return;
  }

  if (cnt_ == NULL) cnt_ = (count_t**)malloc(sizeof(count_t*) * tree_size_);
  if (cnt_ == NULL) {
    printf("cnt_ == NULL");
    return;
  }

  for (int i = 0; i < tree_size_; ++i) {
    int p;
    fread(&p, SIZEOFINT, 1, findex_);

#ifndef OMIT_POS_LIST_
    pos_size_[p] = ScanIntVector(pos_[p]);
    if (pos_size_[p] > tree_width_) tree_width_ = pos_size_[p];
#else
    fread(&self_pos_[p], SIZEOFINT, 1, findex_);
#endif // !OMIT_POS_LIST_

    int n = ScanIntVector(dis_[p]);
    if (n > tree_height_) tree_height_ = n;
    ScanLLVector(cnt_[p]);
  }
  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  printf("Finished loading index. Cost %f seconds\n", std::chrono::duration<double>(dif).count());
  if (findex_ != NULL) {
    fclose(findex_);
    findex_ = NULL;
  }
  printf("tree height: %d\n", tree_height_);
#ifndef OMIT_POS_LIST_
  printf("tree width: %d\n", tree_width_);
#endif  // !OMIT_POS_LIST_
}

inline int TreeQuery::LCAQuery(int _p, int _q) {
  int p = euler_tour_pos_[_p], q = euler_tour_pos_[_q];

  if (p > q) {
    int x = p;
    p = q;
    q = x;
  }
  int len = q - p + 1;

  int i = LOGD[len], k = LOG2[len];

  q = q - i + 1;
  if (tree_node_height_[rmq_index_[k][p]] < tree_node_height_[rmq_index_[k][q]])
    return rmq_index_[k][p];
  else
    return rmq_index_[k][q];
}

int TreeQuery::DistanceQuery(int p, int q) {
  if (p == q) return 0;
  int x = rank_[p], y = rank_[q];
  int lca = LCAQuery(x, y);
  if (lca == x || lca == y) {
    query_count_++;
    if (lca == y) {
      int v = y;
      y = x;
      x = v;
      v = p;
      p = q;
      q = v;
    }
    return dis_[y][pos_[x][pos_size_[x] - 1]];
  } else {
    int res = INF;
    int *dx = dis_[x], *dy = dis_[y], *p2 = pos_[lca];
    int ps = pos_size_[lca];
    for (int i = 0; i < ps; i++) {
      query_count_++;
      int tmp = dx[p2[i]] + dy[p2[i]];
      if (res > tmp) res = tmp;
    }
    return res;
  }
}

int TreeQuery::DistanceCountQuery(int p, int q, count_t& count) {
  if (p == q) {
    count = 1;
    return 0;
  }
  int x = rank_[p], y = rank_[q];
  int lca = LCAQuery(x, y);
  int dis = INF;
  count_t cnt = 0;
 #ifndef OMIT_POS_LIST_
  int ps = pos_size_[lca];
  int position = pos_[lca][pos_size_[lca] - 1];
  #else
  int position = self_pos_[lca];
  #endif
  if (lca == x || lca == y) {
    if (lca == y) {
      int v = y;
      y = x;
      x = v;
      v = p;
      p = q;
      q = v;
    }

#ifndef OMIT_POS_LIST_
    int pos = pos_[x][pos_size_[x] - 1];
#else
    int pos = position;
#endif
    dis = dis_[y][pos];
    count = cnt_[y][pos];
    position--;
  }
  int tmp_dis = INF;

  int *dx = dis_[x], *dy = dis_[y];
  count_t *cx = cnt_[x], *cy = cnt_[y];
 
#ifdef __unix__
  _mm_prefetch(dx, _MM_HINT_T0);
  _mm_prefetch(dy, _MM_HINT_T0);
  _mm_prefetch(cx, _MM_HINT_T0);
  _mm_prefetch(cy, _MM_HINT_T0);
#endif
  // int ps = pos_2Size[lca];

  for (int i = 0; i <= position; i++) {
    query_count_++;
    int tmp = dx[i] + dy[i];
    if (tmp_dis > tmp) {
      tmp_dis = tmp;
      cnt = cx[i] * cy[i];
    } else if (tmp_dis == tmp) {
      cnt += cx[i] * cy[i];
    }
  }

  if (dis > tmp_dis) {
    dis = tmp_dis;
    count = cnt;
  } else if (dis == tmp_dis) {
    count += cnt;
  }

  return dis;
}

void TreeQuery::HopSizeCount(const char* queryfile) {
  FILE* fin = fopen(queryfile, "r");

  int num_queries = 0;  // number of testcases;
  if (fscanf(fin, "%d", &num_queries) != 1)
    printf("wrong input: num_queries\n");
  vector<pair<int, int>> queries;
  for (int i = 0; i < num_queries; ++i) {
    int a, b;
    if (fscanf(fin, "%d %d", &a, &b) != 2)
      printf("wrong input: no enough query data\n");
    queries.push_back({a, b});
  }
  fclose(fin);
  fin = NULL;
  
  long long sum_query = 0;
  long long max_query = 0;
  for (const auto query : queries) {
    count_t count = 0;
    query_count_ = 0;
    DistanceCountQuery(query.first, query.second, count);

    sum_query += query_count_;
    max_query = max_query > query_count_ ? max_query : query_count_;
  }
  printf("Average HopSize: %lld Maximum HopSize: %lld\n",
         sum_query / num_queries, max_query);

}

void TreeQuery::QueryTest(const char* queryfile, const char* resultfile) {
  FILE* fin = fopen(queryfile, "r");
  FILE* fout = fopen(resultfile, "w");

  int num_queries = 0;  // number of testcases;
  if (fscanf(fin, "%d", &num_queries) != 1)
    printf("wrong input: num_queries\n");
  fprintf(fout, "%d\n", num_queries);
  vector<pair<int, int>> queries;
  for (int i = 0; i < num_queries; ++i) {
    int a, b;
    if (fscanf(fin, "%d %d", &a, &b) != 2)
      printf("wrong input: no enough query data\n");
    queries.push_back({a, b});
  }
  fclose(fin);
  fin = NULL;
  vector<count_t> results;
  const auto beg = std::chrono::steady_clock::now();
  for (const auto query : queries) {
    count_t count = 0;
    DistanceCountQuery(query.first, query.second, count);
    results.push_back(count);
  }
  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  printf("Query costs %f micro seconds in average\n",
         std::chrono::duration<double, std::micro>(dif).count() / num_queries);

  for (const auto result : results) {
    fprintf(fout, "%lld\n", result);
  }
  fclose(fout);
  fout = NULL;
}

void TreeQuery::PrintIndex(const char* index_file) {
  FILE* fout = fopen(index_file, "w");
  fprintf(fout, "n_ = %d\n", n_);
  fprintf(fout, "tree_size_ = %d\n", tree_size_);
  
  for (int i = 0; i < tree_size_; ++i) {
    for (int j = 0; j < tree_node_height_[i] - 1; ++j) {
      fprintf(fout, "%d %d;", dis_[i][j], cnt_[i][j]);
    }
    fprintf(fout, "\n");
  }

}

void TreeQuery::ScanIntArray(int* a, int n) { fread(a, SIZEOFINT, n, findex_); }

int TreeQuery::ScanLLVector(count_t*& a) {
  int n;
  fread(&n, SIZEOFINT, 1, findex_);
  a = (count_t*)malloc(sizeof(count_t) * n);
  fread(a, sizeof(count_t), n, findex_);
  return n;
}

int TreeQuery::ScanIntVector(int*& a) {
  int n;
  fread(&n, SIZEOFINT, 1, findex_);
  a = (int*)malloc(sizeof(int) * n);
  ScanIntArray(a, n);
  return n;
}

void TreeQuery::InitLOG() {
  if (LOG2 == NULL) LOG2 = (int*)malloc(sizeof(int) * (n_ * 2 + 10));
  if (LOGD == NULL) LOGD = (int*)malloc(sizeof(int) * (n_ * 2 + 10));

  if (LOG2 == NULL) {
    printf("LOG2 == NULL");
    return;
  }
  if (LOGD == NULL) {
    printf("LOGD == NULL");
    return;
  }

  int k = 0, j = 1;
  for (int i = 0; i < n_ * 2 + 10; i++) {
    if (i > j * 2) {
      j *= 2;
      k++;
    }
    LOG2[i] = k;
    LOGD[i] = j;
  }
}

#define _CRT_SECURE_NO_DEPRECATE

#include <getopt.h>

#include <chrono>
#include <iostream>

#include "Graph.h"
#include "TreeDecomp.h"

using std::string;

void exe(string graph_filename, string index_filename) {
  if (graph_filename == "" || index_filename == "") {
    printf("Error: Missing filename\n");
    return;
  }

  TreeDecomp td(graph_filename.c_str(), index_filename.c_str());
  const auto beg = std::chrono::steady_clock::now();
  td.Reduce();
  const auto reduce_time = std::chrono::steady_clock::now() - beg;
  printf("Reduce is OK, time from begin: %f ms\n",
         std::chrono::duration<double, std::milli>(reduce_time).count());
  td.MakeTree();
  const auto maketree_time = std::chrono::steady_clock::now() - beg;
  printf("MakeTree is OK, time from begin: %f ms\n",
         std::chrono::duration<double, std::milli>(maketree_time).count());
  td.MakeIndex();
  const auto makeindex_time = std::chrono::steady_clock::now() - beg;
  printf("MakeIndex is OK, time from begin: %f ms\n",
         std::chrono::duration<double, std::milli>(makeindex_time).count());

  printf("The index has been written to: %s\n", index_filename.c_str());
}

void help_info() { 
  printf("TC_index -g graph_file -o output_file\n");
}

int main(int argc, char *argv[]) {

  string graph_filename;
  string index_filename;

  const char *optstring = "g:o:h";
  int opt = -1;
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
      case 'g':
        graph_filename = optarg;
        break;
      case 'o':
        index_filename = optarg;
        break;
      case 'h':
        help_info();
        return 0;
    }
  }
  exe(graph_filename, index_filename);
}

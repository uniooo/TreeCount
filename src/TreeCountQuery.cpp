//#define _CRT_SECURE_NO_DEPRECATE

#include <getopt.h>

#include <iostream>

#include "Graph.h"
#include "TreeDecomp.h"

using std::string;

void help_info() {
  printf("TC_query -i index_file -q query_file -o output_file\n");
  printf("-p: print index\n");
}

int main(int argc, char *argv[]) {

  string index_filename;
  string query_filename;
  string answer_filename;
  bool print_index = false;
  bool count_hop_size = false;
  const char *optstring = "i:q:o:phs";
  int opt = -1;
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
      case 'i':
        index_filename = optarg;
        break;
      case 'q':
        query_filename = optarg;
        break;
      case 'o':
        answer_filename = optarg;
        break;
      case 'p':
        print_index = true;
        break;
      case 's':
        count_hop_size = true;
        break;
      case 'h':
        help_info();
        return 0;
    }
  }

  TreeQuery tq(index_filename.c_str());
  if (print_index) {
    string index_printed = index_filename + ".txt";
    tq.PrintIndex(index_printed.c_str());
  } else if (count_hop_size) {
    tq.HopSizeCount(query_filename.c_str());
  } else {
    tq.QueryTest(query_filename.c_str(), answer_filename.c_str());
  }

  return 0;
}

# TreeCount

This is the source code repository for our paper published in VLDB 2022.

## Paper Title: Efficient Shortest Path Counting on Large Road Networks

## Authors: Yu-Xuan Qiu, Dong Wen, Lu Qin, Wentao Li, Ronghua Li, Ying Zhang

## Citation: 

```
@article{DBLP:journals/pvldb/QiuWQLLZ22,
  author    = {Yu{-}Xuan Qiu and
               Dong Wen and
               Lu Qin and
               Wentao Li and
               Ronghua Li and
               Ying Zhang},
  title     = {Efficient Shortest Path Counting on Large Road Networks},
  journal   = {Proc. {VLDB} Endow.},
  volume    = {15},
  number    = {10},
  pages     = {2098--2110},
  year      = {2022},
  url       = {https://www.vldb.org/pvldb/vol15/p2098-qiu.pdf},
  timestamp = {Mon, 26 Sep 2022 17:09:16 +0200},
  biburl    = {https://dblp.org/rec/journals/pvldb/QiuWQLLZ22.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

## Usage

To compile the source, in `src/`:

`make && make clean`

To build index:

`TC_index -g graph_file -o output_file`

To query index:

`TC_query -i index_file -q query_file -o output_file`

`NY/` provides a sample road network, a sample query file, and the corresponding result.
#pragma once
#include<limits>

/*********** control settings ***********/
#define WEIGHTED_GRAPH_

//#define RELEASE_MEM_IMIDIATELY_AFTER_WRITTEN_

//#define OMIT_POS_LIST_

//#define GLOBAL_DIS_ 


#ifdef __GNUC__
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define likely(x) (x)
#define unlikely(x) (x)
#endif


/************ self-defined types & constants ***************/
typedef uint32_t count_t;
typedef int32_t vid_t;
typedef int32_t dis_t;
typedef int32_t treepos_t;

typedef uint32_t ui_t;

typedef std::pair<dis_t, vid_t> disvid_pair_t;

 //#define MAXDIS_ std::numeric_limits<int32_t>::max()
//#define MAXDIS_ = (std::numeric_limits <int32_t>::max() / 2);
#define MAXDIS_ 999999999



/********* debug settings *********/

//#define PRINT_REDUCE_LARGE_DEGREE_
#ifdef PRINT_REDUCE_LARGE_DEGREE_
constexpr auto LARGE_DEGREE_ = 600 ;
#endif // PRINT_REDUCE_LARGE_DEGREE_

//#define DEBUG_




//#include <sys/time.h>
//
//
//double GetTime(void) {
//  struct timeval tv;
//  gettimeofday(&tv, NULL);
//  return tv.tv_sec + tv.tv_usec * 1e-6;
//}
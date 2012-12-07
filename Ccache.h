#ifndef __HOPSCOTCH_HASHMAP__
#define __HOPSCOTCH_HASHMAP__

////////////////////////////////////////////////////////////////////////////////
// ConcurrentHopscotchHashMap Class
//
////////////////////////////////////////////////////////////////////////////////
//TERMS OF USAGE
//------------------------------------------------------------------------------
//
//	Permission to use, copy, modify and distribute this software and
//	its documentation for any purpose is hereby granted without fee,
//	provided that due acknowledgments to the authors are provided and
//	this permission notice appears in all copies of the software.
//	The software is provided "as is". There is no warranty of any kind.
//
//Authors:
//	Maurice Herlihy
//	Brown University
//	and
//	Nir Shavit
//	Tel-Aviv University
//	and
//	Moran Tzafrir
//	Tel-Aviv University
//
//	Date: July 15, 2008.  
//
////////////////////////////////////////////////////////////////////////////////
// Programmer : Moran Tzafrir (MoranTza@gmail.com)
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// INCLUDE DIRECTIVES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Add LRU function to each segment of the hash map
// Add hazard pointer to do memory reclaim
// Programmer: Ma Bo(redshift@outlook.com)
////////////////////////////////////////////////////////////////////////////////


#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include <malloc.h>
//#include <sys/queue.h>
#include <pthread.h>
#include <time.h>
#include <iostream>

#ifdef CCACHE_TEST
#include <time.h>
#define HRTIME_FOREVER  (0x7FFFFFFF)
typedef time_t ink_hrtime;

ink_hrtime ccache_time()
{
   time_t t;
   time(&t);
   return t;
}

ink_hrtime ink_hrtime_from_sec(unsigned int sec)
{
   return sec;
}

#define ink_get_hrtime ccache_time

#else
#include "ink_hrtime.h"
#endif

#define  _NULL_DELTA SHRT_MIN 

/*
 *  * Array-based tail queue functions.
 *   */
#define ARRAYQ_HEAD(name, type)\
    struct name {\
    type * volatile tqh_first;  \
    type * volatile tqh_last;\
}

#define ARRAYQ_ENTRY(type)\
                    struct {\
                    short volatile tqe_next;/* next element */\
                    short volatile tqe_prev;/* address of previous next element */\
                    }

#define ARRAYQ_INIT(head) {\
                    (head)->tqh_first = NULL;\
                    (head)->tqh_last = NULL;\
                    }

#define	ARRAYQ_FIRST(head)	((head)->tqh_first)
#define	ARRAYQ_LAST(head, headname)    ((head)->tqh_last)


#define ARRAYQ_INSERT_HEAD(head, elm, field) {\
                    if ( (head)->tqh_first != NULL){\
                    (elm)->field.tqe_next = (short)((head)->tqh_first - (elm));\
                    (head)->tqh_first->field.tqe_prev = (short)((elm) - (head)->tqh_first); \
                    } \
                    else {\
                    (head)->tqh_last = (elm);\
                    (elm)->field.tqe_next = _NULL_DELTA; \
                    } \
                    (head)->tqh_first = (elm);\
                    (elm)->field.tqe_prev = _NULL_DELTA;\
}

#define ARRAYQ_INSERT_TAIL(head, elm, field) {\
    (elm)->field.tqe_next = _NULL_DELTA;\
    if ( (head)->tqh_last != NULL ) {\
        (elm)->field.tqe_prev = (short)((head)->tqh_last - (elm));  \
        (head)->tqh_last->field.tqe_next = short((elm) - (head)->tqh_last); \
    } \
    else {\
        (elm)->field.tqe_prev = _NULL_DELTA; \
        (head)->tqh_first = (elm);\
    }\
    (head)->tqh_last = (elm);\
}

#define ARRAYQ_REMOVE(head, elm, field) {\
    if (((elm)->field.tqe_next) != _NULL_DELTA) {\
    if (((elm)->field.tqe_prev) == _NULL_DELTA ) {\
        ((elm) + (elm)->field.tqe_next)->field.tqe_prev = _NULL_DELTA; \
        (head)->tqh_first = (elm) + (elm)->field.tqe_next; \
    } \
    else {\
        ((elm) + (elm)->field.tqe_next)->field.tqe_prev += \
        (elm)->field.tqe_prev; \
        ((elm) + (elm)->field.tqe_prev)->field.tqe_next += \
        ((elm)->field.tqe_next); \
    } } \
    else {\
    if (((elm)->field.tqe_prev) == _NULL_DELTA ) \
    { \
        (head)->tqh_last = NULL;\
        (head)->tqh_first = NULL;  \
    } \
    else \
    { \
        (head)->tqh_last = (elm + (elm)->field.tqe_prev); \
        (head)->tqh_last->field.tqe_next = _NULL_DELTA; \
    } } \
}

#define CAS(a_ptr, a_oldVal, a_newVal) __sync_bool_compare_and_swap(a_ptr, a_oldVal, a_newVal)
#define ATOMIC_ADD(a_ptr,value) __sync_fetch_and_add(a_ptr,value)
//////////////////////////////////////////////////////////////////////////
//bit operations
//////////////////////////////////////////////////////////////////////////
#define bit_count(x) __builtin_popcount (x)
#define bit_count64(x) __builtin_popcountll (x)

inline int first_lsb_bit_indx(unsigned int x) {
	if(0==x) 
		return -1;
	return __builtin_ffs(x)-1;
}

inline int first_msb_bit_indx(unsigned int x) {
	if(0==x) 
		return -1;
	return  __builtin_clz(x)-1;
}


class Memory {
	public:
		static const int CACHELINE_SIZE = 64;

		inline static void* byte_malloc(const size_t size)	{return malloc(size);}
		inline static void  byte_free(void* mem)					{free(mem);}

		inline static void* byte_aligned_malloc(const size_t size) {return memalign(CACHELINE_SIZE, size);}
		inline static void* byte_aligned_malloc(const size_t size, const size_t alignment) {return memalign(alignment, size);}
		inline static void  byte_aligned_free(void* mem)	{free(mem);}
		
};

class TTASLock {
public:
    pthread_mutex_t _lock;

	TTASLock() {}
	~TTASLock() {pthread_mutex_destroy(&_lock);}

	inline void init() {pthread_mutex_init(&_lock,NULL);}

	inline void lock() {
	  pthread_mutex_lock(&_lock); 
	}

	inline bool tryLock() {
	   return ( 0 == pthread_mutex_trylock(&_lock));
	}

	inline void unlock() {
		pthread_mutex_unlock(&_lock);
	}
};

//HAZARD POINTER
typedef void (*freecb)(void *data);

template <uint32_t MAX_RETIRENUM>
struct hp_record {
public:
  hp_record() : 
   active(false),
   rcount(0),
   list(NULL)
   {
   }
  ~hp_record() {hpr_destroy();}
  inline bool hpr_isactive(){return active;}
  inline void hpr_setactive(bool _active){active = _active;}
  void hpr_free(void *data,freecb _cb = NULL);
  void hpr_init(void **_g_hps,size_t _size,int _threshold);
  void hpr_destroy();
  int hpr_check() {
    int c=0;
    rnode  *p = list;
    
    while(p)
    {
      c++;
      p = p->next;
    }
    return rcount-c;
  }
private:
  //inner struct
  struct rnode {
    struct rnode *next;
    void *data;
    freecb cb;
  };
  void **g_hps;
  size_t size;
  int threshold;
  //freecb cb;
  
  bool active;
  //retired node
  rnode  *freelist;
  rnode  *list;
  rnode rlist[MAX_RETIRENUM]; 
  int rcount;
  //protect rcount to be modified in destructor
  int _dummy;
};

template < uint32_t MAX_RETIRENUM>
void hp_record<MAX_RETIRENUM>::hpr_init(void **_g_hps,size_t _size,int _threshold)
{
  unsigned int i;
  rcount = 0;
  list = NULL;
  g_hps = _g_hps;
  size = _size;
  threshold = _threshold;

  freelist = &rlist[0];
  for(i=0; i<MAX_RETIRENUM-1; i++)
  {
    rlist[i].data = NULL;
    rlist[i].next = &rlist[i+1];
    rlist[i].cb = NULL;
  }  
  rlist[MAX_RETIRENUM-1].next = NULL;
}


#define HP_BACKOFF_INITIALIZER ((1 << 9) - 1) 
#define HP_BACKOFF_CEILING ((1 << 21) - 1)
typedef volatile unsigned int hp_backoff_t;
 
inline static void
hp_backoff_gb(volatile unsigned int *c)
{
	volatile unsigned int i;
	unsigned int ceiling;

	ceiling = *c;

	for (i = 0; i < ceiling; i++);

	ceiling <<= 1;
	ceiling &= HP_BACKOFF_CEILING;

	*c = ceiling;
	return;
}

template < uint32_t MAX_RETIRENUM>
void hp_record<MAX_RETIRENUM>::hpr_destroy()
{
   hp_backoff_t backoff = HP_BACKOFF_INITIALIZER;

   threshold = 0; 
   //backoff trying until success
   while(list) {
     hpr_free(list->data,list->cb);
     if(list)
       hp_backoff_gb(&backoff);
   }
   assert(rcount == 0);
   assert(list == NULL);
   active = false;
}

template <uint32_t MAX_RETIRENUM>
void hp_record<MAX_RETIRENUM>::hpr_free(void *data,freecb _cb)
{ 
#define MAX_HPNUM 256
  void *mhps[MAX_HPNUM];
  int i;
  rnode *node;
  rnode *prev,*cur;
    
  cur = list;
  prev = NULL;
  while(cur)
  {
    if(data <= cur->data)
      break;
    prev = cur;
    cur = cur->next;
  }
  //ignore dublicate free
  if(cur == NULL || data != cur->data)
  {
    //find a free node
    if(freelist == NULL)
    {
      //no freelist in array,get node from heap
      node = new rnode();
    }
    else
    {
      node = freelist;
      freelist = freelist->next;
    }
    node->data = data;
    node->cb = _cb;
    if(prev)
    {
      prev->next = node;
      node->next = cur;
    }
    else
    {
      node->next = cur;
      list = node;
    } 
    rcount++;
  }  
   
  if(rcount < threshold)
  {
    return; 
  }
  
  assert(list != NULL);

  //use cur as the sentry in qsort to divide the set
  //list is already a sorted list
  
  cur = list;
  int sidx=0,eidx=MAX_HPNUM-1;
  bool find = false;

  //first store global hps to local array
  for(i=0; i < size; i++)
  {
    if(g_hps[i])
    {
      if(g_hps[i] <= cur->data)
      {
        if(g_hps[i] == cur->data) find = true;
        mhps[sidx++] = g_hps[i];
      }
      else
        mhps[eidx--] = g_hps[i];

      if(sidx >= eidx)
        break;        
    }
  }
  
  if(!find)
  {
    void *data = cur->data;
    list = cur->next;
    cur->next = freelist;
    cur->data = NULL;       
    if(cur->cb) cur->cb(data);
    cur->cb = NULL;
    //usually this will not happen
    if(cur < &rlist[0] || cur > &rlist[MAX_RETIRENUM-1]) delete cur;
    else freelist = cur; 
    cur = list;
    prev = NULL;
    rcount--;
  }
  else
  {
    prev = cur;
    cur = cur->next;
  }

  sidx = eidx+1;
  eidx = MAX_HPNUM-1;
  
  //now start second stage to do qsort-similar routine to find cur
  //item larger than cur is put to right half,less equal than cur is put to left half
  //because cur->next must larger than cur,we goto right half then
  //repeat this loop,we can find all the node in list
  void *t;
  
  while(cur)
  {
    find = false;
    
    if(sidx >= eidx)
      break;
    
    t = mhps[sidx];
    while(sidx < eidx)
    {
      while(sidx < eidx && mhps[eidx] > cur->data)
        eidx--;
      mhps[sidx] = mhps[eidx];
      
      while(sidx < eidx && mhps[sidx] <= cur->data)
      {
        if(mhps[sidx] == cur->data) find = true;
        sidx++;
      }  
      mhps[eidx] = mhps[sidx];      
    }
    mhps[sidx] = t;
    if(mhps[sidx] <= cur->data) 
    {
        if(mhps[sidx] == cur->data) find=true; 
        sidx++;
    }
    eidx = MAX_HPNUM-1;
    
    if(!find)
    {
      void *data = cur->data;
      
      if(prev)
        prev->next = cur->next;
      else
        list = cur->next;
      cur->next = freelist;
      cur->data = NULL;
      if(cur->cb) cur->cb(data);
      cur->cb = NULL;
      //usually this will not happen
      if(cur < &rlist[0] || cur > &rlist[MAX_RETIRENUM-1]) delete cur;
      else freelist = cur;
      if(prev)
        cur = prev->next;
      else
        cur = list;
      rcount--;
    }
    else
    {
      prev = cur;
      cur = cur->next;
    }
  }
  return;  
}

template <uint32_t MAX_THREADS, uint32_t HP_NUM_PERTHREAD, uint32_t MAX_RETIRENUM>
class HAZARD_POINTER {
public:
  HAZARD_POINTER() {pool_lock.init();}
  ~HAZARD_POINTER() {}
  void hp_init();

  //return the index of allocated hp record
  int hp_register(int rthreshold=1);
  void hp_unregister(int idx);

  void hp_free(int idx, void *data, freecb cb);
  inline void hp_set(int idx,void *data,int hp_idx=0) 
  {
    g_hps[idx*HP_NUM_PERTHREAD + hp_idx] = data;
  }
  
  inline void hp_clear(int idx)
  {
    int i;
    for(i=0; i<HP_NUM_PERTHREAD; i++)
    {
      g_hps[idx*HP_NUM_PERTHREAD + i] = NULL;
    }
  }

  int hp_check(int idx)
  { 
  	if(idx < 0 || idx >= MAX_THREADS) return 0;
    return hp_record_pool[idx]->hpr_check();
  }
  
private:  
  //hp_record pool for thread
  TTASLock pool_lock;  
  hp_record<MAX_RETIRENUM> *hp_record_pool[MAX_THREADS];
  //global hp
  void *g_hps[MAX_THREADS*HP_NUM_PERTHREAD];
};

template <uint32_t MAX_THREADS, uint32_t HP_NUM_PERTHREAD, uint32_t MAX_RETIRENUM>
int HAZARD_POINTER<MAX_THREADS, HP_NUM_PERTHREAD, MAX_RETIRENUM>::hp_register(int rthreshold)
{
  int i;
  //allocate a hp record
  pool_lock.lock();
  for(i=0; i < MAX_THREADS; i++)
  {
  	if(hp_record_pool[i] == NULL)
	  hp_record_pool[i] = new hp_record<MAX_RETIRENUM>();
	else
     continue;
	
    if(hp_record_pool[i])
    {
      hp_record_pool[i]->hpr_setactive(true);
      break;
    }
  }
  pool_lock.unlock();
  if(i == MAX_THREADS)
    return -1;
  hp_record_pool[i]->hpr_init(&g_hps[0], MAX_THREADS*HP_NUM_PERTHREAD, rthreshold);
  return i;
}

template <uint32_t MAX_THREADS, uint32_t HP_NUM_PERTHREAD, uint32_t MAX_RETIRENUM>
void HAZARD_POINTER<MAX_THREADS, HP_NUM_PERTHREAD, MAX_RETIRENUM>::hp_unregister(int idx)
{
  if(idx < 0 || idx >= MAX_THREADS) return;
  delete hp_record_pool[idx];
  hp_record_pool[idx] = NULL;
}

template <uint32_t MAX_THREADS, uint32_t HP_NUM_PERTHREAD, uint32_t MAX_RETIRENUM>
void HAZARD_POINTER<MAX_THREADS, HP_NUM_PERTHREAD, MAX_RETIRENUM>::hp_free(int idx, void *data, freecb cb)
{
  if(idx < 0 || idx >= MAX_THREADS) return;
  hp_record_pool[idx]->hpr_free(data, cb);  
}


////////////////////////////////////////////////////////////////////////////////
//INNER CLASSES
////////////////////////////////////////////////////////////////////////////////
template <typename _tDATA>
class HASH_INT {
public:
	//you must define the following fields and properties
	static const unsigned int _EMPTY_HASH;
	static const unsigned int _BUSY_HASH;
	static const int _EMPTY_KEY;
	static const _tDATA _EMPTY_DATA;

	inline static unsigned int Calc(int key) {
		key ^= (key << 15) ^ 0xcd7dcd7d;
		key ^= (key >> 10);
		key ^= (key <<  3);
		key ^= (key >>  6);
		key ^= (key <<  2) + (key << 14);
		key ^= (key >> 16);
		return key;
	}

	inline static bool IsEqual(int left_key, int right_key) {
		return left_key == right_key;
	}

	inline static void relocate_key_reference(int volatile& left, const int volatile& right) {
		left = right;
	}

	inline static void relocate_data_reference(_tDATA volatile& left, const _tDATA volatile& right) {
		left = right;
	}
	
	inline static void dump_key(int key,char *title,char *title2)
	{		
		printf("ssl_cache %s %s %d\n",title,title2,key);
	}
};

template <typename _tDATA>
const unsigned int HASH_INT<_tDATA>::_EMPTY_HASH = 0;
template <typename _tDATA>
const unsigned int HASH_INT<_tDATA>::_BUSY_HASH  = 1;
template <typename _tDATA>
const int HASH_INT<_tDATA>::_EMPTY_KEY  = 0;
template <typename _tDATA>
const _tDATA HASH_INT<_tDATA>::_EMPTY_DATA = 0;

//hash str
template <typename _tDATA>
class HASH_STR {
public:
	//you must define the following fields and properties
	static const unsigned int _EMPTY_HASH;
	static const unsigned int _BUSY_HASH;
	static char * _EMPTY_KEY;
	static const _tDATA _EMPTY_DATA;

	inline static unsigned int Calc(char *key) {
	unsigned char *s = (unsigned char *)key;	/* unsigned string */
    unsigned int hval = 0x811c9dc5;
    /*
     * FNV-1a hash each octet in the buffer
     */
    while (*s) {

    	/* xor the bottom with the current octet */
    	hval ^= (unsigned int)*s++;

    	/* multiply by the 32 bit FNV magic prime mod 2^32 */
#if defined(NO_FNV_GCC_OPTIMIZATION)
    	hval *= FNV_32_PRIME;
#else
    	hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
#endif
    }

	 /* return our new hash value */
     return hval;
	}

	inline static bool IsEqual(char *left_key, char *right_key) {
		return !strcmp(left_key,right_key);
	}

	inline static void relocate_key_reference(char * volatile& left, char *volatile& right) {
		left = right;
	}

	inline static void relocate_data_reference(_tDATA volatile& left, const _tDATA volatile& right) {
		left = right;
	}
	
	inline static void dump_key(char *key,char *title,char *title2)
	{		
		printf("%s %s %s\n",title,title2,key);
	}
};

template <typename _tDATA>
const unsigned int HASH_STR<_tDATA>::_EMPTY_HASH = 0;
template <typename _tDATA>
const unsigned int HASH_STR<_tDATA>::_BUSY_HASH  = 1;
template <typename _tDATA>
char * HASH_STR<_tDATA>::_EMPTY_KEY  = NULL;
template <typename _tDATA>
const _tDATA HASH_STR<_tDATA>::_EMPTY_DATA = NULL;

//hash buf
struct s_buf{ 
  unsigned int len;
  void *data;
};

template <typename _tDATA>
class HASH_BUF {
public:
	//you must define the following fields and properties
	static const unsigned int _EMPTY_HASH;
	static const unsigned int _BUSY_HASH;
	static s_buf * _EMPTY_KEY;
	static const _tDATA _EMPTY_DATA;

	inline static unsigned int Calc(void *key) {
	    unsigned int hval = 0x811c9dc5;
		s_buf *buf = (s_buf *)key;
		
	    unsigned char *bp = (unsigned char *)buf->data;	/* start of buffer */
	    unsigned char *be = bp + buf->len;		/* beyond end of buffer */

	    /*
	     * FNV-1a hash each octet in the buffer
	     */
	    while (bp < be) {

		/* xor the bottom with the current octet */
		hval ^= (unsigned int)*bp++;

		/* multiply by the 32 bit FNV magic prime mod 2^32 */
#if defined(NO_FNV_GCC_OPTIMIZATION)
		hval *= FNV_32_PRIME;
#else
		hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
#endif
	    }

	    /* return our new hash value */
	    return hval;
	}

	inline static bool IsEqual(s_buf *left_key, s_buf *right_key) {
		s_buf *buf1 = (s_buf *)left_key;
		s_buf *buf2 = (s_buf *)right_key;
		if(buf1->len != buf2->len) return 0;
		return !memcmp((const void *)buf1->data, (const void *)buf2->data, buf1->len);
	}

	inline static void relocate_key_reference(s_buf * volatile& left, s_buf *volatile& right) {
		left = right;
	}

	inline static void relocate_data_reference(_tDATA volatile& left, const _tDATA volatile& right) {
		left = right;
	}

	inline static void dump_key(s_buf *key,char *title,char *title2)
	{
		char d[256];
		int n=0;
		int off=0;

        for( ; n < key->len; ++n)
		{
		  off += sprintf(d+off,"%02x",((char *)key->data)[n]);
		}
		
		printf("%s %s %s\n",title,title2,d);
	}
};

template <typename _tDATA>
const unsigned int HASH_BUF<_tDATA>::_EMPTY_HASH = 0;
template <typename _tDATA>
const unsigned int HASH_BUF<_tDATA>::_BUSY_HASH  = 1;
template <typename _tDATA>
s_buf *HASH_BUF<_tDATA>::_EMPTY_KEY  = NULL;
template <typename _tDATA>
const _tDATA HASH_BUF<_tDATA>::_EMPTY_DATA = NULL;


////////////////////////////////////////////////////////////////////////////////
// CLASS: NON-BLOCKING Read LRU Queue(concurrent enqueue single dequeue)
// ELEM must be a link list node
// The basic idea is from Google Concurrency Hashmap
// m_theQueue stores every access,and when we get lock(such as write or someting else),
// drain the buffer and reorder the lru list
// The FIFO is implemented following Tsigas Zhang's paper 
// "A Simple, Fast and Scalable Non-Blocking Concurrent FIFO Queue for Shared Memory Multiprocessor Systems"
////////////////////////////////////////////////////////////////////////////////
#define NULL0 (ELEM *)0
#define NULL1 (ELEM *)1

template <typename ELEM, uint32_t BUFF_SIZE>
class NBLRUQueue {
private:
    //max to BUFF_SIZE items
    ELEM *m_theQueue[BUFF_SIZE];

    /// @brief where a new element will be inserted
    volatile uint32_t m_writeIndex;

    /// @brief where the next element where be extracted from
    volatile uint32_t m_readIndex;

    ARRAYQ_HEAD(_lru, ELEM) lru_head;
public:
    volatile int32_t lru_nodenum;
    /// @brief constructor of the class
    NBLRUQueue();
    virtual ~NBLRUQueue();

    uint32_t size();
    
    void init();
    ELEM *lru_new_lock(ELEM* const a_data,volatile int * g_limit);
    bool lru_enqueue(const ELEM *a_data);
    int lru_dequeue(ELEM **a_data_p);
    bool lru_access(const ELEM *a_data);
    void lru_remove_lock(ELEM *a_data,volatile int * g_limit);
    void lru_update_lock();
    void lru_refresh_lock(ELEM* a_data);  
	ELEM *lru_getoldest_lock();
};


template <typename ELEM, uint32_t BUFF_SIZE>
NBLRUQueue<ELEM, BUFF_SIZE>::NBLRUQueue() :
    m_writeIndex(1),
    m_readIndex(0),
    lru_nodenum(0)
{
  int i;
  for(i=0;i<BUFF_SIZE;i++)
    m_theQueue[i] = NULL0;
  m_theQueue[0] = NULL1;
  ARRAYQ_INIT(&lru_head);
}

template <typename ELEM, uint32_t BUFF_SIZE>
void NBLRUQueue<ELEM, BUFF_SIZE>::init()
{
  unsigned int i;
  m_writeIndex = 1; 
  m_readIndex = 0;
  lru_nodenum = 0;

  for(i=0;i<BUFF_SIZE;i++)
    m_theQueue[i] = NULL0;
  m_theQueue[0] = NULL1;
  ARRAYQ_INIT(&lru_head);
}

template <typename ELEM, uint32_t BUFF_SIZE>
NBLRUQueue<ELEM, BUFF_SIZE>::~NBLRUQueue()
{
}

template <typename ELEM, uint32_t BUFF_SIZE>
uint32_t NBLRUQueue<ELEM, BUFF_SIZE>::size()
{
    uint32_t currentWriteIndex = m_writeIndex;
    uint32_t currentReadIndex  = m_readIndex;

    // let's think of a scenario where this function returns bogus data
    // 1. when the statement 'currentWriteIndex = m_writeIndex' is run
    // m_writeIndex is 3 and m_readIndex is 2. Real size is 1
    // 2. afterwards this thread is preemted. While this thread is inactive 2 
    // elements are inserted and removed from the queue, so m_writeIndex is 5
    // m_readIndex 4. Real size is still 1
    // 3. Now the current thread comes back from preemption and reads m_readIndex.
    // currentReadIndex is 4
    // 4. currentReadIndex is bigger than currentWriteIndex, so
    // m_totalSize + currentWriteIndex - currentReadIndex is returned, that is,
    // it returns that the queue is almost full, when it is almost empty
    
    if (currentWriteIndex >= currentReadIndex)
    {
        return (currentWriteIndex - currentReadIndex);
    }
    else
    {
        return (BUFF_SIZE + currentWriteIndex - currentReadIndex);
    }
}

//this function is fully concurrent,means when a enqueue is in processing,
//there may be other enqueue and dequee operations
template <typename ELEM,uint32_t BUFF_SIZE>
bool NBLRUQueue<ELEM, BUFF_SIZE>::lru_enqueue(const ELEM *a_data)
{
    uint32_t te,ate,temp;
    ELEM *tt;
    const ELEM *tnew;    
    int i=0;
    
    do {
      i++;
      te = m_writeIndex;
      ate = te;    
      tt = m_theQueue[ate];

      //the next slot of the tail
      temp = (ate + 1)%BUFF_SIZE;
      
      //we want to find the actual tail
      while(tt != NULL0 && tt != NULL1 )
      {
        //check tail's 
        if(te != m_writeIndex) break;
        //if tail meet head,it is possible that queue is full
        if(temp == m_readIndex) break;
        //now check the next cell
        tt = m_theQueue[temp];
        ate = temp;
        temp = (ate + 1)%BUFF_SIZE;        
      }

      //check the tail's consistency
      if(te != m_writeIndex) continue;
      //check whether queue is full 
      if (temp == m_readIndex)
      {
        ate = (temp + 1) % BUFF_SIZE;
        tt = m_theQueue[ate];
        //the cell after head is occupied
        if(tt != NULL0 && tt != NULL1) {
          
          return false; //queue full
        }
        //help the dequeue to update head
        CAS(&m_readIndex,temp,ate);      
        //try enqueue again
        continue;
      }

      //pointer address is 4-bytes aligned
      if(tt == NULL1)
        tnew = (ELEM *)((uintptr_t)a_data | 0x1);
      else
        tnew = a_data;
      //check the tail consistency
      if(te != m_writeIndex) continue;
      //get the actual tail and enqueue data
      if (CAS(&m_theQueue[ate],tt,tnew))
      {
        if(temp%2==0)
          CAS(&m_writeIndex,te,temp);
        return true;
      }          
    }while(1);
    
    return false;
}

template <typename ELEM,uint32_t BUFF_SIZE>
int NBLRUQueue<ELEM, BUFF_SIZE>::lru_dequeue(ELEM **a_data_p)
{
    uint32_t th,temp;
    ELEM *tt;
    ELEM *tnull;
    int i=0;
    do {
       i++;
      th = m_readIndex; //read the head
      //here is the one we want to dequeue
      temp = (th + 1)%BUFF_SIZE;
      tt = m_theQueue[temp];

      //find the actual head after this loop
      while(tt == NULL0 || tt == NULL1 )
      {
        //check the head's consistency
        if(th != m_readIndex) break;
        //two consecutive NULL means EMPTY return
        if(temp == m_writeIndex) return 1;
        //next call
        temp = (temp + 1)%BUFF_SIZE;   
        tt = m_theQueue[temp];
      }

      //check the head's consistency
      if(th != m_readIndex) continue;
      //check whether queue is empty
      if (temp == m_writeIndex)
      {
        //help the enqueue to update end
        CAS(&m_writeIndex,temp,(temp+1)%BUFF_SIZE);      
        //try dequeue again
        continue;
      }

      //pointer address is 4-bytes aligned
      if((uintptr_t)tt & 0x1)
        tnull = NULL1;
      else
        tnull = NULL0;
      
      //check the head's consistency
      if(th != m_readIndex) continue;
      //get the actual head, null value means empty
      if (CAS(&m_theQueue[temp],tt,tnull))
      {
        if(temp%2==0) CAS(&m_readIndex,th,temp);
        *a_data_p = (ELEM *)((uintptr_t)tt & ~0x1); //return the value
        return 0;
      } 
    }while(1);
    return 0;
}


template <typename ELEM,uint32_t BUFF_SIZE>
bool NBLRUQueue<ELEM, BUFF_SIZE>::lru_access(const ELEM *a_data)
{
  return lru_enqueue(a_data);
}

//this function is partly concurrent,means when a update is in processing,
//there may be enqueue operations on going,BUT no other dequeue can go 
template <typename ELEM, uint32_t BUFF_SIZE>
void NBLRUQueue<ELEM, BUFF_SIZE>::lru_update_lock()
{
  uint32_t th,temp,te;
  ELEM *tt;
  ELEM *tnull;
  
  th = m_readIndex; //read the head
  //here is the one we want to dequeue
  temp = (th + 1)%BUFF_SIZE;
  tt = m_theQueue[temp];

  if(tt == NULL0 || tt == NULL1 )
  {
    if(temp == m_writeIndex) return;
  }

  //check whether queue is empty
  if (temp == m_writeIndex)
  {
    //help the enqueue to update end
    CAS(&m_writeIndex,temp,(temp+1)%BUFF_SIZE);      
  }
  
  te = m_writeIndex;

  //loop 
  for(; temp != te;temp = (temp + 1)%BUFF_SIZE )
  {
    if(m_theQueue[temp] == NULL0 || m_theQueue[temp] == NULL1)
    {
      m_readIndex = temp;
      continue;      
    }
      
    tt = m_theQueue[temp];

    //pointer address is 4-bytes aligned
    if((uintptr_t)tt & 0x1)
      tnull = NULL1;
    else
      tnull = NULL0;
    m_theQueue[temp] = tnull;

    tt = (ELEM *)((uintptr_t)tt & ~0x1);
    //empty key,ignore,the bucket has been deleted
    if(tt->_key != 0)
    {
      ARRAYQ_REMOVE((&lru_head), tt, lru_link);
      ARRAYQ_INSERT_HEAD((&lru_head), tt, lru_link);
    }
    m_readIndex = temp;
  }
}

template <typename ELEM, uint32_t BUFF_SIZE>
void NBLRUQueue<ELEM, BUFF_SIZE>::lru_refresh_lock(ELEM *a_data)
{
    //empty key,ignore,the bucket has been deleted
    if(a_data->_key != 0)
    {
      ARRAYQ_REMOVE((&lru_head), a_data, lru_link);
      ARRAYQ_INSERT_HEAD((&lru_head), a_data, lru_link);
    }
}

template <typename ELEM, uint32_t BUFF_SIZE>
void NBLRUQueue<ELEM, BUFF_SIZE>::lru_remove_lock( ELEM *a_data,volatile int * g_limit)
{
  uint32_t th,temp,te;
  ELEM *tt;
  ELEM *tnull;
  
  th = m_readIndex; //read the head
  //here is the one we want to dequeue
  temp = (th + 1)%BUFF_SIZE;
  tt = m_theQueue[temp];
 
  te = m_writeIndex;

  //loop 
  for(; temp != te;temp = (temp + 1)%BUFF_SIZE )
  {
    if(m_theQueue[temp] == NULL0 || m_theQueue[temp] == NULL1)
      continue;
    
    tt = m_theQueue[temp];

    //pointer address is 4-bytes aligned
    if((uintptr_t)tt & 0x1)
      tnull = NULL1;
    else
      tnull = NULL0;
    
    tt = (ELEM *)((uintptr_t)tt & ~0x1);
    if(tt == a_data)
    {
      m_theQueue[temp] = tnull;
    }
  }
  ARRAYQ_REMOVE(&lru_head, a_data, lru_link);
  a_data->lru_link.tqe_next = _NULL_DELTA;
  a_data->lru_link.tqe_prev = _NULL_DELTA;
      
  lru_nodenum--;
  //assert(lru_nodenum >= 0);
  ATOMIC_ADD(g_limit,1);
}

//add new node to lru list,must under lock
template <typename ELEM, uint32_t BUFF_SIZE>
ELEM *NBLRUQueue<ELEM, BUFF_SIZE>::lru_new_lock(ELEM* const a_data,volatile int * g_limit)
{
  ELEM *t = NULL;

  //if no quota,free a node
  if(*g_limit <= 0)
  {    
    t = ARRAYQ_LAST(&lru_head, _lru);
  }

  ARRAYQ_INSERT_HEAD(&lru_head, a_data, lru_link);
  lru_nodenum++;
  ATOMIC_ADD(g_limit, -1);
  return t;
}

//add new node to lru list,must under lock
template <typename ELEM, uint32_t BUFF_SIZE>
ELEM *NBLRUQueue<ELEM, BUFF_SIZE>::lru_getoldest_lock()
{
  ELEM *t = NULL;

  t = ARRAYQ_LAST(&lru_head, _lru);

  return t;
}

////////////////////////////////////////////////////////////////////////////////
// CLASS: ConcurrentHopscotchHashMap
////////////////////////////////////////////////////////////////////////////////
typedef void *(*timercb)(unsigned int);

template <typename	_tKey, 
          typename	_tData,
			 typename	_tHash,
          typename	_tLock,
			 typename	_tMemory,
			 typename _tHazardPointer,
			 uint32_t BUFFER_SIZE>//actual buffer size is BUFFER_SIZE-2
class HopscotchHashMap {
private:
        // Inner Classes.........................
        struct Bucket {
		short				volatile _first_delta;
		short				volatile _next_delta;
		unsigned int	volatile _hash;
		_tKey				volatile _key;
		_tData			volatile _data;

   	    ARRAYQ_ENTRY(Bucket) lru_link;
        ARRAYQ_ENTRY(Bucket) age_link;//static age link,used for forced timeout
		ink_hrtime _access_time;

		void init() {
			_first_delta	= _NULL_DELTA;
			_next_delta		= _NULL_DELTA;
			_hash			= _tHash::_EMPTY_HASH;
			_key			= _tHash::_EMPTY_KEY;
			_data			= (_tData volatile)_tHash::_EMPTY_DATA;
			_access_time    = 0;
		}
        };

	struct Segment {
		unsigned int volatile	_timestamp;
        void     *evict_timer;
		_tLock	      _lock;
        ink_hrtime    _age_time; //oldest node's access time
        NBLRUQueue<Bucket,BUFFER_SIZE>  _lru;
        ARRAYQ_HEAD(_age, Bucket) _age_head;
    
		void init() {
		_timestamp = 0;
    	_lock.init();
		_lru.init();
        ARRAYQ_INIT(&_age_head);
        _age_time = HRTIME_FOREVER;
		}
	};

	// Fields ...................................................................
	_tHazardPointer * volatile _hazardPointer;
	unsigned int volatile		_segmentShift;
	unsigned int volatile		_segmentMask;
	unsigned int volatile		_bucketMask;
	unsigned int volatile       _segfreeMask;//free bucket border for each segment
	Segment*	volatile	_segments;
    Bucket* volatile	_table;

	freecb data_cb;
	freecb key_cb;
	int _maxNodeNum;
	int volatile _curNodeNum;
	int volatile _timeout;
	bool volatile _forceTimeout;//if true,regardless of node activity
	char * const _name;
	

	const int			_cache_mask;
	const bool			_is_cacheline_alignment;

	// Constants ................................................................
	static const unsigned int _INSERT_RANGE  = 1024*4;
	//static const unsigned int _NUM_SEGMENTS	= 1024*1;
	//static const unsigned int _SEGMENTS_MASK = _NUM_SEGMENTS-1;
	static const unsigned int _RESIZE_FACTOR = 2;

	// Small Utilities ..........................................................
	Bucket* get_start_cacheline_bucket(Bucket* const bucket) {
		return (bucket - ((bucket - _table) & _cache_mask)); //can optimize 
	}

	void remove_key(Segment&			  segment,
                   Bucket* const		  from_bucket,
						 Bucket* const		  key_bucket, 
						 Bucket* const		  prev_key_bucket, 
						 const unsigned int hash) 
	{
		key_bucket->_hash  = _tHash::_EMPTY_HASH;
		key_bucket->_key   = _tHash::_EMPTY_KEY;
		key_bucket->_data  = _tHash::_EMPTY_DATA;

		if(NULL == prev_key_bucket) {
			if (_NULL_DELTA == key_bucket->_next_delta)
				from_bucket->_first_delta = _NULL_DELTA;
			else 
				from_bucket->_first_delta = (from_bucket->_first_delta + key_bucket->_next_delta);
		} else {
			if (_NULL_DELTA == key_bucket->_next_delta)
				prev_key_bucket->_next_delta = _NULL_DELTA;
			else 
				prev_key_bucket->_next_delta = (prev_key_bucket->_next_delta + key_bucket->_next_delta);
		}
		key_bucket->_access_time = 0;

		++(segment._timestamp);
		key_bucket->_next_delta = _NULL_DELTA;
	}
	void add_key_to_begining_of_list(Bucket*	const     keys_bucket, 
										      Bucket*	const		 free_bucket,
												const unsigned int hash,
                                    const _tKey&		 key, 
                                    const _tData& 		 data) 
	{
		free_bucket->_data = data;
		free_bucket->_key  = key;
		free_bucket->_hash = hash;

		if(0 == keys_bucket->_first_delta) {
			if(_NULL_DELTA == keys_bucket->_next_delta)
				free_bucket->_next_delta = _NULL_DELTA;
			else
				free_bucket->_next_delta = (short)((keys_bucket +  keys_bucket->_next_delta) -  free_bucket);
			keys_bucket->_next_delta = (short)(free_bucket - keys_bucket);
		} else {
			if(_NULL_DELTA ==  keys_bucket->_first_delta)
				free_bucket->_next_delta = _NULL_DELTA;
			else
				free_bucket->_next_delta = (short)((keys_bucket +  keys_bucket->_first_delta) -  free_bucket);
			keys_bucket->_first_delta = (short)(free_bucket - keys_bucket);
		}
		free_bucket->_access_time = ink_get_hrtime();
	}

	void add_key_to_end_of_list(Bucket* const      keys_bucket, 
                               Bucket* const		  free_bucket,
                               const unsigned int hash,
                               const _tKey&		  key, 
										 const _tData&		  data,
                               Bucket* const		  last_bucket)
	{
		free_bucket->_data		 = data;
		free_bucket->_key			 = key;
		free_bucket->_hash		 = hash;
		free_bucket->_next_delta = _NULL_DELTA;

		if(NULL == last_bucket)
			keys_bucket->_first_delta = (short)(free_bucket - keys_bucket);
		else 
			last_bucket->_next_delta = (short)(free_bucket - last_bucket);
		free_bucket->_access_time = ink_get_hrtime();
	}

	void optimize_cacheline_use(Segment& segment, Bucket* const free_bucket,int hzd_id) {
		Bucket* const start_cacheline_bucket(get_start_cacheline_bucket(free_bucket));
		Bucket* const end_cacheline_bucket(start_cacheline_bucket + _cache_mask);
		Bucket* opt_bucket(start_cacheline_bucket);

		do {
			if( _NULL_DELTA != opt_bucket->_first_delta) {
				Bucket* relocate_key_last (NULL);
				int curr_delta(opt_bucket->_first_delta);
				Bucket* relocate_key ( opt_bucket + curr_delta);
				do {
					if( curr_delta < 0 || curr_delta > _cache_mask ) {
			            //use local temp to avoid to retire old node
			            volatile int temp=1;
			            //add this node to lru list
			            segment._lru.lru_new_lock(free_bucket,&temp);
                        ARRAYQ_INSERT_TAIL(&segment._age_head, free_bucket, age_link);
                         
						_tHash::relocate_data_reference(free_bucket->_data, relocate_key->_data);
						_tHash::relocate_key_reference(free_bucket->_key, relocate_key->_key);
						free_bucket->_hash  = relocate_key->_hash;
						free_bucket->_access_time = relocate_key->_access_time;

						if(_NULL_DELTA == relocate_key->_next_delta)
							free_bucket->_next_delta = _NULL_DELTA;
						else
							free_bucket->_next_delta = (short)( (relocate_key + relocate_key->_next_delta) - free_bucket );

						if(NULL == relocate_key_last)
							opt_bucket->_first_delta = (short)( free_bucket - opt_bucket );
						else
							relocate_key_last->_next_delta = (short)( free_bucket - relocate_key_last );

						++(segment._timestamp);
						relocate_key->_hash			= _tHash::_EMPTY_HASH;
						_tHash::relocate_key_reference(relocate_key->_key, _tHash::_EMPTY_KEY);
						_tHash::relocate_data_reference(relocate_key->_data, _tHash::_EMPTY_DATA);
			            //remove it from lru list after it's invisible from hash
			            //after new_lock and remove_lock,the node number should keep unchanged
                        ARRAYQ_REMOVE(&segment._age_head, relocate_key, age_link);
                        segment._lru.lru_remove_lock(relocate_key,&temp);
						relocate_key->_next_delta	= _NULL_DELTA;
                        if(_forceTimeout)
                        {
                           opt_bucket = ARRAYQ_FIRST(&segment._age_head);
                           segment._age_time = opt_bucket?opt_bucket->_access_time:HRTIME_FOREVER;
                        }
                        else
                        {
                           opt_bucket = segment._lru.lru_getoldest_lock();
                           segment._age_time =  opt_bucket?opt_bucket->_access_time:HRTIME_FOREVER;
                        }
						return;
					}

					if(_NULL_DELTA == relocate_key->_next_delta)
						break;
					relocate_key_last = relocate_key;
					curr_delta += relocate_key->_next_delta;
					relocate_key += relocate_key->_next_delta;
				} while(true);//for on list
			}
			++opt_bucket;
		} while (opt_bucket <= end_cacheline_bucket);
	}

public:// Ctors ................................................................
	HopscotchHashMap(char *name="",
		_tHazardPointer *hzp = NULL,
        freecb _kcb = NULL,
        freecb _dcb = NULL,
        int maxNodeNum = 32*1024,
				unsigned int inCapacity				= 32*1024,	//init capacity
				unsigned int concurrencyLevel	   = 16,			//num of updating threads
				unsigned int cache_line_size       = 64,			//Cache-line size of machine
				bool is_optimize_cacheline = true
				)		
	:	_cache_mask					( (cache_line_size / sizeof(Bucket)) - 1 ),
		_is_cacheline_alignment	( is_optimize_cacheline ),
		_segmentMask  ( NearestPowerOfTwo(concurrencyLevel) - 1),
		_segmentShift ( CalcDivideShift(NearestPowerOfTwo(concurrencyLevel/(NearestPowerOfTwo(concurrencyLevel)))-1) ),
    _curNodeNum(maxNodeNum),
    _maxNodeNum(maxNodeNum),
    _timeout(0),
    _name(name)
  {
		//ADJUST INPUT ............................
		const unsigned int adjInitCap = NearestPowerOfTwo(inCapacity);
		const unsigned int adjConCURRencyLevel = NearestPowerOfTwo(concurrencyLevel);
		const unsigned int num_buckets( adjInitCap + _INSERT_RANGE + 1);
        _bucketMask = adjInitCap - 1;
		_segmentShift = first_msb_bit_indx(_segmentMask) - first_msb_bit_indx(_bucketMask);
        _segfreeMask = (1<<_segmentShift)-1;
		_hazardPointer = hzp;
        data_cb = _dcb;
        key_cb = _kcb;

		//ALLOCATE THE SEGMENTS ...................
		_segments = (Segment*) _tMemory::byte_aligned_malloc( (_segmentMask + 1) * sizeof(Segment) );
		_table = (Bucket*) _tMemory::byte_aligned_malloc( num_buckets * sizeof(Bucket) );
		Segment* curr_seg = _segments;
		for (unsigned int iSeg = 0; iSeg <= _segmentMask; ++iSeg, ++curr_seg) {
			curr_seg->init();
		}

		Bucket* curr_bucket = _table;
		for (unsigned int iElm=0; iElm < num_buckets; ++iElm, ++curr_bucket) {
			curr_bucket->init();
		}
	}

	~HopscotchHashMap() {
		_tMemory::byte_aligned_free(_table);
		_tMemory::byte_aligned_free(_segments);
	} 

  // Dump segment info ........................................................
  void dumpInfo() {

     for(int i=0;i<=_segmentMask;i++)
     {
	   printf("seg%d: %x %x\n",i,_table + (i<<_segmentShift),_table + (i<<_segmentShift) + _segfreeMask ); 	
       printf("seg%d: %d,",i,_segments[i]._lru.lru_nodenum );
       if((i+1)%4 == 0)
        printf("\n");
     }   
     printf("%d free node in hash map.\n",_curNodeNum);
  }

  //force=true means don't care the activity of an item,evict it when time expires
  //force=false means inactive timeout,since the last access time to current time
  void setTimeOut(int timeout, bool force=false)
  {
  	_timeout = timeout;
	_forceTimeout = force;
  }

  int getTimeOut()
  {
    return _timeout;
  }
	
  // get data Operations .........................................................
  _tData get( const _tKey& key , int hzd_id=-1, ink_hrtime _timecmp=0) {

		//CALCULATE HASH ..........................
		const unsigned int hash( _tHash::Calc(key) );

		//CHECK IF ALREADY CONTAIN ................
		Segment&	segment(_segments[(hash >> _segmentShift) & _segmentMask]);

        //go over the list and look for key
		unsigned int start_timestamp;
		bool retry;
		bool validate = true;
        do {
			retry = false;
			start_timestamp = segment._timestamp;
			Bucket* curr_bucket( &(_table[hash & _bucketMask]) );
			short next_delta( curr_bucket->_first_delta );
			
	        while( _NULL_DELTA != next_delta ) {				
				curr_bucket += next_delta;
			    const _tKey k((_tKey&)(curr_bucket->_key));

				if(key_cb)
		        {
		          _hazardPointer->hp_set(hzd_id, (void *)k, 1);
				}
				
     			//key removed or updated,retry
				if(curr_bucket->_key != k)
				{
				  retry = true;
				  break;
				}
			
				if(hash == curr_bucket->_hash && _tHash::IsEqual(key, k))
		        {
				  const _tData rc((_tData&)(curr_bucket->_data));

				  if(data_cb)
		          {
		            _hazardPointer->hp_set(hzd_id, (void *)rc);
					
		          }

				  //this bucket is expired compared to _timecmp
				  if(curr_bucket->_access_time < _timecmp)
				  	validate = false;
				  if(!_forceTimeout && validate)
				    curr_bucket->_access_time = ink_get_hrtime();
		          
		          if(curr_bucket->_data != rc)
				  {
				  	retry = true;
				  	break;
				  }
		          //avoid ABA problem 
		          if(start_timestamp != segment._timestamp)
		          {
				  	retry = true;
				  	break;
				  }
		          //enqueue access to lru buffer
		          //here curr_bucket may be deleted  
		          if(curr_bucket->_key != _tHash::_EMPTY_KEY && !segment._lru.lru_access(curr_bucket) )
		          {
		            //queue is full,try lock
		            if(segment._lock.tryLock())
		            {
		              segment._lru.lru_update_lock();
		              segment._lru.lru_refresh_lock((Bucket*)curr_bucket);
		              segment._lock.unlock();
		            }
		          }
				  ink_hrtime c = ink_get_hrtime();		   
			 	  if(!validate || (_timeout && c - curr_bucket->_access_time >= ink_hrtime_from_sec(_timeout)))
			 	  {
				  	 //if this item is expired,don't continue to use it,this itme will 
				  	 //be retired in evict timer
				  	 _hazardPointer->hp_set(hzd_id, NULL);
					 _hazardPointer->hp_set(hzd_id, NULL, 1);
				  	 return _tHash::_EMPTY_DATA;
			 	  }		  
		          return rc;
		        }
				next_delta = curr_bucket->_next_delta;
			}
		} while(retry || start_timestamp != segment._timestamp);

		return _tHash::_EMPTY_DATA;
	}

	//modification Operations ...................................................
    //don't try to use the returned data,it is just to indicate key existance
	_tData putIfAbsent(const _tKey& key, const _tData& data,int hzd_id=-1, bool force=false) {
		const unsigned int hash( _tHash::Calc(key) );
        _tData retired = _tHash::_EMPTY_DATA;
		const unsigned int segnum((hash >> _segmentShift) & _segmentMask);
       
		Segment&	segment(_segments[segnum]);

		//go over the list and look for key
		if(force)
		    segment._lock.lock();
        else if(!segment._lock.tryLock())
        {
            //try lock failed,won't do that
            return data;
        }
            
		Bucket* const start_bucket( &(_table[hash & _bucketMask]) );

		Bucket* last_bucket( NULL );
		Bucket* compare_bucket( start_bucket );
		short next_delta( compare_bucket->_first_delta );
		while (_NULL_DELTA != next_delta) {
			compare_bucket += next_delta;
			if( hash == compare_bucket->_hash && _tHash::IsEqual(key, compare_bucket->_key) ) {
				const _tData rc((_tData&)(compare_bucket->_data));
				segment._lock.unlock();
				assert(rc != _tHash::_EMPTY_DATA);
				return rc;
			}
			last_bucket = compare_bucket;
			next_delta = compare_bucket->_next_delta;
		}

		//try to place the key in the same cache-line
		if(_is_cacheline_alignment) {
			Bucket*	free_bucket( start_bucket );
			Bucket*	start_cacheline_bucket(get_start_cacheline_bucket(start_bucket));
			Bucket*	end_cacheline_bucket(start_cacheline_bucket + _cache_mask);
			do {
				if( _tHash::_EMPTY_HASH == free_bucket->_hash ) {
		          segment._lru.lru_update_lock();
		          Bucket *rc = segment._lru.lru_new_lock(free_bucket,&_curNodeNum);
		          if(rc)
		          {
		            retired = remove_direct_lock(rc, segment,hzd_id);
		            assert(retired != _tHash::_EMPTY_DATA);
		          }
                  ARRAYQ_INSERT_TAIL(&segment._age_head, free_bucket, age_link);
				  add_key_to_begining_of_list(start_bucket, free_bucket, hash, key, data);  
                  if(free_bucket->_access_time < segment._age_time)
                     segment._age_time = free_bucket->_access_time;
				  segment._lock.unlock();
		          if(retired != _tHash::_EMPTY_DATA)
		          {
		            if(data_cb)
		            {
		              //_hazardPointer->hp_set(hzd_id, NULL); 
		              _hazardPointer->hp_free(hzd_id, (void *)retired, data_cb);
		            }
		          }
					return _tHash::_EMPTY_DATA;
				}
				++free_bucket;
				if(free_bucket > end_cacheline_bucket)
				  free_bucket = start_cacheline_bucket;
			} while(start_bucket != free_bucket);
		}

		//place key in arbitrary free forward bucket
		Bucket* max_bucket( start_bucket + (SHRT_MAX-1) );
		Bucket* last_table_bucket(_table + (segnum<<_segmentShift) + _segfreeMask );
		if(max_bucket > last_table_bucket)
			max_bucket = last_table_bucket;
		Bucket* free_max_bucket( start_bucket + (_cache_mask + 1) );
		while (free_max_bucket <= max_bucket) {
			if( _tHash::_EMPTY_HASH == free_max_bucket->_hash ) {
        segment._lru.lru_update_lock();
        Bucket *rc = segment._lru.lru_new_lock(free_max_bucket,&_curNodeNum);
        if(rc)
        {
          retired = remove_direct_lock(rc, segment,hzd_id);	
          assert(retired != _tHash::_EMPTY_DATA);
        }
        ARRAYQ_INSERT_TAIL(&segment._age_head, free_max_bucket, age_link);
		//update last_bucket,because we may have just removed the last bucket
		last_bucket = NULL;
        compare_bucket = start_bucket;
		next_delta = compare_bucket->_first_delta;
        while (_NULL_DELTA != next_delta) {
			    compare_bucket += next_delta;
			    last_bucket = compare_bucket;
			    next_delta = compare_bucket->_next_delta;
		}
		
        add_key_to_end_of_list(start_bucket, free_max_bucket, hash, key, data, last_bucket);
        if(free_max_bucket->_access_time < segment._age_time)
           segment._age_time = free_max_bucket->_access_time;
	    segment._lock.unlock();
        if(retired != _tHash::_EMPTY_DATA)
        {
          if(data_cb)
          {
            //_hazardPointer->hp_set(hzd_id, NULL); 
            _hazardPointer->hp_free(hzd_id, (void *)retired, data_cb);
          }
        }
				return _tHash::_EMPTY_DATA;
			}
			++free_max_bucket;
		}

		//place key in arbitrary free backward bucket
		Bucket* min_bucket( start_bucket - (SHRT_MAX-1) );
		if(min_bucket < _table + (segnum<<_segmentShift))
			min_bucket = _table + (segnum<<_segmentShift);
		Bucket* free_min_bucket( start_bucket - (_cache_mask + 1) );
		while (free_min_bucket >= min_bucket) {
			if( _tHash::_EMPTY_HASH == free_min_bucket->_hash ) {
		        segment._lru.lru_update_lock(); 
		        Bucket *rc = segment._lru.lru_new_lock(free_min_bucket,&_curNodeNum);
		        if(rc)
		        {
		          retired = remove_direct_lock(rc, segment,hzd_id);
		          assert(retired != _tHash::_EMPTY_DATA);
		        }
                ARRAYQ_INSERT_TAIL(&segment._age_head, free_min_bucket, age_link);
				//update last_bucket,because we may have just removed the last bucket
				last_bucket = NULL;
		        compare_bucket = start_bucket;
				next_delta = compare_bucket->_first_delta;
		        while (_NULL_DELTA != next_delta) {
					    compare_bucket += next_delta;
					    last_bucket = compare_bucket;
					    next_delta = compare_bucket->_next_delta;
				}
				add_key_to_end_of_list(start_bucket, free_min_bucket, hash, key, data, last_bucket);
                if(free_min_bucket->_access_time < segment._age_time)
                   segment._age_time = free_min_bucket->_access_time;
		        segment._lock.unlock();
		        if(retired != _tHash::_EMPTY_DATA)
		        {
		          if(data_cb)
		          {
		            //_hazardPointer->hp_set(hzd_id, NULL); 
		            _hazardPointer->hp_free(hzd_id, (void *)retired, data_cb);
		          }
		        }
				return _tHash::_EMPTY_DATA;
			}
			--free_min_bucket;
		}

		//NEED TO RESIZE ..........................
		//fprintf(stderr, "ERROR - RESIZE is not implemented - size %u\n", size());
		//exit(1);
		segment._lock.unlock();
	    if(data_cb) data_cb((void *)data);
		if(key_cb) key_cb((void *)key);
		return _tHash::_EMPTY_DATA;
	}

    //put if absent,update if existed
    //don't try to use the returned data,it is just to indicate key existance
	_tData putUpdate(const _tKey& key, const _tData& data,int hzd_id=-1, bool force=false) {
		const unsigned int hash( _tHash::Calc(key) );
        _tData retired = _tHash::_EMPTY_DATA;
	    const unsigned int segnum((hash >> _segmentShift) & _segmentMask);
		Segment&	segment(_segments[segnum]);

		//go over the list and look for key
		if(force)
		    segment._lock.lock();
        else if(!segment._lock.tryLock())
        {
            //try lock failed
            if(data_cb) data_cb((void *)data);
		    if(key_cb) key_cb((void *)key);
            return _tHash::_EMPTY_DATA;
        }
		Bucket* const start_bucket( &(_table[hash & _bucketMask]) );

		Bucket* last_bucket( NULL );
		Bucket* compare_bucket( start_bucket );
		short next_delta( compare_bucket->_first_delta );
		while (_NULL_DELTA != next_delta) {
			compare_bucket += next_delta;
			if( hash == compare_bucket->_hash && _tHash::IsEqual(key, compare_bucket->_key) ) {
				const _tData rc((_tData&)(compare_bucket->_data));
				const _tKey k((_tKey&)(compare_bucket->_key));
				//assert(compare_bucket->_data != data);
				compare_bucket->_key = key;
				compare_bucket->_data = data;
				
				//if(key_cb) key_cb(k);
			    segment._lru.lru_update_lock();
                segment._lru.lru_refresh_lock((Bucket*)compare_bucket);
				segment._lock.unlock();
                if(data_cb)
		          _hazardPointer->hp_free(hzd_id, (void *)rc, data_cb);		        
				if(key_cb)					
				  _hazardPointer->hp_free(hzd_id, (void *)k, key_cb);
				return data;
			}
			last_bucket = compare_bucket;
			next_delta = compare_bucket->_next_delta;
		}

		//try to place the key in the same cache-line
		if(_is_cacheline_alignment) {
			Bucket*	free_bucket( start_bucket );
			Bucket*	start_cacheline_bucket(get_start_cacheline_bucket(start_bucket));
			Bucket*	end_cacheline_bucket(start_cacheline_bucket + _cache_mask);
			do {
				if( _tHash::_EMPTY_HASH == free_bucket->_hash ) {
		          segment._lru.lru_update_lock();
		          Bucket *rc = segment._lru.lru_new_lock(free_bucket,&_curNodeNum);
		          if(rc)
		          {
		            retired = remove_direct_lock(rc, segment,hzd_id);
		            assert(retired != _tHash::_EMPTY_DATA);
		          }
                  ARRAYQ_INSERT_TAIL(&segment._age_head, free_bucket, age_link);
			      add_key_to_begining_of_list(start_bucket, free_bucket, hash, key, data);    
                  if(free_bucket->_access_time < segment._age_time)
                     segment._age_time = free_bucket->_access_time;
                  segment._lock.unlock();
		          if(retired != _tHash::_EMPTY_DATA)
		          {
		            if(data_cb)
		            {
		              //_hazardPointer->hp_set(hzd_id, NULL); 
		              _hazardPointer->hp_free(hzd_id, (void *)retired, data_cb);
		            }
		          }
				  return _tHash::_EMPTY_DATA;
				}
				++free_bucket;
				if(free_bucket > end_cacheline_bucket)
					free_bucket = start_cacheline_bucket;
			} while(start_bucket != free_bucket);
		}

		//place key in arbitrary free forward bucket
		Bucket* max_bucket( start_bucket + (SHRT_MAX-1) );
		Bucket* last_table_bucket(_table + (segnum<<_segmentShift) + _segfreeMask);
		if(max_bucket > last_table_bucket)
			max_bucket = last_table_bucket;
		Bucket* free_max_bucket( start_bucket + (_cache_mask + 1) );
		while (free_max_bucket <= max_bucket) {
			if( _tHash::_EMPTY_HASH == free_max_bucket->_hash ) {
        segment._lru.lru_update_lock();
        Bucket *rc = segment._lru.lru_new_lock(free_max_bucket,&_curNodeNum);
        if(rc)
        {
          retired = remove_direct_lock(rc, segment,hzd_id);	
          assert(retired != _tHash::_EMPTY_DATA);
        }
        ARRAYQ_INSERT_TAIL(&segment._age_head, free_max_bucket, age_link);
		//update last_bucket,because we may have just removed the last bucket
		last_bucket = NULL;
        compare_bucket = start_bucket;
		next_delta = compare_bucket->_first_delta;
        while (_NULL_DELTA != next_delta) {
			    compare_bucket += next_delta;
			    last_bucket = compare_bucket;
			    next_delta = compare_bucket->_next_delta;
		}
        add_key_to_end_of_list(start_bucket, free_max_bucket, hash, key, data, last_bucket);
        if(free_max_bucket->_access_time < segment._age_time)
          segment._age_time = free_max_bucket->_access_time;
        segment._lock.unlock();
        if(retired != _tHash::_EMPTY_DATA)
        {
          if(data_cb)
          {
            //_hazardPointer->hp_set(hzd_id, NULL); 
            _hazardPointer->hp_free(hzd_id, (void *)retired, data_cb);
          }
        }
				return _tHash::_EMPTY_DATA;
			}
			++free_max_bucket;
		}

		//place key in arbitrary free backward bucket
		Bucket* min_bucket( start_bucket - (SHRT_MAX-1) );
		if(min_bucket < _table + (segnum<<_segmentShift))
			min_bucket = _table + (segnum<<_segmentShift);
		Bucket* free_min_bucket( start_bucket - (_cache_mask + 1) );
		while (free_min_bucket >= min_bucket) {
			if( _tHash::_EMPTY_HASH == free_min_bucket->_hash ) {
        segment._lru.lru_update_lock(); 
        Bucket *rc = segment._lru.lru_new_lock(free_min_bucket,&_curNodeNum);
        if(rc)
        {
          retired = remove_direct_lock(rc, segment,hzd_id);
          assert(retired != _tHash::_EMPTY_DATA);
        }
        ARRAYQ_INSERT_TAIL(&segment._age_head, free_min_bucket, age_link);
		//update last_bucket,because we may have just removed the last bucket
		last_bucket = NULL;
        compare_bucket = start_bucket;
		next_delta = compare_bucket->_first_delta;
        while (_NULL_DELTA != next_delta) {
			    compare_bucket += next_delta;
			    last_bucket = compare_bucket;
			    next_delta = compare_bucket->_next_delta;
		}
		add_key_to_end_of_list(start_bucket, free_min_bucket, hash, key, data, last_bucket);
        if(free_min_bucket->_access_time < segment._age_time)
          segment._age_time = free_min_bucket->_access_time;
        segment._lock.unlock();
        if(retired != _tHash::_EMPTY_DATA)
        {
          if(data_cb)
          {
            //_hazardPointer->hp_set(hzd_id, NULL); 
            _hazardPointer->hp_free(hzd_id, (void *)retired, data_cb);
          }
        }
				return _tHash::_EMPTY_DATA;
			}
			--free_min_bucket;
		}

		//NEED TO RESIZE ..........................
		//fprintf(stderr, "ERROR - RESIZE is not implemented - size %u\n", size());
		//exit(1);
		segment._lock.unlock();
		if(data_cb) data_cb((void *)data);
		if(key_cb) key_cb((void *)key);
		return _tHash::_EMPTY_DATA;
	}
	

	bool remove(const _tKey& key, int hzd_id=-1) {
		//CALCULATE HASH ..........................
		const unsigned int hash( _tHash::Calc(key) );

		//CHECK IF ALREADY CONTAIN ................
		Segment&	segment(_segments[(hash >> _segmentShift) & _segmentMask]);
		segment._lock.lock();
		Bucket* const start_bucket( &(_table[hash & _bucketMask]) );
		Bucket* last_bucket( NULL );
		Bucket* curr_bucket( start_bucket );
		short	  next_delta (curr_bucket->_first_delta);
		do {
			if(_NULL_DELTA == next_delta) {
				segment._lock.unlock();
				return false;
			}
			curr_bucket += next_delta;

			if( hash == curr_bucket->_hash && _tHash::IsEqual(key, curr_bucket->_key) ) {
				_tData const rc((_tData&)(curr_bucket->_data));
                _tKey const k((_tKey&)(curr_bucket->_key));
				remove_key(segment, start_bucket, curr_bucket, last_bucket, hash);
                ARRAYQ_REMOVE(&segment._age_head, curr_bucket, age_link);
		        //if(key_cb) key_cb((void *)k);
		        segment._lru.lru_remove_lock(curr_bucket,&_curNodeNum);
		        segment._lru.lru_update_lock();
				if( _is_cacheline_alignment )
					optimize_cacheline_use(segment, curr_bucket, hzd_id);
                if(_forceTimeout)
                {
                   last_bucket = ARRAYQ_FIRST(&segment._age_head);
                   segment._age_time = last_bucket?last_bucket->_access_time:HRTIME_FOREVER;
                }
                else
                {
                   last_bucket = segment._lru.lru_getoldest_lock();
                   segment._age_time =  last_bucket?last_bucket->_access_time:HRTIME_FOREVER;
                }
				segment._lock.unlock();
		        //retire node
		        if(data_cb)       
		          _hazardPointer->hp_free(hzd_id, (void *)rc, data_cb);
				if(key_cb)
				  _hazardPointer->hp_free(hzd_id, (void *)k, key_cb);
        		return true;
			}
			last_bucket = curr_bucket;
			next_delta = curr_bucket->_next_delta;
		} while(true);
		return false;
	}

  //directly remove bucket by pointer,this function must be called under segment lock
  inline _tData remove_direct_lock(Bucket *bucket, Segment &segment,int hzd_id) {
    const unsigned int hash(bucket->_hash);
    Bucket* const start_bucket( &(_table[bucket->_hash & _bucketMask]) );
    Bucket* last_bucket( NULL );
    Bucket *curr_bucket (start_bucket);
    short	  next_delta (curr_bucket->_first_delta);

    do {
      if(_NULL_DELTA == next_delta) {     
	  	const unsigned int th( _tHash::Calc(bucket->_key) );
		assert(th == bucket->_hash);		
        assert(0);
		return _tHash::_EMPTY_DATA;
	  }
      curr_bucket += next_delta;
      if(curr_bucket == bucket)
      {
        _tData const rc((_tData&)(curr_bucket->_data));
        _tKey const k((_tKey&)(curr_bucket->_key));
        assert(_tHash::IsEqual(bucket->_key, curr_bucket->_key));
        remove_key(segment, start_bucket, bucket, last_bucket, hash);
        //if(key_cb) key_cb((void *)k);
        ARRAYQ_REMOVE(&segment._age_head, curr_bucket, age_link);
        segment._lru.lru_remove_lock(curr_bucket,&_curNodeNum);
        
        if(_forceTimeout)
        {
           last_bucket = ARRAYQ_FIRST(&segment._age_head);
           segment._age_time = last_bucket?last_bucket->_access_time:HRTIME_FOREVER;
        }
        else
        {
           last_bucket = segment._lru.lru_getoldest_lock();
           segment._age_time =  last_bucket?last_bucket->_access_time:HRTIME_FOREVER;
        }
        if( _is_cacheline_alignment )
		  optimize_cacheline_use(segment, curr_bucket, hzd_id);
		if(key_cb) 
			_hazardPointer->hp_free(hzd_id, (void *)k, key_cb); 
        return rc;
      }
      last_bucket = curr_bucket;
	  next_delta = curr_bucket->_next_delta;
    } while(true);
    
    return _tHash::_EMPTY_DATA;
  }

//#define DEBUG_EVICT    1
  //call this in your timer
  inline void evict_item(int timeout,unsigned int arg, int hzd_id) {
  #define MAX_EVICT 16
 
     _tData e_data[MAX_EVICT];
	 int idx=0;
	 ink_hrtime c;
	 Bucket *bucket;

     if(arg > _segmentMask)
        return;
	
     c = ink_get_hrtime();

     if(c < _segments[arg]._age_time + ink_hrtime_from_sec(timeout))
     {  
       return;
     }
     
     Segment& segment(_segments[arg]);

	 //try lock failed,this segment have to wait for a whole round complete 	 
	 segment._lock.lock();
     bool t = _forceTimeout;

	 while((bucket=(t?ARRAYQ_FIRST(&segment._age_head):segment._lru.lru_getoldest_lock())) != NULL)
	 {
     	if(c - bucket->_access_time < ink_hrtime_from_sec(timeout))
	      break;
 #ifdef DEBUG_EVICT      
        char title[64];
        sprintf(title,"Evict seg %d item:",arg);
        _tHash::dump_key(bucket->_key, _name, title);
 #endif
	    _tData rc =  remove_direct_lock(bucket, segment, hzd_id);

        assert(rc != _tHash::_EMPTY_DATA);	   
        if(rc != _tHash::_EMPTY_DATA)
        {
		  e_data[idx++] = rc;
		  if(idx >= MAX_EVICT)
		  	break;
        }
	 }
     segment._lock.unlock();
     
     if(data_cb)
     {
	   for(int i=0;i<idx;i++)
	   {
         //_hazardPointer->hp_set(hzd_id, NULL); 
         _hazardPointer->hp_free(hzd_id, (void *)e_data[i], data_cb);
	   }
     }
  }

  void setTimer(timercb fp)
  {
    for(unsigned int i=0; i<=_segmentMask; i++)
    {
      _segments[i].evict_timer = fp(i);
    }    
  }

	//status Operations .........................................................
	unsigned int size() {
		unsigned int counter = 0;
		const unsigned int num_elm( _bucketMask + _INSERT_RANGE );
		for(unsigned int iElm=0; iElm < num_elm; ++iElm) {
			if( _tHash::_EMPTY_HASH != _table[iElm]._hash ) {
				++counter;
			}
		}
		return counter;
	}   

	double percentKeysInCacheline() {
		unsigned int total_in_cache( 0 );
		unsigned int total( 0 );

		Bucket* curr_bucket(_table);
		for(int iElm(0); iElm <= _bucketMask; ++iElm, ++curr_bucket) {

			if(_NULL_DELTA != curr_bucket->_first_delta) {
				Bucket* const startCacheLineBucket( get_start_cacheline_bucket(curr_bucket) );
				Bucket* check_bucket(curr_bucket + curr_bucket->_first_delta);
				int currDist( curr_bucket->_first_delta );
				do {
					++total;
					if( (check_bucket - startCacheLineBucket)  >= 0 && (check_bucket - startCacheLineBucket) <= _cache_mask )
						++total_in_cache;
					if(_NULL_DELTA == check_bucket->_next_delta)
						break;
					currDist += check_bucket->_next_delta;
					check_bucket += check_bucket->_next_delta;
				} while(true);
			}
		}

		//return percent in cache
		return (((double)total_in_cache)/((double)total)*100.0);
	}

private:
	// Private Static Utilities .................................................
	static unsigned int NearestPowerOfTwo(const unsigned int value)	{
		unsigned int rc( 1 );
		while (rc < value) {
			rc <<= 1;
		}
		return rc;
	}

	static unsigned int CalcDivideShift(const unsigned int _value) {
		unsigned int numShift( 0 );
		unsigned int curr( 1 );
		while (curr < _value) {
			curr <<= 1;
			++numShift;
		}
		return numShift;
	}

};

typedef HAZARD_POINTER<256,2,32> WCG_HP;


#ifdef CCACHE_TEST

char *url_key[101] = {
"",
"www.google.com","www.youtube.com","www.yahoo.com","www.wikipedia.org","www.live.com","www.blogspot.com",
"www.amazon.com","www.linkedin.com","www.google.co.in","www.google.de","www.google.com.hk","www.google.co.jp",
"www.google.co.uk","www.bing.com","www.google.fr","www.microsoft.com","www.google.com.br","www.paypal.com",
"www.apple.com","www.google.it","www.tumblr.com","www.google.es","www.google.ru","www.blogger.com","www.fc2.com",
"www.bbc.co.uk","www.ask.com","www.google.com.mx","www.zedo.com","www.google.ca","www.vk.com","www.pinterest.com",
"www.mediafire.com","www.aol.com","www.conduit.com","www.adobe.com","www.google.co.id","www.rakuten.co.jp","www.imgur.com",
"www.wordpress.org","www.google.com.tr","www.godaddy.com","www.pornhub.com","www.google.com.au","www.amazon.de","www.4shared.com",
"www.google.pl","www.amazon.co.jp","www.netflix.com","www.nytimes.com","www.stackoverflow.com","www.alipay.com","www.chinaz.com",
"www.livejournal.com","www.amazon.co.uk","www.dailymotion.com","www.douban.com","www.google.com.sa","www.google.nl","www.addthis.com",
"www.adf.ly","www.google.com.ar","www.badoo.com","www.reddit.com","www.sparkstudios.com","www.bankofamerica.com",
"www.yieldmanager.com","www.vkontakte.ru","www.torrentz.eu","www.google.com.pk","www.aweber.com","www.deviantart.com","www.chase.com",
"www.google.co.th","www.stumbleupon.com","www.spiegel.de","www.rapidshare.com","www.optmd.com","www.blogspot.in","www.google.com.eg",
"www.mozilla.org","www.sourceforge.net","www.myspace.com","www.google.co.za","www.skype.com","www.etsy.com","www.imageshack.us",
"www.reference.com","www.58.com","www.wellsfargo.com","www.taringa.net","www.rediff.com","www.comcast.net","www.liveinternet.ru",
"www.wikimedia.org","www.wikia.com","www.yelp.com","www.google.com.my","www.google.co.ve","www.depositfiles.com"
};

void * mpool[5000];
TTASLock mlock;

void *myallocate(int size)
{   
  mlock.lock();
  for(int i=0;i<5000;i++)
  {
    if(mpool[i] == NULL)
    {
      mpool[i] = malloc(size);
      mlock.unlock();
      return mpool[i];
    }
  }
  mlock.unlock();
  return NULL;
}


void myfree(void *xx)
{
  bool find=false;
  mlock.lock();
  for(int i=0;i<5000;i++)
  {
    if(mpool[i] == xx)
    {
      find = true;
      mpool[i] = NULL;
      break;
    }
  }
  mlock.unlock();
  //printf("freed %x %d\n",xx,find?1:0);
  free(xx); 
}

typedef HAZARD_POINTER<32,2,8> LOCAL_HP;
LOCAL_HP hzp;

HopscotchHashMap<int,char *,HASH_INT<char *>, TTASLock, Memory, LOCAL_HP, 3> hhm("",&hzp,NULL,myfree,32,256);

HopscotchHashMap<char *,char *,HASH_STR<char *>, TTASLock, Memory, LOCAL_HP, 3> hhs("",&hzp,free,myfree,32,256);


pthread_t ntid[16];


struct thread_struct {
  int hp_id;  
};

pthread_key_t my_key;


void *thr_fn(void *arg)
{
  thread_struct* p = NULL;
  int loop = 1000;
  int type = *((int *)arg);

  p = (thread_struct *)pthread_getspecific(my_key);

  if(p == NULL)
  {
    p = (thread_struct *)malloc(sizeof(thread_struct));
    if(p == NULL)
    {
      return NULL;
    }
    p->hp_id = hzp.hp_register(1);
    pthread_setspecific(my_key, p);
  }

  sleep(1);
  
  do {    
    int key = rand()%100+1;
    //int data = rand()%1000+1;
    char *data = (char *)myallocate(32);

    sprintf(data,"%d's string %d",p->hp_id,key);

    if(type)
    {
   
      if( hhs.putUpdate(strdup(url_key[key]),data,p->hp_id) == 0)
      {
         printf("%d put %s ok.\n",p->hp_id,url_key[key]);
      }
      else 
      {
         printf("%d replace %s ok.\n",p->hp_id,url_key[key]);
      } 
    }
    else 
    {
      if( hhm.putUpdate(key,data,p->hp_id) == 0)
      {
         printf("%d put %d ok.\n",p->hp_id,key);
      }
      else 
      {
         printf("%d replace %d ok.\n",p->hp_id,key);
      } 
    }

    for(int i=1;i<=100;i++)
    {
       if(type)
         data = hhs.get(url_key[i],p->hp_id);
       else
         data = hhm.get(i,p->hp_id);
      if(data != 0)
      {
        usleep(3);  
        printf("%d get %s %s ok.\n",p->hp_id,url_key[i],data); 
        break;
      }  
    }

#if 0
    key = rand()%100 + 1;

    if(type)
    {
      if(hhs.remove(url_key[key], p->hp_id))
        printf("%d remove %s ok.\n",p->hp_id,url_key[key]);
    }
    else
    {
      if(hhm.remove(key, p->hp_id))
        printf("%d remove %d ok.\n",p->hp_id,key);
    }
#endif    
   
    usleep(120);
  }while(loop--);


  for(int i=1;i<=100;i++)
  {
     if(type)
       hhs.remove(url_key[i], p->hp_id); 
     else
       hhm.remove(i, p->hp_id);
  }
  hzp.hp_set(p->hp_id, NULL);
  hzp.hp_set(p->hp_id, NULL, 1);
  
  hzp.hp_unregister(p->hp_id);
}

//when we destroy a thread,we must call hp_set(NULL) to releae hp
//then retire all nodes
int main()
{
    int type = 1;
    time_t t;
    srand((unsigned)time(&t));
    mlock.init();

    int nRet = pthread_key_create(&my_key, NULL);

    for(int i=0;i<16;i++)
    {
      if ((nRet = pthread_create(&ntid[i], NULL, thr_fn, (void *)&type)) != 0)
      {
          return -1;
      }      
    }
    
    for(int i=0;i<16;i++)
      pthread_join(ntid[i], NULL);

    printf("check memory leak.....\n");
    for(int i=0;i<5000;i++)
    {
      if(mpool[i] != NULL)
      {
        printf("%x not freeed.\n",mpool[i]);
      }
    }      
}
#endif
#endif

// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"
using namespace std;

// Needed to make frequent large allocations efficient with standard
// malloc implementation.  Otherwise they are allocated directly from
// vm.
#include <malloc.h>
static int __ii =  mallopt(M_MMAP_MAX,0);
static int __jj =  mallopt(M_TRIM_THRESHOLD,-1);

#define newA(__E,__n) (__E*) malloc((__n)*sizeof(__E))

template <class E>
struct identityF { E operator() (const E& x) {return x;}};

template <class E>
struct addF { E operator() (const E& a, const E& b) const {return a+b;}};

template <class E>
struct minF { E operator() (const E& a, const E& b) const {return (a < b) ? a : b;}};

template <class E>
struct maxF { E operator() (const E& a, const E& b) const {return (a>b) ? a : b;}};

#define _SCAN_LOG_BSIZE 10
#define _SCAN_BSIZE (1 << _SCAN_LOG_BSIZE)

template <class T>
struct _seq {
    T* A;
    long n;
    _seq() {A = NULL; n=0;}
    _seq(T* _A, long _n) : A(_A), n(_n) {}
    void del() {free(A);}
};

namespace sequence {
    template <class IntT>
    struct boolGetA {
        bool* A;
        boolGetA(bool* AA) : A(AA) {}
        IntT operator() (IntT i) {return (IntT) A[i];}
    };

    template <class ET, class IntT>
    struct getA {
        ET* A;
        getA(ET* AA) : A(AA) {}
        ET operator() (IntT i) {return A[i];}
    };

    template <class IT, class OT, class IntT, class F>
    struct getAF {
        IT* A;
        F f;
        getAF(IT* AA, F ff) : A(AA), f(ff) {}
        OT operator () (IntT i) {return f(A[i]);}
    };

#define nblocks(_n,_bsize) (1 + ((_n)-1)/(_bsize))

#define blocked_for(_i, _s, _e, _bsize, _body)  {	\
    long _ss = _s;					\
    long _ee = _e;					\
    long _n = _ee-_ss;					\
    long _l = nblocks(_n,_bsize);			\
    parallel_for (long _i = 0; _i < _l; _i++) {		\
      long _s = _ss + _i * (_bsize);			\
      long _e = min(_s + (_bsize), _ee);		\
      _body						\
	}						\
  }

    template <class OT, class IntT, class F, class G>
    OT reduceSerial(IntT s, IntT e, F f, G g) {
        OT r = g(s);
        for (IntT j=s+1; j < e; j++) r = f(r,g(j));
        return r;
    }

    template <class OT, class IntT, class F, class G>
    OT reduce(IntT s, IntT e, F f, G g) {
        IntT l = nblocks(e-s, _SCAN_BSIZE);
        if (l <= 1) return reduceSerial<OT>(s, e, f , g);
        OT *Sums = newA(OT,l);
        blocked_for (i, s, e, _SCAN_BSIZE,
                     Sums[i] = reduceSerial<OT>(s, e, f, g););
        OT r = reduce<OT>((IntT) 0, l, f, getA<OT,IntT>(Sums));
        free(Sums);
        return r;
    }

    template <class OT, class IntT, class F>
    OT reduce(OT* A, IntT n, F f) {
        return reduce<OT>((IntT)0,n,f,getA<OT,IntT>(A));
    }

    template <class OT, class IntT>
    OT plusReduce(OT* A, IntT n) {
        return reduce<OT>((IntT)0,n,addF<OT>(),getA<OT,IntT>(A));
    }

    // g is the map function (applied to each element)
    // f is the reduce function
    // need to specify OT since it is not an argument
    template <class OT, class IT, class IntT, class F, class G>
    OT mapReduce(IT* A, IntT n, F f, G g) {
        return reduce<OT>((IntT) 0,n,f,getAF<IT,OT,IntT,G>(A,g));
    }

    template <class IntT>
    long sum(bool *In, IntT n) {
        return reduce<IntT>((IntT) 0, n, addF<IntT>(), boolGetA<IntT>(In));
    }

    template <class ET, class IntT, class F, class G>
    ET scanSerial(ET* Out, IntT s, IntT e, F f, G g, ET zero, bool inclusive, bool back) {
        ET r = zero;
        if (inclusive) {
            if (back) for (IntT i = e-1; i >= s; i--) Out[i] = r = f(r,g(i));
            else for (IntT i = s; i < e; i++) Out[i] = r = f(r,g(i));
        } else {
            if (back)
                for (IntT i = e-1; i >= s; i--) {
                    ET t = g(i);
                    Out[i] = r;
                    r = f(r,t);
                }
            else
                for (IntT i = s; i < e; i++) {
                    ET t = g(i);
                    Out[i] = r;
                    r = f(r,t);
                }
        }
        return r;
    }

    template <class ET, class IntT, class F>
    ET scanSerial(ET *In, ET* Out, IntT n, F f, ET zero) {
        return scanSerial(Out, (IntT) 0, n, f, getA<ET,IntT>(In), zero, false, false);
    }

    // back indicates it runs in reverse direction
    template <class ET, class IntT, class F, class G>
    ET scan(ET* Out, IntT s, IntT e, F f, G g,  ET zero, bool inclusive, bool back) {
        IntT n = e-s;
        IntT l = nblocks(n,_SCAN_BSIZE);
        if (l <= 2) return scanSerial(Out, s, e, f, g, zero, inclusive, back);
        ET *Sums = newA(ET,nblocks(n,_SCAN_BSIZE));
        blocked_for (i, s, e, _SCAN_BSIZE,
                     Sums[i] = reduceSerial<ET>(s, e, f, g););
        ET total = scan(Sums, (IntT) 0, l, f, getA<ET,IntT>(Sums), zero, false, back);
        blocked_for (i, s, e, _SCAN_BSIZE,
                     scanSerial(Out, s, e, f, g, Sums[i], inclusive, back););
        free(Sums);
        return total;
    }

    template <class ET, class IntT, class F>
    ET scan(ET *In, ET* Out, IntT n, F f, ET zero) {
        return scan(Out, (IntT) 0, n, f, getA<ET,IntT>(In), zero, false, false);}

    template <class ET, class IntT, class F>
    ET scanI(ET *In, ET* Out, IntT n, F f, ET zero) {
        return scan(Out, (IntT) 0, n, f, getA<ET,IntT>(In), zero, true, false);}

    template <class ET, class IntT, class F>
    ET scanBack(ET *In, ET* Out, IntT n, F f, ET zero) {
        return scan(Out, (IntT) 0, n, f, getA<ET,IntT>(In), zero, false, true);}

    template <class ET, class IntT, class F>
    ET scanIBack(ET *In, ET* Out, IntT n, F f, ET zero) {
        return scan(Out, (IntT) 0, n, f, getA<ET,IntT>(In), zero, true, true);}

    template <class ET, class IntT>
    ET plusScan(ET *In, ET* Out, IntT n) {
        return scan(Out, (IntT) 0, n, addF<ET>(), getA<ET,IntT>(In),
                    (ET) 0, false, false);}

#define _F_BSIZE (2*_SCAN_BSIZE)

    // sums a sequence of n boolean flags
    // an optimized version that sums blocks of 4 booleans by treating
    // them as an integer
    // Only optimized when n is a multiple of 512 and Fl is 4byte aligned
    template <class IntT>
    IntT sumFlagsSerial(bool *Fl, IntT n) {
        IntT r = 0;
        if (n >= 128 && (n & 511) == 0 && ((IntT) Fl & 3) == 0) {
            int* IFl = (int*) Fl;
            for (int k = 0; k < (n >> 9); k++) {
                int rr = 0;
                for (int j=0; j < 128; j++) rr += IFl[j];
                r += (rr&255) + ((rr>>8)&255) + ((rr>>16)&255) + ((rr>>24)&255);
                IFl += 128;
            }
        } else for (IntT j=0; j < n; j++) r += Fl[j];
        return r;
    }

    template <class ET, class IntT, class F>
    _seq<ET> packSerial(ET* Out, bool* Fl, IntT s, IntT e, F f) {
        if (Out == NULL) {
            IntT m = sumFlagsSerial(Fl+s, e-s);
            Out = newA(ET,m);
        }
        IntT k = 0;
        for (IntT i=s; i < e; i++) if (Fl[i]) Out[k++] = f(i);
        return _seq<ET>(Out,k);
    }

    template <class ET, class IntT, class F>
    _seq<ET> pack(ET* Out, bool* Fl, IntT s, IntT e, F f) {
        IntT l = nblocks(e-s, _F_BSIZE);
        if (l <= 1) return packSerial(Out, Fl, s, e, f);
        IntT *Sums = newA(IntT,l);
        blocked_for (i, s, e, _F_BSIZE, Sums[i] = sumFlagsSerial(Fl+s, e-s););
        IntT m = plusScan(Sums, Sums, l);
        if (Out == NULL) Out = newA(ET,m);
        blocked_for(i, s, e, _F_BSIZE, packSerial(Out+Sums[i], Fl, s, e, f););
        free(Sums);
        return _seq<ET>(Out,m);
    }

    template <class ET, class IntT>
    IntT pack(ET* In, ET* Out, bool* Fl, IntT n) {
        return pack(Out, Fl, (IntT) 0, n, getA<ET,IntT>(In)).n;}

    template <class IntT>
    _seq<IntT> packIndex(bool* Fl, IntT n) {
        return pack((IntT *) NULL, Fl, (IntT) 0, n, identityF<IntT>());
    }

    template <class ET, class IntT, class PRED>
    IntT filter(ET* In, ET* Out, IntT n, PRED p) {
        bool *Fl = newA(bool,n);
        parallel_for (IntT i=0; i < n; i++) Fl[i] = (bool) p(In[i]);
        IntT  m = pack(In, Out, Fl, n);
        free(Fl);
        return m;
    }
}

template <class ET>
inline bool CAS(ET *ptr, ET oldv, ET newv) {
    if (sizeof(ET) == 4) {
        return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv), *((int*)&newv));
    } else if (sizeof(ET) == 8) {
        return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv), *((long*)&newv));
    }
    else {
        std::cout << "CAS bad length" << std::endl;
        abort();
    }
}

template <class ET>
inline bool writeMin(ET *a, ET b) {
    ET c; bool r=0;
    do c = *a;
    while (c > b && !(r=CAS(a,c,b)));
    return r;
}

template <class ET>
inline void writeAdd(ET *a, ET b) {
    volatile ET newV, oldV;
    do {oldV = *a; newV = oldV + b;}
    while (!CAS(a, oldV, newV));
}

inline uint hash(uint a) {
    a = (a+0x7ed55d16) + (a<<12);
    a = (a^0xc761c23c) ^ (a>>19);
    a = (a+0x165667b1) + (a<<5);
    a = (a+0xd3a2646c) ^ (a<<9);
    a = (a+0xfd7046c5) + (a<<3);
    a = (a^0xb55a4f09) ^ (a>>16);
    return a;
}

inline ulong hash(ulong a) {
    a = (a+0x7ed55d166bef7a1d) + (a<<12);
    a = (a^0xc761c23c510fa2dd) ^ (a>>9);
    a = (a+0x165667b183a9c0e1) + (a<<59);
    a = (a+0xd3a2646cab3487e3) ^ (a<<49);
    a = (a+0xfd7046c5ef9ab54c) + (a<<3);
    a = (a^0xb55a4f090dd4a67b) ^ (a>>32);
    return a;
}

#endif

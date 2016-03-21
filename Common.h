//
// Created by alvis on 16.03.16.
//

#ifndef NETWORKFLOW_COMMON_H
#define NETWORKFLOW_COMMON_H

#include <limits.h>
#include <assert.h>
#include "nheap.h"

#if defined(LONG)
typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#define UINT_T_MAX ULONG_MAX
#else
typedef int intT;
typedef unsigned int uintT;
#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX
#endif

//edges store 32-bit quantities unless EDGELONG is defined
#if defined(EDGELONG)
typedef long intE;
typedef unsigned long uintE;
#define INT_E_MAX LONG_MAX
#define UINT_E_MAX ULONG_MAX
#else
typedef int intE;
typedef unsigned int uintE;
#define INT_E_MAX INT_MAX
#define UINT_E_MAX UINT_MAX
#endif

#define BASH_RED "\033[0;31m"
#define BASH_GREEN "\033[0;32m"
#define BASH_YELLOW "\033[1;33m"
#define BASH_NC "\033[0m"

typedef fHeap<intT> mmHeap;



/// \brief Problem type constants for the \c run() function.
///
/// Enum type containing the problem type constants that can be
/// returned by the \ref run() function of the algorithm.
enum ProblemType {
    /// The problem has no feasible solution (flow).
            INFEASIBLE,
    /// The problem has optimal solution (i.e. it is feasible and
    /// bounded), and the algorithm has found optimal flow and node
    /// potentials (primal and dual solutions).
            OPTIMAL,
    /// The digraph contains an arc of negative cost and infinite
    /// upper bound. It means that the objective function is unbounded
    /// on that arc, however, note that it could actually be bounded
    /// over the feasible flows, but this algroithm cannot handle
    /// these cases.
            UNBOUNDED
};

#endif //NETWORKFLOW_COMMON_H

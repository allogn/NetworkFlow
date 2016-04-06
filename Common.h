//
// Created by alvis on 16.03.16.
//

#ifndef NETWORKFLOW_COMMON_H
#define NETWORKFLOW_COMMON_H

#include <limits.h>
#include <assert.h>
#include "nheap.h"

#define BASH_RED "\033[0;31m"
#define BASH_GREEN "\033[0;32m"
#define BASH_YELLOW "\033[1;33m"
#define BASH_NC "\033[0m"

#define newA(__E,__n) (__E*) malloc((__n)*sizeof(__E))

typedef fHeap<long long> mmHeap;

enum Algorithm {
    ALG_SIA,
    ALG_COST_SCALING,
    ALG_LOCAL_DOMINANT,
    ALG_LEMON_MODIF,
    ALG_SCS,
    ALG_LEMON_ORIG,
    ALG_LSCS
};

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

inline bool isSpace(char c) {
    switch (c)  {
        case '\r':
        case '\t':
        case '\n':
        case 0:
        case ' ' : return true;
        default : return false;
    }
}

#endif //NETWORKFLOW_COMMON_H

//
// Created by alvis on 20.03.16.
//

#ifndef NETWORKFLOW_LEMON_H
#define NETWORKFLOW_LEMON_H


/// \ingroup min_cost_flow_algs
/// \file
/// \brief Cost scaling algorithm for finding a minimum cost flow.

#include <vector>
#include <deque>
#include <limits>

#include <lemon/core.h>
#include <lemon/maps.h>
#include <lemon/math.h>
#include <lemon/static_graph.h>
#include <lemon/circulation.h>
#include <lemon/bellman_ford.h>

namespace lemon {

    /// \brief Default traits class of ModifiedCostScaling algorithm.
    ///
    /// Default traits class of ModifiedCostScaling algorithm.
    /// \tparam GR Digraph type.
    /// \tparam V The number type used for flow amounts, capacity bounds
    /// and supply values. By default it is \c int.
    /// \tparam C The number type used for costs and potentials.
    /// By default it is the same as \c V.
#ifdef DOXYGEN
    template <typename GR, typename V = int, typename C = V>
#else
    template < typename GR, typename V = int, typename C = V,
            bool integer = std::numeric_limits<C>::is_integer >
#endif
    struct ModifiedCostScalingDefaultTraits
    {
        /// The type of the digraph
        typedef GR Digraph;
        /// The type of the flow amounts, capacity bounds and supply values
        typedef V Value;
        /// The type of the arc costs
        typedef C Cost;

        /// \brief The large cost type used for internal computations
        ///
        /// The large cost type used for internal computations.
        /// It is \c long \c long if the \c Cost type is integer,
        /// otherwise it is \c double.
        /// \c Cost must be convertible to \c LargeCost.
        typedef double LargeCost;
    };

    // Default traits class for integer cost types
    template <typename GR, typename V, typename C>
    struct ModifiedCostScalingDefaultTraits<GR, V, C, true>
    {
        typedef GR Digraph;
        typedef V Value;
        typedef C Cost;
#ifdef LEMON_HAVE_LONG_LONG
        typedef long long LargeCost;
#else
        typedef long LargeCost;
#endif
    };


    /// \addtogroup min_cost_flow_algs
    /// @{

    /// \brief Implementation of the Cost Scaling algorithm for
    /// finding a \ref min_cost_flow "minimum cost flow".
    ///
    /// \ref ModifiedCostScaling implements a cost scaling algorithm that performs
    /// push/augment and relabel operations for finding a \ref min_cost_flow
    /// "minimum cost flow" \cite amo93networkflows,
    /// \cite goldberg90approximation,
    /// \cite goldberg97efficient, \cite bunnagel98efficient.
    /// It is a highly efficient primal-dual solution method, which
    /// can be viewed as the generalization of the \ref Preflow
    /// "preflow push-relabel" algorithm for the maximum flow problem.
    /// It is a polynomial algorithm, its running time complexity is
    /// \f$O(n^2m\log(nK))\f$, where <i>K</i> denotes the maximum arc cost.
    ///
    /// In general, \ref NetworkSimplex and \ref ModifiedCostScaling are the fastest
    /// implementations available in LEMON for solving this problem.
    /// (For more information, see \ref min_cost_flow_algs "the module page".)
    ///
    /// Most of the parameters of the problem (except for the digraph)
    /// can be given using separate functions, and the algorithm can be
    /// executed using the \ref run() function. If some parameters are not
    /// specified, then default values will be used.
    ///
    /// \tparam GR The digraph type the algorithm runs on.
    /// \tparam V The number type used for flow amounts, capacity bounds
    /// and supply values in the algorithm. By default, it is \c int.
    /// \tparam C The number type used for costs and potentials in the
    /// algorithm. By default, it is the same as \c V.
    /// \tparam TR The traits class that defines various types used by the
    /// algorithm. By default, it is \ref ModifiedCostScalingDefaultTraits
    /// "ModifiedCostScalingDefaultTraits<GR, V, C>".
    /// In most cases, this parameter should not be set directly,
    /// consider to use the named template parameters instead.
    ///
    /// \warning Both \c V and \c C must be signed number types.
    /// \warning All input data (capacities, supply values, and costs) must
    /// be integer.
    /// \warning This algorithm does not support negative costs for
    /// arcs having infinite upper bound.
    ///
    /// \note %ModifiedCostScaling provides three different internal methods,
    /// from which the most efficient one is used by default.
    /// For more information, see \ref Method.
#ifdef DOXYGEN
    template <typename GR, typename V, typename C, typename TR>
#else
    template < typename GR, typename V = int, typename C = V,
            typename TR = ModifiedCostScalingDefaultTraits<GR, V, C> >
#endif
    class ModifiedCostScaling
    {
    public:

        /// The type of the digraph
        typedef typename TR::Digraph Digraph;
        /// The type of the flow amounts, capacity bounds and supply values
        typedef typename TR::Value Value;
        /// The type of the arc costs
        typedef typename TR::Cost Cost;

        /// \brief The large cost type
        ///
        /// The large cost type used for internal computations.
        /// By default, it is \c long \c long if the \c Cost type is integer,
        /// otherwise it is \c double.
        typedef typename TR::LargeCost LargeCost;

        /// \brief The \ref lemon::ModifiedCostScalingDefaultTraits "traits class"
        /// of the algorithm
        typedef TR Traits;

    public:

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

        /// \brief Constants for selecting the internal method.
        ///
        /// Enum type containing constants for selecting the internal method
        /// for the \ref run() function.
        ///
        /// \ref ModifiedCostScaling provides three internal methods that differ mainly
        /// in their base operations, which are used in conjunction with the
        /// relabel operation.
        /// By default, the so called \ref PARTIAL_AUGMENT
        /// "Partial Augment-Relabel" method is used, which turned out to be
        /// the most efficient and the most robust on various test inputs.
        /// However, the other methods can be selected using the \ref run()
        /// function with the proper parameter.
        enum Method {
            /// Local push operations are used, i.e. flow is moved only on one
            /// admissible arc at once.
                    PUSH,
            /// Augment operations are used, i.e. flow is moved on admissible
            /// paths from a node with excess to a node with deficit.
                    AUGMENT,
            /// Partial augment operations are used, i.e. flow is moved on
            /// admissible paths started from a node with excess, but the
            /// lengths of these paths are limited. This method can be viewed
            /// as a combined version of the previous two operations.
                    PARTIAL_AUGMENT
        };

    private:

        TEMPLATE_DIGRAPH_TYPEDEFS(GR);

        typedef std::vector<int> IntVector;
        typedef std::vector<Value> ValueVector;
        typedef std::vector<Cost> CostVector;
        typedef std::vector<LargeCost> LargeCostVector;
        typedef std::vector<char> BoolVector;
        // Note: vector<char> is used instead of vector<bool>
        // for efficiency reasons

    private:

        template <typename KT, typename VT>
        class StaticVectorMap {
        public:
            typedef KT Key;
            typedef VT Value;

            StaticVectorMap(std::vector<Value>& v) : _v(v) {}

            const Value& operator[](const Key& key) const {
                return _v[StaticDigraph::id(key)];
            }

            Value& operator[](const Key& key) {
                return _v[StaticDigraph::id(key)];
            }

            void set(const Key& key, const Value& val) {
                _v[StaticDigraph::id(key)] = val;
            }

        private:
            std::vector<Value>& _v;
        };

        typedef StaticVectorMap<StaticDigraph::Arc, LargeCost> LargeCostArcMap;

    private:

        // Data related to the underlying digraph
        const GR &_graph;
        int _node_num;
        int _arc_num;
        int _res_node_num;
        int _res_arc_num;
        int _root;

        // Parameters of the problem
        bool _has_lower;
        Value _sum_supply;
        int _sup_node_num;

        // Data structures for storing the digraph
        IntNodeMap _node_id;
        IntArcMap _arc_idf;
        IntArcMap _arc_idb;
        IntVector _first_out;
        BoolVector _forward;
        IntVector _source;
        IntVector _target;
        IntVector _reverse;

        // Node and arc data
        ValueVector _lower;
        ValueVector _upper;
        CostVector _scost;
        ValueVector _supply;

        ValueVector _res_cap;
        LargeCostVector _cost;
        LargeCostVector _pi;
        ValueVector _excess;
        IntVector _next_out;
        std::deque<int> _active_nodes;

        // Data for scaling
        LargeCost _epsilon;
        int _alpha;

        IntVector _buckets;
        IntVector _bucket_next;
        IntVector _bucket_prev;
        IntVector _rank;
        int _max_rank;

    public:

        /// \brief Constant for infinite upper bounds (capacities).
        ///
        /// Constant for infinite upper bounds (capacities).
        /// It is \c std::numeric_limits<Value>::infinity() if available,
        /// \c std::numeric_limits<Value>::max() otherwise.
        const Value INF;

    public:

        /// \name Named Template Parameters
        /// @{

        template <typename T>
        struct SetLargeCostTraits : public Traits {
            typedef T LargeCost;
        };

        /// \brief \ref named-templ-param "Named parameter" for setting
        /// \c LargeCost type.
        ///
        /// \ref named-templ-param "Named parameter" for setting \c LargeCost
        /// type, which is used for internal computations in the algorithm.
        /// \c Cost must be convertible to \c LargeCost.
        template <typename T>
        struct SetLargeCost
                : public ModifiedCostScaling<GR, V, C, SetLargeCostTraits<T> > {
            typedef  ModifiedCostScaling<GR, V, C, SetLargeCostTraits<T> > Create;
        };

        /// @}

    protected:

        ModifiedCostScaling() {}

    public:

        /// \brief Constructor.
        ///
        /// The constructor of the class.
        ///
        /// \param graph The digraph the algorithm runs on.
        ModifiedCostScaling(const GR& graph) :
                _graph(graph), _node_id(graph), _arc_idf(graph), _arc_idb(graph),
                INF(std::numeric_limits<Value>::has_infinity ?
                    std::numeric_limits<Value>::infinity() :
                    std::numeric_limits<Value>::max())
        {
            // Check the number types
            LEMON_ASSERT(std::numeric_limits<Value>::is_signed,
                         "The flow type of ModifiedCostScaling must be signed");
            LEMON_ASSERT(std::numeric_limits<Cost>::is_signed,
                         "The cost type of ModifiedCostScaling must be signed");

            // Reset data structures
            reset();
        }

        /// \name Parameters
        /// The parameters of the algorithm can be specified using these
        /// functions.

        /// @{

        /// \brief Set the lower bounds on the arcs.
        ///
        /// This function sets the lower bounds on the arcs.
        /// If it is not used before calling \ref run(), the lower bounds
        /// will be set to zero on all arcs.
        ///
        /// \param map An arc map storing the lower bounds.
        /// Its \c Value type must be convertible to the \c Value type
        /// of the algorithm.
        ///
        /// \return <tt>(*this)</tt>
        template <typename LowerMap>
        ModifiedCostScaling& lowerMap(const LowerMap& map) {
            _has_lower = true;
            for (ArcIt a(_graph); a != INVALID; ++a) {
                _lower[_arc_idf[a]] = map[a];
            }
            return *this;
        }

        /// \brief Set the upper bounds (capacities) on the arcs.
        ///
        /// This function sets the upper bounds (capacities) on the arcs.
        /// If it is not used before calling \ref run(), the upper bounds
        /// will be set to \ref INF on all arcs (i.e. the flow value will be
        /// unbounded from above).
        ///
        /// \param map An arc map storing the upper bounds.
        /// Its \c Value type must be convertible to the \c Value type
        /// of the algorithm.
        ///
        /// \return <tt>(*this)</tt>
        template<typename UpperMap>
        ModifiedCostScaling& upperMap(const UpperMap& map) {
            for (ArcIt a(_graph); a != INVALID; ++a) {
                _upper[_arc_idf[a]] = map[a];
            }
            return *this;
        }

        /// \brief Set the costs of the arcs.
        ///
        /// This function sets the costs of the arcs.
        /// If it is not used before calling \ref run(), the costs
        /// will be set to \c 1 on all arcs.
        ///
        /// \param map An arc map storing the costs.
        /// Its \c Value type must be convertible to the \c Cost type
        /// of the algorithm.
        ///
        /// \return <tt>(*this)</tt>
        template<typename CostMap>
        ModifiedCostScaling& costMap(const CostMap& map) {
            for (ArcIt a(_graph); a != INVALID; ++a) {
                _scost[_arc_idf[a]] =  map[a];
                _scost[_arc_idb[a]] = -map[a];
            }
            return *this;
        }

        /// \brief Set the supply values of the nodes.
        ///
        /// This function sets the supply values of the nodes.
        /// If neither this function nor \ref stSupply() is used before
        /// calling \ref run(), the supply of each node will be set to zero.
        ///
        /// \param map A node map storing the supply values.
        /// Its \c Value type must be convertible to the \c Value type
        /// of the algorithm.
        ///
        /// \return <tt>(*this)</tt>
        template<typename SupplyMap>
        ModifiedCostScaling& supplyMap(const SupplyMap& map) {
            for (NodeIt n(_graph); n != INVALID; ++n) {
                _supply[_node_id[n]] = map[n];
            }
            return *this;
        }

        /// \brief Set single source and target nodes and a supply value.
        ///
        /// This function sets a single source node and a single target node
        /// and the required flow value.
        /// If neither this function nor \ref supplyMap() is used before
        /// calling \ref run(), the supply of each node will be set to zero.
        ///
        /// Using this function has the same effect as using \ref supplyMap()
        /// with a map in which \c k is assigned to \c s, \c -k is
        /// assigned to \c t and all other nodes have zero supply value.
        ///
        /// \param s The source node.
        /// \param t The target node.
        /// \param k The required amount of flow from node \c s to node \c t
        /// (i.e. the supply of \c s and the demand of \c t).
        ///
        /// \return <tt>(*this)</tt>
        ModifiedCostScaling& stSupply(const Node& s, const Node& t, Value k) {
            for (int i = 0; i != _res_node_num; ++i) {
                _supply[i] = 0;
            }
            _supply[_node_id[s]] =  k;
            _supply[_node_id[t]] = -k;
            return *this;
        }

        /// @}

        /// \name Execution control
        /// The algorithm can be executed using \ref run().

        /// @{

        /// \brief Run the algorithm.
        ///
        /// This function runs the algorithm.
        /// The paramters can be specified using functions \ref lowerMap(),
        /// \ref upperMap(), \ref costMap(), \ref supplyMap(), \ref stSupply().
        /// For example,
        /// \code
        ///   ModifiedCostScaling<ListDigraph> cs(graph);
        ///   cs.lowerMap(lower).upperMap(upper).costMap(cost)
        ///     .supplyMap(sup).run();
        /// \endcode
        ///
        /// This function can be called more than once. All the given parameters
        /// are kept for the next call, unless \ref resetParams() or \ref reset()
        /// is used, thus only the modified parameters have to be set again.
        /// If the underlying digraph was also modified after the construction
        /// of the class (or the last \ref reset() call), then the \ref reset()
        /// function must be called.
        ///
        /// \param method The internal method that will be used in the
        /// algorithm. For more information, see \ref Method.
        /// \param factor The cost scaling factor. It must be at least two.
        ///
        /// \return \c INFEASIBLE if no feasible flow exists,
        /// \n \c OPTIMAL if the problem has optimal solution
        /// (i.e. it is feasible and bounded), and the algorithm has found
        /// optimal flow and node potentials (primal and dual solutions),
        /// \n \c UNBOUNDED if the digraph contains an arc of negative cost
        /// and infinite upper bound. It means that the objective function
        /// is unbounded on that arc, however, note that it could actually be
        /// bounded over the feasible flows, but this algroithm cannot handle
        /// these cases.
        ///
        /// \see ProblemType, Method
        /// \see resetParams(), reset()
        ProblemType run(Method method = PARTIAL_AUGMENT, int factor = 16) {
            LEMON_ASSERT(factor >= 2, "The scaling factor must be at least 2");
            _alpha = factor;
            ProblemType pt = init();
            if (pt != OPTIMAL) return pt;
            start(method);
            return OPTIMAL;
        }

        /// \brief Reset all the parameters that have been given before.
        ///
        /// This function resets all the paramaters that have been given
        /// before using functions \ref lowerMap(), \ref upperMap(),
        /// \ref costMap(), \ref supplyMap(), \ref stSupply().
        ///
        /// It is useful for multiple \ref run() calls. Basically, all the given
        /// parameters are kept for the next \ref run() call, unless
        /// \ref resetParams() or \ref reset() is used.
        /// If the underlying digraph was also modified after the construction
        /// of the class or the last \ref reset() call, then the \ref reset()
        /// function must be used, otherwise \ref resetParams() is sufficient.
        ///
        /// For example,
        /// \code
        ///   ModifiedCostScaling<ListDigraph> cs(graph);
        ///
        ///   // First run
        ///   cs.lowerMap(lower).upperMap(upper).costMap(cost)
        ///     .supplyMap(sup).run();
        ///
        ///   // Run again with modified cost map (resetParams() is not called,
        ///   // so only the cost map have to be set again)
        ///   cost[e] += 100;
        ///   cs.costMap(cost).run();
        ///
        ///   // Run again from scratch using resetParams()
        ///   // (the lower bounds will be set to zero on all arcs)
        ///   cs.resetParams();
        ///   cs.upperMap(capacity).costMap(cost)
        ///     .supplyMap(sup).run();
        /// \endcode
        ///
        /// \return <tt>(*this)</tt>
        ///
        /// \see reset(), run()
        ModifiedCostScaling& resetParams() {
            for (int i = 0; i != _res_node_num; ++i) {
                _supply[i] = 0;
            }
            int limit = _first_out[_root];
            for (int j = 0; j != limit; ++j) {
                _lower[j] = 0;
                _upper[j] = INF;
                _scost[j] = _forward[j] ? 1 : -1;
            }
            for (int j = limit; j != _res_arc_num; ++j) {
                _lower[j] = 0;
                _upper[j] = INF;
                _scost[j] = 0;
                _scost[_reverse[j]] = 0;
            }
            _has_lower = false;
            return *this;
        }

        /// \brief Reset the internal data structures and all the parameters
        /// that have been given before.
        ///
        /// This function resets the internal data structures and all the
        /// paramaters that have been given before using functions \ref lowerMap(),
        /// \ref upperMap(), \ref costMap(), \ref supplyMap(), \ref stSupply().
        ///
        /// It is useful for multiple \ref run() calls. By default, all the given
        /// parameters are kept for the next \ref run() call, unless
        /// \ref resetParams() or \ref reset() is used.
        /// If the underlying digraph was also modified after the construction
        /// of the class or the last \ref reset() call, then the \ref reset()
        /// function must be used, otherwise \ref resetParams() is sufficient.
        ///
        /// See \ref resetParams() for examples.
        ///
        /// \return <tt>(*this)</tt>
        ///
        /// \see resetParams(), run()
        ModifiedCostScaling& reset() {
            // Resize vectors
            _node_num = countNodes(_graph);
            _arc_num = countArcs(_graph);
            _res_node_num = _node_num + 1;
            _res_arc_num = 2 * (_arc_num + _node_num);
            _root = _node_num;

            _first_out.resize(_res_node_num + 1);
            _forward.resize(_res_arc_num);
            _source.resize(_res_arc_num);
            _target.resize(_res_arc_num);
            _reverse.resize(_res_arc_num);

            _lower.resize(_res_arc_num);
            _upper.resize(_res_arc_num);
            _scost.resize(_res_arc_num);
            _supply.resize(_res_node_num);

            _res_cap.resize(_res_arc_num);
            _cost.resize(_res_arc_num);
            _pi.resize(_res_node_num);
            _excess.resize(_res_node_num);
            _next_out.resize(_res_node_num);

            // Copy the graph
            int i = 0, j = 0, k = 2 * _arc_num + _node_num;
            for (NodeIt n(_graph); n != INVALID; ++n, ++i) {
                _node_id[n] = i;
            }
            i = 0;
            for (NodeIt n(_graph); n != INVALID; ++n, ++i) {
                _first_out[i] = j;
                for (OutArcIt a(_graph, n); a != INVALID; ++a, ++j) {
                    _arc_idf[a] = j;
                    _forward[j] = true;
                    _source[j] = i;
                    _target[j] = _node_id[_graph.runningNode(a)];
                }
                for (InArcIt a(_graph, n); a != INVALID; ++a, ++j) {
                    _arc_idb[a] = j;
                    _forward[j] = false;
                    _source[j] = i;
                    _target[j] = _node_id[_graph.runningNode(a)];
                }
                _forward[j] = false;
                _source[j] = i;
                _target[j] = _root;
                _reverse[j] = k;
                _forward[k] = true;
                _source[k] = _root;
                _target[k] = i;
                _reverse[k] = j;
                ++j; ++k;
            }
            _first_out[i] = j;
            _first_out[_res_node_num] = k;
            for (ArcIt a(_graph); a != INVALID; ++a) {
                int fi = _arc_idf[a];
                int bi = _arc_idb[a];
                _reverse[fi] = bi;
                _reverse[bi] = fi;
            }

            // Reset parameters
            resetParams();
            return *this;
        }

        /// @}

        /// \name Query Functions
        /// The results of the algorithm can be obtained using these
        /// functions.\n
        /// The \ref run() function must be called before using them.

        /// @{

        /// \brief Return the total cost of the found flow.
        ///
        /// This function returns the total cost of the found flow.
        /// Its complexity is O(m).
        ///
        /// \note The return type of the function can be specified as a
        /// template parameter. For example,
        /// \code
        ///   cs.totalCost<double>();
        /// \endcode
        /// It is useful if the total cost cannot be stored in the \c Cost
        /// type of the algorithm, which is the default return type of the
        /// function.
        ///
        /// \pre \ref run() must be called before using this function.
        template <typename Number>
        Number totalCost() const {
            Number c = 0;
            for (ArcIt a(_graph); a != INVALID; ++a) {
                int i = _arc_idb[a];
                c += static_cast<Number>(_res_cap[i]) *
                     (-static_cast<Number>(_scost[i]));
            }
            return c;
        }

#ifndef DOXYGEN
        Cost totalCost() const {
            return totalCost<Cost>();
        }
#endif

        /// \brief Return the flow on the given arc.
        ///
        /// This function returns the flow on the given arc.
        ///
        /// \pre \ref run() must be called before using this function.
        Value flow(const Arc& a) const {
            return _res_cap[_arc_idb[a]];
        }

        /// \brief Copy the flow values (the primal solution) into the
        /// given map.
        ///
        /// This function copies the flow value on each arc into the given
        /// map. The \c Value type of the algorithm must be convertible to
        /// the \c Value type of the map.
        ///
        /// \pre \ref run() must be called before using this function.
        template <typename FlowMap>
        void flowMap(FlowMap &map) const {
            for (ArcIt a(_graph); a != INVALID; ++a) {
                map.set(a, _res_cap[_arc_idb[a]]);
            }
        }

        /// \brief Return the potential (dual value) of the given node.
        ///
        /// This function returns the potential (dual value) of the
        /// given node.
        ///
        /// \pre \ref run() must be called before using this function.
        Cost potential(const Node& n) const {
            return static_cast<Cost>(_pi[_node_id[n]]);
        }

        /// \brief Copy the potential values (the dual solution) into the
        /// given map.
        ///
        /// This function copies the potential (dual value) of each node
        /// into the given map.
        /// The \c Cost type of the algorithm must be convertible to the
        /// \c Value type of the map.
        ///
        /// \pre \ref run() must be called before using this function.
        template <typename PotentialMap>
        void potentialMap(PotentialMap &map) const {
            for (NodeIt n(_graph); n != INVALID; ++n) {
                map.set(n, static_cast<Cost>(_pi[_node_id[n]]));
            }
        }

        /// @}

    private:

        // Initialize the algorithm
        ProblemType init() {
            if (_res_node_num <= 1) return INFEASIBLE;

            // Check the sum of supply values
            _sum_supply = 0;
            for (int i = 0; i != _root; ++i) {
                _sum_supply += _supply[i];
            }
            if (_sum_supply > 0) return INFEASIBLE;

            // Check lower and upper bounds
            LEMON_DEBUG(checkBoundMaps(),
                        "Upper bounds must be greater or equal to the lower bounds");


            // Initialize vectors
            for (int i = 0; i != _res_node_num; ++i) {
                _pi[i] = 0;
                _excess[i] = _supply[i];
            }

            // Remove infinite upper bounds and check negative arcs
            const Value MAX = std::numeric_limits<Value>::max();
            int last_out;
            if (_has_lower) {
                for (int i = 0; i != _root; ++i) {
                    last_out = _first_out[i+1];
                    for (int j = _first_out[i]; j != last_out; ++j) {
                        if (_forward[j]) {
                            Value c = _scost[j] < 0 ? _upper[j] : _lower[j];
                            if (c >= MAX) return UNBOUNDED;
                            _excess[i] -= c;
                            _excess[_target[j]] += c;
                        }
                    }
                }
            } else {
                for (int i = 0; i != _root; ++i) {
                    last_out = _first_out[i+1];
                    for (int j = _first_out[i]; j != last_out; ++j) {
                        if (_forward[j] && _scost[j] < 0) {
                            Value c = _upper[j];
                            if (c >= MAX) return UNBOUNDED;
                            _excess[i] -= c;
                            _excess[_target[j]] += c;
                        }
                    }
                }
            }
            Value ex, max_cap = 0;
            for (int i = 0; i != _res_node_num; ++i) {
                ex = _excess[i];
                _excess[i] = 0;
                if (ex < 0) max_cap -= ex;
            }
            for (int j = 0; j != _res_arc_num; ++j) {
                if (_upper[j] >= MAX) _upper[j] = max_cap;
            }

            // Initialize the large cost vector and the epsilon parameter
            _epsilon = 0;
            LargeCost lc;
            for (int i = 0; i != _root; ++i) {
                last_out = _first_out[i+1];
                for (int j = _first_out[i]; j != last_out; ++j) {
                    lc = static_cast<LargeCost>(_scost[j]) * _res_node_num * _alpha; //COST MODIFICATION
                    _cost[j] = lc;
                    if (lc > _epsilon) _epsilon = lc;
                }
            }
            _epsilon /= _alpha;

            // Initialize maps for Circulation and remove non-zero lower bounds
            ConstMap<Arc, Value> low(0);
            typedef typename Digraph::template ArcMap<Value> ValueArcMap;
            typedef typename Digraph::template NodeMap<Value> ValueNodeMap;
            ValueArcMap cap(_graph), flow(_graph);
            ValueNodeMap sup(_graph);
            for (NodeIt n(_graph); n != INVALID; ++n) {
                sup[n] = _supply[_node_id[n]];
            }
            if (_has_lower) {
                for (ArcIt a(_graph); a != INVALID; ++a) {
                    int j = _arc_idf[a];
                    Value c = _lower[j];
                    cap[a] = _upper[j] - c;
                    sup[_graph.source(a)] -= c;
                    sup[_graph.target(a)] += c;
                }
            } else {
                for (ArcIt a(_graph); a != INVALID; ++a) {
                    cap[a] = _upper[_arc_idf[a]];
                }
            }

            _sup_node_num = 0;
            for (NodeIt n(_graph); n != INVALID; ++n) {
                if (sup[n] > 0) ++_sup_node_num;
            }

            // Find a feasible flow using Circulation
            Circulation<Digraph, ConstMap<Arc, Value>, ValueArcMap, ValueNodeMap>
                    circ(_graph, low, cap, sup);
            if (!circ.flowMap(flow).run()) return INFEASIBLE;

            // Set residual capacities and handle GEQ supply type
            assert(_sum_supply == 0);
            for (ArcIt a(_graph); a != INVALID; ++a) {
                Value fa = flow[a];
                _res_cap[_arc_idf[a]] = cap[a] - fa;
                _res_cap[_arc_idb[a]] = fa;
            }
            for (int a = _first_out[_root]; a != _res_arc_num; ++a) {
                int ra = _reverse[a];
                _res_cap[a] = 0;
                _res_cap[ra] = 0;
                _cost[a] = 0;
                _cost[ra] = 0;
            }

            return OPTIMAL;
        }

        // Check if the upper bound is greater than or equal to the lower bound
        // on each forward arc.
        bool checkBoundMaps() {
            for (int j = 0; j != _res_arc_num; ++j) {
                if (_forward[j] && _upper[j] < _lower[j]) return false;
            }
            return true;
        }

        // Execute the algorithm and transform the results
        void start(Method method) {
            const int MAX_PARTIAL_PATH_LENGTH = 4;

            switch (method) {
                case AUGMENT:
                    startAugment(_res_node_num - 1);
                    break;
                case PARTIAL_AUGMENT:
                    startAugment(MAX_PARTIAL_PATH_LENGTH);
                    break;
            }

            // Compute node potentials (dual solution)
            for (int i = 0; i != _res_node_num; ++i) {
                _pi[i] = static_cast<Cost>(_pi[i] / (_res_node_num * _alpha));
            }
//            bool optimal = true;
//            for (int i = 0; optimal && i != _res_node_num; ++i) {
//                LargeCost pi_i = _pi[i];
//                int last_out = _first_out[i+1];
//                for (int j = _first_out[i]; j != last_out; ++j) {
//                    if (_res_cap[j] > 0 && _scost[j] + pi_i - _pi[_target[j]] < 0) {
//                        optimal = false;
//                        break;
//                    }
//                }
//            }
//
//            if (!optimal) {
//                // Compute node potentials for the original costs with BellmanFord
//                // (if it is necessary)
//                typedef std::pair<int, int> IntPair;
//                StaticDigraph sgr;
//                std::vector<IntPair> arc_vec;
//                std::vector<LargeCost> cost_vec;
//                LargeCostArcMap cost_map(cost_vec);
//
//                arc_vec.clear();
//                cost_vec.clear();
//                for (int j = 0; j != _res_arc_num; ++j) {
//                    if (_res_cap[j] > 0) {
//                        int u = _source[j], v = _target[j];
//                        arc_vec.push_back(IntPair(u, v));
//                        cost_vec.push_back(_scost[j] + _pi[u] - _pi[v]);
//                    }
//                }
//                sgr.build(_res_node_num, arc_vec.begin(), arc_vec.end());
//
//                typename BellmanFord<StaticDigraph, LargeCostArcMap>::Create
//                        bf(sgr, cost_map);
//                bf.init(0);
//                bf.start();
//
//                for (int i = 0; i != _res_node_num; ++i) {
//                    _pi[i] += bf.dist(sgr.node(i));
//                }
//            }
//
//            // Shift potentials to meet the requirements of the GEQ type
//            // optimality conditions
//            LargeCost max_pot = _pi[_root];
//            for (int i = 0; i != _res_node_num; ++i) {
//                if (_pi[i] > max_pot) max_pot = _pi[i];
//            }
//            if (max_pot != 0) {
//                for (int i = 0; i != _res_node_num; ++i) {
//                    _pi[i] -= max_pot;
//                }
//            }

            // Handle non-zero lower bounds
            if (_has_lower) {
                int limit = _first_out[_root];
                for (int j = 0; j != limit; ++j) {
                    if (_forward[j]) _res_cap[_reverse[j]] += _lower[j];
                }
            }
        }

        // Initialize a cost scaling phase
        void initPhase() {
            // Saturate arcs not satisfying the optimality condition
            for (int u = 0; u != _res_node_num; ++u) {
                int last_out = _first_out[u+1];
                LargeCost pi_u = _pi[u];
                for (int a = _first_out[u]; a != last_out; ++a) {
                    Value delta = _res_cap[a];
                    if (delta > 0) {
                        int v = _target[a];
                        if (_cost[a] + pi_u - _pi[v] < 0) {
                            _excess[u] -= delta;
                            _excess[v] += delta;
                            _res_cap[a] = 0;
                            _res_cap[_reverse[a]] += delta;
                        }
                    }
                }
            }

            // Find active nodes (i.e. nodes with positive excess)
            for (int u = 0; u != _res_node_num; ++u) {
                if (_excess[u] > 0) _active_nodes.push_back(u);
            }

            // Initialize the next arcs
            for (int u = 0; u != _res_node_num; ++u) {
                _next_out[u] = _first_out[u];
            }
        }

        /// Execute the algorithm performing augment and relabel operations
        void startAugment(int max_length) {
            // Paramters for heuristics
            const int PRICE_REFINEMENT_LIMIT = 2;
            const double GLOBAL_UPDATE_FACTOR = 1.0;
            const int global_update_skip = static_cast<int>(GLOBAL_UPDATE_FACTOR *
                                                            (_res_node_num + _sup_node_num * _sup_node_num));
            int next_global_update_limit = global_update_skip;

            // Perform cost scaling phases
            IntVector path;
            BoolVector path_arc(_res_arc_num, false);
            int relabel_cnt = 0;
            int eps_phase_cnt = 0;
            for ( ; _epsilon >= 1; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                              1 : _epsilon / _alpha )
            {
                ++eps_phase_cnt;

                // Price refinement heuristic
//        if (eps_phase_cnt >= PRICE_REFINEMENT_LIMIT) {
//          if (priceRefinement()) continue;
//        }

                // Initialize current phase
                initPhase();

                // Perform partial augment and relabel operations
                while (true) {
                    // Select an active node (FIFO selection)
                    while (_active_nodes.size() > 0 &&
                           _excess[_active_nodes.front()] <= 0) {
                        _active_nodes.pop_front();
                    }
                    if (_active_nodes.size() == 0) break;
                    int start = _active_nodes.front();

                    // Find an augmenting path from the start node
                    int tip = start;
                    while (int(path.size()) < max_length && _excess[tip] >= 0) {
                        int u;
                        LargeCost rc, min_red_cost = std::numeric_limits<LargeCost>::max();
                        LargeCost pi_tip = _pi[tip];
                        int last_out = _first_out[tip+1];
                        for (int a = _next_out[tip]; a != last_out; ++a) {
                            if (_res_cap[a] > 0) {
                                u = _target[a];
                                rc = _cost[a] + pi_tip - _pi[u];
                                if (rc < 0) {
                                    path.push_back(a);
                                    _next_out[tip] = a;
                                    if (path_arc[a]) {
                                        goto augment;   // a cycle is found, stop path search
                                    }
                                    tip = u;
                                    path_arc[a] = true;
                                    goto next_step;
                                }
                                else if (rc < min_red_cost) {
                                    min_red_cost = rc;
                                }
                            }
                        }

                        // Relabel tip node
                        if (tip != start) {
                            int ra = _reverse[path.back()];
                            min_red_cost =
                                    std::min(min_red_cost, _cost[ra] + pi_tip - _pi[_target[ra]]);
                        }
                        last_out = _next_out[tip];
                        for (int a = _first_out[tip]; a != last_out; ++a) {
                            if (_res_cap[a] > 0) {
                                rc = _cost[a] + pi_tip - _pi[_target[a]];
                                if (rc < min_red_cost) {
                                    min_red_cost = rc;
                                }
                            }
                        }
                        _pi[tip] -= min_red_cost + _epsilon;
                        _next_out[tip] = _first_out[tip];
                        ++relabel_cnt;

                        // Step back
                        if (tip != start) {
                            int pa = path.back();
                            path_arc[pa] = false;
                            tip = _source[pa];
                            path.pop_back();
                        }

                        next_step: ;
                    }

                    // Augment along the found path (as much flow as possible)
                    augment:
                    Value delta;
                    int pa, u, v = start;
                    for (int i = 0; i != int(path.size()); ++i) {
                        pa = path[i];
                        u = v;
                        v = _target[pa];
                        path_arc[pa] = false;
                        delta = std::min(_res_cap[pa], _excess[u]);
                        _res_cap[pa] -= delta;
                        _res_cap[_reverse[pa]] += delta;
                        _excess[u] -= delta;
                        _excess[v] += delta;
                        if (_excess[v] > 0 && _excess[v] <= delta) {
                            _active_nodes.push_back(v);
                        }
                    }
                    path.clear();

                    // Global update heuristic
//          if (relabel_cnt >= next_global_update_limit) {
//            globalUpdate();
                    next_global_update_limit += global_update_skip;
//          }
                }

            }

        }
    }; //class ModifiedCostScaling

    ///@}

} //namespace lemon


#endif //NETWORKFLOW_LEMON_H

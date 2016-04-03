//
// Created by alvis on 20.03.16.
//

#include "LSCS.h"

void LSCS::startAugment(int max_length) {
    IntVector path;
    BoolVector path_arc(_res_arc_num, false);
    int relabel_cnt = 0;
    int eps_phase_cnt = 0;
    for ( ; _epsilon >= 1; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                      1 : _epsilon / _alpha )
    {
        ++eps_phase_cnt;

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

                for (int i = _next_out[tip]; i != _first_out[tip].size(); ++i) {
                    int a = _first_out[tip][i];
                    if (_res_cap[a] > 0) {
                        u = _target[a];
                        rc = _cost[a] + pi_tip - _pi[u];
                        if (rc < 0) {
                            path.push_back(a);
                            _next_out[tip] = i;
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

                //recalculate because some nodes before _next_out could have changed their values
                for (int i = 0; i < _first_out[tip].size(); ++i) {
                    int a = _first_out[tip][i];
                    if (_res_cap[a] > 0) {
                        rc = _cost[a] + pi_tip - _pi[_target[a]];
                        if (rc < min_red_cost) {
                            min_red_cost = rc;
                        }
                    }
                }
                _pi[tip] -= min_red_cost + _epsilon;
                _next_out[tip] = 0;
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
        }
    }
}
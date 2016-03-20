//
// Created by alvis on 17.03.16.
//

#ifndef NETWORKFLOW_TIMERTOOL_H
#define NETWORKFLOW_TIMERTOOL_H

#pragma once

#ifdef WIN32
//under WIN32, we use QueryPerformanceCounter()
#include <windows.h>
#include <ctime>
#endif

#if defined(unix) || defined(__unix__)
//under POSIX system, we use clock_gettime()
//remember we have to use linker option "-lrt"
#include <time.h>
#endif

#include <unordered_map>
#include <iostream>
#include <fstream>

using namespace std;

class Timer
{
public:
    unordered_map<string, double> timings;

    static Timer *instance;
    static double getTime()
    {

#ifdef WIN32
		LONGLONG clock_count;
		QueryPerformanceCounter((LARGE_INTEGER *)&clock_count);
		return (double)clock_count/getInstance()->clock_freq;
#endif

#if defined(unix) || defined(__unix__)
        timespec clock;
        clock_gettime(CLOCK_MONOTONIC, &clock);
        return (double)clock.tv_sec+(double)clock.tv_nsec*1e-9;
#endif

    }

#ifdef WIN32
	LONGLONG clock_freq;
#endif

    static Timer *getInstance()
    {
        if(instance==NULL) instance=new Timer();
        return instance;
    }

    Timer()
    {

#ifdef WIN32
		QueryPerformanceFrequency((LARGE_INTEGER *)&clock_freq);
#endif

    }

    ~Timer()
    {
        if(instance!=NULL) delete instance;
    }

    // create double value outside the class to avoid search in unordered_map
    void save_time(string name, double time) {
        if (timings.find(name) != timings.end()) {
            cout << "Error: Double timing value for " << name << endl;
            exit(1);
        }
        timings.insert({name,time});
    }

    void save_all_timings(string& filename, int experiment_id) {
        ofstream outf(filename,ios::app);

        for ( const auto& n : timings ) {
            outf << experiment_id << "," << n.first << "," << n.second << "\n";
        }
    }

    inline void finish(string title, double start_time) {
        double curtime = this->getTime();
        timings.insert({title, curtime - start_time});
    }
};

//finally our get timer interface
inline double getclock()
{
    return Timer::getTime();
}


#endif //NETWORKFLOW_TIMERTOOL_H

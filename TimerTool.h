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
#include <vector>

using namespace std;

class Timer
{
public:
    unordered_map<string, vector<double>> timings;

    static Timer *instance;
    static inline double getTime()
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


    inline void save_time_total(string title, double total_time) {
        std::unordered_map<std::string,vector<double>>::iterator got = timings.find(title);
        if (got == timings.end()) {
            vector<double> list;
            list.push_back(total_time);
            timings.insert({title, list});
        } else {
            timings[title].push_back(total_time);
        }
    }
    inline void save_time(string title, double start_time) {
        double curtime = this->getTime();
        save_time_total(title, curtime - start_time);
    }

    void output(string filename, int experiment_id) {
        ofstream outf(filename,ios::app);

        for ( auto& n : timings ) {
            vector<double>& curtimes = n.second;
            double sum = 0;
            double max = -1;
            double min = -1;
            for (vector<double>::iterator it = curtimes.begin(); it != curtimes.end(); ++it) {
                sum += *it;
                if (min == -1 || min > *it) min = *it;
                if (max == -1 || max < *it) max = *it;
            }
            sum /= curtimes.size();
            outf << experiment_id << "," << n.first << "," << sum << "\n";
            outf << experiment_id << "," << n.first << " error," << std::max<double>(max - sum, sum - min) << "\n";
        }
        outf.close();
    }

    double get_av_time(string title) {
        std::unordered_map<std::string,vector<double>>::iterator got = timings.find(title);
        if (got == timings.end()) {
            cout << "Error: no title " << title << " is found in the timer" << endl;
            exit(1);
        }
        vector<double>& curtimes = timings[title];
        double sum = 0;
        for (vector<double>::iterator it = curtimes.begin(); it != curtimes.end(); ++it) {
            sum += *it;
        }
        sum /= curtimes.size();
        return sum;
    }

    double get_last_time(string title) {
        std::unordered_map<std::string,vector<double>>::iterator got = timings.find(title);
        if (got == timings.end()) {
            cout << "Error: no title " << title << " is found in the timer" << endl;
            exit(1);
        }
        vector<double>& curtimes = timings[title];
        return curtimes.back();
    }
};

//finally our get timer interface
inline double getclock()
{
    return Timer::getTime();
}


#endif //NETWORKFLOW_TIMERTOOL_H

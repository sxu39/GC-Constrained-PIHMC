//
// Created by Jin Bin on 2022/03/04.
//

#ifndef CPIHMC_RUN_TIME_H
#define CPIHMC_RUN_TIME_H

#include <ctime>
#include <cstdio>

inline void run_time(time_t ts, time_t te){
    time_t delta = te - ts;
    int hour = delta / 3600;
    int minute = (delta % 3600) / 60;
    int second = delta % 60;
    printf("RUN TIME: %02d:%02d:%02d\n", hour, minute, second);
}


#endif //CPIHMC_RUN_TIME_H

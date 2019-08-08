#ifndef __TIMER_H__
#define __TIMER_H__

#include <assert.h>
#include <chrono>

class timer {
public:
  
    timer(const std::string& namee) : name(namee), started(false) {}
    
    void start(){
#ifdef CH_PROFILE
        assert(!started);
        t1 = time::now();
        started = true;
#endif
    }
    
    void stop(){
#ifdef CH_PROFILE
        t2 = time::now();
        accumulated += std::chrono::duration_cast< std::chrono::duration<double> >(t2 - t1);
        started = false;
#endif
    }
    
    void reset(){
#ifdef CH_PROFILE
        accumulated = std::chrono::duration<double>(0.);
        started = false;
#endif
    }
    
    void print() const{
#ifdef CH_PROFILE
        std::cout << "Timer " << name << ": " << accumulated.count() << " [s]" << std::endl;
#endif
    }
    
private:
    bool started;
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> accumulated;
    std::string name;
    using time=std::chrono::high_resolution_clock;
};
#endif

#include <chrono>

class timer {
public:
    timer(const std::string& namee) : name(namee), started(false) {}
    
    void start(){
        assert(!started);
        t1 = high_resolution_clock::now();
        started = true;
    }
    
    void stop(){
        t2 = high_resolution_clock::now();
        accumulated += duration_cast< duration<double> >(t2 - t1);
        started = false;
    }
    
    void reset(){
        accumulated = duration_cast< duration<double> >(0.);
        started = false;
    }
    
    void print() const{
        std::cout << "Timer " << name << ": " << accumulated.count() << " [s]" << std::endl;
    }
    
private:
    bool started;
    high_resolution_clock::time_point t1, t2;
    duration<double> accumulated;
    std::string name;
    std::chrono::high_resolution_clock time;
};
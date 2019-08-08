#include <chrono>

class timer {
public:
    timer(const std::string& namee) : name(namee), started(false) {}
    
    void start(){
        assert(!started);
        t1 = time::now();
        started = true;
    }
    
    void stop(){
        t2 = time::now();
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
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> accumulated;
    std::string name;
    std::chrono::high_resolution_clock time;
};
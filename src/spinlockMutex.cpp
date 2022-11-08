#include "spinlockMutex.h" 

spinlockMutex::spinlockMutex(){}

void spinlockMutex::lock(){while (this->flag.test_and_set(std::memory_order_acquire));}

void spinlockMutex::unlock(){this->flag.clear(std::memory_order_release);}


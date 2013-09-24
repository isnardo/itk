#ifndef _TIMER_H
#define _TIMER_H

#include <sys/time.h>

namespace Linux{
	//Clase para Contar el tiempo en Milisegundos Linux
	class Timer{
		 timeval timer[2];
	  public:
		 timeval start(){
		     gettimeofday(&this->timer[0], NULL);
		     return this->timer[0];
		 }

		 timeval stop(){
		     gettimeofday(&this->timer[1], NULL);
		     return this->timer[1];
		 }

		 int duration() const{
		     int secs(this->timer[1].tv_sec - this->timer[0].tv_sec);
		     int usecs(this->timer[1].tv_usec - this->timer[0].tv_usec);

		     if(usecs < 0){
		         --secs;
		         usecs += 1000000;
		     }
		     return static_cast<int>(secs * 1000 + usecs / 1000.0 + 0.5);
		 }
	};
}

#endif

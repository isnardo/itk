#ifdef __linux
	//Linux Library
	#include <unistd.h>
#elif _WIN32
	//Windows Library
	#include <windows.h>
#else
	//MACOS Libaries
	#include <sys/param.h>
	#include <sys/sysctl.h>
#endif

int getNumberOfCores() {
#ifdef __linux
	//Linux
	return sysconf(_SC_NPROCESSORS_ONLN);
#elif __WIN32
	//Windows 32 bits
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#else
	//MACOS 
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
    return count;    
#endif
}

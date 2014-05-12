/*
  thread.h 2008-01-03

  Modified from CRF++ toolkit http://sourceforge.net/projects/crfpp/
  
  Author taku
*/
#ifndef THREAD_H
#define THREAD_H

#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else			//linux or cygwin
#include <pthread.h>
#endif


#if(defined(_WIN32) && !defined (__CYGWIN__))
#define BEGINTHREAD(src, stack, func, arg, flag, id) \
     (HANDLE)_beginthreadex((void *)(src), (unsigned)(stack), \
                       (unsigned(_stdcall *)(void *))(func), (void *)(arg), \
                       (unsigned)(flag), (unsigned *)(id))
#endif



class thread
{
private:
#ifdef _WIN32
	HANDLE hnd;
#else
    pthread_t hnd;
#endif
public:
    static void* wrapper(void *ptr)
	{
		thread *p = static_cast<thread *>(ptr);
		p->run();
		return 0;
	}

    virtual void run() {};

    void start()
	{
#ifdef _WIN32
      DWORD id;
      hnd= BEGINTHREAD(0, 0, &thread::wrapper, this, 0, &id);
#else 
      pthread_create(&hnd, 0, &thread::wrapper,
                     static_cast<void *>(this));
#endif
    }

    void join()
	{

#ifdef _WIN32
		WaitForSingleObject(hnd, INFINITE);
		CloseHandle(hnd);
#else
		pthread_join(hnd, 0);
#endif
	}
	virtual ~thread() {};
};
#endif

#ifndef __THREADS_H__
#define __THREADS_H__

#include"kmer.h"
#include"read.h"
#include<pthread.h>

void* threadRoutine(void* threadId_p);
void SendSignal(int buffId, int SIG);
int GetEmptyBuff();
void WaitForThreads();
#endif

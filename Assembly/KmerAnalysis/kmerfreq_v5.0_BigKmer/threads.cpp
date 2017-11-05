#include"threads.h"
#include "kmer.h"

extern pthread_mutex_t myMutex;
extern int* threadSIG;
extern int threadNum;
extern int bufferNum;
extern uint64_t** hcBuffer;
extern uint64_t** kmerStat;
extern uint64_t* readNum;
extern uint64_t* kmerNum;
extern int** isBuffFull;
extern int** isBuffEmpty;
extern int mode;
//extern uint64_t* kmerBuffer;

void* threadRoutine(void* threadId_p)
{
	int threadId = *((int*)threadId_p);
	while (1)
	{
		for (int buffId=0; buffId<bufferNum; buffId++)
		{
			if (isBuffFull[buffId][threadId] == 1)
			{
				for (uint64_t i=threadId; i<readNum[buffId]; i+=threadNum)
				{

					Read2Kmer(buffId, i, threadId);
				}

				isBuffFull[buffId][threadId] = 2;

				while (1)
				{
					int j=0;
					for (; j<threadNum; j++)
					{
						if (isBuffFull[buffId][j] != 2)
						{
							break;
						}
					}

					if (j == threadNum)
					{
						break;
					}

					usleep(5);
				}

/*				if (mode == 0)
				{
					for (uint64_t i=0; i<kmerNum[buffId]; i++)
					{
						if (hcBuffer[buffId][i]%threadNum != threadId)
						{
							continue;
						}


						HandleKmer0(buffId, i, threadId);
					}
				}
				else if (mode == 1)
				{
*/
					for (uint64_t i=0; i<kmerNum[buffId]; i++)
					{
						if ((kmerStat[buffId][i] == 0) || (hcBuffer[buffId][i]%threadNum != threadId))
						{
							continue;
						}


						HandleKmer1(buffId, i, threadId);
						kmerStat[buffId][i] = 0;
					}
//				}
				isBuffEmpty[buffId][threadId] = 1;
			}//if
		}//for

		if (threadSIG[threadId] == 1)
		{
			pthread_exit(NULL);
		}

		usleep(1);
	}//while
}

void SendSignal(int buffId, int SIG)
{
	for (int i=0; i<threadNum; i++)
	{
		isBuffFull[buffId][i] = SIG;

		isBuffEmpty[buffId][i] = 0;
	}

}

int GetEmptyBuff()
{
	while (1)
	{
		for (int buffId=0; buffId<bufferNum; buffId++)
		{
			int i=0;

			for (; i<threadNum; i++)
			{
				if (isBuffEmpty[buffId][i] != 1)
				{
					break;
				}
			}

			if (i == threadNum)
			{
				return buffId;
			}
		}

		usleep(10);
	}
}

void WaitForThreads()
{
	while (1)
	{
		int finishedBuffNum = 0;

		for (int buffId=0; buffId<bufferNum; buffId++)
		{
			int i=0;

			for (; i<threadNum; i++)
			{
				if (isBuffEmpty[buffId][i] != 1)
				{
					break;
				}
			}

			if (i == threadNum)
			{
				finishedBuffNum++;
			}
		}

		if (finishedBuffNum == bufferNum)
		{
			break;
		}

		usleep(10);
	}
}

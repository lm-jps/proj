/* The purpose of this module is to demonstrate how to use a multi-threaded alarm.
 * The alarm is implemented as a thread-specific signal in libthreadutil.a.
 *
 *  --Art Amezcua
 */

#include "jsoc_main.h"
#include "tdsignals.h"

char *module_name = "threadalrm";

#define kIn     "in"
#define kNada   "nada"

ModuleArgs_t module_args[] =
{
     {ARG_STRING, kIn, kNada,          "Just for looks"},
     {ARG_END}
};

typedef enum
{
   kThreadSigErr_Success
} ThreadSigErr_t;

static int gLoop = 1;
static pthread_mutex_t mutex;
static struct timeval tv0;

/* Function callback that gets called when alarm 'rings'.  This runs in the signal
 * thread that receives, the alarm signal, not the main thread. */
void shandler(int sig, pthread_mutex_t *mtx)
{
   float elapsed;
   struct timeval tv1;
   gettimeofday(&tv1, NULL);

   elapsed = (float)((tv1.tv_sec * 1000000.0 + tv1.tv_usec -
                      (tv0.tv_sec * 1000000.0 + tv0.tv_usec)) / 1000000.0);

   fprintf(stdout, "Thread '%lld' received alarm signal '%d'.\n", (long long )pthread_self(), sig);
   fprintf(stdout, "Elapsed time is %f seconds.\n", elapsed);

   /* This isn't in the main thread, so put globals in a critical region. */
   pthread_mutex_lock(mtx);
   gLoop = 0;
   pthread_mutex_unlock(mtx);
}

/* DoIt
 *
 *   Module entry point.
 */
int DoIt(void) 
{
   ThreadSigErr_t error = kThreadSigErr_Success;
   td_alarm_t alarm = 0;

   pthread_mutex_init(&mutex, NULL);

   /* Set alarm */
   int seconds = 5;
   gettimeofday(&tv0, NULL);
   if (!td_createalarm(seconds, shandler, &mutex, &alarm))
   {
      /* Do something.  Alarm will cause shandler to be called in <seconds> seconds. */
      while (1)
      {
         /* gLoop is shared by threads, so put access to it in a critical region */
         pthread_mutex_lock(&mutex);
         if (gLoop != 0 && gLoop % 100000 == 0)
         {
            fprintf(stdout, "Iteration # %d\n", gLoop);
         }
         pthread_mutex_unlock(&mutex);

         pthread_mutex_lock(&mutex);
         if (gLoop)
         {
            ++gLoop;
         }
         else
         {
            pthread_mutex_unlock(&mutex);
            break;
         }
         pthread_mutex_unlock(&mutex);
      }

      fprintf(stdout, "Got out of loop - exiting.\n");

      /* You must destory all alarms before calling pthread_mutext_destory() */
      td_destroyalarm(&alarm);
   }

   pthread_mutex_destroy(&mutex);

   return error;
}

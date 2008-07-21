/* The purpose of this module is to demonstrate how to use multi-threaded alarms.
 * Each alarm is implemented as a thread-specific signal in libthreadutil.a.  There
 * are five alarms in this example to demonstrate that any number of alarms
 * can be used, but they can be used sequentially only.  If you create alarms A and B,
 * and A triggers in X seconds, and B triggers in Y seconds, then execution will
 * result in A tigerring in X seconds, and B in X + Y seconds.
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
static int gRearm = 0;
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

void shandler2(int sig, pthread_mutex_t *mtx)
{
   float elapsed;
   struct timeval tv1;
   gettimeofday(&tv1, NULL);

   elapsed = (float)((tv1.tv_sec * 1000000.0 + tv1.tv_usec -
                      (tv0.tv_sec * 1000000.0 + tv0.tv_usec)) / 1000000.0);

   fprintf(stdout, "Thread '%lld' received alarm signal '%d'.\n", (long long )pthread_self(), sig);
   fprintf(stdout, "Elapsed time is %f seconds.\n", elapsed);

   pthread_mutex_lock(mtx);
   gRearm = 1;
   pthread_mutex_unlock(mtx);
}

/* DoIt
 *
 *   Module entry point.
 */
int DoIt(void) 
{
   int err = 0;
   ThreadSigErr_t error = kThreadSigErr_Success;
   td_alarm_t alarm = 0; /* exit loop */
   td_alarm_t alarm2 = 0; /* print something */

   pthread_mutex_init(&mutex, NULL);

   /* Set alarm */
   int seconds = 4;
   int seconds2 = 2;
   int ntimes = 3;
   gettimeofday(&tv0, NULL);

   /* Initially, create two alarms.  The first one is set for two seconds, and triggers in
    * two seconds.  Theh second one is set for four seconds, and triggers in six seconds.  
    * Thereafter, an alarm is set for two seconds three more times.  So, this alarm 
    * triggers at eight, ten, and twelve seconds. */
   err = td_createalarm(seconds2, shandler2, &mutex, &alarm2) || 
     td_createalarm(seconds, shandler, &mutex, &alarm);

   if (!err)
   {
      /* Do something.  Alarm will cause shandler to be called in <seconds> seconds. */
      fprintf(stdout, "Entering loop.\n");

      while (1)
      {
         /* gLoop is shared by threads, so put access to it in a critical region */
         pthread_mutex_lock(&mutex);
         if (gLoop != 0 && gLoop % 1000000 == 0)
         {
            fprintf(stdout, "Iteration # %d\n", gLoop);
         }
         pthread_mutex_unlock(&mutex);

         /* Check second alarm */
         pthread_mutex_lock(&mutex);
         if (gRearm)
         {
            if (ntimes-- > 0)
            {
               pthread_mutex_unlock(&mutex);

               td_destroyalarm(&alarm2, &mutex);

               if (td_createalarm(seconds2, shandler2, &mutex, &alarm2))
               {
                  pthread_mutex_unlock(&mutex);
                  break;
               }

               pthread_mutex_lock(&mutex);
               gRearm = 0;
               pthread_mutex_unlock(&mutex);
            }
            else
            {
               /* We've handled alarm2 ntimes */
               break;
            }
         }
         else
         {
            pthread_mutex_unlock(&mutex);
         }

         pthread_mutex_lock(&mutex);
         if (gLoop)
         {
            ++gLoop;
            pthread_mutex_unlock(&mutex);
         }
         else
         {
            gLoop = 1;
            pthread_mutex_unlock(&mutex);
            td_destroyalarm(&alarm, &mutex);
         }
      }

      fprintf(stdout, "Got out of loop.\n");
   }
   else
   {
      fprintf(stderr, "Failed to create one or more alarms.\n");
   }

   /* You must destroy all alarms before calling pthread_mutext_destory() */
   pthread_mutex_destroy(&mutex);

   return error;
}

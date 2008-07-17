/* The purpose of this module is to demonstrate how to send/handle signals between
 * threads within a DRMS module process.
 *
 *  --Art Amezcua
 */

#include "jsoc_main.h"

char *module_name = "threadsigs";

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

static void *sigthreadfxn(void *arg)
{
   sigset_t bsigs;
   sigset_t rsigs;
   int siggot;
   int err = 0;

   /* block all signals */
   sigfillset(&bsigs);
   pthread_sigmask(SIG_BLOCK, &bsigs, NULL);

   /* create a set of signals that includes only SIGALRM */
   sigemptyset(&rsigs);
   sigaddset(&rsigs, SIGALRM);

   /* Will block until SIGALRM is received. */
   err = sigwait(&rsigs, &siggot);

   if (!err)
   {
      fprintf(stdout, "Got ALRM signal.\n");

      /* gLoop is shared by threads, so put access to it in a critical region */
      pthread_mutex_lock(&mutex);
      gLoop = 0;
      pthread_mutex_unlock(&mutex);
   }

   return arg;
}

/* DoIt
 *
 *   Module entry point.
 */
int DoIt(void) 
{
   int status = 0;
   ThreadSigErr_t error = kThreadSigErr_Success;
   pthread_t sigthread;
   pthread_t mainthread;

   mainthread = pthread_self();

   pthread_mutex_init(&mutex, NULL);

   /* Ensure that the main thread doesn't handle the SIGALRM signal */
   sigset_t rsigs;
   sigemptyset(&rsigs);
   sigaddset(&rsigs, SIGALRM); 
   pthread_sigmask(SIG_BLOCK, &rsigs, NULL);

   /* Create a thread to handle the SIGALRM signal.  It will handle the signal, 
    * and set a global variable that will be used to break out of the while loop
    * below. */
   if((status = pthread_create(&sigthread, NULL, &sigthreadfxn, (void *)&mainthread)))
   {
      fprintf(stderr,"Thread creation failed: %d\n", status);          
      exit(1);
   }

   /* Loop here forever, unless we received a signal telling us to exist */
   while (1)
   {
      /* gLoop is shared by threads, so put access to it in a critical region */
      pthread_mutex_lock(&mutex);
      if (gLoop != 0 && gLoop % 1000 == 0)
      {
         fprintf(stdout, "Iteration # %d\n", gLoop);
      }
      pthread_mutex_unlock(&mutex);

      pthread_mutex_lock(&mutex);
      if (gLoop != 0 && gLoop == 1000)
      {
         /* now send SIGALRM signal */
         pthread_kill(sigthread, SIGALRM);      
      }
      pthread_mutex_unlock(&mutex);

      pthread_mutex_lock(&mutex);
      if (gLoop)
      {
         ++gLoop;
      }
      else
      {
         break;
      }
      pthread_mutex_unlock(&mutex);
   }

   fprintf(stdout, "Got out of loop - exiting.\n");

   /* wait for */
   pthread_join(sigthread, NULL);

   return error;
}

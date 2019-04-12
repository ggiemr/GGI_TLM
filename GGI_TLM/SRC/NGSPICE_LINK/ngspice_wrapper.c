#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER
#include <stdbool.h>
#include <pthread.h>
#else
#define bool int
#define true 1
#define false 0
#define strdup _strdup
#endif
#include <signal.h>
#if defined(__MINGW32__) || defined(_MSC_VER)
#include <windows.h>  // for Sleep
#else
#include <unistd.h>
#include <ctype.h>
#endif
#include "./sharedspice.h"

bool no_bg = true;
int vecgetnumber = 0;
int vecgettime = 0;
double v2dat;
double simtime;

bool initdata_done = false;
int    n_nodes;
int    v_present[100][2];
double v_dat[100];

static bool has_break = false;
void alterp(int sig);
static bool errorflag = false;
static bool errorflagcjs = false;

#ifndef _MSC_VER
pthread_t mainthread;
#endif // _MSC_VER

int
ng_getchar(char* outputreturn, int ident, void* userdata);

int
ng_getstat(char* outputreturn, int ident, void* userdata);

int
ng_exit(int, bool, bool, int ident, void*);

int
ng_thread_runs(bool noruns, int ident, void* userdata);

int
ng_initdata(pvecinfoall intdata, int ident, void* userdata);

int
ng_data(pvecvaluesall vdata, int numvecs, int ident, void* userdata);

int
cieq(register char *p, register char *s);

int
ciprefix(const char *p, const char *s);

int ngspice_wrapper_load_circuit(  bool *error  )
{
    int ret, i;
    char **circarray;
    
#ifndef _MSC_VER
    mainthread = pthread_self();
#endif // _MSC_VER    

    printf("**  ngspice init          **\n");
        
    ret = ngSpice_Init(ng_getchar, ng_getstat, ng_exit,  ng_data, ng_initdata, ng_thread_runs, NULL);
    
    printf("Init thread returned: %d\n", ret);

    printf("**  ngspice load circuit  **\n");
    
    ret = ngSpice_Command("source ./Spice_circuit.cir");
    
#ifndef _MSC_VER
    mainthread = pthread_self();
#endif // _MSC_VER    
            
    *error=errorflagcjs;

    printf("\n** Finished: ngspice loadcircuit. Error flag: %d  **\n",errorflagcjs);
        
    return ret;
}


int ngspice_wrapper_run( double *tout, double *Vout)
{
    int ret, i;
    
/*    char *command;*/   
    
    printf("**  ngspice run  **\n");
     
    ret = ngSpice_Command("bg_run");

    /* wait until simulation stops */
    for (;;) {
        usleep (100000);
        if (no_bg)
            break;
    }
        
    *tout=simtime; 
    
    *Vout=v2dat;
    
    printf("\n** time = %f  voltage= %f **\n\n", *tout, *Vout);
        
    return ret;
}

int ngspice_wrapper_resume( double *tout, double *Vout)
{
    int ret, i;
    
/*    char *command;*/   
    
    printf("**  ngspice resume  **\n");
    
    ret = ngSpice_Command("bg_resume");

    /* wait until simulation stops */
    for (;;) {
        usleep (100000);
        if (no_bg)
            break;
    }
        
    *tout=simtime; 
    
    *Vout=v2dat;
    
    printf("\n** time = %f  voltage= %f **\n\n", *tout, *Vout);
        
    return ret;
}

int ngspice_wrapper_step( double *tout, int *n_nodes_in, int *node_list_in, double *V_array_in)

{
    int ret, i;
    int inode;
    int node_number;
    int vec_get_vnode;
    
    int *pti;
    double *ptd;
        
/*    printf("**  ngspice step  **\n"); */
    
    ret = ngSpice_Command("step 2");
    
    no_bg=false; 

    /* wait until simulation stops */
    for (;;) {
        usleep (10);
        if (no_bg)
            break;
    }
        
    *tout=simtime; 
    
    *n_nodes_in=n_nodes;
    
    pti=node_list_in;
    ptd=V_array_in;
            
    for (inode = 0; inode<n_nodes ; inode++) {
      node_number=v_present[inode][1];
      vec_get_vnode=v_present[inode][2];
          
      *pti=node_number; 
      pti++;
      
      *ptd=v_dat[inode]; 
      ptd++;
    }
    
/*    printf("** time = %f  voltage= %f **\n", *tout, *Vout);  */
        
    return ret;
}

int ngspice_wrapper_run_to_breakpoint( double *tbreak, double *tout, int *n_nodes_in, int *node_list_in, double *V_array_in)

{
    int ret, i;
    int inode;
    int node_number;
    int vec_get_vnode;
    
    bool bret;
    
    int *pti;
    double *ptd;
        
/*    printf("** ngspice_wrapper_run_to_breakpoint  **\n");  */
/*    printf("Set breakpoint at time %e\n", *tbreak);  */
/*    printf("Resume the Ngspice run \n");  */

/*  RESUME COMMAND COULD BE ISSUED IN SUBROUTINE update_ngspice */
     ret = ngSpice_Command("bg_resume"); 
    
    no_bg=false; 
    
    /* wait until simulation stops */
    for (;;) {
         usleep (10);  
       if (no_bg)
            break;
    }
        
    *tout=simtime; 
        
    *n_nodes_in=n_nodes;
    
    pti=node_list_in;
    ptd=V_array_in;
            
    for (inode = 0; inode<n_nodes ; inode++) {
      node_number=v_present[inode][1];
      vec_get_vnode=v_present[inode][2];
          
      *pti=node_number; 
      pti++;
      
      *ptd=v_dat[inode]; 
      ptd++;
    }
    
 /*    printf("** Ngspice halted: time = %e  voltage= %e **\n", *tout, *Vout);   */
        
    return ret;
}

/* Callback function called from bg thread in ngspice to transfer
any string created by printf or puts. Output to stdout in ngspice is
preceded by token stdout, same with stderr.*/
int
ng_getchar(char* outputreturn, int ident, void* userdata)
{
  /*   printf("%s\n", outputreturn);  */
     
    /* setting a flag if an error message occurred */
    
    if (ciprefix("stderr Error:", outputreturn)){
        errorflag = true;
        errorflagcjs = true;}
        
    if (ciprefix("stdout Error", outputreturn))
         errorflagcjs = true;
    return 0;
}


int
ng_getstat(char* outputreturn, int ident, void* userdata)
{
/*    printf("%s\n", outputreturn); */
    return 0;
}

int
ng_thread_runs(bool noruns, int ident, void* userdata)
{
    no_bg = noruns;
/*    if (noruns)
        printf("bg not running\n");
    else
        printf("bg running\n");  
*/
    return 0;
}

/* Callback function called from bg thread in ngspice once per accepted data point */
int
ng_data(pvecvaluesall vdata, int numvecs, int ident, void* userdata)
{
    int *ret;
    int inode;
    int node_number;
    int vec_get_vnode;

    simtime = vdata->vecsa[vecgettime]->creal;
    v2dat = vdata->vecsa[vecgetnumber]->creal;
    
    for (inode = 0; inode<n_nodes ; inode++) {
      node_number=v_present[inode][1];
      vec_get_vnode=v_present[inode][2];
      v_dat[inode]=vdata->vecsa[vec_get_vnode]->creal;
      
/* printf(" inode %d, node_number %d, vec_get_vnode %d, v_dat(inode), %e \n",inode , node_number , vec_get_vnode , v_dat[inode]); */

    }
        
    return 0;
}


/* Callback function called from bg thread in ngspice once upon intialization
   of the simulation vectors)*/
int
ng_initdata(pvecinfoall intdata, int ident, void* userdata)
{
#define MAX_LENGTH_OF_NUMBER 9
    int i;
    int inode;
    int vn = intdata->veccount;
    
    char *br_string;
    char *right_string;
    // number has 10 digits plus \0
    char chnumber[MAX_LENGTH_OF_NUMBER + 1];
    int node_number;
    
    if (initdata_done) return 0;

/*    printf("\n");
    printf("****************************\n");
    printf("**  ngspice initdata  **\n");
    printf("****************************\n"); */

/* Reset the arrays which hold node data to be passed between ngspice and GGI_TLM these are nodes less than 100 */    
    for (inode = 0; inode<100 ; inode++) {
      v_present[inode][1]=0;
      v_present[inode][2]=0;
      v_dat[inode]=0.0;
    }
    
    n_nodes=0;
    
    for (i = 0; i < vn; i++) {

     /*   printf("Vector: %s\n", intdata->vecs[i]->vecname);  */
        
        /* look for a sub-string between brackets - this is a node number  */
        
        br_string  = strchr (intdata->vecs[i]->vecname, '(');
        right_string  = strchr (intdata->vecs[i]->vecname, ')');

        int ich = 0;
        br_string++;

            while (br_string < right_string && ich <= MAX_LENGTH_OF_NUMBER) {
                chnumber[ich] = *br_string;
                ich++;
                br_string++;
            }

        // Add a NULL to the end of the string
        chnumber[ich] = '\0';
        node_number=atoi(chnumber);
    
        if ( node_number > 0 && node_number < 101 ) {
      /*     printf("Found output at node: %d, vector element is %d \n", node_number,i); */
          v_present[n_nodes][1]=node_number;
          v_present[n_nodes][2]=i;
          n_nodes++;
        }
        
        /* find the location of V(1) */
        
        if (cieq(intdata->vecs[i]->vecname, "V(1)")) {
            vecgetnumber = i;
        /*       printf("Found node 1 output, vector element is %d \n", i); */
          }
        /* find the location of the time data */
        
        if (cieq(intdata->vecs[i]->vecname, "time")) {
       /*         printf("Found time output, vector element is %d \n", i); */
            vecgettime = i;
          }
    }
    
    initdata_done = true;
    
    return 0;
}


/* Callback function called from bg thread in ngspice if fcn controlled_exit()
   is hit. Do not exit, but unload ngspice. */
int
ng_exit(int exitstatus, bool immediate, bool quitexit, int ident, void* userdata)
{

    if(quitexit) {
        printf("DNote: Returned form quit with exit status %d\n", exitstatus);
        exit(exitstatus);
    }
    if(immediate) {
        printf("DNote: Unloading ngspice inmmediately is not possible\n");
        printf("DNote: Can we recover?\n");
    }

    else {
        printf("DNote: Unloading ngspice is not possible\n");
        printf("DNote: Can we recover? Send 'quit' command to ngspice.\n");
        errorflag = true;
        ngSpice_Command("quit 5");
//        raise(SIGINT);
    }

    return exitstatus;
}

/* Funcion called from main thread upon receiving signal SIGTERM */
void
alterp(int sig) {
    ngSpice_Command("bg_halt");
}


/* Case insensitive str eq. */
/* Like strcasecmp( ) XXX */
int
cieq(register char *p, register char *s)
{
    while (*p) {
        if ((isupper(*p) ? tolower(*p) : *p) !=
            (isupper(*s) ? tolower(*s) : *s))
            return(false);
        p++;
        s++;
    }
    return (*s ? false : true);
}

/* Case insensitive prefix. */
int
ciprefix(const char *p, const char *s)
{
    while (*p) {
        if ((isupper(*p) ? tolower(*p) : *p) !=
            (isupper(*s) ? tolower(*s) : *s))
            return(false);
        p++;
        s++;
    }
    return (true);
}


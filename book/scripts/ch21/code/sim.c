// Independent simulator of Mira's 17 state beeping Collatz machine.
// Usage: ./sim trace NSTEPS      -> run from blank tape, report beeps and checkpoints
//        ./sim macros NMAX VMAX  -> check the three macro rules from C(n,v,K)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *TM = "1RL0RB_1LC1RB_0LK0LD_1RE1LD_0RF0RB_0LG0LG_1LH1LG_0RA0RI_1RJ1RI_1RF1RF_0RA1LK_0LM0LM_0LN1LM_1RO1LN_0RA0RP_0RQ1RP_1RL1RQ";
static int wr[17][2], mv[17][2], nx[17][2];
#define BEEP 14 /* state O */

static void parse(void){
  int s=0; const char *p=TM;
  for(s=0;s<17;s++){
    for(int r=0;r<2;r++){
      wr[s][r]=p[0]-'0'; mv[s][r]=(p[1]=='R')?1:-1; nx[s][r]=p[2]-'A'; p+=3;
    }
    if(*p=='_') p++;
  }
}

#define OFF 400000
#define TAPE (2*OFF)
static unsigned char tape[TAPE];

int main(int argc,char**argv){
  parse();
  if(argc>=3 && strcmp(argv[1],"trace")==0){
    long long N=atoll(argv[2]);
    memset(tape,0,sizeof tape);
    long long pos=OFF; int st=0; long long beeps=0;
    long long lastbeep=-1;
    // checkpoints: record time of each beep for reconstruction
    for(long long t=0;t<N;t++){
      int r=tape[pos];
      if(st==BEEP){ beeps++; lastbeep=t+1; }
      tape[pos]=wr[st][r]; int ns=nx[st][r]; pos+=mv[st][r]; st=ns;
      if(pos<10||pos>TAPE-10){fprintf(stderr,"tape bound\n");return 1;}
    }
    printf("steps %lld beeps %lld last_beep_time %lld state %c head %lld\n",N,beeps,lastbeep,'A'+st,pos-OFF);
    // decode current configuration if in state A at head 0
    return 0;
  }
  if(argc>=4 && strcmp(argv[1],"macros")==0){
    int NMAX=atoi(argv[2]), VMAX=atoi(argv[3]);
    long long bad=0, cases=0;
    for(int n=1;n<=NMAX;n++) for(int v=1;v<=VMAX;v++){
      for(int K=v; K<=v+3; K++){
        memset(tape,0,sizeof tape);
        for(int c=-n;c<=-2;c++) tape[OFF+c]=1;
        for(int c=0;c<K;c++) tape[OFF+c]=1;
        tape[OFF+v-1]=0;
        long long pos=OFF; int st=0; long long steps=0, beeps=0;
        do{
          int r=tape[pos];
          if(st==BEEP) beeps++;
          tape[pos]=wr[st][r]; int ns=nx[st][r]; pos+=mv[st][r]; st=ns; steps++;
          if(steps>100000000LL){bad++;break;}
        }while(!(st==0 && pos==OFF));
        // expected
        int en, ev, eK; long long esteps, ebeeps;
        if(v==1){ en=n+1; ev=n+1; eK=(K>n+1?K:n+1); esteps=2LL*n*(n+2); ebeeps=n; }
        else if(v%2==0){ int k=v/2; en=n; ev=k; eK=K; esteps=2LL*k*(k+1); ebeeps=0; }
        else { int k=(v-1)/2; en=n; ev=3*k+2; eK=(K>3*k+2?K:3*k+2); esteps=5LL*k*k+9*k+4; ebeeps=0; }
        // check tape equals C(en,ev,eK)
        int ok=1;
        for(int c=-OFF+2;c<OFF-2;c++){
          int want=0;
          if(c>=-en && c<=-2) want=1;
          else if(c>=0 && c<eK && c!=ev-1) want=1;
          if(tape[OFF+c]!=want){ok=0;break;}
        }
        if(!ok||steps!=esteps||beeps!=ebeeps){bad++; if(bad<5) printf("FAIL n=%d v=%d K=%d steps=%lld exp=%lld beeps=%lld\n",n,v,K,steps,esteps,beeps);}
        cases++;
      }
    }
    printf("macro cases %lld failures %lld\n",cases,bad);
    return 0;
  }
  return 1;
}

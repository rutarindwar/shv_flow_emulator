#include <stdio.h>
#include <stdlib.h>

void usage(char *prog){
   printf("usage: %s [start_idx] [end_idx]\n", prog);
   printf("\twhere 0 < start_idx < end_idx\n");
   exit(-1);
}

int min2(int a, int b) {
   return a < b ? a : b;
}

void runit(int start_idx, int end_idx) {
   char fname[128];
   sprintf(fname, "job-%d-%d.qsub", start_idx, end_idx);

   //create the job submit file
   FILE *fp = fopen(fname, "w+");
   fprintf(fp, "#!/bin/sh\n");
   fprintf(fp,"#PBS -N auto5\n");
   fprintf(fp,"#PBS -l nodes=1:ppn=1,walltime=12:00:00\n");
   fprintf(fp,"#PBS -q batch\n");
   fprintf(fp,"#PBS -j oe\n");
   fprintf(fp,"#PBS -o log-$PBS_JOBNAME-$PBS_JOBID\n");
   fprintf(fp,"\n");
   fprintf(fp,"cd $PBS_O_WORKDIR\n");
   fprintf(fp,"for id in `seq %d %d`;\n", start_idx, end_idx);
   fprintf(fp,"do  \n");
   fprintf(fp,"  ./designAutoNr $id;\n");
   fprintf(fp,"done\n");
   fprintf(fp,"\n");
   fclose(fp);


   //submit the job
   char cmd[256];
   sprintf(cmd, "qsub %s", fname);
   system(cmd);

}

int main(int argc, char *argv[]) {
   int i;
   for (i=0; i<argc; i++) {
      printf("argv[%d] = %s\n", i, argv[i]);
   }
   
   if (argc < 3) {
       usage(argv[0]);
   }

   int start_idx = atoi(argv[1]);
   int end_idx = atoi(argv[2]);
   if (end_idx < start_idx || start_idx <= 0) {
       usage(argv[0]);
   }

   int block_size = end_idx - start_idx + 1;
   if (argc >=3 ) {
       block_size = atoi(argv[3]);
   }

   int id = start_idx;
   for (; id<=end_idx; id+=block_size) {
      int block_start_id = id;
      int block_end_id = min2(id + block_size - 1, end_idx);
      runit(block_start_id, block_end_id);
   }

   return 0;
}

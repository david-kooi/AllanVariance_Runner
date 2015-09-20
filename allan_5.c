/* Allan variance, read a data file. Use function contiguous to get the data back in the*/
/* time sequence. First find out how many data points there are, and what the time */
/* interval is for data point. Finally call function allan in a loop and then */
/* write the allan variance out to a file. */
/* The noise fluctuation BW is automatically calculated. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define PI 3.14159265

/*  *ntp is a pointer to the address. Expect nothing back, just pointer to memory locations. */
void read_file1(float * t1, float * tp1, float *to_tp, char *readfile1, int *ntp);
void contiguous(float * t1, float * tp1, float to_tp, int ntp, int intcount, float *ta, float * s, float * x, int *imax);
void allan(float * s, float * x, int imax, float *variance);
void write_file(float * Tint, float * var, float * var_norm, float * rad_norm, float * sigma_norm, char *writefile1, int intcount);

main(int argc, char *argv[])  
{
         //Setup variables

         int i, k, j, ntp, intcount, imax, istart, second, choice;
         float tpsum, * t1, * tp1, * s, * x;
         float to_tp, ta, accum, mean, variance, bfi, Tadev, Tamin, accum1, accum2, accumbeta;
         float beta, betanom, betadenom, betaprim, points, ti, xrange;
         float * Tint, * var, * var_norm, * sigma_norm, * rad_norm;
      	 char* readfile;
         char*  writefile1;

         t1 = (float *) malloc( sizeof(float) * 190000 );
         tp1 = (float *) malloc( sizeof(float) * 190000 );
         s = (float *) malloc( sizeof(float) * 100000 );
         x = (float *) malloc( sizeof(float) * 100000 );
         Tint = (float *) malloc( sizeof(float) * 100000 );
         var = (float *) malloc( sizeof(float) * 100000 );
         var_norm = (float *) malloc( sizeof(float) * 100000 );
         sigma_norm = (float *) malloc( sizeof(float) * 100000 );
         rad_norm = (float *) malloc( sizeof(float) * 100000 );

       //Process Arguments
         //Check for proper usage
         if (argc != 4){
            printf("Usage: \n");
            printf("run <input filename> <Num nodes> <outputfile>\n");
            return 0;
         }


         //Setup input variables
         readfile = argv[1];
         points = (float)atoi(argv[2]); //Convert string argument to int then to float
         writefile1 = argv[3];

         printf("InputFile: %s\nSamples/Second: %lf\nOutputFile: %s\n\n", readfile, points, writefile1);

         //Check if input file exists
         if ( access(readfile, F_OK) == -1){
           printf("Input File Does not exist\n");
           return 0;
         }

           /* readfile = (char *) malloc( sizeof(char) * 50 ); */
           /* writefile1 = (char *) malloc( sizeof(char) * 50 ); */

	 xrange=20.0; /* for allan plot*/


	/*system("clear);*/
        for(i=1;i<=100;i++){
           Tint[i]=0.0, var[i]=0.0, var_norm[i]=0.0, sigma_norm[i]=0.0, rad_norm[i]=0.0;
	}

/*	printf("\nEnter data file name (exmpl: test.dat)");
	printf("\nData structure: Time, Tp (volts)"); 
        printf("\n\nNote: Contiguous data samples are clipped at 5 sigma from mean");
	printf("\n \nFile name? ");
	scanf("%s",&readfile); */           /*read in filename*/
	                                  /*printf("%s \n", readfile1 );*/
	                                  /*printf("%d \t %d \n", i, nfiles);*/

/*	printf("\nNumber of data points / Second (rate data is acquired)? ");
        scanf("%f",&points); */           /*integration points per second*/
        ti=1.0/points;                  /*integration time per sample of data in HP powerhead*/
        printf("\nTheoretical Sample Time: %5.4f mS\n",ti*1000.0);

        /* call function read_file to aquire data in arrays t1, tp1 and return pointers*/
        /* arrays already point to memory, so no pointers or brackets needed */ 
	read_file1(t1, tp1, &to_tp, readfile, &ntp);      
	   printf("\nNumber of data points: %d", ntp);
	    
         /*find the <s(t)>, so that the allan variance can be normalized. */
	 accum=0.0;  
           for(i=1;i<=ntp;i++){
	    /*printf("%d %8.7f %8.7f \n", i, tp1[i], t1[i]);*/
   		accum=accum+tp1[i];
	     }
	 mean=accum/ntp;
//	 printf("\nOriginal data Voltage Mean: s(t)= %8.7f (mV)", mean*1000.0);    		
	   /*print the data*/
	   /* for (i=1;i<=ntp;i++){*/
   	   /* printf("%d %4.3f %3.2f %3.2f %d \n",i,t1[i],tp1[i],to_tp,ntp);*/
	   /*  }	*/

           /*call function allan_var to aquire data in the allen time (ta) and arrays s[] and x[]*/
  	   /*need to pass integration time counter (intcount)*/
	   /*loop until imax =1, no more data points left over*/ 
	   /*error bars increase as less data points are left over*/
	
	intcount=1, imax=2;	/*intcount sets how many contiguous samples*/
	while (imax>1){ 	/*should be 1 in the end I think*/  	 	
	   contiguous(t1, tp1, to_tp, ntp, intcount, &ta, s, x, &imax);
	     /*printf("\nintcount sets the number of contiguous samples in s[],x[]");*/ 
             /*printf("\nimax,intcount,Tint: %d %d %4.3f\n",imax,intcount,ta);*/
	      /*for (i=1;i<=imax;i++){ */
		/*printf("\ni, s[i], x[i]: %d %4.3e %4.3e",i,s[i],x[i]);*/
   		 /*} */ 
	 	  allan(s, x, imax, &variance);	 /*function allan, now calculate the allan variance, var*/
		    /*now normalize*/
		    var[intcount]=variance, Tint[intcount]=ta;  /*ta=to_tp*i, then allan sampling time*/
		    var_norm[intcount]=var[intcount]/pow(mean,2);      /*normalized variance*/
		    sigma_norm[intcount]=pow(var_norm[intcount],0.5);  /*normalized sigma (variance=sigme^2)*/
		  /*printf("\n\nTint, Allan var %5.4f %5.4e ",Tint[intcount],var[intcount]);*/
		  /*printf("\nvar_norm, sigma_norm %5.4e %5.4e \n",var_norm[intcount],sigma_norm[intcount]);*/
	   intcount=intcount+1;         
       /*printf("\n------------------------------------------------------------------------------------------\n");*/

        }	
	   /*calculate the fluctuation (noise bw)*/
	   /*if i=1, the Tint can be zero, thus use i=2, could check later but does not matter really*/
  		/* For the noise we should use the integration time, not the sample time*/
		/* However Tint[] contains the sample time, thus correct by ti/tsample*/
		/* factor 1 needed according to Schieder, S, R half time*/
		bfi=0.0;
		for(i=1;i<=10;i++){
                  bfi=bfi+pow(mean,2)/((var[i]*Tint[i]*(ti/to_tp)));
		  /* printf("\ni=%d, Noise fluxuation BW = %6.5f MHz, Mean = %6.5f mV\n",i, bfi/1E6,mean*1000);*/
                }
	        bfi=bfi/10.0; /*take the mean of 10 bfi*/
  	        /*printf("\nTint[2], Var[2] %6.5e %6.5e\n",Tint[2],var[2]);*/
	
	/*print results and caculate the ideal (mean) normalized radiometer equation values */
	  /*printf("\n i   Tint   Allan Var       Var_norm         Rad_norm       Sigma_norm");*/
	    for(i=1;i<=intcount-2;i++){  /*last data point in NaN*/
		/* Ideal Radiometer variance, corrected for integration (ti) vs sampling (to_tp) time*/ 
	  	rad_norm[i]=1/(1.0*bfi*Tint[i]*(ti/to_tp)); 
	      /*printf("\n %d \t%5.4f \t%5.4e \t%5.4e \t%5.4e \t%5.4e",  i,Tint[i],var[i],var_norm[i],rad_norm[i],sigma_norm[i]);*/
   	    }

        /* calculate when the allan variance deviates more than 20% from radiometer equation*/
        /* note that sigma is root(var), thus 20% variation of sigma is 44% variation in variance*/
	/* bin the number as to average the noise a bit more*/
          for (i=10;i<=intcount-10;i++){
             if (((var_norm[i-3]+var_norm[i-2]+var_norm[i-1]+var_norm[i]+var_norm[i+1]+var_norm[i+2]+var_norm[i+3])/7) > 1.44*rad_norm[i]){
		Tadev=Tint[i];
	        break;
		}
             else
	 	Tadev=0.0, istart=0;
 	  }

	/*caculate the minima Tamin*/
	  for (i=30;i<=intcount-30;i++){
	    accum1=0.0,accum2=0.0;
	      for (k=i-10;k<=i+10;k++){ /* +- 0.5 seconds, this seems best to average out the noise*/
		accum1=accum1+var_norm[k];	 	
		accum2=accum2+var_norm[k+21]; /*add 2*delta+1*/
		/*printf("\ni, k, accum1, accum2: %d %d %e %e", i,k,accum1,accum2);*/
	      }
            if (accum2>accum1){
      		Tamin=Tint[i];
		accum1=0.0,accum2=0.0;
		istart=i; 		/*for caculating beta*/
	  	break;
		}
	      else{
		accum1=0.0,accum2=0.0;
	      }   
           } 
	
	/*calculate beta(T) for T > Tadev. Find mean <beta> */
	   accumbeta=0.0, beta=1.0, betaprim=0.0;
	for(j=1;j<=50;j++){ /*did not use while loop, need to be fail save (<1%)*/
	   for (i=istart+20;i<=istart+120;i++){ /*calculate beta from (Tamin+1) to (Tamin+7) seconds integration time T*/
	      betanom=log(fabs(beta*(var_norm[i]*1.0*bfi*Tamin-(Tamin/Tint[i])))); /*factor 1 needed according to Schieder, S, R half time*/ 
	      betadenom=log(fabs(Tint[i]/Tamin));
	      accumbeta=accumbeta+(betanom/betadenom); /*<need to get <beta>*/
	        /*printf("\ni, T, var_norm, betanom, betadenom %d %f %e %f %f", i, Tint[i], var_norm[i], betanom, betadenom);*/
	   }
	   betaprim=fabs(accumbeta/(intcount-2-istart));
             /*printf("\nbeta, betaprim , delta %4.3f %4.3f %4.3f",beta, betaprim, betaprim-beta);*/
	     beta=betaprim;		
	}
            
		second=(int)1/to_tp;		  /*1 second, force to int*/
      /*exit(1); /*terminates program at this spot*/

   /*now write the data to a file*/
       // printf("\n\nenter file name? (exmpl: allan.dat) :");
       // scanf("%s",&writefile1);
         
        /*call function write_file to save data to disk*/
        write_file(Tint, var, var_norm, rad_norm, sigma_norm, writefile1, intcount);
            

		free(t1);
		free(tp1);
		free(s);
		free(x);
		free(Tint);
		free(var);
		free(var_norm);
		free(sigma_norm);
		free(rad_norm);

       	return 0;
}

void read_file1(float t1[], float tp1[], float *to_tp, char *readfile, int *ntp)
/* *ntp point to the address where this value is located */
/* & is used to assign the address to a variable (not an array)*/
/* variable in function are local, thus can have different names */
/* read data from a file using function read_file*/
/* specify %f, not %6.3f etc, & needed in scanf.*/
{  

        int i, imax, n, itot;	
	char  tmp[50], exten[5];
	float to, tsample,count;                        
	FILE *pread;   /*declare pread a pointer to type-name FILE, so pread is*/
                       /*pointing to FILE, which holds info obtained from fopen*/
	
	
	    /*sprintf(tmp, "%s.%.3d",readfile1, n);  /* this puts everything in string tmp*/
	    /*strcat(readfile1, "asdf");*/   /*this connects strings together*/
            /*sprintf(tmp, "%s.%s",readfile1, exten); */
            /*printf("tmp: %s\n", tmp);*/

	    /* Open file, testing for success, remember that fp is a pointer to FILE */
	    if ((pread = fopen(readfile, "r")) == NULL) { /*fopen returns a pointer (fp)*/
	       printf("Error opening text file for reading\n");   /*which points to FILE*/
	       exit(1);     /*do not use return here but exit(), that contains file info*/  
	   }  


	 i=1;
	 while ((fscanf(pread, "%f %f", &t1[i], &tp1[i]) != EOF)){
	   /* printf("%d %8.3f Seconds %8.3f dBm \n", i, t1[i], tp1[i]);*/ 
	       /*tp1[i]=pow(10.0,(tp1[i]/20.0))/1000.0; /*tp1 in dBm, want Voltage */
               /*tp1[i]=pow((tp1[i]*50.0),0.5); /*in Watt, want Voltage */
               /*t1[i]=t1[i]/1000.0;            /*in mS, should be in seconds*/
           i++;
	}
         imax=i-1;

	/*for(i=1;i<=20;i++){*/
	/*  printf("%d %8.3f Seconds %8.7f V \n", i, t1[i], tp1[i]);*/
        /*  } */

        *ntp=imax;  /* ntp point to the memory location where value is stored */
        *to_tp=(t1[20]-t1[10])/10.0; /*means Sample time*/
        printf("\nActual sample time: %5.4f mS", *to_tp*1000.0);
        /*printf("ntp: %d\n", *ntp);*/
            
        fclose(pread);  /*close FILE*/
}

void contiguous(float t1[], float tp1[], float to_tp, int ntp, int intcount, float *ta, float s[], float x[], int *imax)
{
  /*In this function calculate two contiguous data samples from data array tp1*/
  /*return the contigous data values in arrays s, x. The integration time is ta */
  /*the number of data points in the arrays is imax*/
  /*intcount sets the number of contiguous data samples*/

  int i, k, counter;    
  float saccum, xaccum;  
  counter=1, saccum=0.0, xaccum=0.0;

  for (i=1;i<=ntp;i++){
        /*printf("\nallan (time, tp): %d %4.3f %3.2f",i,t1[i],tp1[i]);*/
          }
	/*printf("\ntime interval, number of points: %3.2f %d\n",to_tp,ntp);*/
	 /*intcount sets how many contiguous samples*/
	for (i=1;i<=ntp-2*intcount;i=i+2*intcount){
	     for (k=1;k<=intcount;k++){
		saccum=saccum+tp1[i+k-1];
               	xaccum=xaccum+tp1[i+k+intcount-1];
	/*printf("\ni, kloop, intcount, s(i+k-1), x(i+k+intcount-1): %d %d %d %d %d",i,k,intcount,i+k-1,i+k+intcount-1);*/
	     }
	s[counter]=saccum/(k-1);
   	x[counter]=xaccum/(k-1);
        /*printf("\niloop, k, intcount,saccum,xaccum: %d %d %d %4.3f %4.3f ",i,k,intcount,saccum,xaccum);*/
	/*printf("\niloop, k, counter,s,x           : %d %d %d %4.3f %4.3f ",i,k,counter,s[counter],x[counter]);*/
	saccum=0, xaccum=0;
	counter=counter+1;
   }
   *imax=counter-1;
   *ta=to_tp*intcount;
   /*printf("\ncounter-1= %d",counter-1);*/
}

void allan(float s[], float x[], int imax, float *variance)
{
  int i;
  float accum_dif, accum_s, accum_x, mean_x, mean_s, sum_x, sum_s;
  float delta, accum, mean_d, sigma_s, sigma_x;

  accum_dif=0.0, accum_s=0.0, accum_x=0.0, sum_x=0.0, sum_s=0.0;

  /* first find mean and sigma for contiguous data samples s[i] and x[i]*/
  /* then braket s[i], and x[i] with mean +- 2Sigma to avoid that the mean*/
  /* of the difference gets blown up by some noise spike in the data set*/
    
    /*calculate the mean*/ 
    for (i=1;i<=imax;i++){
      accum_s=accum_s+s[i];
      accum_x=accum_x+x[i];
     }          
      mean_s=accum_s/imax;
      mean_x=accum_x/imax;

    /*calculate the standard deviation (1sigma=65%, 2sigma=95%)*/
    for (i=1;i<=imax;i++){
      sum_s=sum_s+pow((s[i]-mean_s),2);
      sum_x=sum_x+pow((x[i]-mean_x),2);
    }
      sigma_s=pow((sum_s/(imax-1)),0.5);
      sigma_x=pow((sum_x/(imax-1)),0.5);
   
    /*printf("\nmean_s, mean_x, sigma_s, sigma_x = %6.5e %6.5e %6.5e %6.5e",mean_s,mean_x,sigma_s,sigma_x);*/
    /*printf("\nimax= %d",imax);*/
   
    /* now check if s[i] or x[i] exceed mean +- 5sigma, if so use the mean of s[i] or x[i]*/
    /* 3 sigma is not enough, 5 sigma certainly looks good. so 4-5 sigma is ok*/
    for (i=1;i<=imax;i++){
	/*printf("\ni, s = %d %6.5e",i,s[i]);*/
        if (s[i]>(mean_s+sigma_s*5.0) || s[i]<(mean_s-sigma_s*5.0)){ 
          s[i]=mean_s;
	}
        /*printf("\ni, x = %d %6.5e",i,x[i]);*/
        if (x[i]>(mean_x+sigma_x*5.0) || x[i]<(mean_x-sigma_x*5.0)){
          x[i]=mean_x;
        }
     }
     
    /*calculate the mean of the difference with corrected contiguous data samples*/
       for (i=1;i<=imax;i++){
     	/*printf("\nin allan: i, s[i], x[i]= %d %6.5e %6.5e ",i,s[i],x[i]);*/
	delta=(s[i]-x[i]);
	accum_dif=accum_dif+delta;
	 /*printf("\n in allan: i, s, x, delta, accum_dif = %d %6.5e %6.5e %6.5e %6.5e",i,s[i],x[i],delta,accum_dif);*/
	}
    mean_d=accum_dif/imax;
    /*printf("\nmean_d = %5.4e\n",mean_d);  /*mean of the difference between contiguous data samples*/	

    /*now calculate the allan variance; var_a and return via a pointer*/
    accum=0.0;
      for (i=1;i<=imax;i++){
        delta=(s[i]-x[i]);
        accum=accum+pow((delta-mean_d),2);
	/*printf("\nAllan delta, Allan accum: %6.5e %6.5e", delta, accum);*/ 
       }
      *variance=accum/(2.0*(imax-1));
     /*printf("\nAllan variance = %8.7e",*variance); /*variance is sigma^2*/ 
}

void  write_file(float Tint[], float var[], float var_norm[], float rad_norm[], float sigma_norm[], char *writefile1, int intcount)
{
 int count,sample;
        int i;
        FILE *pwrite;  /*declare pread a pointer to type-name FILE, so pread is*/
                       /*pointing to FILE, which holds info obtained from fopen*/
          
        /* Open file, testing for success, remember that fp is a pointer to FILE */
          
           if ((pwrite = fopen(writefile1, "a")) == NULL) { /*fopen returns a pointer (fp)*/
               printf("Error opening text file for writing\n");   /*which points to FILE*/
               exit(1);     /*do not use return here but exit(), that contains file info*/
           }

           /*fprintf(pwrite, " Tint     Var     Var_norm    Rad_norm  Sigma_norm\n");*/

           for (i=1;i<=intcount-2;i++) { /*loop number of sections. Substract two, NaN*/
            /*fprintf(pwrite, "%5.4f %5.4e %5.4e %5.4e %5.4e \n",  Tint[i], var[i], var_norm[i], pow(rad_norm[i],0.5), sigma_norm[i]);*/
	    fprintf(pwrite, "%5.4f %5.4e %5.4e \n",  Tint[i], pow(rad_norm[i],0.5), sigma_norm[i]);
           }
        fclose(pwrite);  /*close FILE*/

}

void gnuplotlog( float *x, float *y1, float *y2, int length, float xrange) {
  char foo[80];
  int  i;

  FILE* pipeToPlot;
  // Run gnuplot
  pipeToPlot = popen("gnuplot", "w");
  
  // A standard plotting command, with the "-" specifying that
  // the data points will be specified on this command line.
  fprintf(pipeToPlot, "set logscale\n");
  fprintf(pipeToPlot, "plot \"-\" with linespoints\n");
  
  // The data points
  //fprintf(pipeToPlot, "-2 4\n");
  //fprintf(pipeToPlot, "-1 1\n");
  // The data points are terminated by a solitary "e".  Try using gnuplot like
  // this by hand so you have a feel of how it works.

  for (i=1;i<=length-2;i++) { /*loop number of sections. Substract two, NaN*/
    if (x[i]>xrange){
      break;
    }
    fprintf(pipeToPlot, "%5.4e %5.4e \n",  x[i], y1[i]);
  }

fprintf(pipeToPlot, "e\n");

  // This is necessary to force the buffer to be flushed and the pipe to
  // be written.  By default, pipes, files, and other similar entities in unix
  // *may* be buffered in such a way that they only finalize the data transfer
  // when they have to.  This is great for performance in the long run, but when
  // you need to know that data has gone through, you must manually flush the
  // buffer.
  fflush(pipeToPlot);

  // Just wait for the user to hit enter.
  printf("\nHit enter to continue.\n");
  
  fgetc(stdin);
  fgetc(stdin);
 
  // Don't forget to close your pipe like this, or there could be trouble!
  pclose(pipeToPlot);
}

void gnuplotlin( float *x, float *y, int length) {
  char foo[80];
  int  i;

  FILE* pipeToPlot;
  // Run gnuplot
  pipeToPlot = popen("gnuplot", "w");
  
  // A standard plotting command, with the "-" specifying that
  // the data points will be specified on this command line.
  
  fprintf(pipeToPlot, "plot \"-\" with linespoints\n");
  
  //The data points
  //fprintf(pipeToPlot, "-2 4\n");
  //fprintf(pipeToPlot, "-1 1\n");
  //The data points are terminated by a solitary "e".  Try using gnuplot like
  //this by hand so you have a feel of how it works.

  for (i=1;i<=length-2;i++) { /*loop number of sections. Substract two, NaN*/
       fprintf(pipeToPlot, "%5.4e %5.4e \n",  x[i], y[i]);
  }

fprintf(pipeToPlot, "e\n");

  // This is necessary to force the buffer to be flushed and the pipe to
  // be written.  By default, pipes, files, and other similar entities in unix
  // *may* be buffered in such a way that they only finalize the data transfer
  // when they have to.  This is great for performance in the long run, but when
  // you need to know that data has gone through, you must manually flush the
  // buffer.
  fflush(pipeToPlot);

  // Just wait for the user to hit enter.
  printf("\nHit enter to continue.\n");
  
  fgetc(stdin);
  fgetc(stdin);
  
  // Don't forget to close your pipe like this, or there could be trouble!
  pclose(pipeToPlot);
}

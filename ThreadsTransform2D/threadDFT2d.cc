// Threaded two-dimensional Discrete FFT transform
// Marcus Bolden
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>
//#include <cmath>

#include "Complex.h"
#include "InputImage.h"

#include <pthread.h>

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.
int activeThreads, totalThreads, writeflag;

pthread_mutex_t coutMutex;
pthread_mutex_t arrayMutex;
pthread_mutex_t arrayMutex2;
pthread_mutex_t arrayMutex3;
pthread_mutex_t activeMutex;
pthread_cond_t allDoneCondition;
pthread_barrier_t barrier; // barrier synchronization object
pthread_barrier_t barrier2;
pthread_barrier_t barrier3;
pthread_barrier_t barrier4;

unsigned N;
Complex* W;
Complex* h;
Complex* h2;
Complex* h3;
InputImage *p;

using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

//function to precompute W's
void ComputeW(Complex* W){
    for(int n = 0; n < N / 2; n++){
        W[n].real = cos(2*M_PI*n / N);
        W[n].imag = -sin(2*M_PI*n / N);
    }
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier() // Again likely need parameters
{
}



void Transform1D(Complex* h0, Complex* W, int N, int start)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)

    //reorder using bit-reverse (analogous to the size 1 transform)
    Complex temp;
   
    for(int i = 0; i < N ; i++){
      //need to only switch places if i hasn't been switched
      //also, nothing needs to be done if i is the same as reverse i
      if(i != ReverseBits(i) && i < ReverseBits(i)){
        temp = h0[start+i];
        h0[start+i] = h[start+ReverseBits(i)];
        h0[start+ReverseBits(i)] = temp;
      }
    }
   
    int M = N / 2; //number of time each transform is done
    int trans = 2; //starting with 2-pt transform
    Complex t1, t2;
    int start2 = start; //need a new start so I can change the start time and keep the original

    //calculate the FFT

    while(trans <= N){
        for(int i = 0; i < M; i++){ //number of times transform is done
            for(int j = start2; j < start2 + trans/2; j++){ //transform size loop
               //DFT algorithm using two temp variables    
               t1 = h0[j];
               t2 = W[(j-start2)*N/trans] * h0[j + trans/2];    
               h0[j] = t1 + t2; 
               h0[j+trans/2] = t1 - t2;
            }
            //update indices
            start2 = start2 + trans;
        } //number of times transform is done (512) 
        trans = trans * 2;
        M = M / 2;
        start2 = start;
    } //loop until no more transforms needed (2->1024)
    
     /* Complex u,t;
      Complex ww, wm;
      int m = 1;

      for(int s = 1; s <= log2(N); s++){
        m *= 2;
        wm.real = cos(2*M_PI / m);
        wm.imag = -sin(2*M_PI / m);
        for(int k = 0; k <= N - 1; k += m){
          ww.real = 1;
          ww.imag = 0;
          for(int j = start; j <= start + m/2 - 1; j++){
            pthread_mutex_lock(&arrayMutex);
            t = ww * h[k+j+m/2];
            u = h[k + j];
            h[k + j] = u + t;
            h[k+j+m/2] = u - t;
            pthread_mutex_unlock(&arrayMutex);

            if(s==3&&k==0&&m==8){
              pthread_mutex_lock(&coutMutex);
              cout<<"ww: "<<ww<<endl;
              pthread_mutex_unlock(&coutMutex);
            }
             ww = ww * wm;  
          }
        }
      }*/
    /*if(start==0){
      for(int i = 0; i < 4; i++){
        cout << h[i] << " ";
      }
    }*/
}

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number

    unsigned long tn = (unsigned long)v;

  // Calculate 1d DFT for assigned rows
  int begin = (N / totalThreads) * tn * N;
  Complex temp2;

  //int begin = 0;
  for(int i = 0; i < N / totalThreads; i++){
    pthread_mutex_lock(&arrayMutex);
    Transform1D(h, W, N, begin);
    pthread_mutex_unlock(&arrayMutex);

    begin += N;
  } 
   
  // wait for all threads to complete at this point
    pthread_barrier_wait(&barrier);
  //save 1D image
      
      pthread_mutex_lock(&arrayMutex);
    if(tn == 0){     
      (*p).SaveImageData("MyAfter1D.txt", h, N, N);
    }
      pthread_mutex_unlock(&arrayMutex);
  //wait until save finished
    pthread_barrier_wait(&barrier);

  // Calculate 1d DFT for assigned columns

  //transpose full 1D array ///////////////////////
    pthread_mutex_lock(&arrayMutex);
  if(tn==1){
    for(int j = 0; j < N; j++){
        for(int i = 0; i < N; i++){          
            h2[i + j * N] = h[j + i * N];  
        }
    }
  }
  pthread_mutex_unlock(&arrayMutex);
  //wait for transpose to complete
  pthread_barrier_wait(&barrier);

  //calculate DFT
  begin = (N / totalThreads) * tn * N;
  
  for(int i = 0; i < N / totalThreads; i++){
    pthread_mutex_lock(&arrayMutex);
    Transform1D(h2, W, N, begin);
    pthread_mutex_unlock(&arrayMutex);
    begin += N;
  } 
   
  //wait for 1D column calculation to finish
  pthread_barrier_wait(&barrier);
  //transpose again///////////////////////////////////
  pthread_mutex_lock(&arrayMutex);
  if(tn==1){
   for(int j = 0; j < N; j++){
        for(int i = 0; i < N; i++){          
            h3[i + j * N] = h2[j + i * N];        
        }
    } 
  }
  pthread_mutex_unlock(&arrayMutex);
  //wait for transpose to finish
   pthread_barrier_wait(&barrier);
   
  // Decrement active count and signal main if all complete
    pthread_mutex_lock(&activeMutex);
    activeThreads--;
    if (activeThreads == 0)
    {
        pthread_cond_signal(&allDoneCondition);
    }
     pthread_mutex_unlock(&activeMutex);
    return 0;
}

void Transform2D(const char* inputFN)
{ // Do the 2D transform here.
    InputImage image(inputFN);  // Create the helper object for reading the image
  // Create the global pointer to the image array data
    //also need temps used during transpose
    h = image.GetImageData();
    p = &image; //used for saving within thread
    h2 = image.GetImageData();
    h3 = image.GetImageData();
    N = image.GetWidth();
    
    //precompute W's
    W = new Complex[N / 2];
    ComputeW(W);

  // Create 16 threads
    
    //lock mutex used in condition wait(cond_wait automatically unlocks->locks)
    pthread_mutex_lock(&activeMutex);
    
    activeThreads = totalThreads;
    
    //create threads
    for (int i = 0; i < totalThreads; ++i){
        pthread_t t;
        pthread_create(&t, 0, Transform2DTHread,  (void*)i);
    }

  // Wait for all threads complete
    pthread_cond_wait(&allDoneCondition, &activeMutex);

  // Write the transformed data
    //debug loop
  for(int i = 0; i < 4; i++){
    cout << h3[i] << " ";
  }
  cout<<endl;

    image.SaveImageData("MyAfter2D.txt", h3, N, N);

    pthread_mutex_unlock(&activeMutex);
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  // Threads initialization here
   
    pthread_mutex_init(&arrayMutex, 0); 
    //pthread_mutex_init(&arrayMutex2, 0); 
    //pthread_mutex_init(&arrayMutex3, 0); 
    pthread_mutex_init(&coutMutex, 0);
    pthread_mutex_init(&activeMutex, 0);
    pthread_cond_init(&allDoneCondition, 0);
    totalThreads = 16; //change back to 16
    pthread_barrier_init(&barrier, NULL, totalThreads);
    //pthread_barrier_init(&barrier2, NULL, totalThreads);
    //pthread_barrier_init(&barrier3, NULL, totalThreads);
    //pthread_barrier_init(&barrier4, NULL, totalThreads);
    
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  

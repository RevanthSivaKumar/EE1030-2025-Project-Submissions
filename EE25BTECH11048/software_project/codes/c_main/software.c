#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//for dot product
double dotproduct(double *a1,double *a2,int imagescount){
    double total=0; 
    for(int i=0;i<imagescount;i++){
        total+=a1[i]*a2[i];
    }
    return total;
}
//to get unit vector
void unitvector(double *vec,int imagescount){
    double total=sqrt(dotproduct(vec,vec,imagescount));
    if(total==0){
        total=1;
    }
    for(int i=0;i<imagescount;i++){
         vec[i]/=total;
    }
}
//to get Y=AX
void Axmatrix(double *A,double *x,double *y,int rows,int cols){
    for(int i=0;i<rows;i++){
        double sum=0;
        for(int j=0;j<cols;j++){
             sum+=A[i*cols+j]*x[j];
        }
        y[i]=sum;
    }
}
//to get Y=ATX
void ATxmatrix(double *A,double *x,double *y,int rows,int cols){
    for(int j=0;j<cols;j++){
        double sum=0;
        for(int i=0;i<rows;i++){
            sum+=A[i*cols+j]*x[i];
            }
        y[j]=sum;
    }
}
//largest singular value and vectors (rank-1 SVD) of A using 10 power iterations
void rank1(double *A,int rows,int cols,double *u,double *vec,double *total){
    int i,j;
    for(i=0;i<cols;i++){
         vec[i]=1;
    }
    unitvector(vec,cols);

    for(i=0;i<10;i++){
        double *utemp=malloc(rows*sizeof(double));
        double *vtemp=malloc(cols*sizeof(double));

        Axmatrix(A,vec,utemp,rows,cols);
        unitvector(utemp,rows);

        ATxmatrix(A,utemp,vtemp,rows,cols);
        unitvector(vtemp,cols);

        for(j=0;j<rows;j++) {
            u[j]=utemp[j];
        }
        for(j=0;j<cols;j++){
             vec[j]=vtemp[j];
        }
        free(utemp);
        free(vtemp);
    }

    double *Y=malloc(rows*sizeof(double));
    Axmatrix(A,vec,Y,rows,cols);
    *total=sqrt(dotproduct(Y,Y,rows));
    free(Y);
}
//main compression code
void compress(char *inputimage, int kvalue, char *outputimage){
    FILE *fp=fopen(inputimage,"rb");

    fseek(fp,0,SEEK_END);
    long size=ftell(fp);
    fseek(fp,0,SEEK_SET);

    unsigned char *temparory=malloc(size);
    fread(temparory,1,size,fp);
    fclose(fp);

    int width,height,cols;
    unsigned char *img=stbi_load_from_memory(temparory,size,&width,&height,&cols,1);
    free(temparory);

    cols=1;

    printf("k=%d processing...\n", kvalue);

    double *A=malloc(height*width*sizeof(double));
    double *Original=malloc(height*width*sizeof(double));

    int i,j;
    double Frobeniusnorm = 0;//initialise f norm = 0

    for(i=0;i<height;i++)
        for(j=0;j<width;j++){
            double vec = img[i*width+j];
            A[i*width+j] = vec;
            Original[i*width+j] = vec;
            Frobeniusnorm += vec*vec;
        }

    Frobeniusnorm = sqrt(Frobeniusnorm);

    double **U=malloc(kvalue*sizeof(double*));
    double **V=malloc(kvalue*sizeof(double*));
    double *S=malloc(kvalue*sizeof(double));

    for(i=0;i<kvalue;i++){
        U[i]=malloc(height*sizeof(double));
        V[i]=malloc(width*sizeof(double));

        rank1(A,height,width,U[i],V[i],&S[i]);

        int x,y;
        for(x=0;x<height;x++)
            for(y=0;y<width;y++)
                A[x*width+y]-=S[i]*U[i][x]*V[i][y];
    }

    double *Recompressed=calloc(height*width,sizeof(double));
    for(i=0;i<kvalue;i++){
        for(j=0;j<height;j++){
            for(int x=0;x<width;x++){
                Recompressed[j*width+x]+=S[i]*U[i][j]*V[i][x];
            }
        }
    }

    double error = 0;
    for(i=0;i<height;i++)
        for(j=0;j<width;j++){
            double d = Original[i*width+j] - Recompressed[i*width+j];
            error += d*d;
        }

    error = sqrt(error);

    double Relativeerror = (error / Frobeniusnorm) * 100.0;//relative error finding

    printf("Relative error (k=%d) = %.2f%%\n", kvalue, Relativeerror);

    unsigned char *outimg=malloc(height*width);
    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            double vec = Recompressed[i*width+j];
            if(vec<0){
                vec=0;//start from 0 pixel value
            }
            if(vec>255){
                vec=255;//goes till only 255
            }
            outimg[i*width+j]=(unsigned char)vec;
        }
    }

    stbi_write_jpg(outputimage,width,height,1,outimg,90);

    printf("%s saved.\n", outputimage);
}

int main(){
    int imagescount;
    printf("How many images: ");
    scanf("%d",&imagescount);

    int kvals[4]={5,20,50,100};

    for(int i=0;i<imagescount;i++){
        char imgname[100];
        printf("Enter image %d: ", i+1);
        scanf("%s", imgname);

        for(int j=0;j<4;j++){
            char outname[100];
            sprintf(outname, "output%d3k%d.jpg", i+1, kvals[j]);
            compress(imgname, kvals[j], outname);
        }
    }

    printf("done\n");
    return 0;
}







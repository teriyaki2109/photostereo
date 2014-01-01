// Tony TUNG -- copyright 2011-2013
// photometricSVD -- Photometric stereo using SVD technique
// * We assume a Lambertian model, varying lighting conditions, unknown light source directions
// * One camera view
// * Require OpenCV library 
// g++ photostereo.cpp `pkg-config --libs --cflags opencv`
// ./a.out 4

//#include "stdafx.h"
#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

inline float signof(float a) { return (a == 0.0f) ? 0.0f : (a<0.0f ? -1.0f : 1.0f); }

using namespace std;

void colorMap(unsigned char *rgb,float value,float min,float max)
{
  float max3=(max-min)/3;
  value-=min;
  if(value==HUGE_VAL) {rgb[0]=rgb[1]=rgb[2]=255;}
  else if(value<0) {rgb[0]=rgb[1]=rgb[2]=0;}
  else if(value<max3) {rgb[0]=(unsigned char)(255*value/max3);rgb[1]=0;rgb[2]=0;}
  else if(value<2*max3) {rgb[0]=255;rgb[1]=(unsigned char)(255*(value-max3)/max3);rgb[2]=0;}
  else if(value<max) {rgb[0]=255;rgb[1]=255;rgb[2]=(unsigned char)(255*(value-2*max3)/max3);}
  else {rgb[0]=rgb[1]=rgb[2]=255;}
}

class photostereo
{
	int numImg;
	vector<IplImage *> tmpImg;
	int width, height;

	CvMat* U;
	float *dataU;

	IplImage* Mask;
	uchar* dataMask;

	IplImage* Z;
	float *dataZ;

	void computeSVD();
	void integrate();
	void outputMesh(char* fName);

	public:
	photostereo(int aNumImg);
	~photostereo();

};

photostereo::photostereo(int aNumImg)
{
	numImg=aNumImg;
	tmpImg.reserve(numImg);

	//Input files should be 00.jpg...04.jpg
	char fileName[10][64];
	for(int i=0; i<numImg; i++)
	{
		sprintf(fileName[i], "%02d.jpg", i);
		IplImage* aImage=cvLoadImage(fileName[i], 0);
		tmpImg.push_back(aImage);
	}

	//Retrieve image dimension
	height=tmpImg[0]->height;
	width =tmpImg[0]->width;
	cout << "Image size: width=" << width << " height=" << height << "\n";

	computeSVD();
	integrate();
	outputMesh(fileName[0]);

}


photostereo::~photostereo()
{
	for(int i=0; i<numImg; i++) { if(tmpImg[i]) cvReleaseImage(&tmpImg[i]); }
	if(U) cvReleaseMat(&U);
	if(Mask) cvReleaseImage(&Mask);
	if(Z) cvReleaseImage(&Z);
}

void photostereo::computeSVD()
{
	int i,j,k;

	//SVD
	U = cvCreateMat(height*width, numImg, CV_32FC1);
	dataU = U->data.fl;
	CvMat* A = cvCreateMat(height*width, numImg, CV_32FC1);
	CvMat* D = cvCreateMat(numImg, numImg, CV_32FC1);
	float* data = A->data.fl;
	float* dataD = D->data.fl;
	CvMat* V = cvCreateMat(numImg, numImg, CV_32FC1);

	//Populate A
	for(k=0; k<numImg; k++) {
	for(i=0; i<height; i++) {
	for(j=0; j<width; j++) {
		uchar* datak = (uchar*) tmpImg[k]->imageData; 
		data[i*width*numImg+j*numImg+k]=(float) datak[i*width+j];
	} } }
	
	//Output matrix A
	// A = U D V^T
	cvSVD(A, D, U, V, CV_SVD_V_T);
	
	//Output square root eigenvalues of AAT
	for(i=0; i<numImg; i++) { for(j=0; j<numImg; j++) { printf("%2.4f ", dataD[i*numImg+j]); } printf("\n"); }	

	//Memory release
	cvReleaseMat(&A); cvReleaseMat(&V); cvReleaseMat(&D); 

	//U contains the eigenvectors of AAT, which are as well the z,x,y components of the surface normals for each pixel	
	IplImage* Szxy[3]; for(i=0; i<3; i++) { Szxy[i] = cvCreateImage(cvSize(width,height), IPL_DEPTH_8U, 1); }
	uchar* dataSz;	uchar* dataSx; 	uchar* dataSy;
	dataSz = (uchar*) Szxy[0]->imageData;
	dataSx = (uchar*) Szxy[1]->imageData;
	dataSy = (uchar*) Szxy[2]->imageData;

	IplImage* S = cvCreateImage(cvSize(width,height), IPL_DEPTH_8U, 3);
	uchar* dataS = (uchar*) S->imageData;
	
	//Get normalized surface normal components from U
	float Min, Max; Min=dataU[0]; Max=dataU[0];
	for(i=0; i<height*width; i++) {
		for(j=0; j<3; j++) {
			if(dataU[i*numImg+j]<Min) Min=dataU[i*numImg+j];
			if(dataU[i*numImg+j]>Max) Max=dataU[i*numImg+j];			
		}
	}

	float rS, rSxyz, sx, sy, sz;
	for(i=0; i<height; i++) {
	for(j=0; j<width; j++) {		
		rSxyz=1/sqrt(dataU[i*width*numImg+j*numImg]*dataU[i*width*numImg+j*numImg]+
				dataU[i*width*numImg+j*numImg+1]*dataU[i*width*numImg+j*numImg+1]+
					dataU[i*width*numImg+j*numImg+2]*dataU[i*width*numImg+j*numImg+2]);

		dataSz[i*width+j]=128.0f+127.0f*signof(dataU[i*width*numImg+j*numImg])*fabs(dataU[i*width*numImg+j*numImg])*rSxyz;
		dataSx[i*width+j]=128.0f+127.0f*signof(dataU[i*width*numImg+j*numImg+1])*fabs(dataU[i*width*numImg+j*numImg+1])*rSxyz;
		dataSy[i*width+j]=128.0f+127.0f*signof(dataU[i*width*numImg+j*numImg+2])*fabs(dataU[i*width*numImg+j*numImg+2])*rSxyz;

		//printf("%d %d %d \n", dataU[i*width*numImg+j*numImg], dataU[i*width*numImg+j*numImg+1], dataU[i*width*numImg+j*numImg+2]);
		
		dataS[i*width*3+j*3]  =dataSz[i*width+j];
		dataS[i*width*3+j*3+1]=dataSx[i*width+j];
		dataS[i*width*3+j*3+2]=dataSy[i*width+j];
	} } 

	//Output surface images
	cvSaveImage("base0.png", Szxy[0]);
	cvSaveImage("base1.png", Szxy[1]);
	cvSaveImage("base2.png", Szxy[2]);
	cvSaveImage("normalmap.png", S);

	//Memory release
	for(i=0; i<3; i++) { cvReleaseImage(&(Szxy[i])); }
	for(i=0; i<numImg; i++) { cvReleaseImage((IplImage **) &(tmpImg[i])); }

	//Mask computation
	Mask = cvCreateImage(cvSize(width,height), IPL_DEPTH_8U, 1);
	dataMask = (uchar*) Mask->imageData;
	int BG[3];
	BG[0]=255.0f*(dataU[10*width*numImg+10*numImg]-Min)/(Max-Min);
	BG[1]=255.0f*(dataU[10*width*numImg+10*numImg+1]-Min)/(Max-Min);
	BG[2]=255.0f*(dataU[10*width*numImg+10*numImg+2]-Min)/(Max-Min);
	cout << "BG(10,10)=" << BG[0] << "," << BG[1] << "," << BG[2] << "\n";

	for(i=0; i<height; i++) {
	for(j=0; j<width; j++) {
		if(fabs(255.0f*(dataU[i*width*numImg+j*numImg]-Min)/(Max-Min)-BG[0])<2 &&
		   fabs(255.0f*(dataU[i*width*numImg+j*numImg+1]-Min)/(Max-Min)-BG[1])<2 &&
		   fabs(255.0f*(dataU[i*width*numImg+j*numImg+2]-Min)/(Max-Min)-BG[2])<2) { dataMask[i*width+j]=0; }
		else
		{ dataMask[i*width+j]=255;}
	} }
	cvErode(Mask, Mask, NULL, 5);
	//Ouput mask
	cvSaveImage("mask.png", Mask);

	//Memory release
	cvReleaseImage(&S);

	for(i=0; i<height; i++) {
	for(j=0; j<width; j++) {
		rSxyz=2*dataU[i*width*numImg+j*numImg]*dataU[i*width*numImg+j*numImg]+
				dataU[i*width*numImg+j*numImg+1]*dataU[i*width*numImg+j*numImg+1]+
					dataU[i*width*numImg+j*numImg+2]*dataU[i*width*numImg+j*numImg+2];
		if(rSxyz!=0) { rSxyz=1/sqrt(rSxyz); };
		dataU[i*width*numImg+j*numImg]=2*signof(dataU[i*width*numImg+j*numImg])*fabs(dataU[i*width*numImg+j*numImg])*rSxyz;
		dataU[i*width*numImg+j*numImg+1]=signof(dataU[i*width*numImg+j*numImg+1])*fabs(dataU[i*width*numImg+j*numImg+1])*rSxyz;
		dataU[i*width*numImg+j*numImg+2]=signof(dataU[i*width*numImg+j*numImg+2])*fabs(dataU[i*width*numImg+j*numImg+2])*rSxyz;
	} }
	
}

//Normalmap integration to recover depthmap
//zk+1(i,j) = 1/4 * [ zk(i+1,j) + zk(i-1,j) + zk(i,j+1) + zk(i,j-1) + nx(i-1,j) - nx(i,j) + ny(i,j-1) - ny(i,j)]
void photostereo::integrate()
{	
	int i,j,k;	
	Z= cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
	dataZ = (float*) Z->imageData; cvSetZero(Z);

	k=0;
	while(k<1000)
	{
		for(i=1; i<height-1; i++) {
		for(j=1; j<width-1; j++) {
		if(dataMask[i*width+j]>0) {
			//inside the foreground
	 		if(dataMask[(i-1)*width+j]>0 && dataMask[(i+1)*width+j]>0 && dataMask[i*width+j-1]>0 && dataMask[i*width+j+1]>0) {
				dataZ[i*width+j]=0.25*(dataZ[(i+1)*width+j]+dataZ[(i-1)*width+j]+dataZ[i*width+j+1]+dataZ[i*width+j-1]
						+dataU[(i-1)*width*numImg+j*numImg+1]-dataU[i*width*numImg+j*numImg+1]
						+dataU[i*width*numImg+(j-1)*numImg+2]-dataU[i*width*numImg+j*numImg+2]);
			}
			else //top side outside
	 		if(dataMask[(i-1)*width+j]>0 && dataMask[(i+1)*width+j]==0 && dataMask[i*width+j-1]>0 && dataMask[i*width+j+1]>0) {
				dataZ[i*width+j]=0.33*(dataZ[(i-1)*width+j]+dataZ[i*width+j+1]+dataZ[i*width+j-1])
						+0.25*(dataU[(i-1)*width*numImg+j*numImg+1]-dataU[i*width*numImg+j*numImg+1]
						+dataU[i*width*numImg+(j-1)*numImg+2]-dataU[i*width*numImg+j*numImg+2]);
			}
			else //bottom side outside
	 		if(dataMask[(i-1)*width+j]==0 && dataMask[(i+1)*width+j]>0 && dataMask[i*width+j-1]>0 && dataMask[i*width+j+1]>0) {
				dataZ[i*width+j]=0.33*(dataZ[(i+1)*width+j]+dataZ[i*width+j+1]+dataZ[i*width+j-1])
						+0.25*(-dataU[(i+1)*width*numImg+j*numImg+1]+dataU[i*width*numImg+j*numImg+1]
						+dataU[i*width*numImg+(j-1)*numImg+2]-dataU[i*width*numImg+j*numImg+2]);
			}
			else //left side outside
	 		if(dataMask[(i-1)*width+j]>0 && dataMask[(i+1)*width+j]>0 && dataMask[i*width+j-1]==0 && dataMask[i*width+j+1]>0) {
				dataZ[i*width+j]=0.33*(dataZ[(i+1)*width+j]+dataZ[(i-1)*width+j]+dataZ[i*width+j+1])
						+0.25*(dataU[(i-1)*width*numImg+j*numImg+1]-dataU[i*width*numImg+j*numImg+1]
						-dataU[i*width*numImg+(j+1)*numImg+2]+dataU[i*width*numImg+j*numImg+2]);
			}
			else //right side outside
	 		if(dataMask[(i-1)*width+j]>0 && dataMask[(i+1)*width+j]>0 && dataMask[i*width+j-1]>0 && dataMask[i*width+j+1]==0) {
				dataZ[i*width+j]=0.33*(dataZ[(i+1)*width+j]+dataZ[(i-1)*width+j]+dataZ[i*width+j-1])
						+0.25*(dataU[(i-1)*width*numImg+j*numImg+1]-dataU[i*width*numImg+j*numImg+1]
						+dataU[i*width*numImg+(j-1)*numImg+2]-dataU[i*width*numImg+j*numImg+2]);
			}
			else //top left corner outside
	 		if(dataMask[(i-1)*width+j]>0 && dataMask[(i+1)*width+j]==0 && dataMask[i*width+j-1]==0 && dataMask[i*width+j+1]>0) {
				dataZ[i*width+j]=0.5*(dataZ[(i-1)*width+j]+dataZ[i*width+j+1])
						+0.25*(dataU[(i-1)*width*numImg+j*numImg+1]-dataU[i*width*numImg+j*numImg+1]
						-dataU[i*width*numImg+(j+1)*numImg+2]+dataU[i*width*numImg+j*numImg+2]);
			}
			else //top right corner outside
	 		if(dataMask[(i-1)*width+j]>0 && dataMask[(i+1)*width+j]==0 && dataMask[i*width+j-1]>0 && dataMask[i*width+j+1]==0) {
				dataZ[i*width+j]=0.5*(dataZ[(i-1)*width+j]+dataZ[i*width+j-1])
						+0.25*(dataU[(i-1)*width*numImg+j*numImg+1]-dataU[i*width*numImg+j*numImg+1]
						+dataU[i*width*numImg+(j-1)*numImg+2]-dataU[i*width*numImg+j*numImg+2]);
			}
			else //bottom left corner outside
	 		if(dataMask[(i-1)*width+j]==0 && dataMask[(i+1)*width+j]>0 && dataMask[i*width+j-1]==0 && dataMask[i*width+j+1]>0) {
				dataZ[i*width+j]=0.5*(dataZ[(i+1)*width+j]+dataZ[i*width+j+1])
						+0.25*(-dataU[(i+1)*width*numImg+j*numImg+1]+dataU[i*width*numImg+j*numImg+1]
						-dataU[i*width*numImg+(j+1)*numImg+2]+dataU[i*width*numImg+j*numImg+2]);
			}
			else //bottom right corner outside
	 		if(dataMask[(i-1)*width+j]==0 && dataMask[(i+1)*width+j]>0 && dataMask[i*width+j-1]>0 && dataMask[i*width+j+1]==0) {
				dataZ[i*width+j]=0.5*(dataZ[(i+1)*width+j]+dataZ[i*width+j-1])
						+0.25*(-dataU[(i+1)*width*numImg+j*numImg+1]+dataU[i*width*numImg+j*numImg+1]
						+dataU[i*width*numImg+(j-1)*numImg+2]-dataU[i*width*numImg+j*numImg+2]);
			}
		} } }
		k++;
	}

}

//Output 3D mesh in off format. fName is used to retrieve the colors
void photostereo::outputMesh(char* fName)
{
	int i,j,k;

	float minZ=0, maxZ=0;
	for(i=1; i<height-1; i++) {
	for(j=1; j<width-1; j++) {
		if(dataMask[i*width+j]>0) {
			if(dataZ[i*width+j]<minZ) minZ=dataZ[i*width+j];
			if(dataZ[i*width+j]>maxZ) maxZ=dataZ[i*width+j];
		}
	} }

	cout << "minZ=" << minZ << ", maxZ=" << maxZ << "\n";

	IplImage* Zc= cvCreateImage(cvSize(width,height), IPL_DEPTH_8U, 3);
	char* dataZc = (char*) Zc->imageData;

	FILE* pointcloud; pointcloud=fopen("poincloud.off", "w");
	fprintf(pointcloud, "COFF\n"); fprintf(pointcloud, "%d %d 0\n", width*height, (width-1)*(height-1)*2);
	float val;
	IplImage *colorImg=cvLoadImage(fName, 1);
	uchar* dataColor = (uchar*) colorImg->imageData;

	for(i=0; i<height; i++) {
	for(j=0; j<width; j++) {
		if(dataMask[i*width+j]>0) {
			val=1-(dataZ[i*width+j]-minZ)/(maxZ-minZ);
			dataZ[i*width+j]=255*(1-(dataZ[i*width+j]-minZ)/(maxZ-minZ));
		} else { val=0.0f; dataZ[i*width+j]=0; }

		unsigned char color[3]; colorMap(color, val, 0, 1);
		dataZc[i*width*3+j*3]=color[0];
		dataZc[i*width*3+j*3+1]=color[1];
		dataZc[i*width*3+j*3+2]=color[2];
	
		fprintf(pointcloud, "%d %d %f %2.3f %2.3f %2.3f 1.0\n", i, j, 150*val,
			float(dataColor[i*width*3+j*3+2])/255.0f, float(dataColor[i*width*3+j*3+1])/255.0f, float(dataColor[i*width*3+j*3])/255.0f);

	} }

	for(i=0; i<height-1; i++) {
	for(j=0; j<width-1; j++) {
		fprintf(pointcloud, "3 %d %d %d\n", i*width+j, i*width+j+1, (i+1)*width+j);
		fprintf(pointcloud, "3 %d %d %d\n", i*width+j+1, (i+1)*width+j+1, (i+1)*width+j);		
	} }

	fclose(pointcloud);

	cvSaveImage("depthmap.png", Z);
	cvSaveImage("depthmapc.png", Zc);
	
	//Memory release
	cvReleaseImage(&colorImg);
	cvReleaseImage(&Zc);
}


int main (int argc, char **argv)
{
	cout << "Photostereo v0.1b: 3D reconstruction using monitor as light source - Tony Tung 2013 (c)\n";
	cout << "** Have fun and feel free to improve it :) **\n";
	cout << "tonytung.org\n\n";

	if(argc != 2) { cout << "Usage: ./photostereo 5=[nb of input images]\nPlease try again :)\n"; return(0); }
	const int numImg=atoi(argv[1]);//argc-1; //Number of images (=number of different light sources)
	cout << "Input: " << numImg << "images\n";

	photostereo photos(numImg);

	return 1;
}

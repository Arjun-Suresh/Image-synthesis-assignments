// =============================================================================
// VIZA656/CSCE647 at Texas A&M University
// Homework 5
// Texture mapping
// grass.ppm - plane texture; sky.ppm - Infinite sphere texture; stone.ppm - sphere texture; normal.ppm - bump map
// output files are generated in the same folder
// =============================================================================

#include <cstdlib>
#include <iostream>
#include <GL/glut.h>
#include <math.h>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cstring>
#include <malloc.h>
#include "OBJ_Loader.h"


#define WIDTHVALUEID 1
#define HEIGHTVALUEID 2
#define MAXCOLORVALUEID 3
#define COLORFOREGROUND 0
#define COLORBACKGROUND 1

#define REDOFFSET 0
#define GREENOFFSET 1
#define BLUEOFFSET 2


#define maximum(x, y, z) ((x) > (y)? ((x) > (z)? (x) : (z)) : ((y) > (z)? (y) : (z)))
#define minimum(x, y, z) ((x) < (y)? ((x) < (z)? (x) : (z)) : ((y) < (z)? (y) : (z)))
#define SX 100
#define SY 100

#define SLX 100
#define SLY 100

#define SPOTLIGHTMIN 0.707
#define MINANGLE 0.0
#define MAXANGLE 1.0
#define MINREFLECTION 0.9
#define MAXREFLECTION 1.0
#define KD 10
#define KS 10
#define KB 0


using namespace std;
// =============================================================================
// These variables will store the input ppm image's width, height, and color
// =============================================================================
int width[4], height[4], maxColorValue=255, magicNo, widthComputed=750, heightComputed=750;
unsigned char *pixmapOrig[4], *pixmapComputed;

class point
{
  public:
  double x;
  double y;
  double z;
  point(){}
  point(double xVal, double yVal, double zVal)
  {
    x=xVal, y=yVal, z=zVal;
  }
  point(const point& p)
  {
    x = p.x;
    y=p.y;
    z=p.z;
  }
};


class gVector
{  
  public:
  double x;
  double y;
  double z;
  gVector(){}
  gVector(double xVal, double yVal, double zVal)
  {
    x=xVal, y=yVal, z=zVal;
  }

  gVector(const gVector& v)
  {
    x = v.x;
    y = v.y;
    z = v.z;
  }
  
  double length()
  {
    double dist = pow(x,2)+pow(y,2)+pow(z,2);
    return sqrt(dist);
  }
  
  gVector* operator * (const gVector& vObj)
  {
    double xVal = y*vObj.z-z*vObj.y;
    double yVal = z*vObj.x-x*vObj.z;
    double zVal = x*vObj.y-y*vObj.x;
    return new gVector(xVal, yVal, zVal);
  }

  double dotProduct(const gVector& vObj)
  {
    return x*vObj.x+y*vObj.y+z*vObj.z;
  }
  
  void scalarMultiply(double val)
  {
    x=x*val, y=y*val, z=z*val;
  } 
  
  gVector* operator + (const gVector& vObj)
  {
    double xVal = x+vObj.x;
    double yVal = y+vObj.y;
    double zVal = z+vObj.z;
    return new gVector(xVal, yVal, zVal);
  }
  
  gVector* operator - (const gVector& vObj)
  {
    double xVal = x-vObj.x;
    double yVal = y-vObj.y;
    double zVal = z-vObj.z;
    return new gVector(xVal, yVal, zVal);
  }   
};

class color
{
  public:
    double red;
    double green;
    double blue;
    double alpha;
    color(double r, double g, double b, double a)
    {
      red=r, green=g, blue=b, alpha=a;
    }
    color(){}
    void multiply(double t)
    {
      if(t == 0)
      	alpha=red=green=blue=0;
      else
      {
        red = red*t;
        green = green*t;
        blue = blue*t;
      }
    }
    
    void multiplyColor(color& input)
    {
      red = red*input.red;
      green = green*input.green;
      blue = blue*input.blue;
      alpha = alpha*input.alpha;
    }
    
    void addColor(color& input)
    {
      red = red+input.red;
      green = green+input.green;
      blue = blue+input.blue;
      alpha = alpha+input.alpha;
    }
    
    void setColor(double r, double g, double b, double a)
    {
      red=r, green=g, blue=b, alpha=a;
    }

    void getResult(double& rVal, double& gVal, double& bVal)
    {
      rVal=red/alpha;
      gVal=green/alpha;
      bVal=blue/alpha;
      if (gVal>1)
        gVal=1;
      if (rVal>1)
        rVal=1;
      if (bVal>1)
        bVal=1;
    }
};
class matrix
{
  public:
  float a[4][4];
  void rotateAroundX(float theta)
  {
    for(int i=0;i<4;i++)
    {
      a[i][i]=1;
      if(i==1)
      {
        a[i][1]=cos(theta);
        a[i][2]=-sin(theta);
      }
      if(i==2)
      {
        a[i][1]=sin(theta);
        a[i][2]=cos(theta);
      }
    }
  }
  void rotateAroundY(float theta)
  {
    for(int i=0;i<4;i++)
    {
      a[i][i]=1;
      if(i==2)
      {
        a[i][2]=cos(theta);
        a[i][0]=-sin(theta);
      }
      if(i==0)
      {
        a[i][2]=sin(theta);
        a[i][0]=cos(theta);
      }
    }
  }
  void rotateAroundZ(float theta)
  {
    for(int i=0;i<4;i++)
    {
      a[i][i]=1;
      if(i==0)
      {
        a[i][0]=cos(theta);
        a[i][1]=-sin(theta);
      }
      if(i==1)
      {
        a[i][0]=sin(theta);
        a[i][1]=cos(theta);
      }
    }
  }
  void translate(float tx, float ty, float tz)
  {
    for(int i=0;i<4;i++)
    {
      a[i][i]=1;
      if (i==0)
        a[i][3]=tx;
      if (i==1)
        a[i][3]=ty;
      if (i==2)
        a[i][3]=tz;
    }
  }
  matrix* operator * (const matrix& m)
  {
    matrix* result = new matrix();
    for(int i=0;i<4;i++)
    {
      for(int j=0;j<4;j++)
      {
        int sum=0;
        for(int k=0;k<4;k++)
          sum+=a[i][k]*m.a[k][j];
        result->a[i][j]=sum;
      }
    }
    return result;
  }
  template <typename T>
  void operate(T& input, T& output)
  {
    output.x = a[0][0] * input.x + a[0][1] * input.y + a[0][2] * input.z;
    output.y = a[1][0] * input.x + a[1][1] * input.y + a[1][2] * input.z;
    output.z = a[2][0] * input.x + a[2][1] * input.y + a[2][2] * input.z;
  }
};

class triangle
{
  public:
  point* vertices[3];
  gVector* normals[3];
  gVector* triangleNormal;
  float texture[3][2];
  bool colorAdded;

  triangle(point& p1, point& p2, point& p3, gVector& v1, gVector& v2, gVector& v3, float tex1[2] = NULL, float tex2[2] = NULL, float tex3[2] = NULL, bool rotate=false)
  {
    vertices[0] = new point(p1.x, p1.y, p1.z);
    vertices[1] = new point(p2.x, p2.y, p2.z);
    vertices[2] = new point(p3.x, p3.y, p3.z);

    normals[0] = new gVector(v1.x, v1.y, v1.z);
    normals[1] = new gVector(v2.x, v2.y, v2.z);
    normals[2] = new gVector(v3.x, v3.y, v3.z);

    if (tex1 && tex2 && tex3)
    {
      colorAdded=true;
      texture[0][0] = tex1[0];
      texture[0][1] = tex1[1];
      texture[1][0] = tex2[0];
      texture[1][1] = tex2[1];
      texture[2][0] = tex3[0];
      texture[2][1] = tex3[1];
    }
    else
      colorAdded=false;

    if (rotate)
    {
      for(int i=0; i<3;i++)
      {
        matrix m1, m2, m3, m4;
        m1.translate(-vertices[i]->x, -vertices[i]->y, -vertices[i]->z);
        m2.rotateAroundX(3.1416/4.0);
        m3.rotateAroundY(3.1416/4.0);
        m4.translate(vertices[i]->x, vertices[i]->y, vertices[i]->z);
        matrix* res1 = m1 * m2;
        matrix* res2 = *res1 * m3;
        matrix* res3 = *res2 * m4;
        point result;
        res3->operate<point>(*vertices[i], result);
        vertices[i]->x = result.x;
        vertices[i]->y = result.y;
        vertices[i]->z = result.z;
        delete res1;
        delete res2;
        delete res3;

        m1.translate(-normals[i]->x, -normals[i]->y, -normals[i]->z);
        m4.translate(normals[i]->x, normals[i]->y, normals[i]->z);
        matrix* res4 = m1 * m2;
        matrix* res5 = *res4 * m3;
        matrix* res6 = *res5 * m4;
        gVector vResult;
        res6->operate<gVector>(*normals[i], vResult);
        normals[i]->x = vResult.x;
        normals[i]->y = vResult.y;
        normals[i]->z = vResult.z;
        delete res4;
        delete res5;
        delete res6;
      }
    }
    gVector t1(vertices[1]->x-vertices[0]->x, vertices[1]->y-vertices[0]->y, vertices[1]->z-vertices[0]->z);
    gVector t2(vertices[2]->x-vertices[0]->x, vertices[2]->y-vertices[0]->y, vertices[2]->z-vertices[0]->z);
    gVector *t3 = t1 * t2;
    triangleNormal = new gVector(t3->x, t3->y, t3->z);
    triangleNormal->scalarMultiply(1.0/triangleNormal->length());

  }

  triangle(point& p1, point& p2, point& p3, float tex1[2] = NULL, float tex2[2] = NULL, float tex3[2] = NULL, bool rotate=false)
  {
    vertices[0] = new point(p1.x, p1.y, p1.z);
    vertices[1] = new point(p2.x, p2.y, p2.z);
    vertices[2] = new point(p3.x, p3.y, p3.z);

    if (tex1 && tex2 && tex3)
    {
      colorAdded=true;
      texture[0][0] = tex1[0];
      texture[0][1] = tex1[1];
      texture[1][0] = tex2[0];
      texture[1][1] = tex2[1];
      texture[2][0] = tex3[0];
      texture[2][1] = tex3[1];
    }
    else
      colorAdded=false;

    if (rotate)
    {
      for(int i=0; i<3;i++)
      {
        matrix m1, m2, m3, m4;
        m1.translate(-vertices[i]->x, -vertices[i]->y, -vertices[i]->z);
        m2.rotateAroundX(3.1416/4.0);
        m3.rotateAroundY(3.1416/4.0);
        m4.translate(vertices[i]->x, vertices[i]->y, vertices[i]->z);
        matrix* res1 = m1 * m2;
        matrix* res2 = *res1 * m3;
        matrix* res3 = *res2 * m4;
        point result;
        res3->operate<point>(*vertices[i], result);
        vertices[i]->x = result.x;
        vertices[i]->y = result.y;
        vertices[i]->z = result.z;
        delete res1;
        delete res2;
        delete res3;
      }
    }
    
    gVector t1(vertices[1]->x-vertices[0]->x, vertices[1]->y-vertices[0]->y, vertices[1]->z-vertices[0]->z);
    gVector t2(vertices[2]->x-vertices[0]->x, vertices[2]->y-vertices[0]->y, vertices[2]->z-vertices[0]->z);
    gVector *t3 = t1 * t2;
    triangleNormal = new gVector(t3->x, t3->y, t3->z);
    triangleNormal->scalarMultiply(1.0/triangleNormal->length());

    normals[0] = new gVector(triangleNormal->x, triangleNormal->y, triangleNormal->z);
    normals[1] = new gVector(triangleNormal->x, triangleNormal->y, triangleNormal->z);
    normals[2] = new gVector(triangleNormal->x, triangleNormal->y, triangleNormal->z);    
  }

  ~triangle()
  {
    for(int i=0;i<3;i++)
    {
      delete vertices[i];
      delete normals[i];
    }
  }
};

point *pc, *p0, *pe; //Field of view and eye
point* spheres[2], *cylinder, *plane00; //Shapes to be rendered- spheres and plane
double radius[2]; //radius of 2 spheres
double d; //distance between field of view and the eye
gVector *n2, *n1, *n0; //field of view related normal gVectors
gVector *planeN0, *planeN1, *planeN2; //gVectors of an external plane which is rendered

point* lightPosition; //Only for point and spot light
gVector* lightDirection; //Only for spot and direction light
point* lightCorner;
gVector* lightN0, *lightN1; //For area light

point *ppc, *pp0, *ppe; //Field of view and eye
double pd; //distance between field of view and the eye
gVector *pn2, *pn1, *pn0; //field of view related normal gVectors

triangle** tObjs;
int numOfMeshes;

inline double mod(double x) 
{
  return (x > 0 ? (x) : (-1*x));
}

//Function to resize ppm file buffer array
unsigned char* resizeArray(unsigned char* oldArray, long int oldSize, long int& newSize) 
{
    newSize = oldSize * 2;
    unsigned char* newArray = new unsigned char[newSize];
    std::memcpy( newArray, oldArray, oldSize * sizeof(unsigned char) );
    return newArray;
}

//Set red, green and blue in the new pixMap array to be written to the ppm file
void setPixelColor(int y, int x, int red, int green, int blue)
{
  int i = ((heightComputed-y-1) * widthComputed + x) * 3; 
  pixmapComputed[i++] = red;
  pixmapComputed[i++] = green;
  pixmapComputed[i] = blue;
}

void setPixelColorOrig(long int &val, long int& index, unsigned char* fileBuffer, int iter)
{
  pixmapOrig[iter][val++]=fileBuffer[index++];
}
// =============================================================================
// OpenGL Display and Mouse Processing Functions.
//
// You can read up on OpenGL and modify these functions, as well as the commands
// in main(), to perform more sophisticated display or GUI behavior. This code
// will service the bare minimum display needs for most assignments.
// =============================================================================
static void windowResize(int w, int h)
{   
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,(w/2),0,(h/2),0,1); 
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity() ;
}
static void windowDisplay(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  glRasterPos2i(0,0);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  glDrawPixels(widthComputed, heightComputed, GL_RGB, GL_UNSIGNED_BYTE, pixmapComputed);
  glFlush();
}
static void processMouse(int button, int state, int x, int y)
{
  if(state == GLUT_UP)
  exit(0);               // Exit on mouse click.
}
static void init(void)
{
  glClearColor(1,1,1,1); // Set background color.
}




//****************************************************************************************************************************
//**********************Functions to read data from ppm file******************************************************************
//**********************Not used in this project******************************************************************************
//****************************************************************************************************************************
void parseCommentLine(long int& index, unsigned char* fileBuffer, long int numOfCharacters)
{
  while(index<numOfCharacters && fileBuffer[index]!='\n')
  index++;
}

bool parseMagicNumber(long int& index, unsigned char* fileBuffer, long int numOfCharacters, int& magicNumberParsed)
{
  while(index<numOfCharacters)
  {
    if(fileBuffer[index]=='#')
    {
      parseCommentLine(index, fileBuffer, numOfCharacters);
      break;
    }
    if(fileBuffer[index]=='P' || fileBuffer[index]=='p')
    {
      index++;
      if(fileBuffer[index]=='1' || fileBuffer[index]=='2'|| fileBuffer[index]=='3' || fileBuffer[index]=='4' || fileBuffer[index]=='5' || fileBuffer[index]=='6')
      {
        magicNumberParsed=1;
	magicNo=fileBuffer[index]-48;
        break;
      }
      else 
        return false;
    }
    if(isspace(fileBuffer[index]) && fileBuffer[index]!='\n')
      index++;
    else
      return false;
  }
  if(magicNumberParsed)
  {
    index++;
    if(index<numOfCharacters && !(isspace(fileBuffer[index]) || fileBuffer[index]=='#'))
      return false;
    else
    {
      if(fileBuffer[index]=='#')
      {
        parseCommentLine(index, fileBuffer, numOfCharacters);
      }
    }
        
  }
  return true;
}

bool parseValue(long int& index, unsigned char* fileBuffer, long int numOfCharacters, int& parseValue, int valueId, int iter)
{
  int value;
  while(index<numOfCharacters)
  {
    if(fileBuffer[index]=='#')
    {
      parseCommentLine(index, fileBuffer, numOfCharacters);
      break;
    }
    if(isdigit(fileBuffer[index]))
    {
      char valueString[20];
      int k=0;
      valueString[k++]=fileBuffer[index++];

      while(isdigit(fileBuffer[index]) && index<numOfCharacters)
        valueString[k++]=fileBuffer[index++];
      valueString[k]='\0';
      value = atoi(valueString);
      parseValue=1;
      break;
    }
    if(isspace(fileBuffer[index]))
      index++;
    else
      return false;
  }
  if(parseValue)
  {
    if(index<numOfCharacters && !((isspace(fileBuffer[index]) && valueId != MAXCOLORVALUEID) || (valueId == MAXCOLORVALUEID && fileBuffer[index] == '\n') || fileBuffer[index]=='#'))
      return false;
    else
    {
      if(fileBuffer[index]=='#')
      {
        parseCommentLine(index, fileBuffer, numOfCharacters);
      }
    }        
  }
  switch(valueId)
  {
    case WIDTHVALUEID:
      width[iter]=value;
      break;
    case HEIGHTVALUEID:
      height[iter]=value;
      break;
    case MAXCOLORVALUEID:
      maxColorValue=value;
      break;
  }
  return true;
}

bool fillPixelsBin(long int& index, unsigned char* fileBuffer, long int numOfCharacters, int iter)
{
  int wVal,hVal;
  wVal=width[iter];
  hVal=height[iter];
  long int rowVal=0,colVal=0;
  while(rowVal<hVal)
  {
    long int val=(((hVal-rowVal-1)*wVal)+colVal)*3;
    for(int i=0;i<3;i++)
    	setPixelColorOrig(val, index, fileBuffer, iter);
    colVal=(colVal+1)%wVal;
    if(!colVal)
      rowVal++;
  }
  return true;
}

bool readPPMFile(char* filePath, int iter)
{
  std::fstream ppmFile;
  std::ifstream checkFile(filePath);
  long int numOfCharacters=0;
  int magicNumberParsed=0, widthParsed=0, heightParsed=0, maxColorValueParsed=0;
  if(!checkFile.good())
  {
    cout<<"Filename not found\n";
    return false;
  }
  ppmFile.open (filePath, std::ios::in|std::ios::binary);
  ppmFile.seekg(0,std::ios::end);
  ppmFile.seekg(0,std::ios::end);
  numOfCharacters = ppmFile.tellg();
  ppmFile.seekg(0,std::ios::beg);
  unsigned char* fileBuffer = new unsigned char[numOfCharacters+1];
  ppmFile.read((char *)fileBuffer, numOfCharacters);
  fileBuffer[numOfCharacters]='\0';
  ppmFile.close();
  for(long int index=0;index<numOfCharacters;index++)
  {
    if(fileBuffer[index]=='#')
    {
      parseCommentLine(index, fileBuffer, numOfCharacters);
      continue;
    }
    if(!magicNumberParsed)
    {
      if(!parseMagicNumber(index, fileBuffer, numOfCharacters, magicNumberParsed))
      {
        return false;
      }
        
      continue;
    }
    if(magicNumberParsed && !widthParsed)
    {
      if(!parseValue(index, fileBuffer, numOfCharacters, widthParsed, WIDTHVALUEID, iter))
      {
        return false;
      }
      continue;
    }
    if(magicNumberParsed && widthParsed && !heightParsed)
    {
      if(!parseValue(index, fileBuffer, numOfCharacters, heightParsed, HEIGHTVALUEID, iter))
      {
        return false;
      }
      continue;
    }
    if(magicNumberParsed && widthParsed && heightParsed && !maxColorValueParsed)
    {
      if(!parseValue(index, fileBuffer, numOfCharacters, maxColorValueParsed, MAXCOLORVALUEID, iter))
      {
        return false;
      }
      continue;
    }
    if(magicNumberParsed && widthParsed && heightParsed && maxColorValueParsed)
    {
      pixmapOrig[iter] = new unsigned char[width[iter]*height[iter]*3];
      if(!(fillPixelsBin(index, fileBuffer, numOfCharacters, iter)))
        return false;
      delete[] fileBuffer;
      return true;
    }
  }
}





//************************************************************************************************************
//**********************Functions to save from pixelMap to ppm file*******************************************
//************************************************************************************************************
void fillCharacters(unsigned char* fileBuffer, long int& index, char* data)
{
  for(int i=0;i<strlen(data);i++)
    fileBuffer[index++]=data[i];
}

void writeToFile(unsigned char* fileBuffer, long int numberOfCharacters, fstream& ppmFile)
{  
  for(long int i=0; i< numberOfCharacters; i++)
  {
    ppmFile << fileBuffer[i];
  }
}

void writeHeader(unsigned char* fileBuffer, long int& index)
{
  maxColorValue=255;    
  char widthString[5], heightString[5], maxColorString[5];
  sprintf(widthString, "%d", widthComputed);
  sprintf(heightString, "%d", heightComputed);
  sprintf(maxColorString, "%d", maxColorValue);
  char magicNumberString[3];
  strcpy(magicNumberString, "P6");
  fillCharacters(fileBuffer, index, magicNumberString);
  fileBuffer[index++]='\n';
  fillCharacters(fileBuffer, index, widthString);
  fileBuffer[index++]=' ';
  fillCharacters(fileBuffer, index, heightString);
  fileBuffer[index++]='\n';
  fillCharacters(fileBuffer, index, maxColorString);
  fileBuffer[index++]='\n';  
}


unsigned char* writePixelBufferToFileBuffer(unsigned char* fileBuffer, long int& index, long int origSize)
{
  int charCount=0;
  for(int i=0;i<heightComputed;i++)
  {
    for(int j=0;j<widthComputed;j++)
    {
      if(index >= origSize-50)
      {
        unsigned char* tempBuffer = resizeArray(fileBuffer, origSize, origSize);
        delete[] fileBuffer;
        fileBuffer=tempBuffer;
      }
      int k=((heightComputed-i-1)*widthComputed+j)*3;
      for(int x=0;x<3;x++)
      fileBuffer[index++]=pixmapComputed[k++];
    }
  }
  return fileBuffer;
}

void generatePPMFile()
{
  std::fstream ppmFile;
  char fileName[50];
  strcpy(fileName,"outputRasterConversion.ppm");
  ppmFile.open(fileName,std::fstream::out);
  long int index=0;
  unsigned char* fileBuffer = new unsigned char[10000];
  writeHeader(fileBuffer, index);
  fileBuffer=writePixelBufferToFileBuffer(fileBuffer, index, 10000);
  writeToFile(fileBuffer, index, ppmFile);
  ppmFile.close();
}


//*****************************************************************************************************
//*************************************Rasterization functions*****************************************
//*****************************************************************************************************

bool sphereIntersection(point* center, point* pe, double radius, gVector* npe, double& t)
{
  gVector* pce = new gVector(center->x - pe->x, center->y - pe->y, center->z - pe->z);
  double c = pow(pce->length(),2) - pow(radius,2);
  double b = pce->dotProduct(*npe);
  double delta = pow(b,2) - c;
  delete pce;
  if (b <= 0 || delta <=0)
    return false;
  double t1 = b + pow(delta,0.5);
  double t2 = b - pow(delta,0.5);
  t = (t1<t2)?(t1):(t2);
  return true;
}

bool planeIntersection(point* corner, point* pe, gVector* normal, gVector* npe, double& t)
{
  gVector* test = new gVector(pe->x - corner->x, pe->y - corner->y, pe->z - corner->z);
  if (normal->dotProduct(*npe) == 0)
    return false;
  t = -1*((normal->dotProduct(*test))/(normal->dotProduct(*npe)));
  delete test;
  if (t<=0 || t>1000)
    return false;
  return true;
}

void clamp(double& value, double min, double max)
{
  if (value<min)
    value=min;
  if (value>max)
    value=max;
}

void initEnvironment()
{
  d=50;
  n2 = new gVector(0,0,-1);
  pe = new point (250, 250, 150);
  pc = new point (pe->x+(n2->x)*d,pe->y+(n2->y)*d,pe->z+(n2->z)*d);
  gVector* Vref = new gVector(0, 1, 2);
  Vref->scalarMultiply(((double)1)/Vref->length());
  
  n0 = (*n2) * (*Vref);
  n0->scalarMultiply(((double)1)/n0->length());
  n1 = (*n0) * (*n2);
  n1->scalarMultiply(((double)1)/n1->length());
  p0 = new point (pc->x-((SX/2)*(n0->x))-((SY/2)*(n1->x)),pc->y-((SX/2)*(n0->y))-((SY/2)*(n1->y)),pc->z-((SX/2)*(n0->z))-((SY/2)*(n1->z)));

}

void initLight(int option)
{
  
  if(option == 1)
    lightPosition = new point(150,150,0);
  else if(option == 2)
    lightPosition = new point(100,200,-70);
  else
    lightPosition = new point(100,160,-50);
}

void initMeshes(int option)
{
  if (option ==1)
  {
    numOfMeshes=12;
    point p1(182, 251, -150);
    point p2(182, 251, -50);
    point p3(252, 181, -50);
    point p4(252, 181, -150);
    point p5(111, 181, -50);
    point p6(111, 181, -150);
    point p7(182, 110 , -150);
    point p8(182, 110, -50);

    float* dummy = NULL;
    tObjs = (triangle**)malloc(sizeof(triangle*) * numOfMeshes);
    tObjs[0] = new triangle(p1, p2, p3);
    tObjs[1] = new triangle(p1, p3, p4);
    tObjs[2] = new triangle(p5, p6, p7);
    tObjs[3] = new triangle(p5, p7, p8);
    tObjs[4] = new triangle(p2, p5, p8);
    tObjs[5] = new triangle(p2, p8, p3);
    tObjs[6] = new triangle(p4, p7, p6);
    tObjs[7] = new triangle(p4, p6, p1);
    tObjs[8] = new triangle(p2, p1, p6);
    tObjs[9] = new triangle(p2, p6, p5);
    tObjs[10] = new triangle(p8, p7, p4);
    tObjs[11] = new triangle(p8, p4, p3);
  }
  else
  {
    numOfMeshes=4;
    point p1(286, 163, -111);
    point p2(97, 158, -110);
    point p3(184, 360, -110);
    point p4(167, 149, -71);
    float* dummy = NULL;
    tObjs = (triangle**)malloc(sizeof(triangle*) * numOfMeshes);
    tObjs[0] = new triangle(p1, p3, p2);
    tObjs[1] = new triangle(p1, p3, p4);
    tObjs[2] = new triangle(p1, p4, p2);
    tObjs[3] = new triangle(p2, p4, p3);
  }
}

double clamp(double val)
{
  if (val<0)
    return 0;
  if (val>1)
    return 1;
  return val;
}

void initMeshFromFile(int option)
{
  char meshFile[100];
  cout<<"Enter the name of obj file\n";
  cin>>meshFile;
  objl::Loader Loader;
	bool loadout = Loader.LoadFile(meshFile);
	if (!loadout)
	{
    cout<<"Obj file doesn't exist\n";
    exit(0);
  }
	objl::Mesh curMesh = Loader.LoadedMeshes[0];
  numOfMeshes =curMesh.Indices.size()/3;
  tObjs = (triangle**)malloc(sizeof(triangle*) * numOfMeshes);
  int k1=0, k2=0, k3=0;
  
	for (int j = 0; j < curMesh.Indices.size(); j+=3)
	{
	  if (curMesh.Vertices[curMesh.Indices[j]].Position.X < 100)
            k1=200;
          if (curMesh.Vertices[curMesh.Indices[j]].Position.Y < 100)
            k2=200;
          if (curMesh.Vertices[curMesh.Indices[j]].Position.Z > -50)
            k3=-200;
          if (k1 && k2 && k3)
            break;
        }
	for (int j = 0; j < curMesh.Indices.size(); j+=3)
	{
	  point p1(curMesh.Vertices[curMesh.Indices[j]].Position.X+k1, curMesh.Vertices[curMesh.Indices[j]].Position.Y+k2, curMesh.Vertices[curMesh.Indices[j]].Position.Z+k3);
    point p2(curMesh.Vertices[curMesh.Indices[j+1]].Position.X+k1, curMesh.Vertices[curMesh.Indices[j+1]].Position.Y+k2, curMesh.Vertices[curMesh.Indices[j+1]].Position.Z+k3);
    point p3(curMesh.Vertices[curMesh.Indices[j+2]].Position.X+k1, curMesh.Vertices[curMesh.Indices[j+2]].Position.Y+k2, curMesh.Vertices[curMesh.Indices[j+2]].Position.Z+k3);
    
    gVector v1(curMesh.Vertices[curMesh.Indices[j]].Normal.X, curMesh.Vertices[curMesh.Indices[j]].Normal.Y, curMesh.Vertices[curMesh.Indices[j]].Normal.Z);
    gVector v2(curMesh.Vertices[curMesh.Indices[j+1]].Normal.X, curMesh.Vertices[curMesh.Indices[j+1]].Normal.Y, curMesh.Vertices[curMesh.Indices[j+1]].Normal.Z);
    gVector v3(curMesh.Vertices[curMesh.Indices[j+2]].Normal.X, curMesh.Vertices[curMesh.Indices[j+2]].Normal.Y, curMesh.Vertices[curMesh.Indices[j+2]].Normal.Z);
    if (option == 2)
    {
      float t1[2], t2[2], t3[2];
      t1[0] = curMesh.Vertices[curMesh.Indices[j]].TextureCoordinate.X;
      t1[1] = curMesh.Vertices[curMesh.Indices[j]].TextureCoordinate.Y;

      t2[0] = curMesh.Vertices[curMesh.Indices[j+1]].TextureCoordinate.X;
      t2[1] = curMesh.Vertices[curMesh.Indices[j+1]].TextureCoordinate.Y;

      t3[0] = curMesh.Vertices[curMesh.Indices[j+2]].TextureCoordinate.X;
      t3[1] = curMesh.Vertices[curMesh.Indices[j+2]].TextureCoordinate.Y;

      tObjs[j/3] = new triangle(p1, p2, p3, v1, v2, v3, t1, t2, t3);
    }
    else
      tObjs[j/3] = new triangle(p1, p2, p3, v1, v2, v3);
  }
}

float abs(float val)
{
  if (val < 0)
    return val * (-1.0);
  return val;
}
float maxAbs(float a, float b, float c)
{
  if (abs(a)>abs(b))
  {
    if (abs(a)>abs(c))
      return a;
  }
  else
  {
    if (abs(b)>abs(c))
      return b;
  }
  return c;
}

void computeTriangleValues(point& hitpoint, triangle& tObj, float& max, float& val1, float& val2, float& val3)
{
  gVector v1(tObj.vertices[0]->x-hitpoint.x, tObj.vertices[0]->y-hitpoint.y, tObj.vertices[0]->z-hitpoint.z);
  gVector v2(tObj.vertices[1]->x-hitpoint.x, tObj.vertices[1]->y-hitpoint.y, tObj.vertices[1]->z-hitpoint.z);
  gVector v3(tObj.vertices[2]->x-hitpoint.x, tObj.vertices[2]->y-hitpoint.y, tObj.vertices[2]->z-hitpoint.z);

  gVector* a1 = v2 * v3;
  gVector* a2 = v3 * v1;
  gVector* a3 = v1 * v2;
  gVector t1(tObj.vertices[1]->x-tObj.vertices[0]->x, tObj.vertices[1]->y-tObj.vertices[0]->y, tObj.vertices[1]->z-tObj.vertices[0]->z);
  gVector t2(tObj.vertices[2]->x-tObj.vertices[0]->x, tObj.vertices[2]->y-tObj.vertices[0]->y, tObj.vertices[2]->z-tObj.vertices[0]->z);
  gVector *a = t1 * t2;

  max = maxAbs(a->x, a->y, a->z);
  if (max == a->x)
  {
    val1 = a1->x;
    val2 = a2->x;
    val3 = a3->x;
  }
  else if (max == a->y)
  {
    val1 = a1->y;
    val2 = a2->y;
    val3 = a3->y;
  }
  else
  {
    val1 = a1->z;
    val2 = a2->z;
    val3 = a3->z;
  }
  delete a1;
  delete a2;
  delete a3;
  delete a;
}

bool liesInside(point& hitpoint, triangle& tObj)
{
  float max;
  float val1, val2, val3;

  computeTriangleValues(hitpoint, tObj, max, val1, val2, val3);

  if (val1/max > 0 && val2/max > 0 && val3/max > 0 && (val1/max + val2/max + val3/max >= 0.99) && (val1/max + val2/max + val3/max <= 1.01))
    return true;
  return false;
}

void computeParameters(triangle& tObj, point& hitpoint, gVector& normal, float tex[2]=NULL)
{
  float max;
  float val1, val2, val3;

  computeTriangleValues(hitpoint, tObj, max, val1, val2, val3);

  normal.x = tObj.normals[0]->x * (val1/max) + tObj.normals[1]->x * (val2/max) + tObj.normals[2]->x * (val3/max);
  normal.y = tObj.normals[0]->y * (val1/max) + tObj.normals[1]->y * (val2/max) + tObj.normals[2]->y * (val3/max);
  normal.z = tObj.normals[0]->z * (val1/max) + tObj.normals[1]->z * (val2/max) + tObj.normals[2]->z * (val3/max);

  if (tex)
  {
    tex[0] = tObj.texture[0][0] * (val1/max) + tObj.texture[1][0] * (val2/max) + tObj.texture[2][0] * (val3/max);
    tex[1] = tObj.texture[0][1] * (val1/max) + tObj.texture[1][1] * (val2/max) + tObj.texture[2][1] * (val3/max);
  }
}

void getColor(color& surfaceColor, color& lightColor, point* lightPosition, gVector& nh, gVector& npe, point& hitPoint)
{
  double d, s, b;
  
  gVector lightgVector(lightPosition->x-hitPoint.x,lightPosition->y-hitPoint.y,lightPosition->z-hitPoint.z);
  double distance = lightgVector.length();
  lightgVector.scalarMultiply(1.0/distance);
    
  double t;
  bool shadow=false;
  
  for (int k=0; k<numOfMeshes; k++)
  {
    if (planeIntersection(tObjs[k]->vertices[0], &hitPoint, tObjs[k]->triangleNormal, &lightgVector, t))
    {
      point intersect(hitPoint.x+(lightgVector.x*t),hitPoint.y+(lightgVector.y*t),hitPoint.z+(lightgVector.z*t));
      if (liesInside(intersect, *tObjs[k]) && t>0 && t<=distance)
      {
        shadow=true;
        break;
      }
    }
  }
  
  if (shadow)
    return;
  double res = lightgVector.dotProduct(nh);
  if (res<0.2)
    d=0;
  else
  d=clamp((res-MINANGLE)/(MAXANGLE-MINANGLE));

  gVector reflection(-lightgVector.x+2*res*nh.x,-lightgVector.y+2*res*nh.y,-lightgVector.z+2*res*nh.z);
  reflection.scalarMultiply(-1.0/reflection.length());
  double reflectValue = npe.dotProduct(reflection);
  if (reflectValue < 0.95)
    s=0;
  else
    s=clamp((reflectValue-MINREFLECTION)/(MAXREFLECTION-MINREFLECTION));

  double outline = -1.0 * npe.dotProduct(nh);
  b = clamp((outline-MINANGLE)/(MAXANGLE-MINANGLE));
  color diffuse = surfaceColor, specular = surfaceColor, border = surfaceColor;
  diffuse.multiplyColor(lightColor);
  diffuse.multiply(KD*d);
  specular.multiplyColor(lightColor);
  specular.multiply(KS*s);
  border.multiplyColor(lightColor);
  border.multiply(KB*b);
  surfaceColor.addColor(diffuse);
  surfaceColor.addColor(specular);
  surfaceColor.addColor(border);   
}


void applyRasterization()
{
  cout<<"Enter:\n1. Without obj file\n2. With obj file\n";
  int option;
  cin>>option;
  if (option == 1)
  {
    cout << "Enter:\n1. Cube\n2. Tetrahedron\n";
    int value;
    cin>>value;
    initMeshes(value);
    initLight(value);
  }
  else
  {
    cout << "Enter:\n1. without texture\n2. With texture\n";
    int value;
    cin >> value;
    initMeshFromFile(value);
    initLight(3);
  }
  point* testPoint = new point(0,0,0);
  point* lightPos = new point(0,0,0);
  color lightColor(10,10,10,10);
  for(int j=0;j<heightComputed;j++)
  {
    for(int i=0;i<widthComputed;i++)
    {
      double rChannel=0, gChannel=0, bChannel=0;
      double randomValX = ((double)rand() / (double)RAND_MAX)/4;
      double randomValY = ((double)rand() / (double)RAND_MAX)/4;
      for(int m=0;m<4;m++)
      {
        for(int l=0;l<4;l++)
        {
          double x= i + (double)l/4 + randomValX;
	        double y= j + (double)m/4 + randomValY;
      	  testPoint->x = p0->x + (n0->x)*SX*(x/widthComputed) + (n1->x)*SY*(y/heightComputed);
          testPoint->y = p0->y + (n0->y)*SX*(x/widthComputed) + (n1->y)*SY*(y/heightComputed);
          testPoint->z = p0->z + (n0->z)*SX*(x/widthComputed) + (n1->z)*SY*(y/heightComputed);
          gVector* npe = new gVector(testPoint->x-pe->x,testPoint->y-pe->y,testPoint->z-pe->z);
          double fieldLength = npe->length();
          npe->scalarMultiply(((double)1)/fieldLength);
          double tMin = 99999;
          bool changedValue = false;
          double t;
          color ambient;
          int shape=-1;
          lightPos->x=lightPosition->x;
          lightPos->y=lightPosition->y;
          lightPos->z=lightPosition->z;

          for (int k=0; k<numOfMeshes; k++)
          {
            if (planeIntersection(tObjs[k]->vertices[0], pe, tObjs[k]->triangleNormal, npe, t))
            {
              point hitPoint(pe->x+(npe->x*t),pe->y+(npe->y*t),pe->z+(npe->z*t));
              if (liesInside(hitPoint, *tObjs[k]) && (t<tMin || !changedValue) && fieldLength < t)
              {
                shape = k;
                tMin = t;
                changedValue = true;
              }
            }
          }
          if (changedValue)
          {
            int red, green, blue;
            
            point hitPoint(pe->x+(npe->x*tMin),pe->y+(npe->y*tMin),pe->z+(npe->z*tMin));
            gVector normalAtPoint;
            float texture[2];
            computeParameters(*tObjs[shape], hitPoint, normalAtPoint, texture);
            normalAtPoint.scalarMultiply(1.0/normalAtPoint.length());

            if (!tObjs[shape]->colorAdded)
            {
              red = 20, green = 90, blue = 70;
            }
            else
            {
              int xVal = texture[0]*width[0];
              int yVal = texture[1]*height[0];
              int num = ((yVal*width[0])+xVal)*3;
              red = pixmapOrig[0][num++];
              green = pixmapOrig[0][num++];
              blue = pixmapOrig[0][num];
            }
            
            color planeColor(red, green, blue, 255);
            getColor(planeColor, lightColor, lightPos, normalAtPoint, *npe, hitPoint); 
            ambient=planeColor;
            double redVal, greenVal, blueVal;
            ambient.getResult(redVal, greenVal, blueVal);
            rChannel+=redVal;
            gChannel+=greenVal;
            bChannel+=blueVal;
          }
          else
          {
            int red=40, green=30, blue=20;
            rChannel+=red/255.0;
            gChannel+=green/255.0;
            bChannel+=blue/255.0;
          }          
          delete npe;
        }
      }
      int rVal, gVal, bVal;
      rVal = (int)((rChannel/16.0)*255.0);
      gVal = (int)((gChannel/16.0)*255.0);
      bVal = (int)((bChannel/16.0)*255.0);
      setPixelColor(j,i,rVal,gVal,bVal);
    }
  }        
}

// =============================================================================
// main() Program Entry
// =============================================================================
int main(int argc, char *argv[])
{
  char texturePPM[100]="texture.ppm";
  readPPMFile(texturePPM, 0);
  pixmapComputed = new unsigned char[widthComputed * heightComputed * 3];

  initEnvironment();
  applyRasterization();
  generatePPMFile();


  glutInit(&argc, argv);
  glutInitWindowPosition(100, 100); // Where the window will display on-screen.
  glutInitWindowSize(widthComputed, heightComputed);
  glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
  glutCreateWindow("Homework Two");
  init();
  glutReshapeFunc(windowResize);
  glutDisplayFunc(windowDisplay);
  glutMouseFunc(processMouse);
  glutMainLoop();

  return 0; //This line never gets reached. We use it because "main" is type int.
}
     


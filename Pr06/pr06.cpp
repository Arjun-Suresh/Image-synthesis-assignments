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
#define KD 5
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

  point(double xVal, double yVal, double zVal)
  {
    x=xVal, y=yVal, z=zVal;
  }
};


class vector
{  
  public:
  double x;
  double y;
  double z;

  vector(double xVal, double yVal, double zVal)
  {
    x=xVal, y=yVal, z=zVal;
  }
  
  double length()
  {
    double dist = pow(x,2)+pow(y,2)+pow(z,2);
    return sqrt(dist);
  }
  
  vector* operator * (const vector& vObj)
  {
    double xVal = y*vObj.z-z*vObj.y;
    double yVal = z*vObj.x-x*vObj.z;
    double zVal = x*vObj.y-y*vObj.x;
    return new vector(xVal, yVal, zVal);
  }

  double dotProduct(const vector& vObj)
  {
    return x*vObj.x+y*vObj.y+z*vObj.z;
  }
  
  void scalarMultiply(double val)
  {
    x=x*val, y=y*val, z=z*val;
  } 
  
  vector* operator + (const vector& vObj)
  {
    double xVal = x+vObj.x;
    double yVal = y+vObj.y;
    double zVal = z+vObj.z;
    return new vector(xVal, yVal, zVal);
  }
  
  vector* operator - (const vector& vObj)
  {
    double xVal = x-vObj.x;
    double yVal = y-vObj.y;
    double zVal = z-vObj.z;
    return new vector(xVal, yVal, zVal);
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


point *pc, *p0, *pe; //Field of view and eye
point* spheres[2], *cylinder, *plane00; //Shapes to be rendered- spheres and plane
double radius[2]; //radius of 2 spheres
double d; //distance between field of view and the eye
vector *n2, *n1, *n0; //field of view related normal vectors
vector *planeN0, *planeN1, *planeN2; //vectors of an external plane which is rendered

point* lightPosition; //Only for point and spot light
vector* lightDirection; //Only for spot and direction light
point* lightCorner;
vector* lightN0, *lightN1; //For area light

point *ppc, *pp0, *ppe; //Field of view and eye
double pd; //distance between field of view and the eye
vector *pn2, *pn1, *pn0; //field of view related normal vectors

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

bool sphereIntersection(point* center, point* pe, double radius, vector* npe, double& t)
{
  vector* pce = new vector(center->x - pe->x, center->y - pe->y, center->z - pe->z);
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

bool planeIntersection(point* corner, point* pe, vector* normal, vector* npe, double& t)
{
  vector* test = new vector(pe->x - corner->x, pe->y - corner->y, pe->z - corner->z);
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
  n2 = new vector(0,0,-1);
  pe = new point (250, 250, 150);
  pc = new point (pe->x+(n2->x)*d,pe->y+(n2->y)*d,pe->z+(n2->z)*d);
  vector* Vref = new vector(0, 1, 2);
  Vref->scalarMultiply(((double)1)/Vref->length());
  
  n0 = (*n2) * (*Vref);
  n0->scalarMultiply(((double)1)/n0->length());
  n1 = (*n0) * (*n2);
  n1->scalarMultiply(((double)1)/n1->length());
  p0 = new point (pc->x-((SX/2)*(n0->x))-((SY/2)*(n1->x)),pc->y-((SX/2)*(n0->y))-((SY/2)*(n1->y)),pc->z-((SX/2)*(n0->z))-((SY/2)*(n1->z)));

  spheres[0] = new point (325,325,-25);
  spheres[1] = new point (225,280,-155);
  radius[0]=50;
  radius[1]=75;

  planeN2 = new vector(0,-0.996,-0.087);  //The external plane's normal is the normal of view plane by 95 degrees around y axis
  plane00 = new point(50, 500, -100);
  planeN0 = (*planeN2) * (*Vref);
  planeN0->scalarMultiply(1.0/planeN0->length());
  planeN1 = (*planeN0) * (*planeN2);
  planeN1->scalarMultiply(1.0/planeN1->length());


  pd=50;
  pn2 = new vector(-1,0,0);
  ppe = new point (550, 250, -150);
  ppc = new point (ppe->x+(pn2->x)*pd,ppe->y+(pn2->y)*pd,ppe->z+(pn2->z)*pd);
  vector* pVref = new vector(2, 1, 0);
  pVref->scalarMultiply(((double)1)/pVref->length());
  
  pn0 = (*pn2) * (*pVref);
  pn0->scalarMultiply(((double)1)/pn0->length());
  pn1 = (*pn0) * (*pn2);
  pn1->scalarMultiply(((double)1)/pn1->length());
  pp0 = new point (ppc->x-((SX/2)*(pn0->x))-((SY/2)*(pn1->x)),ppc->y-((SX/2)*(pn0->y))-((SY/2)*(pn1->y)),ppc->z-((SX/2)*(pn0->z))-((SY/2)*(pn1->z)));
 
}

void initLight(int option)
{
  if (option == 4)
  {
    lightCorner = new point(150,270,55);
    lightN0 = new vector(1,0,0);
    lightN1 = new vector(0,1,0);
  }
  else
  {
    if(option == 1 || option == 3)
      lightPosition = new point(250,250,-60);
    if(option>1)
    {
      lightDirection = new vector(115,20,-120);
      lightDirection->scalarMultiply(1.0/lightDirection->length());
    }
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

void getColor(color& surfaceColor, color& lightColor, point* lightPosition, vector& nh, vector& npe, point& hitPoint)
{
  double d, s, b;
  
  vector lightVector(lightPosition->x-hitPoint.x,lightPosition->y-hitPoint.y,lightPosition->z-hitPoint.z);
  double distance = lightVector.length();
  lightVector.scalarMultiply(1.0/distance);
    
  double t;
  bool shadow=false;
  if (sphereIntersection(spheres[0], &hitPoint, radius[0], &lightVector, t))
  {
    if (t>0 && t<=distance)
      shadow=true;
  }
  if (!shadow && sphereIntersection(spheres[1], &hitPoint, radius[1], &lightVector, t))
  {
    if (t>0 && t<=distance)
      shadow=true;
  }
  if (!shadow && planeIntersection(plane00, &hitPoint, planeN2, &lightVector, t))
  {
    if (t>0 && t<=distance)
      shadow=true;
  }
  
  if (shadow)
    return;
  double res = lightVector.dotProduct(nh);
  if (res<0.2)
    d=0;
  else
  d=clamp((res-MINANGLE)/(MAXANGLE-MINANGLE));

  vector reflection(-lightVector.x+2*res*nh.x,-lightVector.y+2*res*nh.y,-lightVector.z+2*res*nh.z);
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

/*void getSurfaceColor(point& hitPoint, int shape, vector& npe, int& red, int& green, int& blue, vector* normalAdd)
{
  double xCord, yCord;
  if (shape == 0 || shape == 1)
  {
    vector p0h(hitPoint.x - spheres[shape]->x, hitPoint.y - spheres[shape]->y, hitPoint.z - spheres[shape]->z);
    p0h.scalarMultiply(1.0/(double)p0h.length());
    vector sphereN2(0,-1,0);
    vector sphereN1(0,0,1);
    vector sphereN0(1,0,0);
    double z = sphereN2.dotProduct(p0h);
    double y = sphereN1.dotProduct(p0h);
    double x = sphereN0.dotProduct(p0h);
    
    xCord = acos(z);
    double val = y/sin(xCord);
    if (val < -1)
      val=-1;
    if (val>1)
      val=1;
    yCord = acos(val); 
    if(x<0)
      yCord = 2*3.14159 - yCord;
    xCord = xCord/3.14159;
    yCord = yCord/(2*3.14159);
    int xVal = (int)(xCord*width[3]);
    int yVal = height[3]-1-(int)(yCord*height[3]);
    int num = ((yVal*width[3])+xVal)*3;
    double red = (double)pixmapOrig[3][num++]/255.0;
    double green = (double)pixmapOrig[3][num++]/255.0;
    double blue = (double)pixmapOrig[3][num]/255.0;
    normalAdd->x = (2.0*red - 1.0)*sphereN0.x+(2.0*green - 1.0)*sphereN1.x+(2.0*blue - 1.0)*sphereN2.x;
    normalAdd->y = (2.0*red - 1.0)*sphereN0.y+(2.0*green - 1.0)*sphereN1.y+(2.0*blue - 1.0)*sphereN2.y;
    normalAdd->z = (2.0*red - 1.0)*sphereN0.z+(2.0*green - 1.0)*sphereN1.z+(2.0*blue - 1.0)*sphereN2.z;
    normalAdd->scalarMultiply(1.0/normalAdd->length());
  }
  else if (shape == 2)
  {
    vector p0h(hitPoint.x - floor(hitPoint.x), hitPoint.y - floor(hitPoint.y), hitPoint.z - floor(hitPoint.z));
    p0h.scalarMultiply(1.0/p0h.length());
    xCord = p0h.dotProduct(*planeN0);
    yCord = p0h.dotProduct(*planeN1);
    if (xCord<0)
      xCord+=1;
    if (yCord<0)
      yCord+=1;
    
  }
  else
  {
    vector sphereN2(0,-1,0);
    vector sphereN1(0,0,1);
    vector sphereN0(1,0,0);
    double z = sphereN2.dotProduct(npe);
    double y = sphereN1.dotProduct(npe);
    double x = sphereN0.dotProduct(npe);
    xCord = acos(z);
    double val = y/sin(xCord);
    if (val < -1)
      val=-1;
    if (val>1)
      val=1;
    yCord = acos(val);
    if(x<0)
      yCord = 2*3.14159 - yCord;
    xCord = xCord/3.14159;
    yCord = yCord/(2*3.14159);
  }
  if (shape>0)
    shape--;
  xCord=xCord*width[shape];
  yCord=yCord*height[shape];
  int iCord = floor(xCord+0.5)-1;
  int jCord = floor(yCord+0.5)-1;
  double u = xCord - ((double)(iCord)+0.5);
  double v = yCord - ((double)(jCord)+0.5);
  
  if (iCord < 0)
    iCord+=width[shape];
  if (jCord < 0)
    jCord+=height[shape];
  jCord=height[shape]-1-jCord;
  int num1 = (((jCord%height[shape])*width[shape])+(iCord%width[shape]))*3;
  int red1=pixmapOrig[shape][num1++];
  int green1=pixmapOrig[shape][num1++];
  int blue1=pixmapOrig[shape][num1];

  int num2 = ((((jCord+1)%height[shape])*width[shape])+(iCord%width[shape]))*3;
  int red2=pixmapOrig[shape][num2++];
  int green2=pixmapOrig[shape][num2++];
  int blue2=pixmapOrig[shape][num2];

  int num3 = ((((jCord+1)%height[shape])*width[shape])+((iCord+1)%width[shape]))*3;
  int red3=pixmapOrig[shape][num3++];
  int green3=pixmapOrig[shape][num3++];
  int blue3=pixmapOrig[shape][num3];

  int num4 = (((jCord%height[shape])*width[shape])+((iCord+1)%width[shape]))*3;
  int red4=pixmapOrig[shape][num4++];
  int green4=pixmapOrig[shape][num4++];
  int blue4=pixmapOrig[shape][num4];
    
  red = (1.0-u)*(1.0-v)*red1 + u*(1.0-v)*red4 + (1.0-u)*v*red2 + u*v*red3;
  green = (1.0-u)*(1.0-v)*green1 + u*(1.0-v)*green4 + (1.0-u)*v*green2 + u*v*green3;
  blue = (1.0-u)*(1.0-v)*blue1 + u*(1.0-v)*blue4 + (1.0-u)*v*blue2 + u*v*blue3;
}*/

bool withinBound(point& intersect, int category)
{
  point *corner;
  vector *nx, *ny;
  if (category == 1)
  {
    corner = p0;
    nx = n0;
    ny = n1;
  }
  else
  {
    corner = pp0;
    nx = pn0;
    ny = pn1;
  }
  point* p[4];
  p[0] = new point(corner->x, corner->y, corner->z);
  p[1] = new point(corner->x+ n0->x*SX, corner->y+ n0->y*SX, corner->z+ n0->z*SX);
  p[2] = new point(corner->x+ n1->x*SY, corner->y+ n1->y*SY, corner->z+ n1->z*SY);
  p[3] = new point(corner->x+ n0->x*SX + n1->x*SY, corner->y+ n0->y*SX + n1->y*SY, corner->z+ n0->z*SX + n1->z*SY);
  int xmin=p[0]->x, xmax=p[0]->x, ymin=p[0]->y, ymax=p[0]->y, zmin=p[0]->z, zmax=p[0]->z;
  for (int i=1;i<4;i++)
  {
    if (xmin > p[i]->x)
      xmin = p[i]->x;
    if (xmax < p[i]->x)
      xmax = p[i]->x;
    if (ymin > p[i]->y)
      ymin = p[i]->y;
    if (ymax < p[i]->y)
      ymax = p[i]->y;
    if (zmin > p[i]->z)
      zmin = p[i]->z;
    if (zmax < p[i]->z)
      zmax = p[i]->z;
  }
  if (intersect.x >= xmin && intersect.x <= xmax &&intersect.y >= ymin && intersect.y <= ymax &&intersect.z >= zmin && intersect.z <= zmax)
    return true;
  return false;
}

bool texturize(point* hitpoint, int& red,int& green, int& blue, int option, vector* npe=NULL, bool finite=true)
{
  point intersect(0,0,0);
  if (finite)
  {
    vector ray(0,0,0);
    if (option == 1)
    {
      if (hitpoint->x > ppe->x)
        return false;
      ray.x = hitpoint->x - ppe->x;
      ray.y = hitpoint->y - ppe->y;
      ray.z = hitpoint->z - ppe->z;
      ray.scalarMultiply(1.0/ray.length());
    }
    else if (option == 2)
    {
      ray.x = -1.0 * pn2->x;
      ray.y = -1.0 * pn2->y;
      ray.z = -1.0 * pn2->z;
    }
  
    vector pn2Comp(-1.0 * pn2->x, -1.0 * pn2->y, -1.0 * pn2->z);
    double t=0.0;
    if (planeIntersection(pp0, hitpoint, &pn2Comp, &ray, t))
    {
      intersect.x = hitpoint->x+t*ray.x;
      intersect.y = hitpoint->y+t*ray.y;
      intersect.z = hitpoint->z+t*ray.z;
      if (!withinBound(intersect,1))
        return false;
    }
    else
      return false; 
  }
  else
  {
    double t=0.0;
    if (planeIntersection(p0, pe, n2, npe, t))
    {
      intersect.x = p0->x+t*npe->x;
      intersect.y = p0->y+t*npe->y;
      intersect.z = p0->z+t*npe->z;
      if (!withinBound(intersect,2))
        return false;
    }
    else
      return false; 
  } 
  vector vTemp(intersect.x-pp0->x,intersect.y-pp0->y,intersect.z-pp0->z);
  int y = (int)(((vTemp.dotProduct(*pn1))*(double)height[0])/SY);
  int x = (int)(((vTemp.dotProduct(*pn0))*(double)width[0])/SX);
  int loc = ((y*width[0])+x)*3;
  red = pixmapOrig[0][loc++];
  green = pixmapOrig[0][loc++];
  blue = pixmapOrig[0][loc];
  return true;
} 
  
void applyRasterization()
{
  cout<<"Enter:\n1. Perspective image projection\n2. Parallel image projection\n";
  int option;
  cin>>option;
  initLight(1);
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
          vector* npe = new vector(testPoint->x-pe->x,testPoint->y-pe->y,testPoint->z-pe->z);
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

          if (sphereIntersection(spheres[0], pe, radius[0], npe, t))
          {
            if (t<tMin || !changedValue)
            {
              shape=0;
              tMin=t;
              changedValue=true;
            }
          }
          if (sphereIntersection(spheres[1], pe, radius[1], npe, t))
          {
            if (t<tMin || !changedValue)
            {
              shape=1;
              tMin=t;
              changedValue=true;
            }
          }
          if (planeIntersection(plane00, pe, planeN2, npe, t))
          {
            if ((t<tMin || !changedValue) && fieldLength < t)
            {
     	      shape=2;
              tMin=t;
              changedValue=true;
            }
          }
          if (shape == 0)
          {
            int red, green, blue;
            point hitPoint(pe->x+(npe->x*tMin),pe->y+(npe->y*tMin),pe->z+(npe->z*tMin));
            vector nh(hitPoint.x-spheres[0]->x,hitPoint.y-spheres[0]->y,hitPoint.z-spheres[0]->z);
            nh.scalarMultiply(1.0/nh.length());
            texturize(&hitPoint, red, green, blue, option);                   
            ambient.setColor(red, green, blue, 255);
            getColor(ambient, lightColor, lightPos, nh, *npe, hitPoint);
          }
          if (shape == 1)
          {
            int red, green, blue;
            point hitPoint(pe->x+(npe->x*tMin),pe->y+(npe->y*tMin),pe->z+(npe->z*tMin));
            vector nh(hitPoint.x-spheres[1]->x,hitPoint.y-spheres[1]->y,hitPoint.z-spheres[1]->z);
            nh.scalarMultiply(1.0/nh.length());
            texturize(&hitPoint, red, green, blue, option);    
                      
            ambient.setColor(red, green, blue, 255);
            getColor(ambient, lightColor, lightPos, nh, *npe, hitPoint);
          }
          if (shape == 2)
          {
            int red, green, blue;
            point hitPoint(pe->x+(npe->x*tMin),pe->y+(npe->y*tMin),pe->z+(npe->z*tMin));
            texturize(&hitPoint, red, green, blue, option);    
            color planeColor(red,green,blue,255);
            getColor(planeColor, lightColor, lightPos, *planeN2, *npe, hitPoint); 
            ambient=planeColor;
          }
          if (shape!=-1)
          {
            double red, green, blue;
            ambient.getResult(red, green, blue);
            rChannel+=red;
            gChannel+=green;
            bChannel+=blue;
          }
          else
          {
            int red, green, blue;
            point dummyPoint(0,0,0);
            texturize(NULL, red, green, blue, option, npe, false);    
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
  char spherePPM[100], planePPM[100], infiniteSpherePPM[100], normalMap[100];
  cout<<"Enter the PPM file name for objects texture\n";
  cin>>spherePPM;
  readPPMFile(spherePPM, 0);
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
     


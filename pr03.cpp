// =============================================================================
// VIZA656/CSCE647 at Texas A&M University
// Homework 2
// Basic camera
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

#define SPOTLIGHTMIN 0.707
#define MINANGLE 0.0
#define MAXANGLE 1.0
#define MINREFLECTION 0.707
#define MAXREFLECTION 1.0
#define KD 5
#define KS 5
#define KB 1
using namespace std;
// =============================================================================
// These variables will store the input ppm image's width, height, and color
// =============================================================================
int width, height, maxColorValue=255, magicNo, widthComputed=750, heightComputed=750;
unsigned char *pixmapOrig, *pixmapComputed;

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
    
    void multiply(double t)
    {
      red = red*t;
      green = green*t;
      blue = blue*t;
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
point* spheres[2], *cylinder, *planeCorner; //Shapes to be rendered- spheres and plane
double radius[2]; //radius of 2 spheres
double d; //distance between field of view and the eye
vector *n2, *n1, *n0; //field of view related normal vectors
vector *planeNormal; //Normal of an external plane which is rendered

point* lightPosition; //Only for point and spot light
vector* lightDirection; //Only for spot and direction light

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

void setPixelColorOrig(long int &val, long int& index, unsigned char* fileBuffer)
{
  pixmapOrig[val++]=fileBuffer[index++];
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

bool parseValue(long int& index, unsigned char* fileBuffer, long int numOfCharacters, int& parseValue, int valueId)
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
      width=value;
      break;
    case HEIGHTVALUEID:
      height=value;
      break;
    case MAXCOLORVALUEID:
      maxColorValue=value;
      break;
  }
  return true;
}

bool fillPixelsBin(long int& index, unsigned char* fileBuffer, long int numOfCharacters)
{
  int wVal,hVal;
  wVal=width;
  hVal=height;
  long int rowVal=0,colVal=0;
  while(rowVal<hVal)
  {
    long int val=(((hVal-rowVal-1)*wVal)+colVal)*3;
    for(int i=0;i<3;i++)
    	setPixelColorOrig(val, index, fileBuffer);
    colVal=(colVal+1)%wVal;
    if(!colVal)
      rowVal++;
  }
  return true;
}

bool readPPMFile(char* filePath)
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
      if(!parseValue(index, fileBuffer, numOfCharacters, widthParsed, WIDTHVALUEID))
      {
        return false;
      }
      continue;
    }
    if(magicNumberParsed && widthParsed && !heightParsed)
    {
      if(!parseValue(index, fileBuffer, numOfCharacters, heightParsed, HEIGHTVALUEID))
      {
        return false;
      }
      continue;
    }
    if(magicNumberParsed && widthParsed && heightParsed && !maxColorValueParsed)
    {
      if(!parseValue(index, fileBuffer, numOfCharacters, maxColorValueParsed, MAXCOLORVALUEID))
      {
        return false;
      }
      continue;
    }
    if(magicNumberParsed && widthParsed && heightParsed && maxColorValueParsed)
    {
      pixmapOrig = new unsigned char[width*height*3];
      if(!(fillPixelsBin(index, fileBuffer, numOfCharacters)))
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

bool sphereIntersection(point* center, double radius, vector* npe, double& t)
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

bool planeIntersection(point* corner, vector* normal, vector* npe, double& t)
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
  spheres[1] = new point (225,250,-50);
  radius[0]=50;
  radius[1]=75;

  planeNormal = new vector(-0.996,0,-0.087);  //The external plane's normal is the normal of view plane by 95 degrees around y axis
  planeCorner = new point(325, 325, -100);
}

void initLight(int option)
{
  if(option == 1 || option == 3)
    lightPosition = new point(225,250,50);
  if(option>1)
  {
    lightDirection = new vector(50,30,-25);
    lightDirection->scalarMultiply(1.0/lightDirection->length());
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

void getColor(color& surfaceColor, color& lightColor, vector& nh, vector& npe, point& hitPoint, int option)
{
  double d, s, b;
  if(option == 3)
  {
    vector lightVector(hitPoint.x-lightPosition->x,hitPoint.y-lightPosition->y,hitPoint.z-lightPosition->z);
    lightVector.scalarMultiply(1.0/lightVector.length());
    if (lightVector.dotProduct(*lightDirection)<SPOTLIGHTMIN)
      return;
    else
    {
      lightVector.scalarMultiply(-1);

      double res = lightVector.dotProduct(nh);
      d=clamp((res-MINANGLE)/(MAXANGLE-MINANGLE));

      vector reflection(-lightVector.x+2*res*nh.x,-lightVector.y+2*res*nh.y,-lightVector.z+2*res*nh.z);
      reflection.scalarMultiply(-1.0/reflection.length());
      s=clamp((npe.dotProduct(reflection)-MINANGLE)/(MAXANGLE-MINANGLE));
        
      double outline = -1.0 * npe.dotProduct(nh);
      b = clamp((outline-MINANGLE)/(MAXANGLE-MINANGLE));
    }
  }
  else if (option == 2)
  {
    vector lightVector(-lightDirection->x, -lightDirection->y, -lightDirection->z);
    double res = lightVector.dotProduct(nh);
    d=clamp((res-MINANGLE)/(MAXANGLE-MINANGLE));
 
    vector reflection(-lightVector.x+2*res*nh.x,-lightVector.y+2*res*nh.y,-lightVector.z+2*res*nh.z);
    reflection.scalarMultiply(-1.0/reflection.length());
    s=clamp((npe.dotProduct(reflection)-MINANGLE)/(MAXANGLE-MINANGLE));
      
    double outline = -1.0 * npe.dotProduct(nh);
    b = clamp((outline-MINANGLE)/(MAXANGLE-MINANGLE));
  }
  
  else
  {
    vector lightVector(lightPosition->x-hitPoint.x,lightPosition->y-hitPoint.y,lightPosition->z-hitPoint.z);
    lightVector.scalarMultiply(1.0/lightVector.length());
    
    double res = lightVector.dotProduct(nh);
    d=clamp((res-MINANGLE)/(MAXANGLE-MINANGLE));

    vector reflection(-lightVector.x+2*res*nh.x,-lightVector.y+2*res*nh.y,-lightVector.z+2*res*nh.z);
    reflection.scalarMultiply(-1.0/reflection.length());
    s=clamp((npe.dotProduct(reflection)-MINREFLECTION)/(MAXREFLECTION-MINREFLECTION));

    double outline = -1.0 * npe.dotProduct(nh);
    b = clamp((outline-MINANGLE)/(MAXANGLE-MINANGLE));
  }
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
  int lightOption;
  cout<<"Enter\n1. Point Light\n2. Direction light\n3. Spot light\n";
  cin>>lightOption;
  initLight(lightOption);
  point* testPoint = new point(0,0,0);
  color lightColor(10,10,10,10);
  for(int j=0;j<heightComputed;j++)
  {
    for(int i=0;i<widthComputed;i++)
    {
      int r = pixmapOrig[((((heightComputed-j-1)%height)*width)+(i%width))*3];
      int g = pixmapOrig[(((((heightComputed-j-1)%height)*width)+(i%width))*3)+1];
      int b = pixmapOrig[(((((heightComputed-j-1)%height)*width)+(i%width))*3)+2];
      double rChannel=0, gChannel=0, bChannel=0;
      for(int m=0;m<4;m++)
      {
        for(int l=0;l<4;l++)
        {
          double x= i + (double)l/4 + ((double)rand() / (double)RAND_MAX)/4;
	  double y= j + (double)m/4 + ((double)rand() / (double)RAND_MAX)/4;
      	  testPoint->x = p0->x + (n0->x)*SX*(x/widthComputed) + (n1->x)*SY*(y/heightComputed);
          testPoint->y = p0->y + (n0->y)*SX*(x/widthComputed) + (n1->y)*SY*(y/heightComputed);
          testPoint->z = p0->z + (n0->z)*SX*(x/widthComputed) + (n1->z)*SY*(y/heightComputed);
          vector* npe = new vector(testPoint->x-pe->x,testPoint->y-pe->y,testPoint->z-pe->z);
          double fieldLength = npe->length();
          npe->scalarMultiply(((double)1)/fieldLength);
          double tMin = 99999;
          bool changedValue = false;
          double t;
          color ambient(r,g,b,maximum(r,g,b));
          int shape=-1;
          if (sphereIntersection(spheres[0], radius[0], npe, t))
          {
            if (t<tMin || !changedValue)
            {
              shape=0;
              tMin=t;
              changedValue=true;
            }
          }
          if (sphereIntersection(spheres[1], radius[1], npe, t))
          {
            if (t<tMin || !changedValue)
            {
              shape=1;
              tMin=t;
              changedValue=true;
            }
          }
          if (planeIntersection(planeCorner, planeNormal, npe, t))
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
            point hitPoint(pe->x+(npe->x*tMin),pe->y+(npe->y*tMin),pe->z+(npe->z*tMin));
            vector nh(hitPoint.x-spheres[0]->x,hitPoint.y-spheres[0]->y,hitPoint.z-spheres[0]->z);
            nh.scalarMultiply(1.0/nh.length());
            getColor(ambient, lightColor, nh, *npe, hitPoint, lightOption);
          }
          if (shape == 1)
          {
            point hitPoint(pe->x+(npe->x*tMin),pe->y+(npe->y*tMin),pe->z+(npe->z*tMin));
            vector nh(hitPoint.x-spheres[1]->x,hitPoint.y-spheres[1]->y,hitPoint.z-spheres[1]->z);
            nh.scalarMultiply(1.0/nh.length());
            getColor(ambient, lightColor, nh, *npe, hitPoint, lightOption);
          }
          if (shape == 2)
          {
            color planeColor(20,130,40,255);
            point hitPoint(pe->x+(npe->x*tMin),pe->y+(npe->y*tMin),pe->z+(npe->z*tMin));
            getColor(planeColor, lightColor, *planeNormal, *npe, hitPoint, lightOption); 
            ambient=planeColor;
          }
          delete npe;
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
            rChannel+=150.0/255.0;
            gChannel+=100.0/255.0;
            bChannel+=180.0/255.0;
          }
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
  char inputPPMFile[100];
  cout<<"Enter the goal image ppm file name\n";
  cin>>inputPPMFile;
  readPPMFile(inputPPMFile);
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
     


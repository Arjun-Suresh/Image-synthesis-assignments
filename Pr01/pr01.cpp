// =============================================================================
// VIZA656/CSCE647 at Texas A&M University
// Homework 1
// 3D to 2D Raster conversion
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

using namespace std;
// =============================================================================
// These variables will store the input ppm image's width, height, and color
// =============================================================================
int width=750, height=750, maxColorValue=255, magicNo;
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
    return pow((pow(x,2)+pow(y,2)+pow(z,2)),0.5);
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


point* planeCorner, *planeCorner2;
point* spheres[2], *cylinder;
double radius[2];
vector *n2, *n1, *n0, *planeNormal2;


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
  int i = ((height-y-1) * width + x) * 3; 
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
  glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, pixmapComputed);
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
  sprintf(widthString, "%d", width);
  sprintf(heightString, "%d", height);
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
  for(int i=0;i<height;i++)
  {
    for(int j=0;j<width;j++)
    {
      if(index >= origSize-50)
      {
        unsigned char* tempBuffer = resizeArray(fileBuffer, origSize, origSize);
        delete[] fileBuffer;
        fileBuffer=tempBuffer;
      }
      int k=((height-i-1)*width+j)*3;
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

bool isInside(point* center, double radius, point* testPoint)
{
  if ((pow((testPoint->x - center->x),2) + pow((testPoint->y - center->y),2) + pow((testPoint->z - center->z),2) - pow (radius,2)) <= 0)
  {
    return true;
  }
  return false;
}

bool isInsidePlane(point* center, vector* n2, point* testPoint)
{
  vector* test = new vector(testPoint->x - center->x, testPoint->y - center->y, testPoint->z - center->z);
  if (n2->dotProduct(*test) <= 0)
    return true;
  return false;
}

bool insideCylinder(point* center, point* testPoint)
{
  vector* test = new vector(testPoint->x - center->x, testPoint->y - center->y, testPoint->z - center->z);
  double result = 15* pow((n0->dotProduct(*test)/(double)SX),2)+10* pow((n1->dotProduct(*test)/(double)SY),2) - 1;
  if (result<=0)
    return true;
  return false;
}

void initEnvironment()
{
  n2 = new vector(0,1,1);
  n2->scalarMultiply(((double)1)/n2->length());
  vector* Vref = new vector(0, 10, 5);
  Vref->scalarMultiply(((double)1)/Vref->length());
  n0 = (*n2) * (*Vref);
  n0->scalarMultiply(((double)1)/n0->length());
  n1 = (*n2) * (*n0);
  n1->scalarMultiply(((double)1)/n1->length());
  planeCorner = new point(100,100,100);
  spheres[0] = new point (80,40,160);
  spheres[1] = new point (60,45,160);
  radius[0]=10;
  radius[1]=15;

  cylinder = new point (50, 25, 60);

  planeNormal2 = new vector(0,-1,1);
  planeNormal2->scalarMultiply(((double)1)/planeNormal2->length());
  planeCorner2 = new point(0, 0, 100);
}

void applyRasterization()
{
  initEnvironment();
  point* testPoint = new point(0,0,0);
  for(int j=0;j<height;j++)
  {
    for(int i=0;i<width;i++)
    {
      double redChannel=0, greenChannel=0, blueChannel=0;
      for(int m=0;m<4;m++)
      {
        for(int l=0;l<4;l++)
        {
          double x= i + (double)l/4 + ((double)rand() / (double)RAND_MAX)/4;
	  double y= j + (double)m/4 + ((double)rand() / (double)RAND_MAX)/4;
      	  testPoint->x = planeCorner->x + (n0->x)*SX*(x/width) + (n1->x)*SY*(y/height);
          testPoint->y = planeCorner->y + (n0->y)*SX*(x/width) + (n1->y)*SY*(y/height);
          testPoint->z = planeCorner->z + (n0->z)*SX*(x/width) + (n1->z)*SY*(y/height);
          if (insideCylinder(cylinder, testPoint))
          {
            redChannel += (120.0/(255.0*16.0));
            greenChannel +=(65.0/(255.0*16.0));
            blueChannel +=(25.0/(255.0*16.0));
          }
          else if (isInside(spheres[0], radius[0], testPoint))
          {
            redChannel += (80.0/(255.0*16.0));
            greenChannel +=(90.0/(255.0*16.0));
            blueChannel +=(40.0/(255.0*16.0));
          }
          else if (isInside(spheres[1], radius[1], testPoint))
          {
            redChannel += (65.0/(255.0*16.0));
            greenChannel +=(75.0/(255.0*16.0));
            blueChannel +=(25.0/(255.0*16.0));
          }
          else if (isInsidePlane(planeCorner2, planeNormal2, testPoint))
          {
            redChannel += (20.0/(255.0*16.0));
            greenChannel +=(30.0/(255.0*16.0));
            blueChannel +=(40.0/(255.0*16.0));            
          }
          else
          {
            redChannel += (20.0/(255.0*16.0));
            greenChannel +=(70.0/(255.0*16.0));
            blueChannel +=(50.0/(255.0*16.0));
          }
        }
      }
      setPixelColor(j,i,(int)(redChannel*255.0),(int)(greenChannel*255.0),(int)(blueChannel*255.0));
    }
  }       
}

// =============================================================================
// main() Program Entry
// =============================================================================
int main(int argc, char *argv[])
{
  char inputPPMFile[100], controlPPMFile[100];
  pixmapComputed = new unsigned char[width * height * 3];

  applyRasterization();

  generatePPMFile();


  glutInit(&argc, argv);
  glutInitWindowPosition(100, 100); // Where the window will display on-screen.
  glutInitWindowSize(width, height);
  glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
  glutCreateWindow("Homework Four");
  init();
  glutReshapeFunc(windowResize);
  glutDisplayFunc(windowDisplay);
  glutMouseFunc(processMouse);
  glutMainLoop();

  return 0; //This line never gets reached. We use it because "main" is type int.
}
     


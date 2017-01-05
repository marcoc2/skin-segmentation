/*
  Name:        PixelMap.h
  Copyright:   Version 0.1
  Author:      Rodrigo Luis de Souza da Silva
  Date:        04/07/2009
  Description: Simple class to load BMP images and modify its data.
               You can create an image manually and convert it to a
               pixel matrix and vice-versa.
*/

//////////////////////////////////////////////////////////////////
// Preprocessor directives
#define GLUT_NO_LIB_PRAGMA
#define GLUT_NO_WARNING_DISABLE
#define GLUT_DISABLE_ATEXIT_HACK

#define FILTER_X  1
#define FILTER_Y  2
#define FILTER_XY 3

#define FILTER_3x3 1
#define FILTER_5x5 2

//////////////////////////////////////////////////////////////////
// New types
typedef char Byte;
typedef unsigned char uByte;

//////////////////////////////////////////////////////////////////
// Classes

// Simple pixel Structure
struct pixel{
   uByte value;      // Grayscale value
   uByte R, G, B, A; // Color RGBA values
};

typedef struct _COMPLEX{
	double real;
	double imag;
} COMPLEX, *pCOMPLEX;

// Pixel map class
class PixelMap{
   public:
      // Constructors and Destructors
      PixelMap(const char *fname);
      PixelMap();
      ~PixelMap();

      // Read image from a bitmap file
      void Read(const char *fname);

      // Create image from a 'uByte' array
      void CreateImage(int w, int h, uByte *d);

      // view stored image
      void ViewImage();
      void ViewOldImage();

      // Modify and get pixel information
      uByte  pixel_elem(int x, int y, int elem);
      uByte *pixel_pos(int x, int y);

      // Set/Get image's width and height
      int  GetWidth();
      void SetWidth( int);
      int  GetHeight();
      void SetHeight( int);

      // Set/Get image data
      uByte* GetData();
      uByte* GetOldData();
      void   SetData(uByte*);

      // Convertion methods
      void ConvertPixelsToData(pixel **m);
      void ConvertDataToPixels(pixel **m);
      void ConvertToGrayScale();

      //Transformações implementadas
      void fourierTransform();
      void ifourierTransform();

      //Filtros implementados
      void filtroMedia(int n);
      void filtroMediana(int n);

      void ruidoSalPimenta(double percentual);
      void ruidoPersonalzado(double percentual, double increase);


      void reescreve();
      void removeRuido(int x, int y, int raio);
      void adicionaRuido(int x, int y, int raio);


      int  FFT2D(int inverse);
      int  FFT(int dir,int m,double *x,double *y);
      int  Powerof2(int nx, int *m, int *twopm);

      ///Filtros
      void passaAlta(int R);
      void passaBaixa(int R);

   protected:
      int width, height;
      uByte *data;
      uByte *oldData;

      float *imgRE;
      float *imgIG;

      float *fourierRE;
      float *fourierIG;

      int inverse;

      bool grayScale;
};

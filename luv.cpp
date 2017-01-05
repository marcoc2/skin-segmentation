#include <stdlib.h>
#include <stdio.h>
#include "luv.h"

using namespace std;

int scale = 255;
IplImage* imagem;
IplImage* imagem2;

int RGBToLUV(IplImage* bitmap, int step, int nChannels)
{

   float rgbToxyz[3][3] = {
     { 0.490f, 0.310f, 0.200f },
     { 0.177f, 0.812f, 0.011f },
     { 0.000f, 0.010f, 0.990f }
   };

   float Xn, Yn, Zn, Un, Vn;

   Xn = rgbToxyz[0][0] + rgbToxyz[0][1] + rgbToxyz[0][2];
   Yn = rgbToxyz[1][0] + rgbToxyz[1][1] + rgbToxyz[1][2];
   Zn = rgbToxyz[2][0] + rgbToxyz[2][1] + rgbToxyz[2][2];

   Un = (4.0f * Xn) / (Xn + 15.0f * Yn + 3.0f * Zn);
   Vn = (9.0f * Yn) / (Xn + 15.0f * Yn + 3.0f * Zn);

//   X = 0.412453*R + 0.35758 *G + 0.180423*B
//   Y = 0.212671*R + 0.71516 *G + 0.072169*B
//   Z = 0.019334*R + 0.119193*G + 0.950227*B
//   http://software.intel.com/sites/products/documentation/hpc/ipp/ippi/ippi_ch6/functn_RGBToXYZ.html

//   BMP* channelR = new float[width*height];
//   getChannel(source, 0, channelR);
//   BMP* channelG = new float[width*height];
//   ImageBase::GetChannel(source, 1, channelG);
//   float* channelB = new float[width*height];
//   ImageBase::GetChannel(source, 2, channelB);

   for (unsigned int i = 0; i < bitmap->height; i++)
   {
     for (unsigned int j = 0; j < bitmap->width; j++)
     {

       float s[3] = {
           (float) bitmap->imageData[i*step + j*nChannels+2] / (float) scale,
           (float) bitmap->imageData[i*step + j*nChannels+1] / (float) scale,
           (float) bitmap->imageData[i*step + j*nChannels] / (float) scale
           };

//       dest.SetPixel(s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2]
//       * rgbToxyz[0][2], s[0] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1]
//       + s[2] * rgbToxyz[1][2], s[0] * rgbToxyz[2][0] + s[1]
//       * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2], j, i);

//       cout << " | " << (int) bitmap->pixel_data[i * bitmap->width + j].red << " , "
//       << (int) bitmap->pixel_data[i * bitmap->width + j].green
//       << " , " << (int) bitmap->pixel_data[i * bitmap->width + j].blue;

//       cout << " | " << s[0] << " , " << s[1] << " , " << s[2];

//       bitmap->pixel_data[i * bitmap->width + j].red = s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2] * rgbToxyz[0][2];
//       bitmap->pixel_data[i * bitmap->width + j].blue = s[1] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1] + s[2] * rgbToxyz[1][2];
//       bitmap->pixel_data[i * bitmap->width + j].green = s[2] * rgbToxyz[2][0] + s[1] * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2];

//       cout << " | " << (int) bitmap->pixel_data[i * bitmap->width + j].red << " , "
//       << (int) bitmap->pixel_data[i * bitmap->width + j].green
//       << " , " << (int) bitmap->pixel_data[i * bitmap->width + j].blue;

       // now transform xyz to luv
       float x = s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2] * rgbToxyz[0][2];
       float y = s[0] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1] + s[2] * rgbToxyz[1][2];
       float z = s[0] * rgbToxyz[2][0] + s[1] * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2];

       // compute L*
       float L = y / Yn;

       if (L > 0.008856)
         {
           L = 116.0f * pow(L, 1.0f / 3.0f) - 16.0f;
         }
       else
         {
           L = 903.3f * L;
         }

       //    //compute u* and v*
       float denominator = x + 15.0f * y + 3.0f * z;
       float u_star = 0;
       float v_star = 0;
       if (denominator != 0)
         {
           u_star = 13.0f * L * ((4.0f * x / denominator) - Un);
           v_star = 13.0f * L * ((9.0f * y / denominator) - Vn);

         }

       bitmap->imageData[i*step + j*nChannels+2] = L;
       bitmap->imageData[i*step + j*nChannels+1] = u_star;
       bitmap->imageData[i*step + j*nChannels] = v_star;
       }
    }



   return 0;
 }

int LUVToRGB(IplImage* bitmap, int step, int nChannels){

  float rgbToxyz[3][3] =
  {
    { 0.490f, 0.310f, 0.200f },
    { 0.177f, 0.812f, 0.011f },
    { 0.000f, 0.010f, 0.990f } };

  float Xn, Yn, Zn, Un, Vn;

  Xn = rgbToxyz[0][0] + rgbToxyz[0][1] + rgbToxyz[0][2];
  Yn = rgbToxyz[1][0] + rgbToxyz[1][1] + rgbToxyz[1][2];
  Zn = rgbToxyz[2][0] + rgbToxyz[2][1] + rgbToxyz[2][2];

  Un = (4.0f * Xn) / (Xn + 15.0f * Yn + 3.0f * Zn);
  Vn = (9.0f * Yn) / (Xn + 15.0f * Yn + 3.0f * Zn);

  for (unsigned int i = 0; i < bitmap->height; i++)
    {
      for (unsigned int j = 0; j < bitmap->width; j++)
        {
          float luv[3] = {
               ((float) bitmap->imageData[i*step + j*nChannels+2] / (float) scale) * 100,
               (float) bitmap->imageData[i*step + j*nChannels+1] / (float) scale,
               (float) bitmap->imageData[i*step + j*nChannels] / (float) scale};

          //-------------------------------
          // transform luv to xyz

          float L = luv[0];
          float u_star = luv[1];
          float v_star = luv[2];
          float x = 0;
          float y = 0;
          float z = 0;

          if (L > 7.9996248)
            {
              y = Yn * pow((float) ((L + 16.0) / 116.0), 3);
            }
          else
            {
              y = Yn * (L / 903.3f);
            }

          if (L > 0)
            {
              float u = u_star / (13.0f * L) + Un;
              float v = v_star / (13.0f * L) + Vn;
              z = 3.0f * y / v - 5.0f * y - (3.0f / 4.0f) * u * y / v;
              x = 9.0f * y / v - 15.0f * y - 3.0f * z;
            }

//          bitmap->pixel_data[i * bitmap->width + j].red = x;
//          bitmap->pixel_data[i * bitmap->width + j].blue = y;
//          bitmap->pixel_data[i * bitmap->width + j].green = z;

       float xyzTorgb[3][3] =
     {
       { 2.3649f, -0.8971f, -0.4678f },
       { -0.5156f, 1.4273f, 0.0883f },
       { 0.0052f, -0.0144f, 1.0092f } };

           float s[3] =
             { (float) x, y, z};

           bitmap->imageData[i*step + j*nChannels+2] = ((s[0] * xyzTorgb[0][0] + s[1] * xyzTorgb[0][1] + s[2] * xyzTorgb[0][2])*scale);
           bitmap->imageData[i*step + j*nChannels+1] = ((s[0] * xyzTorgb[1][0] + s[1] * xyzTorgb[1][1] + s[2] * xyzTorgb[1][2])*scale);
           bitmap->imageData[i*step + j*nChannels] = ((s[0] * xyzTorgb[2][0] + s[1]* xyzTorgb[2][1] + s[2] * xyzTorgb[2][2])*scale);
        }
    }

  return 0;

}

int RGBToXYZ(IplImage* bitmap, int step, int nChannels)
{

   float rgbToxyz[3][3] = {
     { 0.490f, 0.310f, 0.200f },
     { 0.177f, 0.812f, 0.011f },
     { 0.000f, 0.010f, 0.990f }
   };

//   X = 0.412453*R + 0.35758 *G + 0.180423*B
//   Y = 0.212671*R + 0.71516 *G + 0.072169*B
//   Z = 0.019334*R + 0.119193*G + 0.950227*B
//   http://software.intel.com/sites/products/documentation/hpc/ipp/ippi/ippi_ch6/functn_RGBToXYZ.html

//   BMP* channelR = new float[width*height];
//   getChannel(source, 0, channelR);
//   BMP* channelG = new float[width*height];
//   ImageBase::GetChannel(source, 1, channelG);
//   float* channelB = new float[width*height];
//   ImageBase::GetChannel(source, 2, channelB);

   for (unsigned int i = 0; i < bitmap->height; i++)
   {
     for (unsigned int j = 0; j < bitmap->width; j++)
     {
       float s[3] = {
       (float) bitmap->imageData[i*step + j*nChannels+2] / (float) scale,
       (float) bitmap->imageData[i*step + j*nChannels+1]/ (float) scale,
       (float) bitmap->imageData[i*step + j*nChannels]/ (float) scale};

//       dest.SetPixel(s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2]
//       * rgbToxyz[0][2], s[0] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1]
//       + s[2] * rgbToxyz[1][2], s[0] * rgbToxyz[2][0] + s[1]
//       * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2], j, i);

//       cout << " | " << (int) bitmap->pixel_data[i * bitmap->width + j].red << " , "
//       << (int) bitmap->pixel_data[i * bitmap->width + j].green
//       << " , " << (int) bitmap->pixel_data[i * bitmap->width + j].blue;

//       cout << " | " << s[0] << " , " << s[1] << " , " << s[2];

//       bitmap->pixel_data[i * bitmap->width + j].red = s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2] * rgbToxyz[0][2];
//       bitmap->pixel_data[i * bitmap->width + j].blue = s[1] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1] + s[2] * rgbToxyz[1][2];
//       bitmap->pixel_data[i * bitmap->width + j].green = s[2] * rgbToxyz[2][0] + s[1] * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2];

//       cout << " | " << (int) bitmap->pixel_data[i * bitmap->width + j].red << " , "
//       << (int) bitmap->pixel_data[i * bitmap->width + j].green
//       << " , " << (int) bitmap->pixel_data[i * bitmap->width + j].blue;

       // now transform xyz to luv

       bitmap->imageData[i*step + j*nChannels+2] = scale * (s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2] * rgbToxyz[0][2]);
       bitmap->imageData[i*step + j*nChannels+1] = scale * (s[0] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1] + s[2] * rgbToxyz[1][2]);
       bitmap->imageData[i*step + j*nChannels] = scale * (s[0] * rgbToxyz[2][0] + s[1] * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2]);

   }
  }

   return 0;
}

int XYZToRGB(IplImage* bitmap, int step, int nChannels)
{
    float xyzTorgb[3][3] =
    {
      { 2.3649f, -0.8971f, -0.4678f },
      { -0.5156f, 1.4273f, 0.0883f },
      { 0.0052f, -0.0144f, 1.0092f } };

//   X = 0.412453*R + 0.35758 *G + 0.180423*B
//   Y = 0.212671*R + 0.71516 *G + 0.072169*B
//   Z = 0.019334*R + 0.119193*G + 0.950227*B
//   http://software.intel.com/sites/products/documentation/hpc/ipp/ippi/ippi_ch6/functn_RGBToXYZ.html

//   BMP* channelR = new float[width*height];
//   getChannel(source, 0, channelR);
//   BMP* channelG = new float[width*height];
//   ImageBase::GetChannel(source, 1, channelG);
//   float* channelB = new float[width*height];
//   ImageBase::GetChannel(source, 2, channelB);

   for (unsigned int i = 0; i < bitmap->height; i++)
   {
     for (unsigned int j = 0; j < bitmap->width; j++)
     {
       float s[3] = {
       (float) bitmap->imageData[i*step + j*nChannels+2] / (float) scale,
       (float) bitmap->imageData[i*step + j*nChannels+1] / (float) scale,
       (float) bitmap->imageData[i*step + j*nChannels  ] / (float) scale};

//       dest.SetPixel(s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2]
//       * rgbToxyz[0][2], s[0] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1]
//       + s[2] * rgbToxyz[1][2], s[0] * rgbToxyz[2][0] + s[1]
//       * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2], j, i);

//       cout << " | " << (int) bitmap->pixel_data[i * bitmap->width + j].red << " , "
//       << (int) bitmap->pixel_data[i * bitmap->width + j].green
//       << " , " << (int) bitmap->pixel_data[i * bitmap->width + j].blue;

//       cout << " | " << s[0] << " , " << s[1] << " , " << s[2];

//       bitmap->pixel_data[i * bitmap->width + j].red = s[0] * rgbToxyz[0][0] + s[1] * rgbToxyz[0][1] + s[2] * rgbToxyz[0][2];
//       bitmap->pixel_data[i * bitmap->width + j].blue = s[1] * rgbToxyz[1][0] + s[1] * rgbToxyz[1][1] + s[2] * rgbToxyz[1][2];
//       bitmap->pixel_data[i * bitmap->width + j].green = s[2] * rgbToxyz[2][0] + s[1] * rgbToxyz[2][1] + s[2] * rgbToxyz[2][2];

//       cout << " | " << (int) bitmap->pixel_data[i * bitmap->width + j].red << " , "
//       << (int) bitmap->pixel_data[i * bitmap->width + j].green
//       << " , " << (int) bitmap->pixel_data[i * bitmap->width + j].blue;

       // now transform xyz to luv

       bitmap->imageData[i*step + j*nChannels+2] = scale * (s[0] * xyzTorgb[0][0] + s[1] * xyzTorgb[0][1] + s[2] * xyzTorgb[0][2]);
       bitmap->imageData[i*step + j*nChannels+1] = scale * (s[0] * xyzTorgb[1][0] + s[1] * xyzTorgb[1][1] + s[2] * xyzTorgb[1][2]);
       bitmap->imageData[i*step + j*nChannels  ] = scale * (s[0] * xyzTorgb[2][0] + s[1] * xyzTorgb[2][1] + s[2] * xyzTorgb[2][2]);

     }
   }

       return 0;
}

void convertion(IplImage* bitmap, char* conv, int step, int nChannels)
{
    short * luv_data = (short*) malloc(sizeof(short)*step*sizeof(short)*(bitmap->height));
    uchar *data = ( uchar* )bitmap->imageData;
    if (!strcmp(conv, "rgb2luv"))
    {
        int luv[3];
        for (unsigned int i = 0; i < bitmap->height; i++)
        {
            for (unsigned int j = 0; j < bitmap->width; j++)
            {
                rgb2luv(data[i*step + j*nChannels], // nao funciona na formula antiga pq o canal 01 tava trocado com o 03
                        data[i*step + j*nChannels+1],
                        data[i*step + j*nChannels+2],
                        luv);

                data[i*step + j*nChannels] = luv[0];   // -> L
                data[i*step + j*nChannels+1] = luv[1]; //ordem trocada -> v
                data[i*step + j*nChannels+2] = luv[2]; //ordem trocada -> u


            }
        }
    } else if (!strcmp(conv, "luv2rgb"))
    {
        for (unsigned int i = 0; i < bitmap->height; i++)
        {
            for (unsigned int j = 0; j < bitmap->width; j++)
            {
                int vec[3];
                int rgb[3];
                vec[0] = data[i*step + j*nChannels+2];
                vec[1] = data[i*step + j*nChannels+1];
                vec[2] = data[i*step + j*nChannels  ];
                luv2rgb(vec , rgb);

                data[i*step + j*nChannels] = rgb[0];
                data[i*step + j*nChannels+1] = rgb[1];
                data[i*step + j*nChannels+2] = rgb[2];


            }
        }
    } else {printf("Escolha um tipo de conversao");}
}


void rgb2luv(int R, int G, int B, int *luv) {
    //http://www.brucelindbloom.com

    double r, g, b, X, Y, Z, yr;
    double L;
    double eps = 216./24389.;
    double k = 24389./27.;

    double Xr = 0.964221;  // reference white D50
    double Yr = 1.0;
    double Zr = 0.825211;


//			double ur=4.*.964221/(.964221+15.+3.*.825211);
//			// vr=9*Yr/(Xr+15*Yr+3*Zr)
//			double vr=9./(.964221+15.+3.*.825211);

//			double dmax=(double)255;
//			double dmaxdiv2=(double)(255./2.+1);

    // multiplication is faster than division
    double div12p92=1./12.92;
    double div1p055=1./1.055;
//			double aThird=1./3.;
    double divmax=1./(double)255;

    // RGB to XYZ

    r = R*divmax; //R 0..1
    g = G*divmax; //G 0..1
    b = B*divmax; //B 0..1

    // assuming sRGB (D65)
    if (r <= 0.04045)
        r *=div12p92;
    else
        r = (float) pow((r+0.055)*div1p055,2.4);

    if (g <= 0.04045)
        g *=div12p92;
    else
        g = (float) pow((g+0.055)*div1p055,2.4);

    if (b <= 0.04045)
        b *=div12p92;
    else
        b = (float) pow((b+0.055)*div1p055,2.4);


//     X =  0.412424f * r + 0.357579f * g + 0.180464f  * b;
//     Y =  0.212656f * r + 0.715158f * g + 0.0721856f * b;
//     Z = 0.0193324f * r + 0.119193f * g + 0.950444f  * b;
    /*
     // chromatic adaptation transform from D65 to D50
      X =  1.047835f * X_ + 0.022897f * Y_ - 0.050147f * Z_;
      Y =  0.029556f * X_ + 0.990481f * Y_ - 0.017056f * Z_;
      Z = -0.009238f * X_ + 0.015050f * Y_ + 0.752034f * Z_;
     */

    X =  0.436052025f*r	+ 0.385081593f*g + 0.143087414f *b;
    Y =  0.222491598f*r	+ 0.71688606f *g + 0.060621486f *b;
    Z =  0.013929122f*r	+ 0.097097002f*g + 0.71418547f  *b;

    // XYZ to Luv

    double u, v, u_, v_, ur_, vr_;

    u_ = 4*X / (X + 15*Y + 3*Z);
    v_ = 9*Y / (X + 15*Y + 3*Z);

    ur_ = 4*Xr / (Xr + 15*Yr + 3*Zr);
    vr_ = 9*Yr / (Xr + 15*Yr + 3*Zr);

    yr = Y/Yr;

    if ( yr > eps )
        L =  (float) (116*pow(yr, 1/3.) - 16);
    else
        L = k * yr;

    u = 13*L*(u_ -ur_);
    v = 13*L*(v_ -vr_);

    luv[0] = (int) (2.55*L + .5);
//    luv[1] = (int) (u + .5) + 84;
//    luv[2] = (int) (v + .5) + 124;
//    luv[0] = (int) L + .5;
    luv[1] = (int) (u + .5) + 84;
    luv[2] = (int) (v + .5) + 124;
}

void luv2rgb(int* vec, int* rgb) {

    double r,g,b;
    double X,Y,Z;
    double L,u,v;
    double a_,b_,d_;
    double tL;

    // all matrices, formulas etc. used here are from www.brucelindbloom.com
    // const double eps=216./24389.;
    // const double k=24389./27.;
    double keps=8; // k*eps

    // XYZ reference white D50:
    double Xr=0.964221;
    // not Yr==1 for D50...
    // const double Yr=1.0;
    double Zr=0.825211;
    double f1=100./255.;

    double u0t13=13.*((4.*Xr)/(Xr+15.+3.*Zr));	// Yr==1
    double v0t13=13.*((9.)/(Xr+15.+3.*Zr));		// Yr==1
    double div3=1./3.;

    // note multiplication is faster than division:
    double div116=1./116.;
    double divk=27./24389.;
    double div2p4=1./2.4;
    double dmax=255;
    //double sub=maxChannelValue()/2+1;

    L=tL=vec[2]*f1;
    u=vec[0];
    v=vec[1];

    // LAB to XYZ:
    if(tL>keps)
        Y=pow((L+16.)*div116,3.);
    else
        Y=L*divk;

    a_=div3*((52.*L/(u+L*u0t13))-1.);
    b_=-5.*Y;
    d_=Y*((39.*L/(v+L*v0t13))-5.);

    X=(d_-b_)/(a_+div3);
    Z=X*a_+b_;

    // chromatic adaption to d65 ref white:

    // X_d65=.955556*X+(-.028302)*Y+.012305*Z;
    // Y_d65=(-.023049)*X+1.009944*Y+(-.020494)*Z;
    // Z_d65=.063197*X+.021018*Y+1.330084*Z;

    // convert to rgb
    // r=3.24071*X_d65+(-1.53726)*Y_d65+(-0.498571)*Z_d65;
    // g=(-0.969258)*X_d65+1.87599*Y_d65+(0.0415557)*Z_d65;
    // b=(0.0556352)*X_d65+(-0.203996)*Y_d65+1.05707*Z_d65;

    // both transforms at once:

    // note this is propably not right:
    // we may better use the inverse of RGB2LAB()
    // because RGB2LAB() followed by LAB2RGB() is not a nop, currently
    // r=3.100603999*X+(-1.654744053)*Y+(-0.591759767)*Z;
    // g=(-0.966793795)*X+1.922950202*Y+(0.004899313)*Z;
    // b=(0.124668106)*X+(-0.185381626)*Y+1.410857179*Z;

    // inverse:
    r=3.134051341*X+(-1.617027711)*Y+(-0.49065221)*Z;
    g=(-0.97876273)*X+1.916142228*Y+(0.033449629)*Z;
    b=(0.071942577 )*X+(-0.22897118)*Y+1.405218305 *Z;

    // there are different RGB color spaces
    // for most of these, we have to compute
    // r=r^(1/G) for the transform, G being the gamma
    // of the RGB colorspace.

    // anyway, we don't know which RGB color space to use (why not?)
    // thus, we assume sRGB which is widely used in the win domain
    // for sRGB, gamma is approx. 2.2 and using the formula above is called
    // 'simplified sRGB', because the real formula is different:
    if(r<=.0031308)
        r*=12.92;
    else
        r=pow(1.055*r,div2p4)-0.055;

    if(g<=.0031308)
        g*=12.92;
    else
        g=pow(1.055*g,div2p4)-0.055;

    if(b<=.0031308)
        b*=12.92;
    else
        b=pow(1.055*b,div2p4)-0.055;

    // transform r,g,b from [0,1] to [0,maxChannelValue()]
    // take care: r, g, b may exceed [0,1] now, this may be true
    // if we do RGB2LUV() followed by LUV2RGB() due to cutting dezimals
    // when converting from double to type T
    // thus:
    rgb[0] =(int)(max(0.,min(255.,dmax*r))+.5);
    rgb[1] =(int)(max(0.,min(255.,dmax*g))+.5);
    rgb[2] =(int)(max(0.,min(255.,dmax*b))+.5);

}

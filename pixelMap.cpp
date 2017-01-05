/*
  Name:        PixelMap.cpp
  Copyright:   Version 0.1
  Author:      Rodrigo Luis de Souza da Silva
  Date:        04/07/2009
  Description: Simple class to load BMP images and modify its data.
               You can create an image manually and convert it to a
               pixel matrix and vice-versa.
*/

#include "pixelMap.h"
#include <fstream>
#include <math.h>
#include <iostream>
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/glut.h>
#include <ctime>

#define FALSE 0
#define TRUE 1


//////////////////////////////////////////////////////////////////
PixelMap::PixelMap(const char *fname)
	: width(0), height(0), data(0)
{
   this->grayScale = false;
	this->Read(fname);
}

//////////////////////////////////////////////////////////////////
PixelMap::PixelMap()
	: width(0), height(0), data(0)
{
   this->grayScale = false;
}

//////////////////////////////////////////////////////////////////
PixelMap::~PixelMap()
{
	if( data )
		delete[] data;
}

//////////////////////////////////////////////////////////////////
void PixelMap::Read(const char *fname)
{
	using namespace std;

	unsigned short planes;	// number of planes in image (must be 1)
	unsigned short bpp;			// number of bits per pixel (must be 24)

	ifstream fin(fname, ios::in | ios::binary);
	if( !fin )
	{
		cerr << "File not found " << fname << '\n';
		exit(1);
	}

	fin.seekg(18, ios::cur);

	fin.read((Byte *)&width, sizeof(int));
	fin.read((Byte *)&height, sizeof(int));
	//cout << "width: " << width << " height: " << height << '\n';

	fin.read((Byte *)&planes, sizeof(unsigned short));
	if( planes != 1 )
	{
		cout << "Planes from " << fname << " is not 1: " << planes << "\n";
		exit(1);
	}

	fin.read((Byte *)&bpp, sizeof(unsigned short));
	if( bpp != 24 )
	{
		cout << "Bpp from " << fname << " is not 24: " << bpp << "\n";
		exit(1);
	}

	fin.seekg(24, ios::cur);

	int size(width * height * 3); // size of the image in bytes (3 is to RGB component).
	data = new uByte[size];


    fourierRE = (float*) malloc(width*height*3*4);
    fourierIG = (float*) malloc(height*width*3*4);

	fin.read((Byte *)data, size);

	uByte tmp;	// temporary color storage for bgr-rgb conversion.

	for( int i(0); i <  size; i += 3 ){
		tmp = data[i];
		data[i] = data[i+2];
		data[i+2] = tmp;
	}

	inverse = 0;

}

//////////////////////////////////////////////////////////////////
void PixelMap::CreateImage(int w, int h, uByte *d)
{
   this->width  = w;
   this->height = h;
   this->data   = d;
}

//////////////////////////////////////////////////////////////////
void PixelMap::ViewImage(){
   glDrawPixels(this->GetWidth(), this->GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, (uByte*) this->GetData());
}

//////////////////////////////////////////////////////////////////
void PixelMap::ViewOldImage(){
   glDrawPixels(this->GetWidth(), this->GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, (uByte*) this->GetOldData());
}

//////////////////////////////////////////////////////////////////
uByte PixelMap::pixel_elem(int x, int y, int elem){

	int pos = (y*width+x) * 3 + elem;
	return data[pos];
}

//////////////////////////////////////////////////////////////////
uByte *PixelMap::pixel_pos(int x, int y)
{
	int pos = (y * width + x) * 3;
	return &data[pos];
}

//////////////////////////////////////////////////////////////////
int PixelMap::GetWidth()
{
   return this->width;
}

//////////////////////////////////////////////////////////////////
void PixelMap::SetWidth(int w)
{
   this->width = w;
}

//////////////////////////////////////////////////////////////////
int PixelMap::GetHeight()
{
   return this->height;
}

//////////////////////////////////////////////////////////////////
void PixelMap::SetHeight(int h)
{
   this->height = h;
}

//////////////////////////////////////////////////////////////////
uByte* PixelMap::GetData()
{
   return this->data;
}

uByte* PixelMap::GetOldData()
{
   return this->oldData;
}

//////////////////////////////////////////////////////////////////
void PixelMap::SetData(uByte *d)
{
   this->data = d;
}

//////////////////////////////////////////////////////////////////
void PixelMap::ConvertPixelsToData(pixel **m)
{
   int k = 0;
   for(int j = 0; j < this->GetHeight(); j++)
      for(int i = 0; i < this->GetWidth(); i++)
      {
         // If grayscale, set all components to value
         if(this->grayScale)
            m[i][j].R = m[i][j].G = m[i][j].B = m[i][j].value;
         data[k]   = m[i][j].R;
         data[k+1] = m[i][j].G;
         data[k+2] = m[i][j].B;
         k+=3;
      }
}

//////////////////////////////////////////////////////////////////
void PixelMap::ConvertDataToPixels(pixel **m)
{
   uByte *dataAux = this->data;
   for(int j = 0; j < this->GetHeight(); j++)
      for(int i = 0; i < this->GetWidth(); i++)
      {
         if(this->grayScale)
            m[i][j].value = *dataAux; // All components are the same in grayscale mode
         m[i][j].R = *dataAux++;
         m[i][j].G = *dataAux++;
         m[i][j].B = *dataAux++;
      }
}

//////////////////////////////////////////////////////////////////
void PixelMap::ConvertToGrayScale()
{
   this->grayScale = true;
   double value;
	int k = 0;
   for(int j = 0; j < this->GetHeight(); j++)
      for(int i = 0; i < this->GetWidth(); i++){
         value =  0.3 * (double) data[k] + 0.59 * (double) data[k+1] + 0.11 * (double) data[k+2];
         data[k] = data[k+1] = data[k+2] = (uByte) value;
         k+=3;
      }
}

//////////////////////////////////////////////////////////////////
void PixelMap::filtroMedia(int n)
{
    double somatorio = 0.0;
    int iteracoes = n/2;

    uByte *newData;
    newData = new uByte[width*height*3];


    for (int m = 0; m<width*height*3; m++)
        newData[m] = (uByte) 255.0;

    for(int i = 0; i < this->GetWidth(); i++){
        for(int j = 0; j < this->GetHeight(); j++){
            for (int c = i-iteracoes; c<=i+iteracoes; c++){
                for (int k = j-iteracoes; k<=j+iteracoes; k++){
                    if (c>=0 && c<this->GetWidth() && k>=0 && k<this->GetHeight())
                         somatorio += (double) pixel_elem(c,k,0);
//                if (i==1 && j==1)
//                    printf("\nA[%d,%d] += B[%d,%d]",i,j,c,k);
                }
            }

        somatorio /= n*n;
        int pos = (j*width+i)*3;
        newData[pos] = newData[pos+1] = newData[pos+2] = (uByte) somatorio;
        somatorio  = 0.0;
        }
    }
    delete[] data;
    data = newData;
}


///////////////////////////////////////////////////////////////////
void PixelMap::filtroMediana(int n)
{
    int iteracoes = n/2;
    int valorMedio = (n*n)/2 + 1;
    //printf("Valor medio = %d",valorMedio);
    double vetorMediana[n*n];

    uByte *newData;
    newData = new uByte[width*height*3];

    for (int m = 0; m<width*height*3; m++)
        newData[m] = (uByte) 0.0;

    int cont_pos_vetor = 0;
    int cont_sem_valores = 0; //bordas

    for(int i = 0; i < this->GetWidth(); i++){
        for(int j = 0; j < this->GetHeight(); j++){
            for (int c = i-iteracoes; c<=i+iteracoes; c++){
                for (int k = j-iteracoes; k<=j+iteracoes; k++){
                    if (c>=0 && c<this->GetWidth() && k>=0 && k<this->GetHeight())
                         vetorMediana[cont_pos_vetor++] = (double) pixel_elem(c,k,0);
                     else{
                         vetorMediana[cont_pos_vetor++] = -999;
                         cont_sem_valores++;
                     }
                }
            }

            //Ordenação ridícula pelo método bolha
            double aux;
            for(int p=0; p<(n*n); p++)
            for(int o=0; o<(n*n)-1; o++)
            if(vetorMediana[o]>vetorMediana[o+1]){
                aux = vetorMediana[o];
                vetorMediana[o] = vetorMediana[o+1];
                vetorMediana[o+1] = aux;
            }


            //Tratando problema de bordas
            //Modo 1:
            vetorMediana[valorMedio] = vetorMediana[valorMedio + (int) (cont_sem_valores/2)];


            //Modo2:
            /*
            int prox = valorMedio;
            do{
                vetorMediana[valorMedio] = vetorMediana[prox++];
            }while (vetorMediana[valorMedio] == -999);*/

            cont_pos_vetor = 0;
            cont_sem_valores = 0;
            int pos = (j*width+i)*3;
            newData[pos] = newData[pos+1] = newData[pos+2] = (uByte) vetorMediana[valorMedio];
        }
    }
    delete[] data;
    data = newData;
}

int randomNumber(int hi){  //the correct random number generator for [0,hi]
       // scale in range [0,1)
       const float scale = rand()/float(RAND_MAX);
       // return range [0,hi]
       return int(scale*hi + 0.5); // implicit cast and truncation in return
}

void PixelMap::ruidoSalPimenta(double perc){

    srand(time(NULL)); //Inicializa a cadeia randomica

    ///Calculamos a quantidade de pixels que iremos modificar por linha primeiro.
    int pmax = (int) (perc*this->GetWidth());

    for (int j = 0; j<this->GetHeight(); j++)//para cada linha..
       for (int i = 0; i<pmax; i++){//..alteramos até pmax pixels
           ///Pega um dos pixels da linha de forma aleatória.
           int x = randomNumber(this->GetWidth());
           int pos = (j * width + x) * 3;

           ///Vamos atribuir 0.0 ou 255.0 para cada pixel de forma aleatória.
           data[pos] = data[pos+1] = data[pos+2] = (uByte) randomNumber(1)*255.0;
       }
}

void PixelMap::ruidoPersonalzado(double perc, double increase){

    srand(time(NULL));//Inicializa a cadeia randomica

    int pmax = (int) (perc*this->GetWidth());
    int size = this->GetWidth()*this->GetHeight()*3;

    for (int j = 0; j<this->GetHeight(); j++)//para cada linha..
       for (int i = 0; i<pmax; i++){//..alteramos até pmax pixels
           ///Pega um dos pixels da linha de forma aleatória.
           int x = randomNumber(this->GetWidth());

           int pos = (j * width + x) * 3;

           ///conterá a nova cor do pixel após operação
           double value;

           ///Vamos somar ou subtrair a cor do pixel[pos] de forma aleatória
           if (randomNumber(1)) value = (double) data[pos] + increase;
           else  value = (double) data[pos] - increase;

           if (value<0.0) value = 0.0;
           else if (value>255.0) value = 255.0;

          data[pos] = data[pos+1] = data[pos+2] = (uByte) value;
       }
}

float max = 9.0;

///////////////////////////////////////////////////////////////////
///Aplica transformada de fourier
/*
void PixelMap::fourierTransform(){
    uByte *newData;
    newData = new uByte[width*height*3];

    ///Alocação de um buffer temporario usado para armazenar imagem enquanto processada
    for (int m = 0; m<width*height*3; m++)
        newData[m] = (uByte) 255.0;

    int N = this->GetWidth();
    int M = this->GetHeight();

    ///Valor do @spectro @maximo da imagem.
    ///Será usado para visualizar a imagem no spaço de fourier
    float maxSpectrum = -99999999999.9;

    for(int u = 0; u < N; u++){
        for(int v = 0; v < M; v++){
          ///Indice da imagem final f(u,v)
          int index = (v * N + u)*3;
          ///Zerando o buffer que contera os valores da parte real e
          ///imaginária do dominio de fourier
          fourierRE[index] = fourierIG[index] = 0.0;

        ///Itera pela imagem de entrada
          for(int x = 0; x < N; x++){
            for(int y = 0; y < M; y++){
            ///Posição do pixel na imagem de entrada
            int pos = (y * N + x)*3;

            ///Desloca o spectro pro centro da imagem
            float u2 = ((float) u + (float) N / 2.0);
            float v2 = ((float) v + (float) M / 2.0);

            ///Transformada cosseno
            float a = (2.0 * PI) * (((float)v2 * (float)y / (float) M) + ((float)u2 * (float)x / (float) N)) ;

            ///Transformada Normal
            float cosA = cos(-a);
            float senA = sin(-a);
            ///Calcula a parte real e a parte imaginaria da imagem resultante
            fourierRE[index] += ((float) data[pos]) * cosA;
            fourierIG[index] += ((float) data[pos]) * senA;
          }
        }
          ///Divide a parte real e imaginaria da imagem transformada,
          ///segundo a transformada discreta de fourier
          fourierRE[index] /= N*M;
          fourierIG[index] /= N*M;

          ///Pega a maior amplitudo da imagem resultante
          float tempAmplitude = sqrt(fourierRE[index]*fourierRE[index] + fourierIG[index]*fourierIG[index]);
          if(tempAmplitude > maxSpectrum) maxSpectrum = tempAmplitude;

        }
        printf("U Espaco-Fourier: %d\n", u);
    }

    for(int i= 0; i < N; i++)
        for(int j= 0; j < M; j++){

            int index = (j * N + i)*3;

            ///Transformação necessária para colocar os valores dentro do intervalo de [0, 255.0]
            float transform = sqrt(fourierRE[index]*fourierRE[index] + fourierIG[index]*fourierIG[index]);
            float c = 255.0/log10f(1.0 + fabs(maxSpectrum));
            transform = c*log10f(1.0+fabs(transform));

             ///Teste
            if (transform < 0.0) transform = 0.0;
            else if (transform > 255.0) transform = 255.0;

            newData[index] = newData[index+1] = newData[index+2] = (uByte)  transform;
        }

    max = maxSpectrum; //Só é utilizado caso queira testar os filtros (as novas imagens geradas no espectro de fourier)

    delete[] data;
    data = newData;
}


/**--------------------------------------**/
///  TRANSFORMADA INVERSA DE FOURIER     ///
/**--------------------------------------**//*
void PixelMap::ifourierTransform(){
    uByte *newData;
    newData = new uByte[width*height*3];

    for (int m = 0; m<width*height*3; m++)
        newData[m] = (uByte) 255.0;

    int N = this->GetWidth();
    int M = this->GetHeight();

    int size(M * N * 3);
    float *tempRE = (float*) malloc(size*4);
    float *tempIG = (float*) malloc(size*4);

    for(int u = 0; u < N; u++){
        for(int v = 0; v < M; v++){

          int index = (v * N + u)*3;
          tempRE[index] = tempIG[index] = 0.0;

          for(int x = 0; x < N; x++){
            for(int y = 0; y < M; y++){

            int pos = (y * N + x)*3;

            float x2 = ((float) x - (float) N / 2.0);
            float y2 = ((float) y - (float) M / 2.0);

            float a = (2.0 * PI) * (((float)y2 * (float)v / (float) M) + ((float)x2 * (float)u / (float) N)) ;

            ///Transformada cosseno
            float cosA = cos(a);
            float senA = sin(a);

            tempRE[index] += fourierRE[pos]*cosA - fourierIG[pos]*senA;
            tempIG[index] += fourierRE[pos]*senA + fourierIG[pos]*cosA;
          }
        }

        }
        printf("U Inversa: %d\n", u);
    }


    for(int i= 0; i < N; i++)
        for(int j= 0; j < M; j++){

            int index = (j * N + i)*3;
            float transform = sqrt(tempRE[index]*tempRE[index] + tempIG[index]*tempIG[index]);
            if (transform < 0.0) transform = 0.0;
            else if (transform > 255.0) transform = 255.0;

            newData[index] = newData[index+1] = newData[index+2] = (uByte)  transform;
        }

    delete[] data;
    data = newData;
}
*/


///Função utilizada para exibir a imagem de frequencia após alguma modificação, como os filtros.
void PixelMap::reescreve(){
    int N = this->GetWidth();
    int M = this->GetHeight();

    float maxSpectrum = max;

    for(int i= 0; i < N; i++)
    for(int j= 0; j < M; j++){

        int index = (j * N + i)*3;

        double transform = sqrt(fourierRE[index]*fourierRE[index] + fourierIG[index]*fourierIG[index]);

        float c = 255.0/log10f(1.0 + fabs(maxSpectrum));
        transform = c*log10f(1.0+fabs(transform));

        ///Teste
        if (transform < 0.0) transform = 0.0;
        else if (transform > 255.0) transform = 255.0;

        data[index] = data[index+1] = data[index+2] = (uByte)  transform;
    }

}

int PixelMap::FFT2D(int inverse){

    uByte *newData;
    newData = new uByte[width*height*3];

    pixel **pixels;
    COMPLEX **c;

    int altura = this->GetHeight();
    int largura = this->GetWidth();

    pixels = (pixel**)malloc(altura*sizeof(pixel*));
    for (int i = 0; i < altura; i++)
          pixels[i] = (pixel*)malloc(largura*sizeof(pixel));

    this->ConvertDataToPixels(pixels);

    c = (COMPLEX**)malloc(altura*sizeof(COMPLEX*));
    for (int i = 0; i < altura; i++)
        c[i] = (COMPLEX*)malloc(largura*sizeof(COMPLEX));

    if(inverse == 1){//FFT
        for (int i = 0; i < altura; i++)
            for (int j = 0; j < largura; j++) {
                c[i][j].real = pow(-1, i+j)*pixels[i][j].value;
                c[i][j].imag = 0.0;
            }
    }

    if(inverse == 0){//iFFT
        for (int i = 0; i < largura; i++)
            for (int j = 0; j < altura; j++) {
                int ind = (j * this->GetWidth() + i)*3;
                c[i][j].real = fourierRE[ind];
                c[i][j].imag = fourierIG[ind];
            }
    }

    int nx = altura;
    int ny = largura;
    int dir = inverse;
    int i, j;
    int m, twopm;
    double *real,*imag;

    /* Transform the rows */
    real = (double *)malloc(nx * sizeof(double));
    imag = (double *)malloc(nx * sizeof(double));
    if (real == NULL || imag == NULL)
        return(FALSE);

    if (!Powerof2(nx,&m,&twopm) || twopm != nx)
        return(FALSE);

    for (j = 0 ; j < ny; j++){
        for (i = 0; i < nx; i++){
            real[i] = c[i][j].real;
            imag[i] = c[i][j].imag;
        }

        FFT(dir, m, real, imag);
        for (i = 0; i < nx; i++) {
            c[i][j].real = real[i];
            c[i][j].imag = imag[i];
        }
    }
    free(real);
    free(imag);

    /* Transform the columns */
    real = (double *)malloc(ny * sizeof(double));
    imag = (double *)malloc(ny * sizeof(double));

    if (real == NULL || imag == NULL)
        return(FALSE);

    if (!Powerof2(ny,&m,&twopm) || twopm != ny)
        return(FALSE);

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            real[j] = c[i][j].real;
            imag[j] = c[i][j].imag;
        }
        FFT(dir,m,real,imag);
        for (j = 0; j < ny; j++) {
            c[i][j].real = real[j];
            c[i][j].imag = imag[j];
        }
    }
    free(real);
    free(imag);

    float maxSpectrum = 0.0;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++){
            int index = (j * this->GetWidth() + i)*3;

            fourierRE[index] = c[i][j].real;
            fourierIG[index] = c[i][j].imag;

            float amplitude = sqrt(fourierIG[index]*fourierIG[index] + fourierRE[index]*fourierRE[index]);
            if(amplitude > maxSpectrum) maxSpectrum = amplitude;
        }

      for(int i= 0; i < largura; i++)
        for(int j= 0; j < altura; j++){
            int index = (j * largura + i)*3;
            ///Transformação necessária para colocar os valores dentro do intervalo de [0, 255.0]
            double transform = sqrt(fourierRE[index]*fourierRE[index] + fourierIG[index]*fourierIG[index]);
//            double transform = (fourierRE[index]*fourierRE[index] + fourierIG[index]*fourierIG[index]);
//            double transform = 1/(atan(fourierIG[index]/fourierRE[index]));

            if(inverse == 1){
                float c = 255.0/log10f(1.0 + fabs(maxSpectrum));
                transform = c*log10f(1.0+fabs(transform));
            }

            ///Teste
            if (transform < 0.0) transform = 0.0;
            else if (transform > 255.0) transform = 255.0;

            newData[index] = newData[index+1] = newData[index+2] = (uByte)  transform;
        }

    max = maxSpectrum; //Só é utilizado caso queira testar os filtros (as novas imagens geradas no espectro de fourier)

    delete[] data;
    data = newData;
}

int PixelMap::FFT(int dir, int m, double *x, double *y) {
    long nn, i, i1, j, k, i2, l, l1, l2;
    double c1, c2, tx, ty, t1, t2, u1, u2, z;

    /* Calculate the number of points */
    nn = 1;
    for (i = 0; i < m; i++)
        nn *= 2;

    /* Do the bit reversal */
    i2 = nn >> 1;
    j = 0;
    for (i = 0; i < nn - 1; i++){
        if (i < j) {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    /* Compute the FFT */
    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;
    for (l = 0; l < m; l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j = 0; j < l1; j++) {
            for (i=j;i<nn;i+=l2) {
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }

    for (i = 0; i < nn; i++){
        x[i] /= (double)sqrt(nn);
        y[i] /= (double)sqrt(nn);
    }

    return(TRUE);
}

int PixelMap::Powerof2(int nx, int *m, int *twopm) {
    int pwr;
	*m = 0;
	for (pwr = 1; pwr < nx; pwr = pwr*2) {
		*m = *m + 1;
	}
	*twopm = pwr;
	return(TRUE);
}


void PixelMap::removeRuido(int x, int y, int raio){
    int M = this->GetWidth();
    int N = this->GetHeight();

    int newy = M - y;
    int newx = N - x;

    for(int u = 0; u < M; u++)
        for(int v = 0; v < N; v++){

            float u2 = (float) u - x;
            float v2 = (float) v - y;

            float u3 = (float) u - newx;
            float v3 = (float) v - newy;

            if (u2*u2+v2*v2<=(raio*raio) || u3*u3+v3*v3<=(raio*raio) ){
                int index = (v * M + u)*3;
                fourierRE[index] = 0.0;
                fourierIG[index] = 0.0;
            }
        }

    reescreve();
}

void PixelMap::adicionaRuido(int x, int y, int raio){
    int M = this->GetWidth();
    int N = this->GetHeight();
    //printf ("Adicionando ruido.\n");

   int newy = M - y;
   int newx = N - x;

    for(int u = 0; u < M; u++)
        for(int v = 0; v < N; v++){

            float u2 = (float) u - x;
            float v2 = (float) v - y;

            float u3 = (float) u - newx;
            float v3 = (float) v - newy;

            double dist1 = u2*u2+v2*v2;
            double dist2 = u3*u3+v3*v3;

            if (dist1<=(raio*raio) || dist2<=(raio*raio)){
                int index = (v * M + u)*3;

                double value = 255.0;
                double c = 255.0/raio;
                double dist;
                float U,V;


                if (dist1<=(raio*raio)) {
                    dist = dist1;
                    U = u2;
                    V = v2;
                }
                else{
                    dist = dist2;
                    U = u3;
                    V = v3;
                }

                if (V==0 && U==0) value = 900.0;
                else if (V==0 && U<0) value = (int)((1/-U)*value)%255 + value - sqrt(dist)*c;
                else if (V==0) value = (int)((1/U)*value)%255 + value - sqrt(dist)*c;
                else if (U==0 && V<0) value = (int)((1/-V)*value)%255 + value - sqrt(dist)*c;
                else if (U==0) value = (int)((1/V)*value)%255 + value - sqrt(dist)*c;
              /*  else if (U>0 && V>0) value = (int)((1/U*2+1/V*2)*value)%255 + value - sqrt(dist)*c;
                else if (U>0 && V<0) value = (int)((1/U*2+1/-V*2)*value)%255 + value - sqrt(dist)*c;
                else if (U<0 && V>0) value = (int)((1/-U*2+1/V*2)*value)%255 + value - sqrt(dist)*c;
                else if (U<0 && V<0) value = (int)((1/-U*2+1/-V*2)*value)%255 + value - sqrt(dist)*c;
*/
               else value = (int)((1/sqrt(U*U+V*V))*value)%255 + value - sqrt(dist)*c;
            /*  Desenha só a circunferencia em dégradé
                if (u == x || v == y) value = ((value - sqrt(dist1)*c)+value)/2.0;
                else value = value - sqrt(dist1)*c;*/



                if (fourierRE[index]>0) fourierRE[index]= value;
                else fourierRE[index] = -value;

                if (fourierIG[index]>0) fourierIG[index]= value;
                else fourierIG[index]= -value;

//                fourierRE[index] = 255.0;
//                fourierIG[index] = 255.0;

            }
        }
    reescreve();
}

///FILTROS Passa alta, Passa Baixa: nao estao sendo usados, mas a implementação é como se segue
void PixelMap::passaAlta(int R){
    int M = this->GetWidth();
    int N = this->GetHeight();

    for(int u = 0; u < M; u++)
        for(int v = 0; v < N; v++){

            float u2 = (float) u - (float)M/2.0;
            float v2 = (float) v - (float) N/2.0;

            if (u2*u2+v2*v2<(R*R)){
                int index = (v * M + u)*3;
                fourierRE[index] = 0.0;
                fourierIG[index] = 0.0;
            }
        }
    reescreve();
}

void PixelMap::passaBaixa(int R){
    int M = this->GetWidth();
    int N = this->GetHeight();

    for(int u = 0; u < M; u++)
        for(int v = 0; v < N; v++){

            float u2 = (float) u - (float)M/2.0;
            float v2 = (float) v - (float)N/2.0;

            if (u2*u2+v2*v2>=(R*R)){
                int index = (v * M + u)*3;
                    fourierRE[index] = 0.0;
                    fourierIG[index] = 0.0;
            }
        }
    reescreve();
}




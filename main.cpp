#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>

#include "luv.h"

int imgLoad;
int window1, window2;
IplImage* img = 0;
IplImage* imgAux = 0;
int step = 0;
int nChannels = 0;
int mouseX, mouseY;
bool* vflags; // Matriz de flags do tamanho da imagem;
unsigned char R, G, B;
bool luv_activated = false;

CvHaarClassifierCascade* cascade_f = 0; // para detecção de face
CvHaarClassifierCascade* cascade_e = 0; // para detecção de olhos
CvMemStorage* storage;
int desl_x, desl_y, init_x, init_y; // posição da seleção
int buttonGlobal, stateGlobal; // para acessar os estados do mouse em outras funções
IplImage* tempData = 0; // imagem temp
IplImage* img_gt = 0; // ground truth
IplImage* img_rgb = 0; // imagem ORIGINAL guardada, e não a renderizada em RGB
IplImage* img_luv = 0; // imagem renderizada na janela LUV
CvRect ROI; // Região de Interesse que é o retangulo marcado
int maxL, minL, maxU, minU, maxV, minV; //limites dos canais no retangulo selecionado
CvRect* eye;  // retangulo do olho
int janID = 0; // id da janela do LUV para redisplay

// Inverte imagem na vertical para mostrar corretamente na tela e a ordem dos canais
void revert( IplImage* src, IplImage* dst )
{
    for( int i = 0; i < src->height; i++ )
    {
        for( int j = 0; j < src->width; j++ )
        {
            src->imageData[ i * step + j * nChannels +
                            2 ] = dst->imageData[ ( dst->height - i ) * step + j * nChannels ];
            src->imageData[ i * step + j * nChannels +
                            1 ] = dst->imageData[ ( dst->height - i ) * step + j * nChannels + 1 ];
            src->imageData[ i * step + j *
                            nChannels ] = dst->imageData[ ( dst->height - i ) * step + j * nChannels + 2 ];
        }
    }
}


// Retorna posição da imagem para salvar
void revertBack()
{
    for( int i = 0; i < img->height; i++ )
    {
        for( int j = 0; j < img->width; j++ )
        {
            imgAux->imageData[ i * step + j * nChannels +
                               2 ] = img->imageData[ ( imgAux->height - i ) * step + j * nChannels ];
            imgAux->imageData[ i * step + j * nChannels +
                               1 ] = img->imageData[ ( imgAux->height - i ) * step + j * nChannels + 1 ];
            imgAux->imageData[ i * step + j *
                               nChannels ] = img->imageData[ ( imgAux->height - i ) * step + j * nChannels + 2 ];
        }
    }
}


// Matriz de flags para controle do crescimento de região
void fillFlag()
{
    vflags = ( bool* )malloc( sizeof( bool ) * step * ( img->width ) + ( img->height ) );
    int total, total2;
    for( int i = 0; i < img->height; i++ )
    {
        for( int j = 0; j < img->width; j++ )
        {
            vflags[ i * step + j ] = false;
            //total = i*step + j;
            //total2 = step*(img->width) + img->height;
        }
    }
    printf( "Total vflags: %d\n", total );
    printf( "Total2 vflags: %d\n", total2 );
}


// Standart display
static void display( void )
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //Desenha Imagem
    glDrawPixels( img->width, img->height, GL_RGB, GL_UNSIGNED_BYTE, ( uchar* ) img->imageData );


    glutSwapBuffers();
    glutReshapeWindow( img->width, img->height );
}


// Ground truth display
static void display_gt( void )
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //Desenha Imagem
    glDrawPixels( img_gt->width, img_gt->height, GL_RGB, GL_UNSIGNED_BYTE, ( uchar* ) img_gt->imageData );


    glutSwapBuffers();
    glutReshapeWindow( img_gt->width, img_gt->height );
}


// LUV display
static void display_luv( void )
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //Desenha Imagem
    glDrawPixels( img_luv->width, img_luv->height, GL_RGB, GL_UNSIGNED_BYTE, ( uchar* ) img_luv->imageData );


    glutSwapBuffers();
    glutReshapeWindow( img_luv->width, img_luv->height );
}


// Region growing
void evaluatePixel( int x, int y )
{
    if( !( vflags[ ( y * step ) + x ] ) )
    {
        R = img->imageData[ ( y * step ) + x * nChannels ];
        G = img->imageData[ ( y * step ) + x * nChannels + 1 ];
        B = img->imageData[ ( y * step ) + x * nChannels + 2 ];
        vflags[ ( y * step ) + x ] = true;
        if( !luv_activated )
        {
            if( !( ( R > 95 ) && ( G > 40 ) && ( B > 20 ) &&
                   ( ( MAX( MAX( R, G ),
                            B ) - MIN( MIN( R, G ), B ) ) > 15 ) && ( abs( R - G ) > 15 ) && ( R > G ) && ( R > B ) ) )
            {
                img->imageData[ ( y * step ) + x * nChannels ] = 0;
                img->imageData[ ( y * step ) + x * nChannels + 1 ] = 255;
                img->imageData[ ( y * step ) + x * nChannels + 2 ] = 0;
                return;
            }
        }
        else
        {
            if( !( ( R > 30 ) && ( ( ( G - 84 ) < 00 ) || ( ( G - 84 ) > -20 ) ) && ( ( -( B - 125 ) ) > 10 ) &&
                   ( ( -( B - 125 ) ) < 60 ) ) )
            {
                img->imageData[ ( y * step ) + x * nChannels ] = 0;
                img->imageData[ ( y * step ) + x * nChannels + 1 ] = 255;
                img->imageData[ ( y * step ) + x * nChannels + 2 ] = 0;
                return;
            }
        }
        if( ( ( y < img->height ) && ( x < img->width ) ) && ( ( y > 0 ) && ( x > 0 ) ) )
        {
            evaluatePixel( x, y + 1 );
            evaluatePixel( x, y - 1 );
            evaluatePixel( x + 1, y );
            evaluatePixel( x - 1, y );
        }
    }
    return;
}


// Calculate min and max from channels inside selected area
void calculateMinMax( IplImage* src )
{
    char* pixel;
    int t;
    int total = 0;
    pixel = src->imageData;
    maxL = maxU = maxV = -INT_MAX;
    minL = minU = minV = INT_MAX;

    // +1 pra não pegar borda do retangulo desenhado
    for( int i = ROI.x + 1; i < ROI.x + ROI.width; i++ )
    {
        for( int j = ROI.y + 1; j < ROI.y + ROI.height; j++ )
        {
            t = ( i * step ) + j * nChannels;
            if( ( unsigned char ) pixel[ t ] > maxL )
            {
                maxL = ( unsigned char ) pixel[ t ];
            }
            if( ( unsigned char ) pixel[ t + 1 ] > maxU )
            {
                maxU = ( unsigned char ) pixel[ t + 1 ];
            }
            if( ( unsigned char ) pixel[ t + 2 ] > maxV )
            {
                maxV = ( unsigned char ) pixel[ t + 2 ];
            }
            if( ( unsigned char ) pixel[ t ] < minL )
            {
                minL = ( unsigned char ) pixel[ t ];
            }
            if( ( unsigned char ) pixel[ t + 1 ] < minU )
            {
                minU = ( unsigned char ) pixel[ t + 1 ];
            }
            if( ( unsigned char ) pixel[ t + 2 ] < minV )
            {
                minV = ( unsigned char ) pixel[ t + 2 ];
            }
            //printf("total = %d\n", total);
            //pixel[t] = 255;
            total++;
        }
    }
    printf( "maxL: %d, maxU: %d, maxV: %d, minL: %d, minU: %d, minV: %d\n",
            maxL, maxU, maxV, minL, minU, minV );
    printf( "total = %d\n", total );
}


void calculateMean( IplImage* src )
{
    char* pixel;
    int total = 0;
    int t;
    pixel = src->imageData;

    int soma_L = 0;
    int soma_U = 0;
    int soma_V = 0;

    for( int i = ROI.x + 1; i < ROI.x + ROI.width; i++ )   //+1 pra não pegar borda do retangulo desenhado
    {
        for( int j = ROI.y + 1; j < ROI.y + ROI.height; j++ )
        {
            t = ( i * step ) + j * nChannels;

            soma_L += ( unsigned char ) pixel[ t ];
            soma_U += ( unsigned char ) pixel[ t + 1 ];
            soma_V += ( unsigned char ) pixel[ t + 2 ];

            total++;
        }
    }

    //int threshold = 15;

    int thold_L, thold_U, thold_V;

    int media_L = ( int ) ( soma_L / total );
    int media_U = ( int ) ( soma_U / total );
    int media_V = ( int ) ( soma_V / total );

    long total_desvioL = 0;
    long total_desvioU = 0;
    long total_desvioV = 0;

    for( int i = ROI.x + 1; i < ROI.x + ROI.width; i++ )   //+1 pra não pegar borda do retangulo desenhado
    {
        for( int j = ROI.y + 1; j < ROI.y + ROI.height; j++ )
        {
            t = ( i * step ) + j * nChannels;
            total_desvioL += ( pixel[ t ] - media_L ) * ( pixel[ t ] - media_L );
            total_desvioU += ( pixel[ t + 1 ] - media_U ) * ( pixel[ t + 1 ] - media_U );
            total_desvioV += ( pixel[ t + 2 ] - media_V ) * ( pixel[ t + 2 ] - media_L );
        }
    }

    thold_L = sqrt( ( float ) total_desvioL / ( float )( ( ROI.height * ROI.width ) - 1 ) );
    thold_U = sqrt( ( float ) total_desvioU / ( float )( ( ROI.height * ROI.width ) - 1 ) );
    thold_V = sqrt( ( float ) total_desvioV / ( float )( ( ROI.height * ROI.width ) - 1 ) );

    minL = media_L - thold_L;
    maxL = media_L + thold_L;

    minU = media_U - thold_U;
    maxU = media_U + thold_U;

    minV = media_V - thold_V;
    maxV = media_V + thold_V;

    printf( "maxL: %d, maxU: %d, maxV: %d, minL: %d, minU: %d, minV: %d\n", maxL, maxU, maxV, minL, minU, minV );
    printf( "thold_L: %d, thold_U: %d, thold_V: %d", thold_L, thold_U, thold_V );
    printf( "total = %d\n", total );
}


// Apply limits
void applyFormula( IplImage* src )
{
    for( int i = 0; i < src->width; i++ )
    {
        for( int j = 0; j < src->height; j++ )
        {
            R = src->imageData[ ( j * step ) + i * nChannels ];
            G = src->imageData[ ( j * step ) + i * nChannels + 1 ];
            B = src->imageData[ ( j * step ) + i * nChannels + 2 ];
            if( !( ( R > minL && R < maxL ) && ( G > minU && G < maxU ) && ( B > minV && B < maxV ) ) )
            {
                src->imageData[ j * step + i * nChannels + 2 ] = 0;
                src->imageData[ j * step + i * nChannels + 1 ] = 0;
                src->imageData[ j * step + i * nChannels ] = 0;
            }
            else
            {
                src->imageData[ j * step + i * nChannels + 2 ] = 255;
                src->imageData[ j * step + i * nChannels + 1 ] = 255;
                src->imageData[ j * step + i * nChannels ] = 255;
            }
        }
    }
}


// Calculate chance of success and failure
void compareGroundtruth()
{
    int verd_pos, verd_neg, falso_pos, falso_neg;
    verd_pos = verd_neg = falso_pos = falso_neg = 0;
    int sub = 0;
    int n_pele = 0;
    int n_naopele = 0;

    for( int i = 0; i < img_gt->width; i++ )
    {
        for( int j = 0; j < img_gt->height; j++ )
        {
            //printf("cor aqui: %d", (unsigned char) img_gt->imageData[j*step + i*nChannels]);
            if( ( unsigned char ) img_gt->imageData[ j * step + i * nChannels ] == 255 )
            {
                n_pele++;
            }
            else
            {
                n_naopele++;
            }
        }
    }

    for( int i = 0; i < img->width; i++ )
    {
        for( int j = 0; j < img->height; j++ )
        {
            sub = ( unsigned char )img_luv->imageData[ j * step + i * nChannels ] -
                  ( unsigned char )img_gt->imageData[ j * step + i * nChannels ];
            //img->imageData[j*step + i*nChannels] = (unsigned char)img->imageData[j*step + i*nChannels] - (unsigned char)img_gt->imageData[(j)*step + i*nChannels];
            switch( sub )
            {
                case 0:
                    if( ( unsigned char ) img_luv->imageData[ j * step + i * nChannels ] == 255 )
                    {
                        verd_pos++;
                    }
                    else
                    {
                        verd_neg++;
                    }
                    break;

                case -255:
                    falso_neg++;
                    break;

                case 255:
                    falso_pos++;
                    break;
            }
        }
    }

    printf( "Total: %d\n", img_gt->width * img_gt->height );
    printf( "Total de pixels com pele: %d\n", n_pele );
    printf( "Total de pixels sem pele: %d\n", n_naopele );
    printf( "Total de verd_pos: %d\n", verd_pos );
    printf( "Total de verd_neg: %d\n", verd_neg );
    printf( "Total de falso_pos: %d\n", falso_pos );
    printf( "Total de falso_neg: %d\n", falso_neg );

    printf( "Porcentagem de verdadeiro positivo: %f\n",
            100 * ( ( float )verd_pos / ( float )( img_gt->width * img_gt->height ) ) );
    printf( "Porcentagem de verdadeiro negativo: %f\n",
            100 * ( ( float )verd_neg / ( float )( img_gt->width * img_gt->height ) ) );
    printf( "Porcentagem de falso positivo: %f\n",
            100 * ( ( float )falso_pos / ( float )( img_gt->width * img_gt->height ) ) );
    printf( "Porcentagem de falso negativo: %f\n",
            100 * ( ( float )falso_neg / ( float )( img_gt->width * img_gt->height ) ) );
}


// Detect face and eyes
void detectEyes( IplImage* src, std::string spaceColor, IplImage* src_luv )
{
    //int i;
    //int fcenter_x = 0;
    //int fcenter_y = 0;

    /* detect faces */
    CvSeq* faces = cvHaarDetectObjects(
        src, cascade_f, storage,
        1.1, 3, 0, cvSize( 40, 40 ) );

    /* return if not found */
    if( faces->total == 0 )
    {
        return;
    }

    /* draw a rectangle */
    CvRect* face = ( CvRect* )cvGetSeqElem( faces, 0 );
//	cvRectangle(src,
//				cvPoint(face->x, face->y),
//				cvPoint(face->x + face->width, face->y + face->height),
//				CV_RGB(255, 0, 0), 1, 8, 0);

    /* reset buffer for the next object detection */
    cvClearMemStorage( storage );

//    /* Set the Region of Interest: estimate the eyes' position */
//    cvSetImageROI(src, cvRect(face->x, face->y + (face->height/5.5), face->width, face->height/3.0));
//
//    /* detect eyes */
//	CvSeq* eyes = cvHaarDetectObjects(
//        src, cascade_e, storage,
//		1.15, 3, 0, cvSize(25, 15));
//
//    /* draw a rectangle for each eye found */
//	for( i = 0; i < (eyes ? eyes->total : 0); i++ ) {
//		eye = (CvRect*)cvGetSeqElem( eyes, i );
//		cvRectangle(src,
//					cvPoint(eye->x, eye->y),
//					cvPoint(eye->x + eye->width, eye->y + eye->height),
//					CV_RGB(255, 0, 0), 1, 8, 0);
//        //ROI = cvRect(face->x + (face->x/4), face->y + (face->y/4), face->width/2, (face->height/6));
//        ROI = cvRect(face->x + (face->x/4), face->y + (face->y/4), face->width/2, (face->height/6));
//        calcula_min_max(src);
//	}
    //ROI = cvRect(face->x + (face->x/4), face->y + (face->y/4), face->width/2, (face->height/6));
    //ROI = cvRect(face->y + (face->y/0.8), face->x + (face->x/6), (face->height/5),face->width/1.5 );
    ROI = cvRect( face->y + 2.5 * ( face->width / 5 ), face->x + face->width / 5, face->height / 4, face->width / 1.7 );
    if( spaceColor.compare( "luv" ) == 0 )
    {
        calculateMinMax( src_luv );
        applyFormula( src_luv );
    }
    else
    {
        calculateMinMax( src );
        applyFormula( src );
    }
    cvResetImageROI( src );
}


void drawRectangle( int init_x, int init_y, int x, int y )
{
    cvResetImageROI( img );
    cvCopyImage( tempData, img );
    //img->imageData = tempData;
    cvRectangle(
        img,
        cvPoint( init_x, init_y ),
        cvPoint(
            init_x + desl_x,
            init_y + desl_y
            ),
        CV_RGB( 255, 0, 0 ),
        1, 8, 0
        );

    ROI = cvRect( init_y + desl_y, init_x, abs( desl_y ), abs( desl_x ) );

    cvSetImageROI( img, ROI );
}


void mouse( int button, int state, int x, int y )
{
    mouseX = ( int ) x;
    mouseY = ( int ) ( img->height - y );
    desl_x = init_x = mouseX;
    desl_y = init_y = mouseY;
    buttonGlobal = button;
    stateGlobal = state;

    //printf("mousex: %d, mousey: %d", mouseX, mouseY);

    if( button == GLUT_LEFT_BUTTON )
    {
        if( state == GLUT_DOWN )
        {
            //img->imageData[(x*step) + y*nChannels] = 255; img->imageData[(x*step) + (y*nChannels)+1] = 255; img->imageData[(x*step) + (y*nChannels)+2] = 255;
            //img->imageData[((int)mouseY*step) + (int)mouseX*nChannels] = 255; img->imageData[((int)mouseY*step) + ((int)mouseX*nChannels)+1] = 255; img->imageData[((int)mouseY*step) + ((int)mouseX*nChannels)+2] = 255;
            printf( "Valores do Pixel - R: %d - G: %d - B: %d\n",
                    ( uchar ) img->imageData[ ( mouseY * step ) + mouseX * nChannels ],
                    ( uchar ) img->imageData[ ( mouseY * step ) + mouseX * nChannels + 1 ],
                    ( uchar ) img->imageData[ ( mouseY * step ) + mouseX * nChannels + 2 ] );
            printf( "Valores do Pixel - L: %d - U: %d - V: %d\n",
                    ( uchar ) img->imageData[ ( mouseY * step ) + mouseX * nChannels ],
                    ( uchar ) img->imageData[ ( mouseY * step ) + mouseX * nChannels + 1 ] - 84,
                    ( uchar ) img->imageData[ ( mouseY * step ) + mouseX * nChannels + 2 ] - 125 );
            //avalia_pixel(mouseX, mouseY);
//             printf("Valores do Pixel - R: %d - G: %d - B: %d\n", luv_data[(mouseY)*step + mouseX*nChannels],
//                luv_data[(mouseY)*step + mouseX*nChannels+1],
//                luv_data[(mouseY)*step + mouseX*nChannels+2]);
            printf( "mouseX: %d - mouseY: %d\n", mouseX, mouseY );
        }
    }

    glutPostRedisplay();
    glutPostWindowRedisplay( janID );
}


void mouse_gt( int button, int state, int x, int y )
{
    int mouseX_gt, mouseY_gt;
    mouseX_gt = ( int ) x;
    mouseY_gt = ( int ) ( img_gt->height - y );
//    desl_x = init_x = mouseX;
//    desl_y = init_y = mouseY;
//    buttonGlobal = button;
//    stateGlobal = state;

    //printf("mousex: %d, mousey: %d", mouseX, mouseY);

    if( button == GLUT_LEFT_BUTTON )
    {
        if( state == GLUT_DOWN )
        {
            printf( "Ground Truth: Valores do Pixel - R: %d - G: %d - B: %d\n",
                    ( uchar ) img_gt->imageData[ ( mouseY_gt * step ) + mouseX_gt * nChannels ],
                    ( uchar ) img_gt->imageData[ ( mouseY_gt * step ) + mouseX_gt * nChannels + 1 ],
                    ( uchar ) img_gt->imageData[ ( mouseY_gt * step ) + mouseX_gt * nChannels + 2 ] );
        }
    }

    glutPostRedisplay();
}


void idle()
{
    glutPostRedisplay();
}


//MouseMotion callback
void mouseMotion( int x, int y )
{
    mouseX = ( float ) x;
    mouseY = ( float ) ( img->height - y );
    desl_x = mouseX - init_x;
    desl_y = mouseY - init_y;

    //printf("mousex: %d, mousey: %d", mouseX, mouseY);
    if( buttonGlobal == GLUT_LEFT_BUTTON )
    {
        if( stateGlobal == GLUT_DOWN )
        {
            drawRectangle( init_x, init_y, desl_x, desl_y );
        }
    }

    printf( "init_x: %d, init_y: %d, desl_x: %d, desl_y: %d\n", init_x, init_y, desl_x, desl_y ),

    glutPostRedisplay();
}


void passiveMotion( int x, int y )
{
    mouseX = ( float ) x;
    mouseY = ( float ) ( img->height - y );

    glutPostRedisplay();
}


static void key( unsigned char key, int x, int y )
{
    switch( key )
    {
        case 's':
            int p[ 3 ];
            p[ 0 ] = CV_IMWRITE_PNG_COMPRESSION;
            p[ 1 ] = 100;
            p[ 2 ] = 0;
            revertBack();
            cvSaveImage( "teste.png", imgAux, p );
            break;

        case 'r':
            //avalia_pixel(50, 50);
            //preenche_flag();
            break;

        case 't':
            printf( "tamanho do velor em bytes: %d", sizeof( bool ) * ( img->width ) * ( img->height ) );
            printf( "numero de pixels %d", ( img->width ) * ( img->height ) );
            break;

        case '1':
            img_luv = cvCreateImage( cvSize( img->width, img->height ), IPL_DEPTH_8U, 3 );
            cvCopyImage( img, img_luv );
            convertion( img_luv, std::string("rgb2luv"), step, nChannels );
            luv_activated = true;
            //cvCopyImage(img, tempData);
            //cvCopyImage(img, img_rgb);
            //LUVToRGB(img, step, nChannels);
            break;

        case '2':
            convertion( img, "luv2rgb", step, nChannels );
            cvCopyImage( img, tempData );
            //LUVToRGB(img, step, nChannels);
            break;

        case 'a':
            for( int i = 0; i < img->width; i++ )
            {
                for( int j = 0; j < img->height; j++ )
                {
                    evaluatePixel( i, j );
                }
            }
            break;

        case 'd':
            cvCopyImage( img_rgb, tempData );
            revert( imgAux, img_luv );

            // detecta pelo img_rgb (imagem original) e modifica na que esta na tela
            detectEyes( tempData, std::string( "luv" ), imgAux );
            revert( img_luv, imgAux );

            cvCopyImage( img_rgb, tempData );

            // detecta pelo img_rgb (imagem original) e modifica na que esta na tela
            detectEyes( tempData, std::string( "rgb" ), NULL );
            revert( img, tempData );

            break;

        case 'f':
            calculateMinMax( img );
            applyFormula( img );

            calculateMinMax( img_luv );
            applyFormula( img_luv );
            break;

        case 'm':
            calculateMean( img );
            applyFormula( img );

            calculateMean( img_luv );
            applyFormula( img_luv );
            break;

        case 'c':
            compareGroundtruth();
            break;
    }

    glutPostRedisplay();
    glutPostWindowRedisplay( janID );
}


void init( void )
{
    ///Inicializar sistema de viz
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho( 0.0, img->width, 0.0, img->height, 0.0, 2.0 );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    mouseX = img->width;
    mouseY = img->height;
}


int main( int argc, char* argv[] )
{
    glutInit( &argc, argv );

    std::string file1 ( "haarcascade_frontalface_alt.xml" );
    std::string file2 ( "haarcascade_eye.xml" );

    /* load the face classifier */
    cascade_f = ( CvHaarClassifierCascade* )cvLoad( file1.c_str(), 0, 0, 0 );

    /* load the eye classifier */
    cascade_e = ( CvHaarClassifierCascade* )cvLoad( file2.c_str(), 0, 0, 0 );

    /* setup memory storage, needed by the object detector */
    storage = cvCreateMemStorage( 0 );

    img = cvLoadImage( "images/soldado.bmp", CV_LOAD_IMAGE_COLOR );
    img_gt = cvLoadImage( "images/soldado.bmp", CV_LOAD_IMAGE_COLOR );

    /* always check */
    assert( cascade_f && cascade_e && storage && img );

    //img = cvLoadImage("rgb.png", CV_LOAD_IMAGE_COLOR );
    //img = cvLoadImage("teste.png", CV_LOAD_IMAGE_COLOR );

    imgAux = cvCreateImage( cvSize( img->width, img->height ), IPL_DEPTH_8U, 3 );

    /* print some properties */
    printf( "Filename:    %s\n", argv[ 1 ] );
    printf( "# channels:  %d\n", img->nChannels );
    printf( "Pixel depth: %d bits\n", img->depth );
    printf( "width:       %d pixels\n", img->width );
    printf( "height:      %d pixels\n", img->height );
    printf( "Image size:  %d bytes\n", img->imageSize );
    printf( "Width step:  %d bytes\n", img->widthStep );
    step = img->widthStep;
    nChannels = img->nChannels;

    img_rgb = cvCreateImage( cvSize( img->width, img->height ), IPL_DEPTH_8U, 3 );
    cvCopyImage( img, img_rgb );

    cvCopyImage( img, imgAux );
    revert( img, imgAux );
    cvCopyImage( img_gt, imgAux );
    revert( img_gt, imgAux );

    //preenche_flag();

    //tempData = (char*) malloc(img->imageSize);
    tempData = cvCreateImage( cvSize( img->width, img->height ), IPL_DEPTH_8U, 3 );
    cvCopyImage( img, tempData );

    img_luv = cvCreateImage( cvSize( img->width, img->height ), IPL_DEPTH_8U, 3 );
    cvCopyImage( img, img_luv );
    convertion( img_luv, "rgb2luv", step, nChannels );
    luv_activated = true;

    //convertion(img, "rgb2luv", step, nChannels);
    //RGBToLUV(img, step, nChannels);

    printf( " _____________________________________\n" );
    printf( "Dynamic Segmentation of Human Skin Using LUV Color Space\n" );
    printf( "Press 'd' to automatic skin segmentation using face detection\n" );
    printf( "Press 's' to save image\n" );

    glutInitWindowSize( img->width, img->height );
    glutInitWindowPosition( 500, 100 );
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );

    glutCreateWindow( "Select skin area" );

    glutKeyboardFunc( key );
    glutMouseFunc( mouse );
    glutMotionFunc( mouseMotion );
    glutPassiveMotionFunc( passiveMotion );
    glutDisplayFunc( display );

    //### First window (ground truth) ###
    glutInitWindowPosition( 100, 100 );
    glutCreateWindow( "Ground Truth" );
    glutDisplayFunc( display_gt );
    glutMouseFunc( mouse_gt );
    //######################

    //### Third window (LUV) ###
    glutInitWindowPosition( 100, 500 );
    janID = glutCreateWindow( "LUV" );
    glutDisplayFunc( display_luv );
    glutIdleFunc( idle );
    //######################

    init();

    glutMainLoop();

    return 1;
}



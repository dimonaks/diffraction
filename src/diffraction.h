#include <iostream>//
#include <string>//
#include <cstring> 
#include <algorithm>
#include <fstream>
#include <ctime>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <cmath>
#include <new>
#include <complex>
#include <sstream>

#define type double
#define db double
#define MAX_LATTICE_SIZE 50
#define MAX_NUM_SADP 50     //maximum number of patterns in one run
#define NUM_ATOMIC_NUMBERS 24
#define MAX_NUM_REFLEX 10000 //maximum number of reflexes
using namespace std;

//Обявление структур и типов.
class vector {
    public:
    double x,y,z;
    public:
    vector(double vX=0, double vY=0, double vZ=0) { x=vX; y=vY; z=vZ; } // конструктор
    // ~vector();
    friend vector operator+(vector a, vector b); //сложение векторов
    friend vector operator-(vector a, vector b); //вычитание векторов
    friend vector operator*(vector a, double scalar); //умножение вектора на скаляр
    friend double operator*(vector a, vector b); //скалярное произведение векторов
    friend vector operator%(vector a, vector b); //векторное произведение векторов
    friend vector operator/(vector a, double V);//деление компонент вектора на скаляр
    friend double abs(vector a); //модуль вектора
    void print();//вывод вектора в терминале
    };

vector operator+(vector a, vector b) {
    vector temp;
    temp.x = a.x + b.x;
    temp.y = a.y + b.y;
    temp.z = a.z + b.z;
    return temp;
    }

vector operator-(vector a, vector b) {
    vector temp;
    temp.x = a.x - b.x;
    temp.y = a.y - b.y;
    temp.z = a.z - b.z;
    return temp;
    }

vector operator*(vector a, double scalar) {
    vector temp;
    temp.x = a.x * scalar;
    temp.y = a.y * scalar;
    temp.z = a.z * scalar;
    return temp;
    }

double operator*(vector a, vector b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
    }

vector operator%(vector a, vector b) {
    vector d;
    d.x=(a.y*b.z - a.z*b.y);
    d.y=(a.z*b.x - a.x*b.z);
    d.z=(a.x*b.y - a.y*b.x);
    return d;
    }

double abs(vector a) {
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    }

vector operator/(vector a, double V) {
    vector d;
    d.x=a.x/V;d.y=a.y/V;d.z=a.z/V;
    return d;
    }

void vector::print() {
    printf("(%f,%f,%f) ",x,y,z);
    }

struct vector_int {
    int h;int k;int l; double length();
    };

double vector_int::length(){return(sqrt(h*h+k*k+l*l));}

struct vector_dec {
    double x;double y;double z; double length();
    };

double vector_dec::length(){return(sqrt(x*x+y*y+z*z));}

typedef complex<double> compln;



class LayoutClass {
    
    public:
    db uvw[MAX_NUM_SADP][3]; //Miller indices of axis zone
    int number_of_stuctures, n_of_coloumns, n_of_rows;
    string struct_label[MAX_NUM_SADP];
    string struct_filename[MAX_NUM_SADP];
    string config_filename[MAX_NUM_SADP];

    public:
    // LayoutClass();
    void readin(char* layout_filename);
};



class SADPClass {


//For SADP parameters:     
    int ni, nj, nk; // lattice sizes

    int cartesian_plane;//если 1, то для задания сечения испольуется вектор n_decart[3] в дек. коорд.; он переводится в нормаль для справки, всегда используется внутри программы n_decart[3]  ;
    //если 0, то используется n и вычисляется n_decart для последующего поворота декартовой системы координат Т.е. плоскость в обратном пространстве задается вектором в обратном пространстве;

    //если 2, вводится направленине [UVW], которое переводиться в декартовы координаты.здесь плоскость в обратном задается вектором в прямом;
    //Т.е. это ось зоны

    db ewald_thickness, F_structure_critery, zoom, Fsum_need, r_of_electrongram;

    int scale_index[3];

    db reflex_fontsize; //size of font for reflex names

    db red;


    int n[3], p[3], UVW[3];
    db n_decart[3];


    // int max_num_reflex;//Максимальное кол-во рефлексов на картинке
    int N_of_arrays;//Для динамического определения массивов(количество проверяемых узлов) =(ni+1)*(nj+1)*(nk+1)*8 нет смысла делать


//For crystal structure:
    double b[3][3],a[3][3];
    int b_int[3][3];
    vector_dec basis[10000];
    int nbasis;
    double x_basis[10000],y_basis[10000],z_basis[10000];
    int typat[10000];
    
    //For atomic factors:
    int num_of_data;
    db s_deltaK[NUM_ATOMIC_NUMBERS], F[2][NUM_ATOMIC_NUMBERS];  


    // public:

    public:
    db read_config(string name_of_config_file, db *ind);
    void read_crystal_structure(string name_of_srtucture_file);
    void read_atomic_factors(string atomic_factors_filename);





};














inline string postscript_header(db *boundingbox, db reflex_fontsize) {
    ostringstream strs;
    db conv = 2.834645669;
    strs << "\n\
%!PS-Adobe-3.0 EPSF-3.0\n\
%%Creator: diffraction\n\
%%Title: electronogramm\n\
%%CreationDate: [date the file was created]\n\
%%DocumentData: Clean7Bit\n\
%%Origin: 0 0\n\
%%BoundingBox: " << boundingbox[0] * conv << " " << boundingbox[1] * conv << " " \
                 << boundingbox[2] * conv << " " << boundingbox[3] * conv << " \n\
%%LanguageLevel: 2 \n\
%%Pages: 1\n\
%%Page: 1 1\n\
72 25.4 div                                   % 1 мм = 72/25.4 пунктов\n\
dup                                           % дублировать значение на вершине стека\n\
scale                                         % растянуть в это количество раз по обеим координатам \n\
/fontsize " << reflex_fontsize << " def\n\n\
/k_of_text_pos 1.2 def \n\
/k_horiz_pos k_of_text_pos fontsize mul def   %множ. для горизонтальный сдвига названия\n\
/k_vertical_pos 1.4 def                       %множ. для вертик. сдвига названия\n\
/Times-Roman findfont                         % взять шрифт Times-Roman\n\
fontsize scalefont                            % растянуть до размера fontsize (у нас единица измерения - мм!)\n\
setfont                                       % установить выбранный шрифт\n\
\n\
/reflex {\n\
/r exch def                                   %записать радиус рефлекса в r\n\
r 0 360 arc                                   %нарисовать кружочек в нужном месте радиуса r, в стеке остается только название рефлекса\n\
gsave\n\
fill\n\
grestore                                      %нужно для независимой закраски кружков и возврата в нужную позицию\n\
/l exch def                                   %считывание индексов и их знаков\n\
/sl exch def\n\
/k exch def\n\
/sk exch def\n\
/h exch def\n\
/sh exch def\n\
\n\
0 r sub k_horiz_pos sub                       % высчитать гор. смещение названия (отрицательное)\n\
r k_vertical_pos mul                          % высчитать верт. смещение названия\n\
rmoveto                                       %перейти по смещениям\n\
\n\
\n\
\n\
/drawminus{\n\
(n) eq                                        %равен ли знак(sh,sk,sl) n, т.е. отрицателeн ли индекс?-результат в стек\n\
{\n\
/xdigsize fontsize 3 div def %size of digit\n\
gsave\n\
0 \n\
fontsize fontsize 4 div sub\n\
rmoveto\n\
\n\
\n\
\n\
(/) search\n\
{xdigsize 3 mul 0 rlineto}\n\
{xdigsize 0 rlineto} ifelse\n\
\n\
0.15 setlinewidth                             %толщина минуса \n\
stroke\n\
grestore\n\
} if\n\
\n\
\n\
\n\
} def\n\
\n\
\n\
h sh drawminus\n\
h show\n\
k sk drawminus\n\
k show\n\
l sl drawminus\n\
l show\n\
\n\
\n\
fill} def\n\
";

    return strs.str();

};



inline string postscript_layout(int i_st, int ic, int ir, string SADP, LayoutClass layout, db *boundingbox, db r_of_electrongram ) {

    string SADP_name = layout.struct_label[i_st];

    db xshift_on_figure = boundingbox[2] / layout.n_of_coloumns * ic + boundingbox[2] / layout.n_of_coloumns / 2. ;
    db yshift_on_figure = boundingbox[3] / layout.n_of_rows     * ir + boundingbox[3] / layout.n_of_rows     / 2. ;
    yshift_on_figure += (boundingbox[3] / 50.); //чуть чуть приподнять ряды, чтобы влезли буквы
    
    ostringstream strs;

    //1. Calculate position for the current SADP and prepare Fonts
    strs << xshift_on_figure << " " << yshift_on_figure <<             " translate   % установить начало координат\n" << \
    0 - r_of_electrongram / 7 << " " << 0 - r_of_electrongram - 9.1 << " moveto\n" << \
    "/Times-Roman findfont % взять шрифт Times-Roman\n" << \
    "8 scalefont\n" << \
    "setfont               % установить выбранный шрифт\n" << \
    "(" << SADP_name << ")" << " show fill\n" << \
    "/Times-Roman findfont % взять шрифт Times-Roman\n" << \
    "fontsize scalefont\n" << \
    "setfont               % установить выбранный шрифт\n";


    strs << SADP; // 2. Add coordinates and intensity of reflexes
    

    // 3. Return to the origin
    strs << 0 - xshift_on_figure << " " << 0 - yshift_on_figure <<      " translate   % вернуться на прежнее место\n";

    return strs.str();
    
}








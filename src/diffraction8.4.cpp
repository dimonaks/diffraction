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
#define type double
#define db double
#define macro_size 50
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

//Объявление функций
void convert_b_to_integer(void);
double calculate_structure_factor(double,vector_int);
void read_input(char*,char*, db *ind);//nbasis,basis,typat,num_of_data,s_deltaK,F
int plane(int,int,int);
int plane2(db,db,db);
void calculate_hkl();
void convert_normal();
void calculation_and_output();
extern "C" {void spline_(int&,double*,double*,double*,double*,double*);}
extern "C" {void sevfp_(int&,double&,double*,double*,double*,double*,double*,double&,double&);}
//Конец объявления функций


//Global var
double b[3][3],a[3][3];
int b_int[3][3];
vector_dec basis[10000];







//read from file
db reflex_fontsize; //size of font for reflex names
int ni=macro_size,nj=macro_size,nk=macro_size;
int nbasis;
double x_basis[10000],y_basis[10000],z_basis[10000];
int typat[10000];
int num_of_data=24;
double s_deltaK[24],F[2][24];
int cartesian_plane=1;//если 1, то для задания сечения испольуется вектор n_decart[3] в дек. коорд.; он переводится в нормаль для справки, всегда используется внутри программы n_decart[3]  ;
//если 0, то используется n и вычисляется n_decart для последующего поворота декартовой системы координат Т.е. плоскость в обратном пространстве задается вектором в обратном пространстве;

//если 2, вводится направленине [UVW], которое переводиться в декартовы координаты.здесь плоскость в обратном задается вектором в прямом;
//Т.е. это ось зоны

db ewald_thickness=1e-3;
db F_structure_critery=0.01;

int n[3]={1,1,-1}, p[3]={0,0,0};
db n_decart[3]={1,1,1};
db UVW[3];


db zoom=100,red=10, Fsum_need = 1000;

int scale_index[3];
db r_of_electrongram;
//
//N=100000;//ni*nj*nk*8;

int N=10000;//Максимальное кол-во рефлексов на картинке
int N_of_arrays;//Для динамического определения массивов(количество проверяемых узлов) =(ni+1)*(nj+1)*(nk+1)*8 нет смысла делать



//char names_of_stucture_files[100][100], names_of_config_files[100][100];


db boundingbox[4]={0,0,210,150};//
int n_of_coloumns=4;//количество столбцов и строк электронограм на выводимом рисунке порядок нумерации - направо потом вверх
int n_of_rows=2;


ofstream out;

//Используется кристаллографическое определение векторов обратной решетки.
main(int argc, char *argv[])
{
if(argc!=2){cout<<"Not enouth input data"<< endl;return 0;}
cout<<"\n";
cout <<"Name of program="<<argv[0]<<endl;
cout << "Name of input file="<<argv[1]<<endl;
//local var
//nameofoutps=argv[1];
char nameofoutps[100];
strcpy(nameofoutps,"output/");
strcat(nameofoutps,argv[1]);
//strcat(nameofoutps,"[");
//strcat(nameofoutps,argv[2]);
//strcat(nameofoutps,argv[3]);
//strcat(nameofoutps,argv[4]);
//strcat(nameofoutps,"]");
strcat(nameofoutps,".eps");
cout<<"Name of output file!!!="<<nameofoutps<<endl;




//Считываем имена файлов со структурами и имена кофигурационных файлов к ним из файла списка

int i_structure;
int number_of_stuctures;
string names_of_structures[100];
number_of_stuctures=100;

db index1[number_of_stuctures],index2[number_of_stuctures],index3[number_of_stuctures];
db ind[3];//вспомогательные переменные.

ifstream in; // input

in.open(argv[1], ios::in);
in >>number_of_stuctures;
cout<<"number_of_stuctures= "<<number_of_stuctures<<endl;
in >>n_of_coloumns;
cout<<"n_of_coloumns= "<<n_of_coloumns<<endl;
in >>n_of_rows;
cout<<"n_of_rows= "<<n_of_rows<<endl;

char **names_of_stucture_files=new char*[number_of_stuctures]; for(int i=0; i < number_of_stuctures; i++)\
names_of_stucture_files[i]=new char[100];
char **names_of_config_files=new char *[number_of_stuctures]; for(int i=0; i < number_of_stuctures; i++)\
names_of_config_files[i]=new char[100];

for (int i=0;i<number_of_stuctures;i++){
in >>index1[i];
in >>index2[i];
in >>index3[i];
in >>names_of_stucture_files[i];
in >>names_of_config_files[i];
in >>names_of_structures[i];
cout<<names_of_stucture_files[i]<<endl;
}
in.close();


out.open(nameofoutps, ios::out); // открываем файл для записи

//Выводим postscript
out << "\n\
%!PS-Adobe-3.0 EPSF-3.0\n\
%%Creator: diffraction\n\
%%Title: electronogramm\n\
%%CreationDate: [date the file was created]\n\
%%DocumentData: Clean7Bit\n\
%%Origin: 0 0\n\
%%BoundingBox: ";
out<<boundingbox[0]*2.834645669<<" "<<boundingbox[1]*2.834645669<<\
" "<<boundingbox[2]*2.834645669<<" "<<boundingbox[3]*2.834645669<<" \n\
%%LanguageLevel: 2 \n\
%%Pages: 1\n\
%%Page: 1 1\n\
72 25.4 div           % 1 мм = 72/25.4 пунктов\n\
dup                   % дублировать значение на вершине стека\n\
scale                 % растянуть в это количество раз по обеим координатам\n"; 



int ic=0,ir=0;

db xshift_on_figure, yshift_on_figure;
//1. Считываем информацию о структуре из файлов и трансформируем вектора


for (i_structure=0;i_structure < number_of_stuctures;i_structure++){
if(ic>=n_of_coloumns){ic=0;ir++;}
//if(ir>=n_of_rows){ir=0;}

//просто записываем заданные в файле индексы направления для конкретной структуры.
ind[0] = index1[i_structure];
ind[1] = index2[i_structure];
ind[2] = index3[i_structure];


cout <<"\n\n\n\n!!!!!Name of structure   "<<names_of_stucture_files[i_structure]<<endl;
cout <<"-----------------------------------------------------------------------------\n";
read_input(names_of_config_files[i_structure],names_of_stucture_files[i_structure],ind);//nbasis,basis,typat,num_of_data,s_deltaK,F






if(i_structure==0){
out<<"/fontsize "<< reflex_fontsize<<" def\n";
out<<"\n\
/k_of_text_pos 1.2 def \n\
/k_horiz_pos k_of_text_pos fontsize mul def   %множ. для горизонтальный сдвига названия\n\
/k_vertical_pos 1.4 def %множ. для вертик. сдвига названия\n\
/Times-Roman findfont % взять шрифт Times-Roman\n\
fontsize scalefont          % растянуть до размера fontsize (у нас единица измерения - мм!)\n\
setfont               % установить выбранный шрифт\n\
\n\
/reflex {\n\
/r exch def %записать радиус рефлекса в r\n\
r 0 360 arc %нарисовать кружочек в нужном месте радиуса r, в стеке остается только название рефлекса\n\
gsave\n\
fill\n\
grestore %нужно для независимой закраски кружков и возврата в нужную позицию\n\
/l exch def %считывание индексов и их знаков\n\
/sl exch def\n\
/k exch def\n\
/sk exch def\n\
/h exch def\n\
/sh exch def\n\
\n\
0 r sub k_horiz_pos sub% высчитать гор. смещение названия (отрицательное)\n\
r k_vertical_pos mul % высчитать верт. смещение названия\n\
rmoveto %перейти по смещениям\n\
\n\
\n\
\n\
/drawminus{\n\
(n) eq %равен ли знак(sh,sk,sl) n, т.е. отрицателeн ли индекс?-результат в стек\n\
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
0.15 setlinewidth %толщина минуса \n\
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
}










xshift_on_figure=boundingbox[2]/n_of_coloumns*ic+boundingbox[2]/n_of_coloumns/2;
yshift_on_figure=boundingbox[3]/n_of_rows*ir+boundingbox[3]/n_of_rows/2. ;
//if(ir == 0) 
yshift_on_figure += (boundingbox[3]/50.);//чуть чуть приподнять ряды, чтобы влезли буквы
ic++;


//xshift_on_figure=210./number_of_stuctures*(i+1)-210./number_of_stuctures/2;
//xshift_on_figure=boundingbox[2]/number_of_stuctures*(i+1)-boundingbox[2]/number_of_stuctures/2;
//if(i==0){ yshift_on_figure=75; xshift_on_figure=210./number_of_stuctures/2;}
//else {xshift_on_figure=boundingbox[2]/number_of_stuctures; yshift_on_figure=0;}


//это смещения для надписей !
out<<xshift_on_figure<<" "<<yshift_on_figure<<" translate              % установить начало координат\n";
out<<0-r_of_electrongram/7<<" "<<0-r_of_electrongram-9.1<<" moveto"<<endl; 
out<<"/Times-Roman findfont % взять шрифт Times-Roman\n";
out<<"8 scalefont"<<endl;
out<<"setfont               % установить выбранный шрифт\n";
out<<"("<<names_of_structures[i_structure]<<")"<<" show fill"<<endl;
out<<"/Times-Roman findfont % взять шрифт Times-Roman\n";
out<<"fontsize scalefont"<<endl;
out<<"setfont               % установить выбранный шрифт\n";

calculation_and_output();
out<<0-xshift_on_figure<<" "<<0-yshift_on_figure<<" translate              % вернуться на прежнее место\n";

}
out<<"showpage              % вывести страницу\n%%EOF\n";
out.close();
}











///////////////////////Функции:
void calculation_and_output(){
int h[N],k[N],l[N];
double x[N],y[N],z[N];
double F_structure[N];
//этот вектор от s можно сделать
vector_dec dec_new[N];





cout<<"\n\n\n";
//Нормализуем вектора обратной решетки
//convert_b_to_integer();//работает не всегда, поэтому нигде не используется, только как доп. информация
cout<<"\n\n\n";
//Получаем нормаль в декартовых координатах, если она задана в векторах обратной решетки или
//Переводим нормаль в декартовых координатах n_decart(соответсвтует кубической обратной решетке) в нормль n в векторах обратной решетки
convert_normal();
cout<<"\n\n\n";
//Преобразуем толщину слоя в параметр dt*(n_decart[0]+n_decart[1]+n_decart[2])
//dt-параметр, 2*dt*sqrt(3)-толщина плоскости, когда условие в виде:
//fabs(n_decart[0]*(x-p[0])+n_decart[1]*(y-p[1])+n_decart[2]*(z-p[2]))<=dt*(n_decart[0]+n_decart[1]+n_decart[2])
//Преобразуем толщину слоя в параметр dt*(n_decart[0]+n_decart[1]+n_decart[2])
db n_decart_length=sqrt(n_decart[0]*n_decart[0]+n_decart[1]*n_decart[1]+n_decart[2]*n_decart[2]);//это позволяет уйти от зависимости когда вектор
ewald_thickness=ewald_thickness/2.*(n_decart[0]*n_decart[0]+n_decart[1]*n_decart[1]+n_decart[2]*n_decart[2])/n_decart_length;
//нормали просто увеличен или уменьшен
cout<<"ewald_thickness after transform "<<ewald_thickness<<endl;


//cout <<"N =" <<N;
vector_int tmp, recip;
vector_dec tmp2;
db d,w10,Fsum = 0;
int i,j,ki,s=0;
for(i=-ni;i<ni+1;i++)
	for(j=-nj;j<nj+1;j++)
		for(ki=-nk;ki<nk+1;ki++)
{
//cout << "s= "<<s << "\n";

//tmp.h=i*b_int[0][0]+j*b_int[1][0]+ki*b_int[2][0];
//tmp.k=i*b_int[0][1]+j*b_int[1][1]+ki*b_int[2][1];
//tmp.l=i*b_int[0][2]+j*b_int[1][2]+ki*b_int[2][2];
//tmp2.x=i*b[0][0]+j*b[1][0]+ki*b[2][0];
//tmp2.y=i*b[0][1]+j*b[1][1]+ki*b[2][1];
//tmp2.z=i*b[0][2]+j*b[1][2]+ki*b[2][2];
x[s]=i*b[0][0]+j*b[1][0]+ki*b[2][0];
y[s]=i*b[0][1]+j*b[1][1]+ki*b[2][1];
z[s]=i*b[0][2]+j*b[1][2]+ki*b[2][2];
//if(tmp2.x==0&&tmp2.y==0&&tmp2.z==0){
//cout<<tmp2.x<<" "<<tmp2.y<<" "<<tmp2.z<<endl;
//cout<<plane2(tmp2.x,tmp2.y,tmp2.z)<<"\n";
//}
if(plane2(x[s],y[s],z[s])==1){
//cout<<tmp2.x<<" "<<tmp2.y<<" "<<tmp2.z<<endl;
h[s]=i;
k[s]=j;
l[s]=ki;
recip.h=i;
recip.k=j;
recip.l=ki;
//printf( "%i %i %i \n",h[s],k[s],l[s] );
//K_length=0;
d=sqrt(x[s]*x[s]+y[s]*y[s]+z[s]*z[s]);//верно, только когда p[3]={0,0,0}
//K_length=d[s];
F_structure[s]=calculate_structure_factor(d,recip);//уже возведенный в квадрат структурный фактор по модулю!!!
//cout << "F_structure[s]="<< F_structure[s] << "\n";
//cout << "K_length="<< K_length << "\n";
//cout << "recip[s].h="<< recip[s].h<< "\n";
if(F_structure[s]>=F_structure_critery&&d<=r_of_electrongram/zoom)//доп парметр-радиус, делает электроннограмму круглой, обрезая рефлексы
{
w10 = n_decart[0]*x[s]+n_decart[1]*y[s]+n_decart[2]*z[s];
cout << "F_structure[s]="<< F_structure[s] << ", distance from plane = "<<w10<<"\n";
Fsum += F_structure[s];
s++;
}
}
if(s>=N)cout<<"Error, s>=N, the size of arrays are two small";
}
//Поворачиваем систему координат, чтобы получить координаты точек в плоскости, происходит совмещение нормали с осью Z.
//vector_dec coord_new[N];
//находим угол между осью z и n_decart.
//double n_decart_length=sqrt(n_decart[0]*n_decart[0]+n_decart[1]*n_decart[1]+n_decart[2]*n_decart[2]);
double cosfi,fi=acos(n_decart[2]/n_decart_length);
//cout << "n_decart[2]= " << n_decart[2]<<"\n";
//cout << "n_decart_length= " << n_decart_length<<"\n";
//fi=0;fi=fi*M_PI/180;
cout << "fi= " << fi*180/M_PI<<"\n";
cosfi=cos(fi);
double sinfi=sin(fi);
double omc=1-cosfi;
vector_dec rotate_axis;


rotate_axis.x=n_decart[1];
rotate_axis.y=-n_decart[0];
rotate_axis.z=0;
db rh=rotate_axis.x/rotate_axis.length();
db rk=rotate_axis.y/rotate_axis.length();
db rl=0;
if(rotate_axis.length()==0){rh=0;rk=0;}
cout << "rh,rk= " << rh<<" "<<rk<<"\n";
cout << "rotate_axis.length= " << rotate_axis.length()<<"\n";
for(int i=0;i<s;i++){
dec_new[i].x=x[i]*(cosfi+omc*rh*rh)+y[i]*omc*rh*rk+z[i]*sinfi*rk;
dec_new[i].y=x[i]*omc*rh*rk+y[i]*(cosfi+omc*rk*rk)-z[i]*sinfi*rh;
dec_new[i].z=-x[i]*(sinfi*rk)+y[i]*sinfi*rh+z[i]*cosfi;
cout << "dec_new[i].z=" <<dec_new[i].z <<"\n";
//cout << "z[i]=" <<z[i] <<"\n";
//printf( "%i %i %i \n",h[i],k[i],l[i] );
}

red = Fsum_need / Fsum;
//db max_radius = (F_structure_critery*red)*(F_structure_critery*red)*M_PI;


db R[N];
for(int i=0;i<s;i++)
{
R[i]=sqrt(F_structure[i]/M_PI * red);

if(R[i] < F_structure_critery) continue;
if(h[i]==0&&k[i]==0&&l[i]==0)R[i]=R[i]/7;
cout<<"R= "<<R[i]<<endl;

}
//Позволяет разделить все индексы на scale_index[n] и отобразить в PS  все правильно с дробями
int temphkl[3];
char sign;
for(int i=0;i<s;i++)
{
if(R[i] < F_structure_critery) continue;
temphkl[0]=h[i];temphkl[1]=k[i];temphkl[2]=l[i];
for(int n=0;n<3;n++){
if( ( temphkl[n]/scale_index[n] )*scale_index[n]==temphkl[n]){

if(temphkl[n]<0) out<<"(n) ";
else out<<"(p) ";

out<<"("<<abs(temphkl[n])/scale_index[n]<<") ";




}
else{
int max_div;
//Определение максимального делителя для индекса и  для числа scale_index уменьшающего индексы
for(int j=1;j<=macro_size;j++)
if((temphkl[n]/j)*j==temphkl[n]&&(scale_index[n]/j)*j==scale_index[n]) max_div=j;
//

if(temphkl[n]<0) out<<"(n) ";
else out<<"(p) ";

out<<"("<<abs(temphkl[n])/max_div<<"/"<<scale_index[n]/max_div<<") ";
}

}








out<<dec_new[i].x*zoom<<" "\
<<dec_new[i].y*zoom<<" "<<R[i]<<" reflex"<<endl<<endl;
}
int frametype = 2;
//Draw frame arround pattern
out<<"0.5 setlinewidth "<<endl;
if(frametype == 1) //circle
out<<"0 0 "<<r_of_electrongram<<" 0 360 arc stroke"<<endl;

if(frametype == 2){ //square
db r_of_frame = boundingbox[2]/n_of_coloumns/2;
db r = r_of_frame - (r_of_frame / 20);
cout <<r<<" - r of frame (mm)!\n";
cout <<r/zoom*2<<" - side of the square in reciprocical space (Bohr^-1)!\n";
db nr = 0 - r;
out<<"/r "<<r<<" def"<<endl;
out<<"/nr "<<nr<<" def"<<endl;
out<<"gsave\n\
newpath\n\
nr nr moveto\n\
/nr 2 nr mul def\n\
/r 2 r mul def\n\
r 0 rlineto 0 r rlineto nr 0 rlineto 0 nr rlineto closepath stroke %draw box\n\
grestore"<<endl;





}

//out<<"newpath\n"<<-r_of_electrongram/4<<" "<<-r_of_electrongram-8<<" moveto "<<endl; 
//out<<"("<<names_of_structures[i_structure] <<")"<<" show fill"<<endl;
//out<<r_of_electrongram/4<<" "<<r_of_electrongram+8<<" translate "<<endl; 


//out << dec_new[i].x<<" "<<dec_new[i].y<<" "<<floor(dec_new[i].z*1000)/1000 <<"\n";
//cout << dec_new[i].x<<" "<<dec_new[i].y<< dec_new[i].z<<"\n";
/*for(int i=0;i<s;i++)
{
out << x[i]<<" "<<y[i]<<" "<<z[i]<<" "<<"\n";
//cout << x[i]<<y[i]<<z[i]<<"\n";
}*/



//fprintf(fp, "%e %e %e \n", x[i],y[i],z[i]);
//plane_dec(x[s],y[s],z[s])==1


}































void read_input(char *name_of_config_file, char *name_of_srtucture_file, db *ind){

//считываем основные параметры из конфигурационного файла diffraction.in
ifstream in; // input
in.open(name_of_config_file, ios::in);
in >>ni>>nj>>nk;
printf("ni=%i nj=%i nk=%i\n",ni,nj,nk);
in>>cartesian_plane;
cout<<"cartesian_plane="<<cartesian_plane<<endl;
in>>ewald_thickness;
cout<<"ewald_thickness="<<ewald_thickness<<endl;//Толщина слоя в обратных единицах, что и вектора прямой решетки (Бор-1)
in>>F_structure_critery;

cout<<"F_structure_critery="<<F_structure_critery<<endl;//(мм)вводится радиус, который перводиться в структурный критерий, рефлексы, размер которых меньше этого радиуса отображаться не будут
in>>zoom;
cout<<"zoom="<<zoom<<endl;//увеличение масштаба в zoom раз
in>>Fsum_need ;
cout<<"Fsum_need ="<<Fsum_need <<endl;//Необходимая площадь всех рефлексов (в единицах структурного фактора)
//F_structure_critery=(F_structure_critery*red)*(F_structure_critery*red)*M_PI;

in>>r_of_electrongram;
cout<<"r_of_electrongram="<<r_of_electrongram<<endl;//примерно 85(мм) круглая область отсекающая рефлексы
//r_of_electrongram=r_of_electrongram/zoom;//для перевод в размерность обратного пространства разделить на zoom
in>>scale_index[0]>>scale_index[1]>>scale_index[2];//регулирует размер индексов
cout<<"scale_index 1 2 3= "<<scale_index[0]<<" "<<scale_index[1]<<" "<<scale_index[2]<<endl;
in>>reflex_fontsize;
cout<<"reflex_fontsize= "<<reflex_fontsize<<endl;
in.close();

N_of_arrays==(ni+1)*(nj+1)*(nk+1)*8;

cout<<"\n\n\n";


if(cartesian_plane==1){
n_decart[0]=ind[0];
n_decart[1]=ind[1];
n_decart[2]=ind[2];}
if(cartesian_plane==0){
n[0]=(int)ind[0];
n[1]=(int)ind[1];
n[2]=(int)ind[2];
}
if(cartesian_plane==2){
UVW[0]=(int)ind[0];
UVW[1]=(int)ind[1];
UVW[2]=(int)ind[2];
}

//ifstream in; // input
in.open(name_of_srtucture_file, ios::in);
cout<<"Real vectors:\n";
for(int i=0;i<=2; i++)
{in >> a[i][0] >> a[i][1] >> a[i][2] ;
cout << a[i][0]<<" " << a[i][1]<<" " << a[i][2] <<" \n";
//printf("%.30f %.30f %.30f\n",a[i][0],a[i][1],a[i][2]);
}
cout<<"Number of atoms:\n";
in >> nbasis;
cout << nbasis<<"\n";

cout<<"Reduced coordinates of atoms in rprimd(from abinit) :\n";
for(int i=0;i<nbasis; i++){
in >> x_basis[i]>>y_basis[i]>>z_basis[i];
//cout << x_basis[i]<<" " << y_basis[i]<<" " << z_basis[i] <<" \n";
basis[i].x=x_basis[i];
basis[i].y=y_basis[i];
basis[i].z=z_basis[i];}
cout<<"Types of atoms, note!!! 0-carbon, 1-Titanium:\n";
for(int i=0;i<nbasis; i++){
in >> typat[i];//0-углерод, 1-титан
//cout << "typat ="<<typat[i]<< " ";
}
cout << "\n";
in.close();


in.open("s_c_ti.in", ios::in);
for(int i=0;i<24; i++)
{in >> s_deltaK[i];
}
for(int i=0;i<24; i++)
{in >> F[0][i];
}
for(int i=0;i<24; i++)
{in >> F[1][i];

//cout << s_deltaK[i] <<" "<<F[0][i]<<" "<<F[1][i] <<" \n";//0-углерод, 1-титан
}
in.close();

//Трансформируем вектора прямой решетки в вектора обратной решетки!

vector b1,b2,b3;
vector a1(a[0][0],a[0][1],a[0][2]);
vector a2(a[1][0],a[1][1],a[1][2]);
vector a3(a[2][0],a[2][1],a[2][2]);
double Vcell;
//vec_mul=(a2%a3);
Vcell=a1*(a2%a3);
b1=a2%a3;
b2=a3%a1;
b3=a1%a2;
b1=b1/Vcell;
b2=b2/Vcell;
b3=b3/Vcell;
cout<<"Reciprocal vectors b1,b2,b3:";
b1.print();b2.print();b3.print();cout<<endl;
b[0][0]=b1.x;b[0][1]=b1.y;b[0][2]=b1.z;
b[1][0]=b2.x;b[1][1]=b2.y;b[1][2]=b2.z;
b[2][0]=b3.x;b[2][1]=b3.y;b[2][2]=b3.z;

}










double calculate_structure_factor(double K_length,vector_int recip){
double mul;
compln Fhkl,im_unit(0,1);
double b_coef[24],c_coef[24],d_coef[24],F_temp[24];
double F_atomic,dF_atomic;
//table data are for 200 KeV, we introduce energy_factor
double electron_energy=200.;//energy in KeV of electrons
double energy_factor=(1.+electron_energy/511.)/(1.+200./511.);
Fhkl=(0,0);
//Внимание К_length здесь в Бор-1, если вектора заданы в Бор; в таблице дана величина S=dk/(4*pi) в A-1, где dk - один из векторов 
//обратной решетки с множителем 2*pi(физический выбор). в моей программе вектора обратной решетки 
//без множителя 2*pi(кристаллографический выбор)
//т.е. чтобы получить правильный атомный фактор нужно перевести вектор след. образом: К_length=К_length/2/0.529177
//
K_length=K_length/2./0.529177;

for(int i=0;i<nbasis;i++){
//cout << "b_coef= "<<b_coef[0] <<"\n";


spline_(num_of_data,s_deltaK,F[typat[i]],b_coef,c_coef,d_coef);
sevfp_(num_of_data,K_length,s_deltaK,F[typat[i]],b_coef,c_coef,d_coef,F_atomic,dF_atomic);
//cout << "b_coef= "<<b_coef[0] <<"\n";
//cout<<"typat[i]="<<typat[i]<<", K_length(s)="<<K_length<<", F_atomic="<<F_atomic<<"\n";
F_atomic=F_atomic*energy_factor;
//cout << "F_atomic= "<<F_atomic <<"\n";
mul=basis[i].x*recip.h+basis[i].y*recip.k+basis[i].z*recip.l;
//cout << "mul= " << mul << "\n";
//cout << "basis[i].x= " << basis[i].x << "\n";
//cout << im_unit << "\n";
Fhkl+=F_atomic*exp(2*M_PI*im_unit*mul);
//cout << "Fhkl="<< Fhkl << "\n";
}
//cout << "Fhkl[s]^2="<< Fhkl[s] << "\n";
Fhkl=Fhkl.real()*Fhkl.real()+Fhkl.imag()*Fhkl.imag();
//cout << "Fhkl="<< Fhkl << "\n";
//cout << "Fhkl[s]^2="<< Fhkl[s] << "\n";
return(Fhkl.real());

}

















//Функция 1
void convert_b_to_integer(void)
{
int n_el=9;//частный случай 
int k=0;
double b_tmp[n_el];
bool all_integer=0;
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
{b_tmp[k]=b[i][j];
k++;}


/*//2делаем целыми
while(all_integer!=1)
{all_integer=1;
for(int i=0;i<n_el;i++)
if(fmod(b_tmp[i],1)!=0){
for(int i=0;i<n_el;i++){ b_tmp[i]=b_tmp[i]*10;
//printf("%.30f \n",b_tmp[i]);
}
all_integer=0;}
}*/
for(int i=0;i<n_el;i++)
b_tmp[i]=b_tmp[i]*1e7;


long int b_tmp_int[n_el];
for(int i=0;i<n_el;i++)
{b_tmp_int[i]=b_tmp[i];
cout<<b_tmp_int[i]<<" ";}
cout<<"b converted to integer"<<endl;


//3.Находим наибольший общий множитель
int max_mul,min_val;
for(int i=0;i<n_el;i++){
if(b_tmp_int[i]!=0)
min_val=abs(b_tmp_int[0]);}
for(int i=0;i<n_el;i++){
if(b_tmp_int[i]!=0)
abs(b_tmp_int[i])<min_val?min_val=abs(b_tmp_int[i]):i=i;}
cout<<"min_val= " <<min_val<<"\n";

all_integer=0;
for(int k=1;k<=min_val;k++){
for(int i=0;i<n_el;i++){
if(b_tmp_int[i]%k!=0){all_integer=0;break; }
all_integer=1;}
all_integer==1?max_mul=k:k=k;}

cout<<"max_mul= " <<max_mul<<"\n";
for(int i=0;i<n_el;i++)
{b_tmp_int[i]=b_tmp_int[i]/max_mul;
cout<<b_tmp_int[i]<<" ";}
cout<<" Checking output of b_int in convert_b_to_integer(void) function\n"<<endl;

//4 Заполняем массив b_int[i][j];
k=0;
for(int i=0;i<3;i++)
for(int j=0;j<3;j++){
b_int[i][j]=b_tmp_int[k];
k++;}
}











void convert_normal(){
if(cartesian_plane==0){
for(int i=0;i<3;i++)
n_decart[i]=b[0][i]*n[0]+b[1][i]*n[1]+b[2][i]*n[2];
for(int i=0;i<3;i++)
cout << n[i]<<" ";
cout<<"-normal n"<<endl;
for(int i=0;i<3;i++)
cout << n_decart[i]<<" ";
cout<<"-normal n_decart"<<endl;}

if(cartesian_plane==1){
for(int i=0;i<3;i++)
n[i]=0;
calculate_hkl();
for(int i=0;i<3;i++)
cout << n[i]<<" ";
cout<<"-normal n"<<endl;
for(int i=0;i<3;i++)
cout << n_decart[i]<<" ";
cout<<"-normal n_decart"<<endl;}
//for(int i=0;i<3;i++)
//for(int j=0;j<3;j++)
//cout<<b_int[i][j]<<" ";
//cout<<" -b_int"<<endl;//в некоторых случаях удачно нормализованные вектора обратной решетки.


if(cartesian_plane==2){
for(int i=0;i<3;i++)
n_decart[i]=a[0][i]*UVW[0]+a[1][i]*UVW[1]+a[2][i]*UVW[2];

for(int i=0;i<3;i++)
cout << UVW[i]<<" ";
cout<<"-normal UVW"<<endl;
for(int i=0;i<3;i++)
cout << n[i]<<" ";
cout<<"-normal n"<<endl;
for(int i=0;i<3;i++)
cout << n_decart[i]<<" ";
cout<<"-normal n_decart"<<endl;}





}





void calculate_hkl(){

db a1=-b[0][1]/1./b[0][0];
if(b[0][1]==0)a1=0;
db c=-b[0][2]/1./b[0][1];
if(b[0][2]==0)c=0;
db a11=a1*b[1][0]+b[1][1];
db a12=a1*b[2][0]+b[2][1];
db c11=c*b[1][1]+b[1][2];
db c12=c*b[2][1]+b[2][2];
db A=-c11/a11;
if(c11==0)A=0;
db n12=a1*n_decart[0]+n_decart[1];
db n23=c*n_decart[1]+n_decart[2];
db n_db[3];
 n_db[2]=(n12*A+n23)/(A*a12+c12);
 n_db[1]=(n12-n_db[2]*a12)/a11;
 n_db[0]=(n_decart[0]-n_db[1]*b[1][0]-n_db[2]*b[2][0])/b[0][0];
cout<<"n_db= "<<n_db[0]<<" "<<n_db[1]<<" "<<n_db[2];
cout<<" normal n,double"<<endl;
for(int i=0;i<3;i++)
n[i]=n_db[i]*10;
}











int plane(int x,int y,int z){
if((n_decart[0]*(x-p[0])+n_decart[1]*(y-p[1])+n_decart[2]*(z-p[2]))==0)return 1;
return 0;
}


int plane2(db x,db y,db z){
//cout<<"fabs(n_decart[0]*x+n_decart[1]*y+n_decart[2]*z)= "<<fabs(n_decart[0]*x+n_decart[1]*y+n_decart[2]*z)<<" YO \n";

if(fabs(n_decart[0]*x+n_decart[1]*y+n_decart[2]*z)<=ewald_thickness){

return 1;

}
else
return 0;
}










#include "diffraction.h"

extern "C" {void spline_(int&,double*,double*,double*,double*,double*);}
extern "C" {void sevfp_(int&,double&,double*,double*,double*,double*,double*,double&,double&);}

//Global var
bool g_debug = 1;

//Используется кристаллографическое определение векторов обратной решетки.
main(int argc, char *argv[]) {
    
    if (argc != 2){ cout << "Not enouth input data" << endl; return 0; }
    
    string layout_filename = string(argv[1]);
    string ps_filename = "output/" + layout_filename + ".eps"; //output filename

    cout << "\nProgram " << argv[0] << " have been started ...\n";
    cout << "Name of layout file = " << layout_filename << endl;
    cout << "Name of output file = " << ps_filename     << endl;


    db boundingbox[4] = {0, 0, 210, 150};//
    int ic = 0, ir = 0;
    db reflex_fontsize; //size of font for reflex names
    LayoutClass layout;

    //0. Read layout data
    layout.readin(argv[1]); // Считываем имена файлов со структурами и имена кофигурационных файлов к ним из файла списка


    //Make layout and calculate diffraction patterns
    cout << "\nStart main cycle ...\n  ";
    
    ofstream out;
    out.open(ps_filename.c_str(), ios::out); // открываем файл для записи
    
    for (int i_st = 0; i_st < layout.number_of_stuctures; i_st++) {        
        cout << "Name of structure   " << layout.struct_label[i_st] << endl;


        //1. Read configuration data and crystal structure
        SADPClass SADP;
        SADP.read_config(layout.config_filename[i_st], layout.uvw[i_st]); //Use reflex_fontsize from first config file for all structures
        SADP.read_crystal_structure(layout.struct_filename[i_st]);
        SADP.read_atomic_factors("s_c_ti.in");

        //2. Calculate selected area diffraction patterns (SADP)
        SADP.construct_reiprocal_lattice();
        SADP.rotate();
        SADP.make_SADP_frame(boundingbox, layout.n_of_coloumns);

        // cout << SADP.frame_string;
        // 3. Form postscript output
        if (i_st == 0) 
            out << postscript_header(boundingbox,  SADP.reflex_fontsize      ); //Выводим postscript header
        out     << postscript_layout(i_st, ic, ir, SADP, layout, boundingbox );
        

        ic++; if (ic >= layout.n_of_coloumns) {ic = 0; ir++; } // determine numbers of columns and rows

        // cout << "-----------------------------------------------------------------------------\n\n\n\n\n";
    }

    out << "showpage % вывести страницу\n%%EOF\n";
    
    out.close();
    
}





void LayoutClass::readin(char* layout_filename) {


    ifstream in; 

    cout << "Start reading layout file ..." << endl;

    in.open(layout_filename, ios::in); // open file

    in >> number_of_stuctures;
    in >> n_of_coloumns;
    in >> n_of_rows;

    cout << "number_of_stuctures = " << number_of_stuctures << endl;
    cout << "n_of_coloumns       = " << n_of_coloumns       << endl;
    cout << "n_of_rows           = " << n_of_rows           << endl;

    for (int i = 0; i < number_of_stuctures; i++) {
        in >> uvw[i][0];
        in >> uvw[i][1];
        in >> uvw[i][2];
        in >> struct_filename[i];
        in >> config_filename[i];
        in >> struct_label[i];
        
        cout << struct_filename[i] << endl;
    }
    

    in.close();


}



db SADPClass::read_config(string name_of_config_file, db *ind) {

    ni = MAX_LATTICE_SIZE; nj = MAX_LATTICE_SIZE; nk = MAX_LATTICE_SIZE;

    ifstream in;
    
    cout << "Starting to read config file " << name_of_config_file << endl;    

    in.open(name_of_config_file.c_str(), ios::in);
    
    in >> ni >> nj >> nk;
    in>>cartesian_plane;
    in>>ewald_thickness;
    in>>F_structure_critery;
    in>>zoom;
    in>>Fsum_need ;
    in>>r_of_electrongram;
    in>>scale_index[0]>>scale_index[1]>>scale_index[2];
    in>>reflex_fontsize;

    printf("ni=%i nj=%i nk=%i\n",ni,nj,nk);
    cout<<"cartesian_plane="<<cartesian_plane<<endl;
    cout<<"ewald_thickness="<<ewald_thickness<<endl;//Толщина слоя в обратных единицах, что и вектора прямой решетки (Бор-1)
    cout<<"F_structure_critery="<<F_structure_critery<<endl;//(мм)вводится радиус, который перводиться в структурный критерий, рефлексы, размер которых меньше этого радиуса отображаться не будут
    cout<<"zoom="<<zoom<<endl;//увеличение масштаба в zoom раз
    cout<<"Fsum_need ="<<Fsum_need <<endl;//Необходимая площадь всех рефлексов (в единицах структурного фактора)
    //F_structure_critery=(F_structure_critery*red)*(F_structure_critery*red)*M_PI;
    cout<<"r_of_electrongram="<<r_of_electrongram<<endl;//примерно 85(мм) круглая область отсекающая рефлексы
    //r_of_electrongram=r_of_electrongram/zoom;//для перевод в размерность обратного пространства разделить на zoom
    cout<<"scale_index 1 2 3= "<<scale_index[0]<<" "<<scale_index[1]<<" "<<scale_index[2]<<endl; //регулирует размер индексов
    cout<<"reflex_fontsize= "<<reflex_fontsize<<endl;

    in.close();

    // N_of_arrays = (ni+1) * (nj+1) * (nk+1) * 8;

    cout<<"\n\n\n";


    if (cartesian_plane == 0) {
        n[0]=(int)ind[0];
        n[1]=(int)ind[1];
        n[2]=(int)ind[2];
    }

    if (cartesian_plane == 1) {
        n_decart[0]=ind[0];
        n_decart[1]=ind[1];
        n_decart[2]=ind[2];
    }
    
    
    if (cartesian_plane == 2) {
        UVW[0]=(int)ind[0];
        UVW[1]=(int)ind[1];
        UVW[2]=(int)ind[2];
    }

    return(reflex_fontsize);
}


void SADPClass::read_crystal_structure (string name_of_srtucture_file) {
    
    cout << "Starting to read crystal structure ... " ;
    ifstream in;
    
    in.open(name_of_srtucture_file.c_str(), ios::in);
    
    cout<<"Real vectors:\n";
    
    for(int i=0;i<=2; i++) {
        in >> a[i][0] >> a[i][1] >> a[i][2];
        cout << a[i][0]<<" " << a[i][1]<<" " << a[i][2] <<" \n";
        //printf("%.30f %.30f %.30f\n",a[i][0],a[i][1],a[i][2]);
    }
    
    cout << "Number of atoms:\n";
    in >> nbasis;
    cout << nbasis<<"\n";

    cout << "Reduced coordinates of atoms in rprimd(from abinit) :\n";
    for(int i=0;i<nbasis; i++) {
        in >> x_basis[i]>>y_basis[i]>>z_basis[i];
        //cout << x_basis[i]<<" " << y_basis[i]<<" " << z_basis[i] <<" \n";
        basis[i].x=x_basis[i];
        basis[i].y=y_basis[i];
        basis[i].z=z_basis[i];
    }
    cout<<"Types of atoms, note!!! 0-carbon, 1-Titanium:\n";
    for(int i=0;i<nbasis; i++) {
        in >> typat[i];//0-углерод, 1-титан
        //cout << "typat ="<<typat[i]<< " ";
    }

    cout << "\n";

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




void SADPClass::read_atomic_factors(string atomic_factors_filename) {

    cout << "Starting to read atomic factors ... " ;

    ifstream in;
    num_of_data = NUM_ATOMIC_NUMBERS;

    in.open(atomic_factors_filename.c_str(), ios::in);
    
    for(int i = 0; i < num_of_data; i++)
        in >> s_deltaK[i];
    
    for(int i = 0; i < num_of_data; i++)
        in >> F[0][i];
    
    for(int i = 0;i < num_of_data; i++)
        in >> F[1][i];

    //cout << s_deltaK[i] <<" "<<F[0][i]<<" "<<F[1][i] <<" \n";//0-углерод, 1-титан

    in.close();

}





void SADPClass::construct_reiprocal_lattice() {
    
    cout << "Constructing reciprocal lattice ... \n" ;


    convert_normal();

    //Преобразуем толщину слоя в параметр dt*(n_decart[0]+n_decart[1]+n_decart[2])
    //dt-параметр, 2*dt*sqrt(3)-толщина плоскости, когда условие в виде:
    //fabs(n_decart[0]*(x-p[0])+n_decart[1]*(y-p[1])+n_decart[2]*(z-p[2]))<=dt*(n_decart[0]+n_decart[1]+n_decart[2])
    //это позволяет уйти от зависимости когда вектор нормали просто увеличен или уменьшен
    db n_decart_powered = n_decart[0] * n_decart[0] + n_decart[1] * n_decart[1] + n_decart[2] * n_decart[2];
    n_decart_length  = sqrt(n_decart_powered);
    ewald_thickness     = 0.5 * ewald_thickness * n_decart_powered / n_decart_length;

    cout << "ewald_thickness after transformation " << ewald_thickness << endl;


    //cout <<"MAX_NUM_REFLEX =" <<MAX_NUM_REFLEX;
    vector_int tmp, recip;
    vector_dec tmp2;
    db d, w10;
    Fsum = 0;
    int i, j, ki, s=0;
    for(i=-ni;i<ni+1;i++)
        for(j=-nj;j<nj+1;j++)
            for(ki=-nk;ki<nk+1;ki++) {
                //cout << "s= "<<s << "\n";

                x[s]=i*b[0][0]+j*b[1][0]+ki*b[2][0];
                y[s]=i*b[0][1]+j*b[1][1]+ki*b[2][1];
                z[s]=i*b[0][2]+j*b[1][2]+ki*b[2][2];

                if (plane2(x[s], y[s], z[s], ewald_thickness, n_decart) == 1) {
                    h[s]    = i;
                    k[s]    = j;
                    l[s]    = ki;
                    recip.h = i;
                    recip.k = j;
                    recip.l = ki;
                    //printf( "%i %i %i \n",h[s],k[s],l[s] );
                    d = sqrt(x[s]*x[s]+y[s]*y[s]+z[s]*z[s]);//верно, только когда p[3]={0,0,0}
                    // cout << "d = "<< d << "\n";
                    F_structure[s] = calculate_structure_factor(d, recip);//уже возведенный в квадрат структурный фактор по модулю!!!
                    //cout << "F_structure[s]="<< F_structure[s] << "\n";
                    //cout << "K_length="<< K_length << "\n";
                    //cout << "recip[s].h="<< recip[s].h<< "\n";
                    if( F_structure[s] >= F_structure_critery && d <= r_of_electrongram / zoom ) { //доп парметр-радиус, делает электроннограмму круглой, обрезая рефлексы
                        w10 = n_decart[0] * x[s] + n_decart[1] * y[s] + n_decart[2] * z[s];
                        cout << "F_structure[s]="<< F_structure[s] << ", distance from plane = "<<w10<<"\n";
                        Fsum += F_structure[s];
                        s++;
                    }
                    else {
                        // if (g_debug == 1) cout << "F_structure " << F_structure[s] << " is smaller than F_structure_critery; or d larger than r_of_electrongram /zoom "<<r_of_electrongram / zoom <<"\n";

                    }


                }
                
                if (s >= MAX_NUM_REFLEX) cout << "Error, s>=MAX_NUM_REFLEX, the size of arrays are two small";
                
            }

    N_nodes = s;
};


void SADPClass::rotate() {

    cout << "Rotation of lattice... \n" ;

    //Поворачиваем систему координат, чтобы получить координаты точек в плоскости, происходит совмещение нормали с осью Z.
    //находим угол между осью z и n_decart.
    double fi = acos(n_decart[2] / n_decart_length);
    double cosfi = cos(fi);
    double sinfi = sin(fi);
    double omc   = 1 - cosfi;
    vector_dec rotate_axis;

    //cout << "n_decart[2]= " << n_decart[2]<<"\n";
    //cout << "n_decart_length= " << n_decart_length<<"\n";
    cout << "fi= " << fi * 180 / M_PI << endl;


    rotate_axis.x =  n_decart[1];
    rotate_axis.y = -n_decart[0];
    rotate_axis.z =  0;

    db rh = rotate_axis.x / rotate_axis.length();
    db rk = rotate_axis.y / rotate_axis.length();
    db rl = 0;

    if(rotate_axis.length() == 0) { rh = 0; rk = 0; }

    cout << "rh,rk= " << rh<<" "<<rk<<"\n";
    cout << "rotate_axis.length= " << rotate_axis.length()<<"\n";
    
    for(int i = 0; i < N_nodes; i++) {
        dec_new[i].x=x[i]*(cosfi+omc*rh*rh)+y[i]*omc*rh*rk+z[i]*sinfi*rk;
        dec_new[i].y=x[i]*omc*rh*rk+y[i]*(cosfi+omc*rk*rk)-z[i]*sinfi*rh;
        dec_new[i].z=-x[i]*(sinfi*rk)+y[i]*sinfi*rh+z[i]*cosfi;
        cout << "dec_new[i].z=" <<dec_new[i].z <<"\n";
        //cout << "z[i]=" <<z[i] <<"\n";
        //printf( "%i %i %i \n",h[i],k[i],l[i] );
    }

}


void SADPClass::make_SADP_frame(db *boundingbox, int n_of_coloumns) {

    cout << "Making SADP frame ... \n" ;


    ostringstream out;
    //1. Reduce size of reflexes
    db red = Fsum_need / Fsum; // reduction value
    db R[MAX_NUM_REFLEX];
    for(int i = 0; i < N_nodes; i++) {
        R[i]=sqrt(F_structure[i]/M_PI * red);

        if (R[i] < F_structure_critery) continue;
        if (h[i] == 0 && k[i] == 0 && l[i] == 0)  R[i] = R[i] / 7; //Reduce central reflex
        
        cout << "R = " << R[i] << endl;

    }
    

    //2. Main cycle over reflexes
    int temphkl[3];
    char sign;
    
    for (int i = 0; i < N_nodes; i++) {
        
        if (R[i] < F_structure_critery) continue;
        
        //! Разделить все индексы на scale_index[n] и отобразить в PS все правильно с дробями

        temphkl[0] = h[i]; temphkl[1] = k[i]; temphkl[2] = l[i];
        
        for (int n=0;n<3;n++) {
            
            if ( ( temphkl[n]/scale_index[n] )*scale_index[n]==temphkl[n]){

                if (temphkl[n] < 0) out << "(n) ";
                
                else out << "(p) ";

                out << "(" << abs(temphkl[n]) / scale_index[n] <<") ";
            }
            else {
                int max_div;
                //Определение максимального делителя для индекса и  для числа scale_index уменьшающего индексы
                for (int j=1;j<=MAX_LATTICE_SIZE;j++)
                    if ((temphkl[n]     / j) * j == temphkl[n]              && \
                        (scale_index[n] / j) * j == scale_index[n]) max_div = j;

                if (temphkl[n] < 0) out <<  "(n) ";
                else out << "(p) ";

                out << "(" << abs(temphkl[n]) / max_div << "/" << scale_index[n] / max_div << ") ";
            }

        }

    //Write scaled coordinates and intensities of reflexes
        out << dec_new[i].x * zoom << " " << dec_new[i].y * zoom<<" "<< R[i] <<" reflex" << endl << endl;
    }
    
 
    // 4. Draw frame arround pattern
    int frametype = 2; // Choose type of frame; 1 - circle; 2 - rectangular
    
    out << "0.5 setlinewidth " << endl;
    
    if(frametype == 1)   //circle
        out<<"0 0 "<<r_of_electrongram<<" 0 360 arc stroke"<<endl;

    if(frametype == 2) { // rectangular
        db r_of_frame = boundingbox[2] / n_of_coloumns / 2;
        db r          = r_of_frame - (r_of_frame / 20);
        
        cout << r <<            " - r of frame (mm)!\n";
        cout << r / zoom * 2 << " - side of the square in reciprocical space (Bohr^-1)!\n";
        
        db nr = 0 - r;
        
        out << "/r "  << r  << " def\n";
        out << "/nr " << nr << " def\n";
        
        out << \
"gsave\n\
newpath\n\
nr nr moveto\n\
/nr 2 nr mul def\n\
/r 2 r mul def\n\
r 0 rlineto 0 r rlineto nr 0 rlineto 0 nr rlineto closepath stroke %draw box\n\
grestore\n";

    }


    cout << "Frame is constructed!\n";
    // cout << out.str();
    frame_string = out.str();
    // return (out.str());

}





//Функции:


db SADPClass::calculate_structure_factor(double K_length, vector_int recip) {

    // if (g_debug == 1) cout << "\ncalculate_structure_factor() started ...\n\n";
    int num_of_data = NUM_ATOMIC_NUMBERS;
    double mul;
    compln Fhkl,im_unit(0,1);
    double b_coef[NUM_ATOMIC_NUMBERS],c_coef[NUM_ATOMIC_NUMBERS],d_coef[NUM_ATOMIC_NUMBERS],F_temp[NUM_ATOMIC_NUMBERS];
    double F_atomic, dF_atomic;
    //table data are for 200 KeV, we introduce energy_factor
    double electron_energy=200.;//energy in KeV of electrons
    double energy_factor=(1.+electron_energy/511.)/(1.+200./511.);
    Fhkl=(0,0);
    //Внимание К_length здесь в Бор-1, если вектора заданы в Бор; в таблице дана величина S=dk/(4*pi) в A-1, где dk - один из векторов 
    //обратной решетки с множителем 2*pi(физический выбор). в моей программе вектора обратной решетки 
    //без множителя 2*pi(кристаллографический выбор)
    //т.е. чтобы получить правильный атомный фактор нужно перевести вектор след. образом: К_length=К_length/2/0.529177
    //
    K_length = K_length / 2. / 0.529177;

    for(int i=0; i < nbasis; i++) {
        //cout << "b_coef= "<<b_coef[0] <<"\n";

        spline_(num_of_data,           s_deltaK, F[typat[i]], b_coef, c_coef, d_coef);
        sevfp_( num_of_data, K_length, s_deltaK, F[typat[i]], b_coef, c_coef, d_coef, F_atomic, dF_atomic);
        //cout << "b_coef= "<<b_coef[0] <<"\n";
        //cout<<"typat[i]="<<typat[i]<<", K_length(s)="<<K_length<<", F_atomic="<<F_atomic<<"\n";
        F_atomic = F_atomic * energy_factor;
        //cout << "F_atomic= "<<F_atomic <<"\n";
        mul = basis[i].x * recip.h + basis[i].y * recip.k + basis[i].z * recip.l;
        //cout << "mul= " << mul << "\n";
        //cout << "basis[i].x= " << basis[i].x << "\n";
        //cout << im_unit << "\n";
        Fhkl += F_atomic * exp(2 * M_PI * im_unit * mul);
        //cout << "Fhkl="<< Fhkl << "\n";
        }

    //cout << "Fhkl[s]^2="<< Fhkl[s] << "\n";
    Fhkl = Fhkl.real() * Fhkl.real() + Fhkl.imag() * Fhkl.imag();
    //cout << "Fhkl="<< Fhkl << "\n";
    //cout << "Fhkl[s]^2="<< Fhkl[s] << "\n";
    return(Fhkl.real());

}





void SADPClass::convert_normal(){
    //Получаем нормаль в декартовых координатах, если она задана в векторах обратной решетки или
    //Переводим нормаль в декартовых координатах n_decart(соответсвтует кубической обратной решетке) в нормль n в векторах обратной решетки

    cout << "Converting normal ... \n" ;


    if(cartesian_plane == 0) {
        for(int i=0;i<3;i++)
            n_decart[i]=b[0][i]*n[0]+b[1][i]*n[1]+b[2][i]*n[2];
        for(int i=0;i<3;i++)
            cout << n[i]<<" ";
        cout<<"-normal n"<<endl;
        for(int i=0;i<3;i++)
            cout << n_decart[i]<<" ";
        cout<<"-normal n_decart"<<endl;
    }

    if(cartesian_plane == 1) {
        for(int i=0;i<3;i++)
            n[i]=0;
        calculate_hkl();
        for(int i=0;i<3;i++)
            cout << n[i]<<" ";
        cout<<"-normal n"<<endl;
        for(int i=0;i<3;i++)
            cout << n_decart[i]<<" ";
        cout<<"-normal n_decart"<<endl;
    }

    if(cartesian_plane == 2) {
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
            cout<<"-normal n_decart"<<endl;
    }

}





void SADPClass::calculate_hkl() {

    db a1=-b[0][1]/1./b[0][0];
    
    if(b[0][1]==0) a1=0;
    
    db c=-b[0][2]/1./b[0][1];
    
    if(b[0][2]==0) c=0;
    
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















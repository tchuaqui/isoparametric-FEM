template <typename COORD, typename LINKS> class Mesh{
double Lx, Ly;
int nelements_x, nelements_y;
public: 
typedef Eigen::Matrix<double,4,2> COORD_ELE;
Mesh(double x,double y,int nex,int ney);
double dx();
double dy();
COORD coord();
LINKS links();
COORD_ELE coord_ele(int e);
};

template <typename COORD, typename LINKS> Mesh<COORD,LINKS>::Mesh(double x, double y,int nex,int ney) : Lx(x), Ly(y), nelements_x(nex), nelements_y(ney){};
template <typename COORD, typename LINKS> double Mesh<COORD,LINKS>::dx(){return Lx/nelements_x;};
template <typename COORD, typename LINKS> double Mesh<COORD,LINKS>::dy(){return Ly/nelements_y;};
template <typename COORD, typename LINKS> COORD Mesh<COORD,LINKS>::coord(){
COORD C;
 for (int i=0; i<(nelements_y+1); ++i){
  for (int j=0; j<(nelements_x+1); ++j){
   C((nelements_x+1)*i+j,0)=j*dx();
   C((nelements_x+1)*i+j,1)=i*dy();
  };
 };
return C;};
template <typename COORD, typename LINKS> LINKS Mesh<COORD,LINKS>::links(){
 LINKS L;
 for (int i=0; i<nelements_y; ++i){
  for (int j=0; j<nelements_x; ++j){
   L(0,j+nelements_x*i)=j+i*(nelements_x+1);
   L(1,j+nelements_x*i)=1+j+i*(nelements_x+1);
   L(2,j+nelements_x*i)=nelements_x+2+j+i*(nelements_x+1);
   L(3,j+nelements_x*i)=nelements_x+1+j+i*(nelements_x+1); 
  };
 };
return L;};
template <typename COORD, typename LINKS> typename Mesh<COORD,LINKS>::COORD_ELE Mesh<COORD,LINKS>::coord_ele(int e){
 COORD C=coord();
 LINKS L=links();
 COORD_ELE CE; 
 for (int i=0;i<4;++i){
  CE(i,0)=C(L(i,e),0);
  CE(i,1)=C(L(i,e),1);
 };
 return CE;};


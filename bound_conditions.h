template<typename BOUND> class Bound{
int nelements_x, nelements_y;
char* b_string;
public:
Bound(int nex,int ney,char* s);
int nb_nodes();
BOUND bound_matrix();
};

template<typename BOUND> Bound<BOUND>::Bound(int nex,int ney,char* s) : nelements_x(nex), nelements_y(ney), b_string(s){};
template<typename BOUND> int Bound<BOUND>::nb_nodes(){return 2*(nelements_x+nelements_y);};
template<typename BOUND> BOUND Bound<BOUND>::bound_matrix(){
 BOUND A;
 for (int i=0; i<(nelements_x+1)*(nelements_y+1); ++i){
 A(i,0)=0;
 A(i,1)=0;
 };
 //side 1
 for (int i=0; i<nelements_x+1; ++i){
  if (b_string[0]=='c'){
   A(i,0)=1;
   A(i,1)=1;} 
  };
 //side 2
 for (int i=0; i<nelements_y+1; ++i){
  if (b_string[1]=='c'){
   A((nelements_x+1)*i+nelements_x,0)=1;
   A((nelements_x+1)*i+nelements_x,1)=1;}
 };
 //side 3
 for (int i=((nelements_x+1)*(nelements_y+1)-(nelements_x+1)); i<((nelements_x+1)*(nelements_y+1)); ++i){
  if (b_string[2]=='c'){
   A(i,0)=1;
   A(i,1)=1;}
 };
 //side 4
 for (int i=0; i<nelements_y+1; ++i){
  if (b_string[3]=='c'){
   A((nelements_x+1)*i,0)=1;
   A((nelements_x+1)*i,1)=1;}
 };
 return A;
};



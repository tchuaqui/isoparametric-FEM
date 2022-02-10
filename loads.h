template<typename LOAD_NODES, typename LOAD_NODES_T, typename F_NODES> class Load{
int nelements_x, nelements_y;
public:
Load(int nex,int ney);
LOAD_NODES load_matrix(double side1,double side2, double side3, double side4);
LOAD_NODES_T load_matrix_t(double side1,double side2, double side3, double side4);
};
template <typename LOAD_NODES, typename LOAD_NODES_T, typename F_NODES> Load<LOAD_NODES,LOAD_NODES_T,F_NODES>::Load(int nex, int ney) : nelements_x(nex), nelements_y(ney){};
template<typename LOAD_NODES, typename LOAD_NODES_T, typename F_NODES> LOAD_NODES Load<LOAD_NODES,LOAD_NODES_T,F_NODES>::load_matrix(double side1,double side2, double side3, double side4){
 LOAD_NODES A;
 for (int i=0; i<(nelements_x+1)*(nelements_y+1); ++i){
 A(i,0)=0;
 A(i,1)=0;
 };
 //side 1
 for (int i=0; i<nelements_x+1; ++i){
   A(i,1)=side1;
 };
 //side 2
 for (int i=0; i<nelements_y+1; ++i){
   A((nelements_x+1)*i+nelements_x,0)=side2;
 };
 //side 3
 for (int i=((nelements_x+1)*(nelements_y+1)-(nelements_x+1)); i<((nelements_x+1)*(nelements_y+1)); ++i){
   A(i,1)=side3;
 };
 //side 4
 for (int i=0; i<nelements_y+1; ++i){
   A((nelements_x+1)*i,0)=side4;
 };
 return A;
}; 
template <typename LOAD_NODES, typename LOAD_NODES_T, typename F_NODES> LOAD_NODES_T Load<LOAD_NODES,LOAD_NODES_T,F_NODES>::load_matrix_t(double side1,double side2,double side3,double side4){
LOAD_NODES M=load_matrix(side1,side2,side3,side4);
LOAD_NODES_T F_matrix=M.transpose();
return F_matrix;
};

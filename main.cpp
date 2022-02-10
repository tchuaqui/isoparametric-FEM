#include<iostream>
#include<fstream>
#include<math.h>
#include"Eigen/Eigen/Dense"
#include"Eigen/Eigen/Sparse"
#include<vector>

#include"mesh_grid.h"
#include"bound_conditions.h"
#include"loads.h"
#include"elasticity.h"
#include"shape_functions.h"
#include"integration.h"
#include"local_operator.h"
#include"assembly.h"

int main(){
double length=12e-2;
double height=0.5e-2;
int const nelements_x=80;
int const nelements_y=10;
//Elastic properties
double E=210.0e9;
double v=0.3;
//Boundary conditions (side numbering of structure)
/*      3
    ---------
   |         |
  4|         |2
   |    1    |
    ---------

*/
char b_string[]={'l','c','l','c'}; //clamps (c) or free (l) each side 
//Number of Gauss points
int ngauss=1;
//Define Loads
double l_el_x=0;  //define load on element in x direction
double l_el_y=0;  //define load on element in y direction
double side1_y=0;     //define load on sides
double side2_x=0;       //define load on sides
double side3_y=-10000; //-10000
double side4_x=0;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int const nnodes_ele=4;
int const ndof_ele=2;
int row,col;

typedef Eigen::Matrix<double,ndof_ele,1> LOAD_ELE; 
LOAD_ELE load_ele; load_ele(0)=l_el_x; load_ele(1)=l_el_y;

typedef Eigen::Matrix<double,(nelements_x+1)*(nelements_y+1),ndof_ele> COORD;
typedef Eigen::Matrix<int,nnodes_ele,(nelements_x*nelements_y)> LINKS;
typedef Eigen::Matrix<int,(nelements_x+1)*(nelements_y+1),2> BOUND;
typedef Eigen::Matrix<double,3,3> ELASTIC;
typedef Eigen::Matrix<double,(nelements_x+1)*(nelements_y+1)*ndof_ele,1> F_NODES;
typedef Eigen::Matrix<double,(nelements_x+1)*(nelements_y+1),ndof_ele> LOAD_NODES;
typedef Eigen::Matrix<double,ndof_ele,(nelements_x+1)*(nelements_y+1)> LOAD_NODES_T;
typedef Eigen::Matrix<int,(nelements_x+1)*(nelements_y+1)*ndof_ele,1> ID_NODES_VEC;
typedef Eigen::Matrix<int,ndof_ele,(nelements_x+1)*(nelements_y+1)> ID_NODES;
typedef Eigen::Matrix<int,ndof_ele*nnodes_ele,nelements_x*nelements_y> ID_ELEMENTS;
typedef Eigen::Matrix<int,ndof_ele*nnodes_ele,nelements_x*nelements_y>::Index MAX_INDEX;
typedef Eigen::Matrix<int,nelements_x*nelements_y,1> DIFF;


Mesh<COORD,LINKS> mesh(length,height,nelements_x,nelements_y);
Bound<BOUND> bound(nelements_x,nelements_y,b_string);
Elastic<ELASTIC> elastic(E,v);
Load<LOAD_NODES,LOAD_NODES_T,F_NODES> load_nodes(nelements_x,nelements_y);
Shape_func shape_func;
Integration integration(ngauss);
Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC> local_operator;
LINKS link=mesh.links();

LOAD_NODES_T Aux=load_nodes.load_matrix_t(side1_y,side2_x,side3_y,side4_x);
F_NODES F_nodes(Eigen::Map<F_NODES>(Aux.data(), Aux.cols()*Aux.rows()));

typedef typename Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::K_LOCAL K_LOCAL;
typedef typename Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::F_LOCAL F_LOCAL;


K_LOCAL K_ele; 
F_LOCAL F_ele;

Assemble<ID_NODES,ID_ELEMENTS,BOUND,COORD,LINKS,MAX_INDEX,DIFF> assemble;
ID_NODES Id_Nodes;
ID_ELEMENTS Id_Elements; 
Id_Elements=assemble.id_elements(nnodes_ele,nelements_x,nelements_y,ndof_ele,mesh,bound);
Id_Nodes=assemble.id_nodes(nelements_x,nelements_y,ndof_ele,bound);
int n_adof=assemble.nadof(nelements_x,nelements_y,ndof_ele,bound); //number of active DOF   
int band_width=assemble.bandwidth(nelements_x,nelements_y,ndof_ele,nnodes_ele,mesh,bound); //Bandwidth - estimation of non zero entries in sparse matrix

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Make global load vector from node loads - F_global with F_nodes
Eigen::SparseVector<double> F_global(n_adof), X_global(n_adof);
Eigen::SparseMatrix<int> Index_Id_nodes(2,n_adof); //stores index of active dof's in 1st row and number of global equation in 2nd row
int count=0;
ID_NODES_VEC Id_nodes_vector(Eigen::Map<ID_NODES_VEC>(Id_Nodes.data(), Id_Nodes.cols()*Id_Nodes.rows())); //transform Id_Nodes matrix into vector form

for (int i=0;i<(nelements_x+1)*(nelements_y+1)*ndof_ele;++i){
 if(Id_nodes_vector(i)>-1){count++; Index_Id_nodes.coeffRef(0,count-1)=i; Index_Id_nodes.coeffRef(1,count-1)=Id_nodes_vector(i);}; 
};
for (int i=0;i<n_adof;++i){
F_global.coeffRef(Index_Id_nodes.coeffRef(1,i))=F_nodes(Index_Id_nodes.coeffRef(0,i));
};
////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//LOOP WITH BOUNDARY CONDITIONS AND SPARSE MATRICES - FORM GLOBAL STIFFNESS AND LOAD VECTOR FROM ELEMENTS
typedef Eigen::Triplet<double> T_K;
std::vector<T_K> t_K_global;
t_K_global.reserve(band_width);
Eigen::SparseMatrix<double> K_global(n_adof,n_adof);

double val_K=0;
double val_F=0;

Eigen::Matrix2Xi dynam_mat; dynam_mat.resize(2,0); //Similar to Index_Id_nodes but acts on element level - changes size for each element so dynamic (1st row is index of active dof in element and 2nd row is corresponding number of global equation)
int dummy;
for (int el=0; el<(nelements_x*nelements_y); ++el){
 dummy=0; dynam_mat.resize(2,0); 
 K_ele=local_operator.K_local(el,ngauss,shape_func,integration,mesh,elastic);
 F_ele=local_operator.F_local(el,ngauss,load_ele,shape_func,integration,mesh);
 for (int i=0;i<nnodes_ele*ndof_ele;++i){
 if (Id_Elements(i,el)>-1){
   dummy++; dynam_mat.conservativeResize(2,dummy); dynam_mat(0,dummy-1)=i; dynam_mat(1,dummy-1)=Id_Elements(i,el);};
 };
 for (int i=0;i<dummy;++i){
  val_F=F_ele(dynam_mat(0,i));    //Loads applied on each element
  F_global.coeffRef(dynam_mat(1,i))=F_global.coeffRef(dynam_mat(1,i))+val_F; //Add element loads to F_global(already includes F_nodes)
  for (int j=0;j<dummy;++j){
  val_K=K_ele(dynam_mat(0,i),dynam_mat(0,j)); 
  t_K_global.push_back(T_K(dynam_mat(1,i),dynam_mat(1,j),val_K));
  };
 };
};
K_global.setFromTriplets(t_K_global.begin(), t_K_global.end());
K_global.makeCompressed();

/////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//SOLVE LINEAR PROBLEM
//Conjugate gradient
//Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
//cg.compute(K_global);
//X_global=cg.solve(F_global);                    //Obtain displacements of active DOF's
//std::cout<<std::endl<<"iterations:      "<<cg.iterations()<<std::endl;
//std::cout<<"estimated error:       "<<cg.error()<<std::endl;

//Direct solver
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldlt;
ldlt.compute(K_global);
X_global=ldlt.solve(F_global);                    //Obtain displacements of active DOF's

////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//Final coordinates of nodes
COORD coordinates=mesh.coord(); LOAD_NODES_T coordinates_transp=coordinates.transpose(); //coordinate matrix of nodes and transpose
F_NODES coordinates_vec(Eigen::Map<F_NODES>(coordinates_transp.data(), coordinates_transp.cols()*coordinates_transp.rows())); //convert coordinate matrix to vector of coordinates for all DOF's
F_NODES final_coordinates_vec=coordinates_vec;    //final vector of coordinates (after displacement) of all DOF's initialization 
COORD coordinates_final; //final coordinate matrix of nodes

for (int i=0;i<n_adof;++i){
 final_coordinates_vec(Index_Id_nodes.coeffRef(0,i))+=X_global.coeffRef(i);
};
for (int i=0;i<(nelements_x+1)*(nelements_y+1);i++){
 for (int j=0;j<ndof_ele;++j){
  coordinates_final(i,j)=final_coordinates_vec(i*ndof_ele+j);
 };
};

std::ofstream f_c;
f_c.open ("final_coordinates.txt");
f_c << coordinates_final;
f_c.close();

std::ofstream f_i;
f_i.open ("initial_coordinates.txt");
f_i << coordinates;
f_i.close();
//std::cout<<std::endl<<coordinates_final<<std::endl;
//std::cout<<std::endl<<X_global<<std::endl;
Eigen::VectorXd X_global_dense; X_global_dense=Eigen::VectorXd(X_global);
std::cout<<std::endl<<"Max displacement: " << X_global_dense.maxCoeff()<<std::endl;
std::cout<<std::endl<<"Min displacement: " << X_global_dense.minCoeff()<<std::endl;




///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//GENERATE OUTPUTS

//std::ofstream K_g;
//K_g.open ("K_global.txt");
//K_g << K_global;
//K_g.close();

//std::ofstream F_g;
//F_g.open ("F_global.txt");
//F_g << F_global;
//F_g.close();

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//DISPLAY INTERMEDIATE VARIABLES

//DISPLAY COORD
//COORD A=mesh.coord();
//for (int i=0;i<(nelements_x+1)*(nelements_y+1);++i){
//std::cout<<std::endl;
// for (int j=0;j<2;++j){
//std::cout<<A(i,j)<<"  ";
// };
//};


//DISPLAY LINKS
//LINKS B=mesh.links();
//for (int i=0;i<4;++i){
//std::cout<<std::endl; 
// for (int j=0;j<(nelements_x*nelements_y);++j){
//  std::cout<<B(i,j)<<"  ";
// };/
//};

//DISPLAY COORD_ELE(ELEMENT)
//int e;
//typedef typename Mesh<COORD,LINKS>::COORD_ELE COORD_ELE;
//COORD_ELE CE=mesh.coord_ele(e);
//for (int i=0;i<4;++i){
// std::cout<<std::endl;
// for (int j=0;j<2;++j){
//  std::cout<<CE(i,j)<<"  ";
// };
//};

//DISPLAY BOUND
//BOUND A=bound.bound_matrix();
//for (int i=0;i<(nelements_x+1)*(nelements_y+1);++i){
// std::cout<<std::endl;
// for (int j=0;j<2;++j){
//  std::cout<<A(i,j)<<"  ";
// };
//};

//DISPLAY ID_ELEMENTS
//for (int i=0;i<ndof_ele*nnodes_ele;++i){
// std::cout<<std::endl;
// for (int j=0;j<nelements_x*nelements_y;++j){
//  std::cout<<Id_Elements(i,j)<<"    ";
// };
//};

}

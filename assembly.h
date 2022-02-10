template <typename ID_NODES,typename ID_ELEMENTS, typename BOUND, typename COORD, typename LINKS, typename MAX_INDEX, typename DIFF> class Assemble{
public:
ID_NODES id_nodes(int nelements_x,int nelements_y,int ndof_ele,Bound<BOUND> bound);
ID_ELEMENTS id_elements(int nnodes_ele,int nelements_x,int nelements_y,int ndof_ele,Mesh<COORD,LINKS> mesh,Bound<BOUND> bound);
int const nadof(int nelements_x,int nelements_y,int ndof_ele,Bound<BOUND> bound);
int bandwidth(int nelements_x,int nelements_y,int ndof_ele,int nnodes_ele,Mesh<COORD,LINKS> mesh,Bound<BOUND> bound); //bandwith for sparse matrix K_global
};


template <typename ID_NODES, typename ID_ELEMENTS, typename BOUND, typename COORD, typename LINKS, typename MAX_INDEX, typename DIFF> ID_NODES Assemble<ID_NODES,ID_ELEMENTS,BOUND,COORD,LINKS,MAX_INDEX,DIFF>::id_nodes(int nelements_x,int nelements_y,int ndof_ele,Bound<BOUND> bound){
int var=0; //number of equations to eliminate
int num_equation;
BOUND B_matrix=bound.bound_matrix();
ID_NODES idn;
 for (int i=0;i<(nelements_x+1)*(nelements_y+1);++i){
  for (int j=0;j<ndof_ele;++j){
   if (B_matrix(i,j)==1){
    var++;
    idn(j,i)=-1;
   }
   else{
    num_equation=i*ndof_ele+j-var; //number of global equation
    idn(j,i)=num_equation;
   }; 
  };
 };
return idn;
};

template <typename ID_NODES, typename ID_ELEMENTS, typename BOUND, typename COORD, typename LINKS, typename MAX_INDEX, typename DIFF> ID_ELEMENTS Assemble<ID_NODES,ID_ELEMENTS,BOUND,COORD,LINKS,MAX_INDEX,DIFF>::id_elements(int nnodes_ele,int nelements_x,int nelements_y,int ndof_ele,Mesh<COORD,LINKS> mesh,Bound<BOUND> bound){
ID_ELEMENTS ide;
LINKS link=mesh.links();
ID_NODES idnodes=id_nodes(nelements_x,nelements_y,ndof_ele,bound);
int idof_ele;
 for (int i=0;i<nelements_x*nelements_y;++i){
  for (int j=0;j<nnodes_ele;++j){
   for (int idof=0;idof<ndof_ele;++idof){
    idof_ele=j*ndof_ele+idof;
    ide(idof_ele,i)=idnodes(idof,link(j,i));
   };
  };
 };
return ide;
};

template <typename ID_NODES, typename ID_ELEMENTS, typename BOUND, typename COORD, typename LINKS, typename MAX_INDEX, typename DIFF> int const Assemble<ID_NODES,ID_ELEMENTS,BOUND,COORD,LINKS,MAX_INDEX,DIFF>::nadof(int nelements_x,int nelements_y,int ndof_ele,Bound<BOUND> bound){
int var=0; //number of equations to eliminate
int n;  //number of active degrees of freedom
BOUND B_matrix=bound.bound_matrix();
 for (int i=0;i<(nelements_x+1)*(nelements_y+1);++i){
  for (int j=0;j<ndof_ele;++j){
   if (B_matrix(i,j)==1){
    var++;
   };
  };
 };
n=(nelements_x+1)*(nelements_y+1)*ndof_ele-var;
return n;
};

template <typename ID_NODES, typename ID_ELEMENTS, typename BOUND, typename COORD, typename LINKS, typename MAX_INDEX, typename DIFF> int Assemble<ID_NODES,ID_ELEMENTS,BOUND,COORD,LINKS,MAX_INDEX,DIFF>::bandwidth(int nelements_x,int nelements_y,int ndof_ele,int nnodes_ele,Mesh<COORD,LINKS> mesh,Bound<BOUND> bound){
ID_ELEMENTS ide=id_elements(nnodes_ele,nelements_x,nelements_y,ndof_ele,mesh,bound);
MAX_INDEX max_index;
DIFF diff;
double min;
int min_index,bnw,sum,bandw;
 for (int i=0;i<nelements_x*nelements_y;++i){
  min=ide(0,i); min_index=0;
  ide.col(i).maxCoeff(&max_index);
  for (int j=0;j<nnodes_ele*ndof_ele;++j){
   if (ide(j,i)<min && ide(j,i)>-1){min=ide(j,i); min_index=j;};
  };
 diff(i)=max_index-min_index;
 };
bnw=diff.maxCoeff()+1;
sum=0;
for (int i=0;i<bnw;++i){sum+=i;};
bandw=nadof(nelements_x,nelements_y,ndof_ele,bound)+2*(nadof(nelements_x,nelements_y,ndof_ele,bound)*(bnw-1)-sum);
return bandw;
};



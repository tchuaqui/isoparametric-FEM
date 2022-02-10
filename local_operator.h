template <typename LOAD_ELE, typename COORD, typename LINKS, typename ELASTIC> class Local_operator{
public:
typedef typename Integration::GAUSS GAUSS;
typedef typename Mesh<COORD,LINKS>::COORD_ELE COORD_ELE; 
typedef typename Shape_func::SHAPE_F SHAPE_F;
typedef typename Shape_func::DSHAPE_F DSHAPE_F;
typedef Eigen::Matrix<double,8,8> K_LOCAL;
typedef Eigen::Matrix<double,2,2> JACOB;
typedef Eigen::Matrix<double,3,8> BMATRIX;
typedef Eigen::Matrix<double,8,3> BMATRIX_T;
typedef Eigen::Matrix<double,8,1> F_LOCAL;

JACOB jacobian(int e,DSHAPE_F dNS,COORD_ELE Ce);
BMATRIX Bmatrix(DSHAPE_F dNsdxy);
K_LOCAL K_local(int e,int ngauss,Shape_func shape_func,Integration integration,Mesh<COORD,LINKS> mesh,Elastic<ELASTIC> elastic); 
F_LOCAL F_local(int e, int ngauss, LOAD_ELE load, Shape_func shape_fun, Integration integration,Mesh<COORD,LINKS> mesh);
};

template <typename LOAD_ELE,typename COORD,typename LINKS,typename ELASTIC> typename Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::K_LOCAL Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::K_local(int e,int ngauss,Shape_func shape_func,Integration integration,Mesh<COORD,LINKS> mesh,Elastic<ELASTIC> elastic){
GAUSS p=integration.Gauss_coord();
GAUSS w=integration.Gauss_weight();
COORD_ELE Ce=mesh.coord_ele(e);
ELASTIC C=elastic.C();
K_LOCAL K; K.setZero();

 for (int i=0;i<ngauss;++i){
  for (int j=0;j<ngauss;++j){ 
   DSHAPE_F dNS=shape_func.dshape_f_local(p(i),p(j));
   DSHAPE_F dNsdxy=jacobian(e,dNS,Ce).inverse()*dNS;
   BMATRIX B=Bmatrix(dNsdxy);
   BMATRIX_T Dt=B.transpose()*C;
   K_LOCAL D=Dt*B;
   K=K+D*w(i)*w(j)*jacobian(e,dNS,Ce).determinant();
  };
 };
return K;
};

template <typename LOAD_ELE, typename COORD, typename LINKS, typename ELASTIC> typename Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::F_LOCAL Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::F_local(int e,int ngauss,LOAD_ELE load,Shape_func shape_func,Integration integration,Mesh<COORD,LINKS> mesh){
GAUSS p=integration.Gauss_coord();
GAUSS w=integration.Gauss_weight();
COORD_ELE Ce=mesh.coord_ele(e);
F_LOCAL F; F.setZero();

 for (int i=0;i<ngauss;++i){
  for (int j=0;j<ngauss;++j){
   SHAPE_F Ns=shape_func.shape_f(p(i),p(j));
   DSHAPE_F dNS=shape_func.dshape_f_local(p(i),p(j));
   F=F+Ns.transpose()*load*jacobian(e,dNS,Ce).determinant()*w(i)*w(j);
  };
 };
return F;
};

template <typename LOAD_ELE, typename COORD,typename LINKS,typename ELASTIC> typename Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::JACOB Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::jacobian(int e,DSHAPE_F dNS,COORD_ELE Ce){
JACOB jac;
jac=dNS*Ce;
return jac;
};
template <typename LOAD_ELE, typename COORD, typename LINKS, typename ELASTIC> typename Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::BMATRIX Local_operator<LOAD_ELE,COORD,LINKS,ELASTIC>::Bmatrix(DSHAPE_F dNsdxy){
BMATRIX B;
B(0,0)=dNsdxy(0,0); B(0,1)=0; B(0,2)=dNsdxy(0,1); B(0,3)=0; B(0,4)=dNsdxy(0,2); B(0,5)=0; B(0,6)=dNsdxy(0,3); B(0,7)=0;
B(1,0)=0; B(1,1)=dNsdxy(1,0); B(1,2)=0; B(1,3)=dNsdxy(1,1); B(1,4)=0; B(1,5)=dNsdxy(1,2); B(1,6)=0; B(1,7)=dNsdxy(1,3);
B(2,0)=dNsdxy(1,0); B(2,1)=dNsdxy(0,0); B(2,2)=dNsdxy(1,1); B(2,3)=dNsdxy(0,1); B(2,4)=dNsdxy(1,2); B(2,5)=dNsdxy(0,2); B(2,6)=dNsdxy(1,3); B(2,7)=dNsdxy(0,3);
return B;
};





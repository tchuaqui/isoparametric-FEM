class Shape_func{
public:
typedef Eigen::Matrix<double,2,8> SHAPE_F;
typedef Eigen::Matrix<double,2,4> DSHAPE_F;
SHAPE_F shape_f(double xi, double eta);
DSHAPE_F dshape_f_local(double xi, double eta);
};
typename Shape_func::SHAPE_F Shape_func::shape_f(double xi,double eta){
 SHAPE_F S; 
 S(0,0)=(1.0-xi)*(1.0-eta)/4.0; S(0,1)=0.0; S(1,0)=0.0; S(1,1)=(1.0-xi)*(1.0-eta)/4.0;
 S(0,2)=(1.0+xi)*(1.0-eta)/4.0; S(0,3)=0.0; S(1,2)=0.0; S(1,3)=(1.0+xi)*(1.0-eta)/4.0;
 S(0,4)=(1.0+xi)*(1.0+eta)/4.0; S(0,5)=0.0; S(1,4)=0.0; S(1,5)=(1.0+xi)*(1.0+eta)/4.0;
 S(0,6)=(1.0-xi)*(1.0+eta)/4.0; S(0,7)=0.0; S(1,6)=0.0; S(1,7)=(1.0-xi)*(1.0+eta)/4.0;
 return S;
};
 typename Shape_func::DSHAPE_F Shape_func::dshape_f_local(double xi, double eta){
 DSHAPE_F DS;
 DS(0,0)=-(1.0-eta)/4.0; DS(0,1)=(1.0-eta)/4.0; DS(0,2)=(1.0+eta)/4.0; DS(0,3)=-(1.0+eta)/4.0;
 DS(1,0)=-(1.0-xi)/4.0; DS(1,1)=-(1.0+xi)/4.0; DS(1,2)=(1.0+xi)/4.0; DS(1,3)=(1.0-xi)/4.0;
 return DS;
}; 

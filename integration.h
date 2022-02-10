class Integration{
int gauss_points;
public:
typedef Eigen::Matrix<double,3,1> GAUSS;
Integration(int n);
GAUSS Gauss_coord();
GAUSS Gauss_weight();
};

Integration::Integration(int n) : gauss_points(n){};
typename Integration::GAUSS Integration::Gauss_coord(){
 GAUSS Coor;
 switch (gauss_points){
 case 2:
 Coor(0,0)=-1.0/sqrt(3.0);
 Coor(1,0)=1.0/sqrt(3.0);
 break;
 case 3:
 Coor(0,0)=-sqrt(0.6);
 Coor(2,0)=sqrt(0.6); 
 break;
 return Coor; 
 };
};
typename Integration::GAUSS Integration::Gauss_weight(){
 GAUSS Weight;
 switch (gauss_points){
 case 1:
 Weight(0,0)=2.0;
 break;
 case 2:
 Weight(0,0)=1.0;
 Weight(1,0)=1.0;
 break;
 case 3:
 Weight(0,0)=5.0/9.0;
 Weight(1,0)=8.0/9.0;
 Weight(2,0)=5.0/9.0; 
 break;
 return Weight; 
 };
};

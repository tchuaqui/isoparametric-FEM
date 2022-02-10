template <typename ELASTIC> class Elastic{
double E, v;
public:
Elastic(double Y,double nu);
ELASTIC C();
};

template<typename ELASTIC> Elastic<ELASTIC>::Elastic(double Y, double nu) : E(Y), v(nu){};
template<typename ELASTIC> ELASTIC Elastic<ELASTIC>::C(){
 ELASTIC A;
 A(0,0)=E/(1.0-v*v); A(0,1)=E*v/(1.0-v*v); A(0,2)=0.0;
 A(1,0)=E*v/(1.0-v*v); A(1,1)=E/(1.0-v*v); A(1,2)=0.0;
 A(2,0)=0.0, A(2,1)=0.0; A(2,2)=E/(2.0*(1.0+v));
 return A;
};

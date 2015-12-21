#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/src/LU/FullPivLU.h>
#include <gsl/gsl_sf_bessel.h>

// mathematica->c conversion
#define Power pow
#define BesselJ gsl_sf_bessel_Jn
#define Sqrt sqrt
#define Cos cos
#define Sin sin
#define Abs abs
#define ArcTan atan2

using namespace Eigen;
using namespace std;


//SELECT THE SYSTEM TO INTEGRATE
// #define freeparticle
// #define elastic
// #define elastic_phasespace
#define guidingcenter
// #define guidingcenter3D // (implicit scheme 2)

double h =800; //timestep
const int MAX_T = 1E6; //number of integration steps
ofstream out("out_gen.txt");  //output file
const int NEWTON = 0; //number of implicit newton iterations
const int PRINT_MULTIPLE = 8; //print out every n timesteps
const int PRINT_OFFSET = 1; //start printing from a specific step
const int PRINT_PRECISION = 15; //number of significative digits
const int INIT_STEPS = 0;  //init with auxiliary integrator. 0-> Lagrangian init, ~inf -> Hamiltonian init
const double ORBIT_NORMALIZE = 50;  //normalize the timesteps 

const double hx = 1.E-5;  //step for numerical derivative
// #define richardson   //use richardson extrapolation (default: central difference)


//FIRST GUESS INTEGRATOR
// #define step_guess step_euler
// #define step_guess step_RK4
#define step_guess step_explicit_guidingcenter_momentum
#define explicit // set this if step_explicit_guidingcenter_momentum is set

//implicit integrator, alpha method, verlet or semi-explicit method (qin)
#define discrete_lagrangian discrete_lagrangian_alpha
// #define discrete_lagrangian discrete_lagrangian_verlet
// #define discrete_lagrangian discrete_lagrangian_qin
const double alpha=0.5; //0.5 == Lagrangian midpoint, 0,1 == lagrangian euler


const int DEBUG_TIMESTEP_MULT = 10000; // 0=DISABLE, print to screen the timestep
const bool DEBUG_IMPLICIT_F = 0; //print the implicit function to be converged 
const bool EXIT_ON_ERROR = 0; //exit if the error is too high
const double ERROR_THRESHOLD = 0.1; //error threshold
const int TIME_OFFSET = 0; //start the integrator from a specific timestep

#define LYAPUNOV //compute lyapunov coefficients, not available in implicit 2 (FIXME)
ofstream lyap("out_lyapunov.txt"); //lyapunov output file
double hlyap = 1.E-7; //trajectories separation


#if defined(guidingcenter) || defined(guidingcenter3D)
    //guiding center configuration 

    #define tokamak //(configuration B)
     // #define forcefree //(configuration C)
    // #define twodfield //(configuration A)

    // EXPLICIT PARAMETERS
    const int explicit_scheme = 4; // 1-2-3-4
    const bool DERIVATIVE_2B = 1; //EXCLUDE/INCLUDE SECOND DERIVATIVE OF B IN EXPLICIT 4
    const bool DERIVATIVE_2U = 1; //EXCLUDE/INCLUDE SECOND DERIVATIVE OF U IN EXPLICIT 4

    //field parameters    
    const double B0 = 1.;
    const double R0 = 1.;
    const double q = 2.;        //Safety Factor
    const double kt = 1.;       //toroidal magnetic field multiplicative factor. Set to a low value to stabilize unstable integrators
    const double mu = 2.25E-6;

    // FORCE-FREE
    const double a = 3.;
    const double kff = 0.01;    //perturbation magnitude
    const bool forcefree_pert=0;  //set/unset perturbation

#endif


#ifdef guidingcenter
    // guidingcenter 4D configuration
    const int DIM = 4;

    void momentum(Matrix<double,DIM,1> q,Matrix<double,DIM,1> &p);

    void initial_conditions(Matrix<double,DIM,1> &q0,Matrix<double,DIM,1> &p0){
      q0 << 0.050000000000000, 0.00000000000000 ,0.000000000000000 ,0.000390000000000;
      momentum(q0,p0);
    }

#endif

#ifdef guidingcenter3D
    // guidingcenter 3D configuration
    const int DIM = 3;
    #define step3D

    Vector3d guiding_B(Vector3d x);
    Vector3d guiding_A(Vector3d x);

    void initial_conditions(Matrix<double,DIM,1> &q0,Matrix<double,DIM,1> &p0){

      q0 << 0.050000000000000, 0.00000000000000 ,0.000000000000000;
      double u0 = 0.00039;
      p0 = u0*guiding_B(q0).normalized() + guiding_A(q0);
    }

#endif



#ifdef freeparticle
    const int DIM = 1;

    const double m=1.;
    const double k=1.;

    double energy(Matrix<double,DIM,1> q,Matrix<double,DIM,1> p){return (p(0)*p(0)/(2.*m)-k*q(0));}
    double lagrangian(Matrix<double,DIM,1> q,Matrix<double,DIM,1> v){return (0.5*m*v(0)*v(0) + k*q(0));}
    void f_eq_motion(Matrix<double,2*DIM,1> z,Matrix<double,2*DIM,1> &f){
      f(0)=z(1)/m;
      f(1)=k;
    }
#endif


#ifdef elastic
    const int DIM = 1;
    
    const double m=1.;
    const double k=1.;
    
    void initial_conditions(Matrix<double,DIM,1> &q0, Matrix<double,DIM,1> &p0){
      q0 << 1;
      p0 << 1;
    }
    
    double energy(Matrix<double,DIM,1> q,Matrix<double,DIM,1> p){return (p(0)*p(0)/(2.*m)+0.5*k*q(0)*q(0));}
    double lagrangian(Matrix<double,DIM,1> q,Matrix<double,DIM,1> v){return (0.5*m*v(0)*v(0)-0.5*k*q(0)*q(0));}
    void f_eq_motion(Matrix<double,2*DIM,1> z,Matrix<double,2*DIM,1> &f){
      f(0)=z(1)/m;
      f(1)=-k*z(0);
    }
#endif


#ifdef elastic_phasespace
    const int DIM = 2;

    const double m=1.;
    const double k=1.;
    
    void initial_conditions(Matrix<double,DIM,1> &q0, Matrix<double,DIM,1> &p0){
      q0 << 2,1;
      p0 << 1,0;
    }

    double energy(Matrix<double,DIM,1> q,Matrix<double,DIM,1> p){
      return (p(0)*p(0)/(2.*m)+0.5*k*q(0)*q(0));
    }
    double lagrangian(Matrix<double,DIM,1> q,Matrix<double,DIM,1> v){
      return (v(0)*q(1)-q(1)*q(1)/(2.*m)-0.5*k*q(0)*q(0));
    }
    void f_eq_motion(Matrix<double,2*DIM,1> z,Matrix<double,2*DIM,1> &f){
      f.setZero();
      f(0)=z(1)/m;
      f(1)=-k*z(0);
    }
    void momentum(Matrix<double,DIM,1> q,Matrix<double,DIM,1> &p){
      p(0) = q(1);
      p(1) = 0.;
    }
      
#endif




//GUIDING CENTER SYSTEM DEFINITIONS
#if defined(guidingcenter) || defined(guidingcenter3D)

    #ifdef tokamak
    // TOKAMAK
    Vector3d guiding_A(Vector3d x){
      Vector3d ret;
      ret(0) = 0;
      ret(1) = -B0*R0*log((R0+x(0))/R0)*kt;
      ret(2) = B0/(2.*q*(R0+x(0)))*( 2.*R0*(R0+x(0))*log((R0+x(0))/R0) - 2.*R0*x(0) -2.*x(0)*x(0) - x(1)*x(1) );
      return ret;
    }
    Vector3d guiding_B(Vector3d x){
      Vector3d ret;
      ret(0) = -B0*x(1)/(q*(R0+x(0)));
      ret(1) = -((B0 *( R0 *x(0) *(-2.) -2.*x(0)*x(0) + x(1)*x(1)))/(2. * q *(R0 + x(0))*(R0 + x(0))));
      ret(2) = -B0*R0/(R0+x(0))*kt;
      return ret;
    } 
    #endif


    #ifdef twodfield
    // // 2D FIELD
    Vector3d guiding_B(Vector3d x){
      Vector3d ret(0.,0.,0.);
      ret(2)=1.+0.05*(x(0)*x(0)/4. + x(1)*x(1));
      return ret;
    }
    Vector3d guiding_A(Vector3d x){
      Vector3d ret(0.,0.,0.);
      ret(0)=(-1.)*0.05/3. * x(1)*x(1)*x(1);
      ret(1)=0.05/3. * x(0)*x(0)*x(0)/4. +x(0);
      ret(2)=0;
      return ret;
    }
    #endif


    // FORCE-FREE
    #ifdef forcefree
    Vector3d calc_A_nm(Vector3d x,int n,int m){
      return Vector3d(((sqrt((a - n)*(a + n))*x(1)*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))))/sqrt(x(0)*x(0) + x(1)*x(1)) - 
        (m*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(1)*cos(n*x(2) + m*atan2(x(1),x(0))) + x(0)*sin(n*x(2) + m*atan2(x(1),x(0)))))/(x(0)*x(0) + x(1)*x(1)))/a, 
      (-((sqrt((a - n)*(a + n))*x(0)*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))))/sqrt(x(0)*x(0) + x(1)*x(1))) + 
        (m*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(0)*cos(n*x(2) + m*atan2(x(1),x(0))) - x(1)*sin(n*x(2) + m*atan2(x(1),x(0)))))/(x(0)*x(0) + x(1)*x(1)))/a, 
      gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))));
    }

      
    Vector3d calc_B_nm(Vector3d x,int m,int n){
      return Vector3d
      ((1/(a*(x(0)*x(0) + x(1)*x(1))))*((-m)*(a - n)*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(1)*cos(n*x(2) + m*atan2(x(1),x(0))) + x(0)*sin(n*x(2) + m*atan2(x(1),x(0)))) + 
        sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1))*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(a*x(1)*cos(n*x(2) + m*atan2(x(1),x(0))) - n*x(0)*sin(n*x(2) + m*atan2(x(1),x(0))))), 
      (1/(a*(x(0)*x(0) + x(1)*x(1))))*(m*(a - n)*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(x(0)*cos(n*x(2) + m*atan2(x(1),x(0))) - x(1)*sin(n*x(2) + m*atan2(x(1),x(0)))) - 
        sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1))*gsl_sf_bessel_Jn(-1 + m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*(a*x(0)*cos(n*x(2) + m*atan2(x(1),x(0))) + n*x(1)*sin(n*x(2) + m*atan2(x(1),x(0))))), 
      ((a - n)*(a + n)*gsl_sf_bessel_Jn(m, sqrt(a*a - n*n)*sqrt(x(0)*x(0) + x(1)*x(1)))*cos(n*x(2) + m*atan2(x(1),x(0))))/a);
    }

    Vector3d guiding_B(Vector3d x){
      if (forcefree_pert) return 1.*calc_B_nm(x,0,0)+kff*calc_B_nm(x,1,1)+kff*calc_B_nm(x,1,2);
      else return calc_B_nm(x,0,0);
    } 

    Vector3d guiding_A(Vector3d x){
      if (forcefree_pert) return 1.*calc_A_nm(x,0,0)+kff*calc_A_nm(x,1,1)+kff*calc_A_nm(x,1,2);
      else return calc_A_nm(x,0,0);
    }
    #endif


    //magnetic field computation
    Vector3d Bgrad(Vector3d x){
      Vector3d dx(hx,hx,hx),x1,x0; //dx := (dx,dy,dz)
      Vector3d ret;


      //COMPUTE GRADIENT(B),GRADIENT(phi),JAC(A_dagger)
      for (int j=0;j<3;j++){
        x0 = x1 = x;
        x0(j)-=dx(j);
        x1(j)+=dx(j);
        ret(j) = 0.5*(guiding_B(x1).norm() - guiding_B(x0).norm())/dx(j);
      }

      return ret;
    }


    Matrix3d B_Hessian(Vector3d x){
      Vector3d dx(hx,hx,hx),x1,x0; //dx := (dx,dy,dz)
      Matrix3d ret;


      //COMPUTE GRADIENT(B),GRADIENT(phi),JAC(A_dagger)
      for (int j=0;j<3;j++){
        x0 = x1 = x;
        x0(j)-=dx(j);
        x1(j)+=dx(j);
        ret.col(j) = 0.5*(Bgrad(x1) - Bgrad(x0))/dx(j);
      }

      return ret;
    }

    void guiding_EM_field(Vector3d &x,double &u,Vector3d &B,double &B_norm,Vector3d &b,Vector3d &B_grad,
    Vector3d &B_dag,Vector3d &A_dag,Matrix3d &A_dag_jac,double &phi,Vector3d &phi_grad){
      Matrix3d A_jac;
      Vector3d A;
      A = guiding_A(x);
      Matrix3d b_jac;
      B = guiding_B(x);
      B_norm = B.norm();
      b = B.normalized();
      A_dag = A + u*b;
      
      Vector3d B0,x0,x1,B1;
      
      //COMPUTE SPACESTEP
      Vector3d dx(hx,hx,hx); //dx := (dx,dy,dz)

      //COMPUTE GRADIENT(B),GRADIENT(phi),JAC(A_dagger)
      for (int j=0;j<3;j++){
        x0 = x1 = x;
        x0(j)-=dx(j);
        x1(j)+=dx(j);
        B0 = guiding_B(x0);
        B1 = guiding_B(x1);
        
        A_dag_jac.col(j) = 0.5*(guiding_A(x1) + u*B1.normalized() - guiding_A(x0)-u*B0.normalized())/dx(j);
        A_jac.col(j) = 0.5*(guiding_A(x1) - guiding_A(x0))/dx(j);
        B_grad(j) = 0.5*(B1.norm() - B0.norm())/dx(j);
        b_jac.col(j) = 0.5*(B1.normalized() - B0.normalized())/dx(j);
      }
      
      //COMPUTE B_dagger
      B_dag = B;
      B_dag(0) += u*(b_jac(2,1)-b_jac(1,2));
      B_dag(1) += u*(b_jac(0,2)-b_jac(2,0));
      B_dag(2) += u*(b_jac(1,0)-b_jac(0,1));
      
    }

    void legendre_left_inverse_explicit(double h,Vector4d q0,Vector4d p0, Vector4d &q1){
    // void legendre_left_inverse_explicit(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> p0, Matrix<double,DIM,1> &q1){
      Vector3d B, b, B_grad, B_dag, A_dag,phi_grad,E_dag,x0;
      double B_norm,phi,u0;
      Matrix3d A_dag_jac;
      
      x0 = q0.head(3);
      u0 = q0(3);
      
      guiding_EM_field(x0,u0,B,B_norm,b,B_grad,B_dag,A_dag,A_dag_jac,phi,phi_grad);
      
      Matrix4d M;
      Vector4d W,Q;
      
      //~ //BUILD M
      M(0,1) = B_dag(2);
      M(0,2) = -B_dag(1);
      M(1,0) = -B_dag(2);
      M(1,2) = B_dag(0);
      M(2,0) = B_dag(1);  
      M(2,1) = -B_dag(0);   
      M(0,3) = -b(0);
      M(1,3) = -b(1);
      M(2,3) = -b(2);
      M(3,0)=b(0);
      M(3,1)=b(1);
      M(3,2)=b(2);
      M(0,0) = M(1,1) = M(2,2) = M(3,3) = 0;
      if (explicit_scheme==2) M(3,3) = -h; //QIN version
      else if (explicit_scheme==3){
        //qin modified version
        M(3,3) = -h;
        M(0,3) = M(1,3) = M(2,3) = 0;
        for (int i=0;i<=2;i++) for (int j=0;j<=2;j++) M(i,j) -= 2.*b(i)*b(j)/h;
      }
      else if (explicit_scheme==4){
        //second order hamiltonian
        Matrix4d grad2h = Matrix4d::Zero();
        if (DERIVATIVE_2B) grad2h.block<3,3>(0,0) = mu*B_Hessian(x0);
        if (DERIVATIVE_2U) grad2h(3,3) = 1.;
        M -= h/2.*grad2h;
      }

      M/=2.;

      //BUILD W
      W.head(3) = h/2.*(mu*B_grad) +A_dag -p0.head(3);
      W(3) = (h/2.*u0 -p0(3));
      
      Q = M.inverse()*W;
      
      q1.head(3) = (x0+Q.head(3));
      q1(3) = (u0+Q(3));
    }

    void legendre_right_explicit(double h,Vector4d q0,Vector4d q1, Vector4d &p1){
    // void legendre_right_explicit(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1, Matrix<double,DIM,1> &p1){
      Vector3d B, b, B_grad, B_dag, A_dag,phi_grad,E_dag,x0,x1;
      double B_norm,phi,u0,u1;
      Matrix3d A_dag_jac;
      
      x0 = q0.head(3);
      u0 = q0(3);
      x1 = q1.head(3);
      u1 = q1(3);
      
      guiding_EM_field(x1,u1,B,B_norm,b,B_grad,B_dag,A_dag,A_dag_jac,phi,phi_grad);
      
      Matrix4d M;
      Vector4d W,Q;
      
      //~ //BUILD M
      M(0,1) = B_dag(2);
      M(0,2) = -B_dag(1);
      M(1,0) = -B_dag(2);
      M(1,2) = B_dag(0);
      M(2,0) = B_dag(1);  
      M(2,1) = -B_dag(0);   
      M(0,3) = -b(0);
      M(1,3) = -b(1);
      M(2,3) = -b(2);
      M(3,0)=b(0);
      M(3,1)=b(1);
      M(3,2)=b(2);
      M(0,0) = M(1,1) = M(2,2) = M(3,3) = 0;
      if (explicit_scheme==2) M(3,3) = h; //qin version
      else if (explicit_scheme==3){
        //qin modified version
        M(3,3) = h;
        M(0,3) = M(1,3) = M(2,3) = 0;
        for (int i=0;i<=2;i++) for (int j=0;j<=2;j++) M(i,j) += 2.*b(i)*b(j)/h;
      }
      else if (explicit_scheme==4){
        //second order hamiltonian
        Matrix4d grad2h = Matrix4d::Zero();
        if (DERIVATIVE_2B) grad2h.block<3,3>(0,0) = mu*B_Hessian(x1);
        if (DERIVATIVE_2U) grad2h(3,3) = 1.;
        M += h/2.*grad2h;
      }
      M/=2.;

      
      Vector4d dq;
      dq.head(3) = x1-x0;
      dq(3) = u1-u0;
      Q = M*dq;

      p1.head(3) = Q.head(3) +A_dag -h/2.*mu*B_grad;
      p1(3) = Q(3) -h/2.*u1;
    }
    
    void step_explicit_guidingcenter_momentum(double h,Vector4d q0,Vector4d p0, Vector4d &q1,Vector4d &p1){
    // void step_explicit_guidingcenter_momentum(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> p0, Matrix<double,DIM,1> &q1,Matrix<double,DIM,1> &p1){
      legendre_left_inverse_explicit(h,q0,p0,q1);
      legendre_right_explicit(h,q0,q1,p1);
    }
#endif


#ifdef guidingcenter
    #define compute_momentum

    double energy(Matrix<double,DIM,1> q,Matrix<double,DIM,1> p){
      return (0.5*q(3)*q(3)+mu*guiding_B(q.head(3)).norm());
    }
    double lagrangian(Matrix<double,DIM,1> q,Matrix<double,DIM,1> v){
      Vector3d x = q.head(3);
      Vector3d vx = v.head(3);
      double u = q(3);
      Vector3d B = guiding_B(x);
      Vector3d Adag = guiding_A(x)+u*B.normalized();

      return (Adag.dot(vx)-0.5*u*u-mu*B.norm());
    }


    void f_eq_motion(Matrix<double,2*DIM,1> z,Matrix<double,2*DIM,1> &f){
      Vector3d B, b, B_grad, B_dag, A_dag,phi_grad,E_dag,x;
      double B_norm,phi,B_dag_par,u;
      Matrix3d A_dag_jac;
      
      f.setZero();
      x = z.head(3);
      u = z(3);
      
      guiding_EM_field(x,u,B,B_norm,b,B_grad,B_dag,A_dag,A_dag_jac,phi,phi_grad);
      B_dag_par = B_dag.dot(b);
      E_dag = -mu*B_grad;
      
      f.head(3) = (u* B_dag - b.cross(E_dag))/B_dag_par;
      f(3) = B_dag.dot(E_dag)/B_dag_par;
    }
    void momentum(Matrix<double,DIM,1> q,Matrix<double,DIM,1> &p){
      p(3) = 0.;
      Vector3d x = q.head(3);
      p.head(3) = guiding_A(x)+q(3)*guiding_B(x).normalized();
    }
#endif


#ifdef guidingcenter3D

    double energy(Matrix<double,DIM,1> q,Matrix<double,DIM,1> p){
      double u = (p- guiding_A(q)).norm();
      return (0.5*u*u+mu*guiding_B(q).norm());
    }
    double lagrangian(Matrix<double,DIM,1> q,Matrix<double,DIM,1> v){
      Vector3d B = guiding_B(q);
      Vector3d b = B.normalized();
      double u = b.dot(v);
      Vector3d Adag = guiding_A(q)+u*b;

      #ifdef timestep0
        // FIXME!!
        return 0;
      #else
        return (Adag.dot(v)-0.5*u*u-mu*B.norm());
      #endif
      
    }

    void f_eq_motion(Matrix<double,6,1> z,Matrix<double,6,1> &f){}

    void f_eq_motion_3D(Matrix<double,8,1> z,Matrix<double,8,1> &f){
      Vector3d B, b, B_grad, B_dag, A_dag,phi_grad,E_dag,x;
      double B_norm,phi,B_dag_par,u;
      Matrix3d A_dag_jac;
      
      f.setZero();
      x = z.head(3);
      u = z(3);
      
      guiding_EM_field(x,u,B,B_norm,b,B_grad,B_dag,A_dag,A_dag_jac,phi,phi_grad);
      B_dag_par = B_dag.dot(b);
      E_dag = -mu*B_grad;
      
      f.head(3) = (u* B_dag - b.cross(E_dag))/B_dag_par;
      f(3) = B_dag.dot(E_dag)/B_dag_par;
    }

    void step_RK4_3D(double h,Vector4d q0,Vector4d p0, Vector4d &q1,Vector4d &p1){
      Matrix<double,8,1> k1,k2,k3,k4,z0,z1;
      
      z0.head(4) = q0;
      z0.tail(4) = p0;
      
      f_eq_motion_3D(z0,k1);
      f_eq_motion_3D(z0+0.5*h*k1,k2);
      f_eq_motion_3D(z0+0.5*h*k2,k3);
      f_eq_motion_3D(z0+h*k3,k4);
      
      z1 = z0 + 1./6.*h*(k1+2.*k2+2.*k3+k4);
      
      q1 = z1.head(4);
      p1 = z1.tail(4);
      
      #ifdef init_momentum
        momentum(q0,p0);
      #endif
    }

#endif






double discrete_lagrangian_verlet(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1){
  return 0.5*h*lagrangian(q0,(q1-q0)/h) + 0.5*h*lagrangian(q1,(q1-q0)/h);
}
double discrete_lagrangian_alpha(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1){
  #ifdef timestep0
    return lagrangian((1.-alpha)*q0 + alpha * q1,(q1-q0));
  #else
  return h*lagrangian((1.-alpha)*q0 + alpha * q1,(q1-q0)/h);
  #endif
}


#ifdef guidingcenter
  double discrete_lagrangian_qin(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1){
    Vector3d x0=q0.head(3),x1=q1.head(3);
    double u0 = q0(3),u1=q1(3);
    Vector3d B0 = guiding_B(x0),B1=guiding_B(x1);
    Vector3d Adag0 = guiding_A(x0) + u0*B0.normalized();
    Vector3d Adag1 = guiding_A(x1) + u1*B1.normalized();
    return (0.5*(Adag1+Adag0).dot(x1-x0) -h*mu*B0.norm()-h*(u0*u1)/2.);
  }
#endif


Matrix<double,DIM,1> discrete_lagrangian_central_derivative(double hx,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1,bool dervar){
  //compute central derivative of discrete lagrangian
  //dervar=0 --> derivative of first variable
  //dervar=1 --> derivative of second variable

  Matrix<double,DIM,1> q1_d0,q1_d1,q0_d0,q0_d1,ret;
  
  for (int i=0;i<DIM;i++){
    q1_d0=q1_d1=q1;
    q0_d0=q0_d1=q0;

    if (dervar==1){
      q1_d1(i)+=hx;
      q1_d0(i)-=hx;
    }
    else{ 
      q0_d1(i)+=hx;
      q0_d0(i)-=hx;
    }

    ret(i) = 0.5/hx * (discrete_lagrangian(h,q0_d1,q1_d1)-discrete_lagrangian(h,q0_d0,q1_d0));
  }
  return ret;
}

void legendre_right(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1, Matrix<double,DIM,1> &p1){
  #ifndef richardson
    p1 = discrete_lagrangian_central_derivative(hx,q0,q1,1);
  #else
    p1 = (4.*discrete_lagrangian_central_derivative(hx,q0,q1,1)-discrete_lagrangian_central_derivative(2.*hx,q0,q1,1))/3.;
  #endif
}
void legendre_left(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1, Matrix<double,DIM,1> &p0){
  #ifndef richardson
    p0 = -discrete_lagrangian_central_derivative(hx,q0,q1,0);
  #else
    p0 = -(4.*discrete_lagrangian_central_derivative(hx,q0,q1,0)-discrete_lagrangian_central_derivative(2.*hx,q0,q1,0))/3.;
  #endif
}


void legendre_left_inverse_newton(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> p0, Matrix<double,DIM,1> &q1){
  Matrix<double,DIM,1> Q,q1_d0,q1_d1,f_d0,f_d1,f,p;
  Matrix<double,DIM,DIM> Jf;
  
  legendre_left(h,q0,q1,p);
  f = p0-p;
  
  for (int j=0;j<DIM;j++){
    q1_d0 = q1_d1 = q1;
    q1_d0(j)-=hx;
    q1_d1(j)+=hx;
      
    legendre_left(h,q0,q1_d1,p);
    f_d1 = p0-p;
    legendre_left(h,q0,q1_d0,p);
    f_d0 = p0-p;
    Jf.col(j) = 0.5*(f_d1 - f_d0)/hx;
  }
  q1 = q1 - Jf.inverse()*f;


  if (DEBUG_IMPLICIT_F) cout << "f: " << f.transpose() << endl;
}

void legendre_right_inverse_newton(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> p1, Matrix<double,DIM,1> &q1){
  Matrix<double,DIM,1> Q,q1_d0,q1_d1,f_d0,f_d1,f,p;
  Matrix<double,DIM,DIM> Jf;
  
  legendre_right(h,q0,q1,p);
  f = p1-p;
  
  for (int j=0;j<DIM;j++){
    q1_d0 = q1_d1 = q1;
    q1_d0(j)-=hx;
    q1_d1(j)+=hx;
      
    legendre_right(h,q0,q1_d1,p);
    f_d1 = p1-p;
    legendre_right(h,q0,q1_d0,p);
    f_d0 = p1-p;
    Jf.col(j) = 0.5*(f_d1 - f_d0)/hx;
  }
  q1 = q1 - Jf.inverse()*f;
}


//NON SYMPLECTIC STEPS
void step_euler(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> p0, Matrix<double,DIM,1> &q1,Matrix<double,DIM,1> &p1){
  Matrix<double,2*DIM,1> k1,k2,k3,k4,z0,z1;
  
  z0.head(DIM) = q0;
  z0.tail(DIM) = p0;
  
  f_eq_motion(z0,k1);
  
  z1 = z0 + h*k1;
  
  q1 = z1.head(DIM);
  p1 = z1.tail(DIM);
  
}

void step_RK4(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> p0, Matrix<double,DIM,1> &q1,Matrix<double,DIM,1> &p1){
  Matrix<double,2*DIM,1> k1,k2,k3,k4,z0,z1;
  
  z0.head(DIM) = q0;
  z0.tail(DIM) = p0;
  
  f_eq_motion(z0,k1);
  f_eq_motion(z0+0.5*h*k1,k2);
  f_eq_motion(z0+0.5*h*k2,k3);
  f_eq_motion(z0+h*k3,k4);
  
  z1 = z0 + 1./6.*h*(k1+2.*k2+2.*k3+k4);
  
  q1 = z1.head(DIM);
  p1 = z1.tail(DIM);
  
}


void step_symplectic_lagrangian_momentum(double h,Matrix<double,DIM,1> q0,Matrix<double,DIM,1> p0, Matrix<double,DIM,1> &q1,Matrix<double,DIM,1> &p1){
  legendre_left_inverse_newton(h,q0,p0,q1);
  legendre_right(h,q0,q1,p1);
}

#ifdef LYAPUNOV
#ifndef guidingcenter3D
    //COMPUTE LYAPUNOV COEFFICIENTS AND THE MAXIMUM LYAPUNOV COEFFICIENT (eigenvalues method)
    void compute_lyapunov_exponents(Matrix<double,DIM,1> q0,Matrix<double,DIM,1> q1, Matrix<double,DIM,1> p0, Matrix<double,DIM,1> p1, Matrix<double,2*DIM,2*DIM> &Jproduct,int t){
        Matrix<double,2*DIM,2*DIM> dz1;
        Matrix<double,2*DIM,1> lambda_lyap;
        double maximum_lambda;

        for (int i=0;i<2*DIM;i++) {
          //COMPUTE THE JACOBIAN
          Matrix<double,DIM,1> q0l,p0l,q1l,p1l;

          q0l = q0;
          p0l = p0;
          if (i<DIM) q0l(i)+=hlyap;
          else p0l(i-DIM)+=hlyap;

          step_guess(h,q0l,p0l,q1l,p1l);
          for (int j=0;j<NEWTON;j++) step_symplectic_lagrangian_momentum(h,q0l,p0l,q1l,p1l);

          dz1.block(0,i,DIM,1) = (q1l-q1)/hlyap;
          dz1.block(DIM,i,DIM,1) = (p1l-p1)/hlyap;

        }

        //COMPUTE PRODUCT OF JACOBIANS (PHI)
        if (t==(TIME_OFFSET+1)) Jproduct = dz1;
        else{
          //avoid aliasing
          Matrix<double,2*DIM,2*DIM> temp = dz1*Jproduct;
          Jproduct = temp;
        }

        //COMPUTE PHI' * PHI
        MatrixXd Jtot = Jproduct.transpose()*Jproduct;
        EigenSolver<MatrixXd> es(Jtot);

        for (int i=0;i<(2*DIM);i++) lambda_lyap(i) = log(abs(es.eigenvalues()[i]))/(2.*t);
        maximum_lambda = lambda_lyap.maxCoeff();

        if (((t+PRINT_OFFSET)%PRINT_MULTIPLE)==0){
          lyap << (t-1) << " " << (t-1)/ORBIT_NORMALIZE << " " << maximum_lambda << " " << lambda_lyap.transpose() << endl;
        }
    }

#endif
#endif




int main(int argc, char* argv[]){
  
  cout.setf(std::ios::scientific);
  cout.precision(PRINT_PRECISION);
  out.setf(std::ios::fixed);
  out.precision(PRINT_PRECISION);
  out.setf(ios::showpoint);
  
  Matrix<double,DIM,1> q0,p0,q1,p1,adag;
  double E0; //initial energy

  initial_conditions(q0,p0);

  //initialization with an auxiliary integrator (~lagrangian init)
  #ifndef step3D
    Matrix<double,DIM,1> qt0=q0,pt0=p0;
    for (int i=0;i<INIT_STEPS;i++){
      if (INIT_STEPS!=0) step_RK4(h/double(INIT_STEPS),qt0,pt0,q1,p1);
      qt0 = q1;
      pt0 = p1;
    }
     if (INIT_STEPS>0) legendre_left(h,q0,q1,p0);
     #ifdef explicit
      if ((INIT_STEPS>0) && (NEWTON==0)){
        legendre_right_explicit(h,q0,q1,p1);
        q0 = q1;
        p0 = p1;
      } 
     #endif

  #else
    //3D (implicit2) init

    if (NEWTON==0){
        cout << "Please don't use guiding center implicit scheme 2 without newton loops! Use a 4D explicit scheme instead!" << endl;
        return 0;
    } 

    Vector4d q3t0,p3t0,q3t1,p3t1;
    q3t0.head(3) = q0;
    q3t0(3) = (p0- guiding_A(q0)).norm();

    for (int i=0;i<INIT_STEPS;i++){
      if (INIT_STEPS!=0) step_RK4_3D(h/double(INIT_STEPS),q3t0,p3t0,q3t1,p3t1);
      q3t0 = q3t1;
      p3t0 = p3t1;
    }
    q1 = q3t0.head(3);

    if (INIT_STEPS>0) legendre_left(h,q0,q1,p0);
  #endif


  //compute initial energy
  E0 = energy(q0,p0);
  
  //initialize Lyapunov exponents
  #ifdef LYAPUNOV
    Matrix<double,2*DIM,2*DIM> lyap_Jproduct;
  #endif


  cout << "time step: " << h << endl;
  cout << "Initialization: " << endl;
  cout << "q0:\t" << q0.transpose() << endl;
  cout << "p0:\t" << p0.transpose() << endl;

  

  // ******
  // ******
  //MAIN LOOP
  // ******
  // ******
 
  for (int t=TIME_OFFSET + 1;t<=MAX_T+TIME_OFFSET;t++){

    if ((DEBUG_TIMESTEP_MULT>0) && (t%DEBUG_TIMESTEP_MULT==0)) cout << "Timestep " << t << endl;

    #ifndef step3D
      step_guess(h,q0,p0,q1,p1);
    #else 

      //3D step for guiding center - explicit (implicit scheme 2)
      // (we have to adapt the explicit 4D schemes to our 3D system)
      double u = (p0- guiding_A(q0)).norm();
      Vector4d q30, q31,p30, p31;
      q30.head(3) = q0;
      q30(3) = u;
      p30.head(3) = guiding_A(q0) + u*guiding_B(q0).normalized();
      p30(3) = 0;
      step_explicit_guidingcenter_momentum(h,q30,p30,q31,p31);
      q1 = q31.head(3);
    #endif

    //implicit step
    for (int i=0;i<NEWTON;i++) step_symplectic_lagrangian_momentum(h,q0,p0,q1,p1);


    #ifdef LYAPUNOV
    #ifndef guidingcenter3D
      compute_lyapunov_exponents(q0,q1,p0,p1,lyap_Jproduct,t);
    #endif
    #endif

    if (t==1) {
      cout << "q1:\t" << q1.transpose() << endl;
      cout << "p1:\t" << p1.transpose() << endl;
    }
    
    #ifdef compute_momentum
      momentum(q0,adag);
    #endif
    
    //EXIT IF THE ERROR IS TOO HIGH
    double energy_error = (energy(q0,p0)-E0)/E0;
    if ((EXIT_ON_ERROR) && (abs(energy_error)>ERROR_THRESHOLD)){
      cout << "Timestep: " << t << "\tError too high! Exiting." << endl;
      break;
    }
    
    //PRINT TO FILE
    if (((t+PRINT_OFFSET)%PRINT_MULTIPLE)==0){
      out << (t-1) << " " << (t-1)/ORBIT_NORMALIZE << " " << energy_error
      << " " << q0.transpose() << " " << p0.transpose() 
      << " " << (p0-adag).transpose();
      out << endl;
    }

    //SHIFT
    q0 = q1;
    p0 = p1;

  }
  
  return 0;
}


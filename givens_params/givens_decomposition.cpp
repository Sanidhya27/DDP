#include "givens_decomposition.hpp"

#define I std::complex<double>(0,1)

using namespace itpp;
using std::cout;
using std::endl;

cmat given_matrix(double theta,int a,int b,int s) {
	cmat G = eye_c(s);
	G(a,a) = cos(theta);
	G(a,b) = -sin(theta);
	G(b,a) = sin(theta);
	G(b,b) = cos(theta);
	return G;
}
cmat givens_reconstruction(vec phis, vec thetas,int t,int n){
	cmat M = eye_c(t);
	int q=0,w=0,L=0;
	vec theta,phi;
	cmat G;
	for(int j=0;j<n;j++){
		phi = phis(q,q+t-j-1);
		phi = concat(zeros(j),phi);
		M=M*diag(exp(I*phi));
		q=q+t-j;

		theta = thetas(w,w+t-j-2);
		L=theta.size();
		w=w+t-j-1;
		for(int l=0;l<L;l++){
			G = given_matrix(theta(l),t-l-2,t-l-1,t);
			M=M*G;
		}
	}
	cmat I_ = eye_c(n);
	for(int i=0;i<t-n;i++){
		I_.append_col(zeros_c(n));
	}
	M=M*transpose(I_);
	return M;
}
GIVENSPARAMS givens_decomposition(cmat V){
	vec thetas;
	vec phis;
	cmat D = diag(exp(-1*I*angle(V.get_row(0))));
	// V = V*D;
	int t = V.rows();
	int n = V.cols();
	vec phi;
	double theta;
	cmat D_,G;
	for(int j=0;j<n;j++){
		phi = angle(V.get_col(j)(j,t-1));
		phis=concat(phis,phi);
		phi = concat(zeros(j),phi);
		D_ = diag(exp(-1*I*phi));
		V = D_*V;
		double c,s;
		for(int k=0;k<t-j-1;k++){
			givens(real(V(t-k-2,j)),real(V(t-k-1,j)),c,s);
			theta = atan2(-s,c);
			thetas = concat(thetas,theta);
			G = eye_c(t);
			G(t-k-2,t-k-2) = c;
			G(t-k-2,t-k-1) = -s;
			G(t-k-1,t-k-2) = s;
			G(t-k-1,t-k-1) = c;
			V = G*V;
		}

	}
	GIVENSPARAMS par;
	par.t=thetas;
	par.p=phis;
	return par;

}
// int main(int argc, char const *argv[])
// {
// 	int t=4,r=4;
// 	cmat H = randn_c(r,t);
// 	cmat U,V;

// 	vec s;
// 	bool check = svd(H,U,s,V);
// 	int rank = s.size();
// 	V = V(0,t-1,0,rank-1); //economy mode svd
// 	struct params par = givens_decomposition(V);
// 	cmat M = givens_reconstruction(par.p,par.t,t,rank);
// 	cout<<V<<endl;
// 	cout<<M<<endl;
// 	return 0;
// }

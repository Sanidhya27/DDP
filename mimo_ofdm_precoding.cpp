#include <itpp/itcomm.h>
#include <bits/stdc++.h>
using namespace itpp;
#include <fstream>
using std::cout;
using std::cin;
using std::endl;
#include <vector>
#define I std::complex<double>(0,1)
struct params{
	vec t;
	vec p;
};
cmat given_matrix(double theta,int a,int b,int s){
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
params givens_decomposition(cmat V){
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
	struct params par;
	par.t=thetas;
	par.p=phis;
	return par;

}
cmat generate_channel(int L,int t,int r){
  cmat h;
  for(int i=0;i<t*r;i++){
    h.append_row(randn_c(L));
  }
  return h;
}
cmat take_fft(cmat h,int IFFT_size){
  cmat H;
  for(int i=0;i<h.rows();i++){
    H.append_row(fft(h.get_row(i),IFFT_size));
  }
  return H;
}
cmat get_precoder(int t){
	return eye_c(t);
}
int main()
{
  std::ofstream f;
  f.open("mimo_ofdm.txt");

  double N0;
  int N = 64000; 

  int L = 5;

  int t,r;
  cout<<"No. of transmitter"<<endl;
  cin>>t;
  cout<<"No. of receiver"<<endl;
  cin>>r;
  RNG_randomize();

  bmat bits(t,N),dec_bits(t,N);
  QPSK qpsk;  
  BERC berc,berc2;
  OFDM ofdm;

  int CP_size=6;
  int IFFT_size = 64;

  // N=4;
  // L  = 1;
  // CP_size = 2;
  // IFFT_size = 2;

  ofdm.set_parameters(IFFT_size,CP_size,1);  

  vec snrs_dB="0:2:20";
  vec snr_linear = inv_dB(snrs_dB);

  bits = randb(t,N); // each transmitter send N bits; 
  cvec temp;
  cmat h;
  cmat H;
  cmat V;
  // bits = "1 0 0 0 ;1 1 0 0";
  for (int l = 0; l < snrs_dB.length(); l++) {
    berc.clear();
    N0 = 1 / sqrt(snr_linear(l));     
    cmat symbols,ofdm_out,rec;
    dec_bits = zeros_b(t,N);
    // cout<<"bits"<<bits<<endl;
    for(int i=0;i<t;i++){
      qpsk.modulate_bits(bits.get_row(i), temp);
      symbols.append_row(temp);
    }
    
    V = get_precoder(t);
    symbols = V*symbols;
    for(int i=0;i<t;i++){
      ofdm_out.append_row(ofdm.modulate(symbols.get_row(i)));
    }
    // cout<<"symbols"<<symbols<<endl;
    // cout<<"ofdm_out"<<ofdm_out<<endl;
    int z = IFFT_size+CP_size;
    int x = ofdm_out.cols()/z;    
    // h = ones_c(r*t,L);
    // h.set_row(1, 1*ones_c(L));
    // h.set_row(2, 3*ones_c(L));


    // cout<<"h"<<h<<endl;
    for(int k=0;k<x;k++){
      	// h = "0.449311+1.0786i; 0.0145245+0.010007i";
  		// h="0.449311+1.0786i; 0.0145245+0.010007i; -0.337605+0.778781i;-0.136933+0.0569123i";
      h = generate_channel(L,t,r);
      rec.del_rows(0,rec.rows()-1);
      for(int i=0;i<r;i++){
          temp = zeros_c(z);
          for(int j=0;j<t;j++){
            temp+=filter(h.get_row(i*t + j),1,ofdm_out.get_row(j)(k*z,(k+1)*z -1));
          }
          rec.append_row(temp+1*sqrt(N0) * randn_c(temp.length()));
      }

      H = take_fft(h,IFFT_size); 
      cmat A; 
      cmat dec;
      for(int i=0;i<rec.rows();i++){
        rec.set_row(i,ofdm.demodulate(rec.get_row(i)));
      }
      
      rec.del_cols(IFFT_size,IFFT_size+CP_size-1);
      for(int i=0;i<rec.cols();i++){
        A = transpose(reshape(H.get_col(i),t,r));
        if(r>=t)
          dec.append_col(inv(hermitian_transpose(A)*A)*hermitian_transpose(A)*rec.get_col(i));
        else
          dec.append_col(hermitian_transpose(A)*inv(A*hermitian_transpose(A))*rec.get_col(i));
      } 
      bvec temp2;
      int len;
     for(int i=0;i<dec.rows();i++){
      qpsk.demodulate_bits(dec.get_row(i), temp2);
      len = temp2.length();
      dec_bits.set_submatrix(i,k*len,reshape(temp2,1,len));
     }
    }
    berc.clear();
    for(int i=0;i<bits.rows();i++){
      berc.count(bits.get_row(i), dec_bits.get_row(i));
    }
  
    cout << "There were " << berc.get_errors() << " received bits in error." << endl;
    cout << "There were " << berc.get_corrects() << " correctly received bits." << endl;
    cout << "The error probability was " << berc.get_errorrate() << endl;
    f<<berc.get_errorrate() << "\t" << snrs_dB[l]<<endl;
  } 
  f.close();
  
   return 0;
}


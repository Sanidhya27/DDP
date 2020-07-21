#include <fstream>
#include "givens_decomposition.hpp"

using namespace itpp;
using std::cout;
using std::cin;
using std::endl;

cmat generate_channel(int L,int t,int r) {
  cmat h;
  for(int i=0;i<t*r;i++){
    h.append_row(randn_c(L));
  }
  return h;
}

cmat take_fft(const cmat &h,int IFFT_size) {
  cmat H;
  for(int i=0;i<h.rows();i++){
    H.append_row(fft(h.get_row(i),IFFT_size));
  }
  return H;
}

cmat get_precoder(int t) {
	return eye_c(t);
}

int main(int argc, char *argv[])
{
  std::ofstream f;
  f.open("mimo_ofdm.txt");

  double N0;
  int N = 64000;

  int L = 5;

  int t,r;
  cout<<"No. of transmitters: "<<endl;
  cin>>t;
  cout<<"No. of receivers: "<<endl;
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
    int z = IFFT_size+CP_size;
    int x = ofdm_out.cols()/z;
    for(int k=0;k<x;k++){
      h = generate_channel(L,t,r);
      rec.del_rows(0,rec.rows()-1);
      for(int i=0;i<r;i++){
          temp = zeros_c(z);
          for(int j=0;j<t;j++){
            temp+=filter(h.get_row(i*t + j),1,ofdm_out.get_row(j)(k*z,(k+1)*z -1));
          }
          rec.append_row(temp+1*sqrt(N0) * randn_c(temp.length()));
      }

      H = take_fft(h, IFFT_size);
      cmat A;
      cmat dec;
      for(int i=0;i<rec.rows();i++){
        rec.set_row(i,ofdm.demodulate(rec.get_row(i)));
      }

      rec.del_cols(IFFT_size,IFFT_size+CP_size-1);
      for(int i=0;i<rec.cols();i++){
        A = transpose(reshape(H.get_col(i),t,r));
        // Look at Octave pinv function
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

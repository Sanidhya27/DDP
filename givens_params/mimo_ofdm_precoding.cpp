#include <fstream>
#include "givens_decomposition.hpp"
#define I std::complex<double>(0,1)

// h -> array of length L, each element is a r x t matrix
// H -> array of length NFFT, each element is a r x t matrix
// To get H, we use the formular
// H[k] = sum (n = 0 to (L - 1)) h[n] exp(-1j * 2 * pi / NFFT * n * k)
// Then, precoder for k-th subcarrier is obtained using svd of H[k]

using namespace itpp;
using std::cout;
using std::cin;
using std::endl;

void generate_channel(cmat *h, int L,int t,int r) {
    for(int i=0;i<L;i++){
        h[i] = randn_c(r, t);
    }
}

void take_fft(cmat *H, const cmat *h, int L,int IFFT_size) {
    for (int i = 0; i < IFFT_size; i++) {
        H[i] = zeros_c(h[0].rows(), h[0].cols());
        for (int j = 0; j < L; ++j) {
            H[i] += h[j] * exp(-1*I * 2 * M_PI * j*i / IFFT_size );
        }
        H[i]/=IFFT_size;
    }
}

cmat get_svd_precoder(cmat &A) {
    cmat U,V;
    vec s;
    svd(A,U,s,V);
    // return eye_c(V.rows());
    return V;
}

cmat get_svd_postcoder(cmat &A) {
    cmat U,V;
    vec s;
    svd(A,U,s,V);
    return eye_c(U.rows());
    // return U;
}

void filt(cmat &rec, cmat *h, cmat &ofdm_syms, int L, int ofdm_sym_size, int t, int r) {
    cvec temp;
    for(int i = 0; i < ofdm_sym_size; i++) {
        temp = zeros_c(r);
        for(int j = 0; j < L; j++){
            if( i-j >= 0)
                temp+= h[j] * ofdm_syms.get_col(i-j);
        }
        rec.set_col(i,temp);
    }
}

cmat pinv(cmat &A){
    
    return inv(hermitian_transpose(A)*A)*hermitian_transpose(A);
}

int main(int argc, char const *argv[])
{
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <t> <r>\n";
        return 1;
    }
    int t,r;
    t = atoi(argv[1]);
    r = atoi(argv[2]);
    cout<<"No. of transmitters: " << t << "\n";
    cout<<"No. of receivers: " << r << "\n";

    std::ofstream f;
    f.open("mimo_ofdm.txt");
    
    int N_ofdm_syms = 1000;
    int N_taps = 5;
    int IFFT_size = 64;
    int CP_size = 6;
    int constellation_size = 2;

    QPSK qpsk;
    BERC berc;
    OFDM ofdm;

    ofdm.set_parameters(IFFT_size,CP_size,1);

    vec snrs_dB="0:2:20";
    vec snr_linear = inv_dB(snrs_dB);

    bmat bits[N_ofdm_syms];
    cmat qpsk_syms = zeros_c(t,IFFT_size);
    cmat precoded_syms = zeros_c(t,IFFT_size);
    cmat ofdm_syms = zeros_c(t,IFFT_size + CP_size);
    cmat rec_syms = zeros_c(r,IFFT_size + CP_size);
    cmat dec_ofdm_syms = zeros_c(r,IFFT_size);
    cmat dec_syms = zeros_c(t,IFFT_size);
    bmat dec_bits[N_ofdm_syms];

    cmat h[N_taps];
    cmat H[IFFT_size];

    cmat V,U;
    cvec temp;
    bvec temp2;

    for(int snr_i = 0; snr_i < snrs_dB.length(); snr_i++) {
        
        berc.clear();

        for(int sym_i = 0; sym_i < N_ofdm_syms; sym_i++) {
            
            bits[sym_i] = randb(t,constellation_size*IFFT_size);

            for(int i = 0; i < t; i++) {
                qpsk.modulate_bits(bits[sym_i].get_row(i), temp);
                qpsk_syms.set_row(i, temp);
            }
            generate_channel(h,N_taps,t,r);
            take_fft(H,h,N_taps,IFFT_size);

            for(int i = 0; i < IFFT_size; i++) {
                V = get_svd_precoder(H[i]);
                precoded_syms.set_col(i,V*qpsk_syms.get_col(i));
            }

            for(int i = 0; i < t; i++) {
                ofdm_syms.set_row(i, ofdm.modulate(precoded_syms.get_row(i)));
            }

            filt(rec_syms, h,ofdm_syms,N_taps,IFFT_size + CP_size,t,r);

            rec_syms += 1*randn_c(r,IFFT_size+CP_size)/sqrt(snr_linear(snr_i));

            for(int i = 0; i < r; i++) {
                dec_ofdm_syms.set_row(i, ofdm.demodulate(rec_syms.get_row(i)));
            } 

            for(int i = 0; i < IFFT_size; i++) {
                vec s;
                svd(H[i],U,s,V);
                cmat A = hermitian_transpose(U)*H[i]*V;
                dec_syms.set_col(i,pinv(A)*hermitian_transpose(U)*dec_ofdm_syms.get_col(i));
            }

            dec_bits[sym_i] = zeros_b(t,IFFT_size*constellation_size);

            for(int i = 0; i < t; i++) {
                qpsk.demodulate_bits(dec_syms.get_row(i),temp2);
                dec_bits[sym_i].set_row(i, temp2);
                berc.count(bits[sym_i].get_row(i),temp2);
            }
        }
        cout << "There were " << berc.get_errors() << " received bits in error." << endl;
        cout << "There were " << berc.get_corrects() << " correctly received bits." << endl;
        cout << "The error probability was " << berc.get_errorrate() << endl;
        f<<berc.get_errorrate() << "\t" << snrs_dB[snr_i]<<endl;

    }
    f.close();
    return 0;
}
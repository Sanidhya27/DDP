#include <fstream>
#include <iostream>
#include <string>
#include "givens_decomposition.hpp"
#include "utils.hpp"
#define I std::complex<double>(0,1)

using namespace itpp;
using std::cout;
using std::cin;
using std::endl;
using std::string;

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

void quantise_phi(vec &p, int n_bits){

    // For n bits divide -3.14 to 3.14 in 2^n regions
    vec reg = linspace(-1*M_PI, M_PI, pow2(n_bits));
    int index;
    for(int i = 0; i < p.length(); i++){
        min(abs(reg-p(i)), index);
        p(i) = reg(index);
    }
}

void quantise_theta(vec &theta,int n_bits,int t){

    it_file ff;
    ff.open(std::to_string(n_bits) + "bits.it");
    mat codes;
    ff >> Name("c") >> codes;
    int index;
    int j=0;
    for(int i=0;i<t-1;i++){
        for(int l=0;l<t-1-i;l++){
            min(abs(codes.get_col(l)-theta(j)), index);
            theta(j) = codes.get_col(l)(index);
            j++;
        }
    }
}

GIVENSPARAMS quantise_params(GIVENSPARAMS par,int n_bits,int t){

    quantise_phi(par.p, n_bits);
    quantise_theta(par.t,n_bits,t);
    return par;
}

cmat get_quantised_precoder(const cmat *H, int i, int n_bits, cmat *Vmat, int Vsize, int freq_inter=0){

    if(Vmat[i].rows()!=0){
        return Vmat[i];
    }

    if(freq_inter && i%(freq_inter+1)!=0 ){
        int j = i - i%(freq_inter+1) + freq_inter+1;
        if(j==(Vsize))
            Vmat[i] = Vmat[i-1];

        else {
            int j = i - i%(freq_inter+1) + freq_inter+1;
            Vmat[ j ] = Vmat[j].rows()!=0 ? Vmat[j] : get_quantised_precoder(H,j,n_bits,Vmat,Vsize);
            cmat Sip1, Sim1;
            logm(Vmat[j], Sip1);
            j = i - i%(freq_inter+1); 
            logm(Vmat[j], Sim1);
            j = freq_inter+1;
            expm(((j-i%j)*Sim1 + i%(j)*Sip1)/j, Vmat[i]);
        }

        return Vmat[i];
    }

    else{

        int t = H[i].cols();
        cmat U,V;
        vec s;

        svd(H[i],U,s,V);
        GIVENSPARAMS par = givens_decomposition(V);
        par = quantise_params(par,n_bits,t);
        Vmat[i] = givens_reconstruction(par.p,par.t,t,t);
        return Vmat[i];
    }

}

double capacity(mat &F, int t){
    double sum=0;
    for(int k=0;k<t;k++){
        sum+=log2(1+1/(F(k,k)));
    }
    return sum;
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
    return U;
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

int main(int argc, char const *argv[]) {

    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <t> <r> --optional <quant> <bits> <interpolate>\n";
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

    cmat V,U,A;
    mat F;
    cvec temp;
    bvec temp2;
    double cap;

    bool quantise = false;
    int freq_inter = 0;
    int n_bits=0;

    if (argc == 4) {
        quantise = atoi(argv[3]);
        if(quantise){
            cout << "Usage: " << argv[0] << " <t> <r> --optional <quant> <bits> <interpolate>\n";
            cout<<"Enter number of bits to quantise to and do you need interpolation\n";
            return 1;
        }
    }

    else if (argc == 6) {
        quantise = atoi(argv[3]);
        n_bits = atoi(argv[4]);
        freq_inter = atoi(argv[5]);
    }
    std::ofstream f1;
    f1.open("capacity.txt");

    RNG_randomize();

    for(int snr_i = 0; snr_i < snrs_dB.length(); snr_i++) {

        berc.clear();
        cap = 0;
        
        for(int sym_i = 0; sym_i < N_ofdm_syms; sym_i++) {

            cmat Vmat[IFFT_size];

            bits[sym_i] = randb(t,constellation_size*IFFT_size);

            for(int i = 0; i < t; i++) {
                qpsk.modulate_bits(bits[sym_i].get_row(i), temp);
                qpsk_syms.set_row(i, temp);
            }

            generate_channel(h,N_taps,t,r);
            take_fft(H,h,N_taps,IFFT_size);

            for(int i = 0; i < IFFT_size; i++) {

                if(quantise)
                    V = get_quantised_precoder(H, i, n_bits, Vmat, IFFT_size, freq_inter);
                else
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
                U = get_svd_postcoder(H[i]);
                if(quantise)
                    V = Vmat[i];
                    // V = get_quantised_precoder(H, i, n_bits, Vmat, IFFT_size, freq_inter);
                else
                    V = get_svd_precoder(H[i]);
                A = hermitian_transpose(U)*H[i]*V;
                F = real(inv(snr_linear(snr_i)*hermitian_transpose(A)*A));
                cap+=capacity(F,t);
                dec_syms.set_col(i,pinv(A)*hermitian_transpose(U)*dec_ofdm_syms.get_col(i));
            }

            dec_bits[sym_i] = zeros_b(t,IFFT_size*constellation_size);

            for(int i = 0; i < t; i++) {
                qpsk.demodulate_bits(dec_syms.get_row(i),temp2);
                dec_bits[sym_i].set_row(i, temp2);
                berc.count(bits[sym_i].get_row(i),temp2);
            }
        }
        cap/=N_ofdm_syms;
        cout<<"capacity "<<cap<<endl;
        cout << "There were " << berc.get_errors() << " received bits in error." << endl;
        cout << "There were " << berc.get_corrects() << " correctly received bits." << endl;
        cout << "The error probability was " << berc.get_errorrate() << endl;
        f<<berc.get_errorrate() << "\t" << snrs_dB[snr_i]<<endl;
        f1<<cap<<endl;
        it_file ff;
    ff.open("V.it");
   
    cmat Vmat3[IFFT_size];
    // ff<<Name("IIV0") << Vmat[0];
    // ff<<Name("IIV1") << Vmat[1];
    // ff<<Name("IIV2") << Vmat[2];
    get_quantised_precoder(H, 0,n_bits,Vmat3, IFFT_size, true);
    get_quantised_precoder(H, 1, n_bits,Vmat3, IFFT_size, true);
    get_quantised_precoder(H, 2, n_bits,Vmat3, IFFT_size, true);
    ff<<Name("IV0") << Vmat3[0];
    ff<<Name("IV1") << Vmat3[1];
    ff<<Name("IV2") << Vmat3[2];
    cmat Vmat2[IFFT_size];
    get_quantised_precoder(H, 0,n_bits,Vmat2, IFFT_size, false);
    get_quantised_precoder(H, 1, n_bits,Vmat2, IFFT_size, false);
    get_quantised_precoder(H, 2, n_bits,Vmat2, IFFT_size, false);
    ff<<Name("QV0") << Vmat2[0];
    ff<<Name("QV1") << Vmat2[1];
    ff<<Name("QV2") << Vmat2[2];
    ff<<Name("V0") <<get_svd_precoder(H[0]);
    ff<<Name("V1") <<get_svd_precoder(H[1]);
    ff<<Name("V2") <<get_svd_precoder(H[2]);
    ff.flush();
    ff.close();
    }
    
    f.close();
    return 0;
}

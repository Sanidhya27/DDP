#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include "givens_decomposition.hpp"
#include "mimo-tdl-channel.hpp"
#include "utils.hpp"
#include <assert.h> 

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
void delta_quantise(vec &p,vec &quant, vec &state, vec &delta){

    for(int i =0;i<p.length();i++){

        if(quant(i) < p(i)){
            if(state(i)==1)
                delta(i) = delta(i)*2;
            else
                delta(i)=delta(i)/2;
            state(i) = 1;
            // cout<<"yo"<<quant(i)<<" "<<delta(i)<<endl;
            quant(i) = quant(i) + delta(i);
            // cout<<"yo2"<<quant(i)<<endl;
        }
        else if(quant(i) > p(i)){
            if(state(i)==0)
                delta(i) = delta(i)*2;
            else
                delta(i) = delta(i)/2;
            state(i) = 0;   
            quant(i) = quant(i) - delta(i);    
        }
    }
}
cmat get_delta_quantised_precoder(const cmat *H, int i, cmat *Vmat, GIVENSPARAMS &quant, vec &p_delta,vec &p_state, vec &t_delta, vec &t_state){

    int t = H[i].cols();
    cmat U,V;
    vec s;

    svd(H[i],U,s,V);
    GIVENSPARAMS par = givens_decomposition(V);

    if(quant.p.length()==0){
        quant.p = zeros(par.p.length());
        p_delta = 0.25*ones(par.p.length());
        p_state = zeros(par.p.length());
    }
    if(quant.t.length()==0){
        quant.t = zeros(par.t.length());
        t_delta = 0.25*ones(par.t.length());
        t_state = zeros(par.t.length());
    }

    delta_quantise(par.p,quant.p,p_state,p_delta);
    delta_quantise(par.t,quant.t,t_state,t_delta);
    // cout<<"paramss"<<endl;
    // cout<<par.p<<endl;
    // cout<<quant.p<<endl;
    // cout<<par.t<<endl;
    // cout<<quant.t<<endl;
    Vmat[i] = givens_reconstruction(quant.p,quant.t,t,t);
    return Vmat[i];
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
    int CP_size =6;
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

        MIMO_TDL_Channel mimo_tdl_channel(r,t);
        // mimo_tdl_channel.set_channel_profile_uniform (4);
        mimo_tdl_channel.set_channel_profile(itpp::ITU_Pedestrian_A, 5e-8);
        mimo_tdl_channel.set_norm_doppler(1e-4);
        N_taps = mimo_tdl_channel.taps();
        assert(N_taps < CP_size);

        GIVENSPARAMS quant[IFFT_size];
        vec p_delta[IFFT_size],p_state[IFFT_size],t_state[IFFT_size],t_delta[IFFT_size];

        cmat h[N_taps];
        cmat H[IFFT_size];
        berc.clear();
        cap = 0;
        
        for(int sym_i = 0; sym_i < N_ofdm_syms; sym_i++) {

            cmat Vmat[IFFT_size];

            bits[sym_i] = randb(t,constellation_size*IFFT_size);

            for(int i = 0; i < t; i++) {
                qpsk.modulate_bits(bits[sym_i].get_row(i), temp);
                qpsk_syms.set_row(i, temp);
            }

            mimo_tdl_channel.generate(h);
            take_fft(H,h,N_taps,IFFT_size);
            // cout<<N_taps<<" "<<h[0]<<" "<<H[0]<<endl;
            for(int i = 0; i < IFFT_size; i++) {

                if(quantise)
                    V = get_delta_quantised_precoder(H,i,Vmat,quant[i],p_delta[i],p_state[i],t_delta[i],t_state[i]);
                else
                    V = get_svd_precoder(H[i]);

                precoded_syms.set_col(i,V*qpsk_syms.get_col(i));
                // cout<<quantise<<i<<" "<<V<<"\n"<<get_svd_precoder(H[i])<<endl;
            }
            // cout<<"test"<<get_svd_precoder(H[IFFT_size-1])(0,0)<<" "<<V(0,0)<<endl;

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
    }
    
    f.close();
    return 0;
}

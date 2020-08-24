#include <itpp/itcomm.h>

class MIMO_TDL_Channel
{
public:
    MIMO_TDL_Channel(size_t Nt, size_t Nr, const itpp::vec &avg_power_dB="0", const itpp::ivec &delay_prof="0");
    ~MIMO_TDL_Channel();
    void set_channel_profile(const itpp::vec &avg_power_dB, const itpp::ivec &delay_prof);
    void set_channel_profile (const itpp::Channel_Specification &channel_spec, double sampling_time);
    void set_channel_profile_uniform (int no_taps);
    void set_channel_profile_exponential(int no_taps);
    void generate(itpp::cmat *channel_coeffs);
    void set_norm_doppler (double norm_doppler);
    void get_norm_doppler ();
    int taps() const {
        return tdl_channels[0].taps();
    }
private:
    itpp::TDL_Channel *tdl_channels;
    size_t Nt, Nr;
};

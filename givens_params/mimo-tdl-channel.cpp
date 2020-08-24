#include <itpp/itcomm.h>
#include <algorithm>
#include "mimo-tdl-channel.hpp"

MIMO_TDL_Channel::MIMO_TDL_Channel(size_t Nt,
                                   size_t Nr,
                                   const itpp::vec &avg_power_dB,
                                   const itpp::ivec &delay_prof)
{
    tdl_channels = new itpp::TDL_Channel[Nt * Nr];
    for (size_t i = 0; i < Nt * Nr; ++i)
        tdl_channels[i] = itpp::TDL_Channel();
    this->Nt = Nt;
    this->Nr = Nr;
}

MIMO_TDL_Channel::~MIMO_TDL_Channel()
{
    delete [] tdl_channels;
}

void
MIMO_TDL_Channel::set_channel_profile(const itpp::vec &avg_power_dB, const itpp::ivec &delay_prof)
{
    for (size_t i = 0; i < Nt * Nr; ++i)
        tdl_channels[i].set_channel_profile(avg_power_dB, delay_prof);
}

void
MIMO_TDL_Channel::set_channel_profile_uniform (int no_taps)
{
    for (size_t i = 0; i < Nt * Nr; ++i)
        tdl_channels[i].set_channel_profile_uniform(no_taps);
}

void
MIMO_TDL_Channel::set_channel_profile_exponential(int no_taps)
{
    for (size_t i = 0; i < Nt * Nr; ++i)
        tdl_channels[i].set_channel_profile_exponential(no_taps);
}

void
MIMO_TDL_Channel::generate(itpp::cmat *channel_coeff)
{
    itpp::cmat channel_coeff_one;
    int n_taps = tdl_channels[0].taps();
    for (int l = 0; l < n_taps; ++l) {
        channel_coeff[l].set_size(Nr, Nt);
    }

    for (size_t i = 0; i < Nr; ++i) {
        for (size_t j = 0; j < Nt; ++j) {
            tdl_channels[i * Nt + j].generate(1, channel_coeff_one);
            for (int l = 0; l < n_taps; ++l)
                channel_coeff[l](i, j) = channel_coeff_one(0, l);
        }
    }
}

void
MIMO_TDL_Channel::set_norm_doppler (double norm_doppler)
{
    std::for_each(tdl_channels, tdl_channels + Nt * Nr, [&norm_doppler](itpp::TDL_Channel &c) {
        c.set_norm_doppler(norm_doppler);
    });
}
void
MIMO_TDL_Channel::get_norm_doppler ()
{
    std::for_each(tdl_channels, tdl_channels + Nt * Nr, [](itpp::TDL_Channel &c) {
        std::cout<<c.get_norm_doppler()<<"\n";
    });
}
void
MIMO_TDL_Channel::set_channel_profile (const itpp::Channel_Specification &channel_spec, double sampling_time)
{
    std::for_each(tdl_channels, tdl_channels + Nt * Nr, [&channel_spec, &sampling_time](itpp::TDL_Channel &c) {
        c.set_channel_profile(channel_spec, sampling_time);
    });
}

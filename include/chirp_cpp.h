#ifndef CHIRP_H_
#define CHIRP_H_

#include <stdint.h>
class CHIRP
{
public:
    CHIRP();
    ~CHIRP();

    void init(float _f0, float _f1, uint32_t _N, float _Ts);
    void resetCount();
    float optimizeEndFrequency(const float &_f0, const float &_f1, const float &_t1, const float &_k0, const float &_k1);
    void resetSignals();
    bool update();

    float fchirp();
    float sinarg();
    float exc();

private:
    float m_f0, m_f1, m_t1, m_Ts, m_beta, m_k0, m_k1;
    uint32_t m_count, m_N;
    float m_exc, m_fchirp, m_sinarg;

    /*
    float   calc_beta(float _f0, float _f1, float _t1);
    */
};

#endif // CHIRP_H_
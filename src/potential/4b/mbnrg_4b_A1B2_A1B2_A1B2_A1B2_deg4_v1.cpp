#include "mbnrg_4b_A1B2_A1B2_A1B2_A1B2_deg4_v1.h"

////////////////////////////////////////////////////////////////////////////////

namespace mbnrg_A1B2_A1B2_A1B2_A1B2_deg4 {

mbnrg_A1B2_A1B2_A1B2_A1B2_deg4_v1::mbnrg_A1B2_A1B2_A1B2_A1B2_deg4_v1(const std::string mon1, const std::string mon2, const std::string mon3, const std::string mon4) {

    // =====>> BEGIN SECTION CONSTRUCTOR <<=====
    // =>> PASTE RIGHT BELOW THIS LINE <==

    if (mon1 == "h2o_2c3c4c" and mon2 == "h2o_2c3c4c" and mon3 == "h2o_2c3c4c" and mon4 == "h2o_2c3c4c") {
        coefficients = std::vector<double> {
             1.293592463246415e-02, // 0
            -2.527590670824149e-01, // 1
            -2.803575316329154e-03, // 2
             5.804689413975923e-02, // 3
             2.847408575171950e-01, // 4
            -4.796176341509499e-03, // 5
            -9.637042113523635e-01, // 6
             9.644754625330777e-03, // 7
             4.161496181286579e-02, // 8
             1.601366893925178e-01, // 9
             2.851748446433199e-01, // 10
            -5.268535826545701e-03, // 11
             2.450746069658527e-02, // 12
             9.341768497153155e-02, // 13
             1.084375051321213e-01, // 14
            -4.774604283678897e-02, // 15
            -3.529339195058244e-03, // 16
             1.157155453746819e+01, // 17
             2.916701950365417e-01, // 18
             1.117855386799011e-06, // 19
             7.604375827973731e-03, // 20
            -2.784844256972470e-03, // 21
            -1.890304228747910e+01, // 22
            -1.932832583688824e-01, // 23
            -6.229366460207162e-03, // 24
            -2.122542022932476e-03, // 25
             1.037285719394933e-01, // 26
             4.971006036902755e+00, // 27
            -8.368657484209907e-02, // 28
             1.472859459908757e-01, // 29
            -1.484002864903285e-03, // 30
             5.229591608994186e-01, // 31
             5.106466354200501e-02, // 32
            -2.452781917014893e-03, // 33
            -3.584637218090634e-01, // 34
            -7.677343910193245e-04, // 35
            -1.116827173633696e-02, // 36
            -2.578952045636038e-02, // 37
             5.576972573021179e-02, // 38
             6.109879460804083e-02, // 39
             4.631178597743674e-01, // 40
            -2.313860744430851e-02, // 41
            -1.451661598671804e-05, // 42
             3.916265182297148e-01, // 43
            -6.691834022270447e-02, // 44
             1.708843990744590e+00, // 45
            -1.103087060903730e-01, // 46
             5.364470456655808e-01, // 47
            -1.350827177803829e+00, // 48
            -5.973517100668391e+00, // 49
            -1.598993958592549e+01, // 50
             7.484931881019395e-02, // 51
             2.000326388522788e-02, // 52
            -2.639944295894682e+00, // 53
             1.521865652275861e+01, // 54
             3.548932918839463e-02, // 55
            -5.365655893842659e-02, // 56
             4.916365520620140e-02, // 57
             1.366500718780711e-02, // 58
             1.578615586114169e-03, // 59
            -5.897259000430088e-02, // 60
            -8.752676820146151e-02, // 61
            -5.327985550908109e-02, // 62
             2.206154981151537e-02, // 63
            -2.888503327146351e-04, // 64
            -2.898726652647888e-02, // 65
            -7.233573400005681e-02, // 66
             4.872317917538553e-02, // 67
            -1.968012851369423e-03, // 68
            -1.553018511180274e-01, // 69
            -5.712366672642913e-01, // 70
             1.422476991713441e-03, // 71
             1.606066596063603e-01, // 72
            -2.363771780616139e-01, // 73
            -2.617237945428449e-01, // 74
             1.282797804230772e-03, // 75
             1.116684549281102e+01, // 76
             2.121824543395118e-01, // 77
             7.572688650179683e-02, // 78
             9.638049100377916e-02, // 79
            -3.033294998426339e-03, // 80
            -4.829243977746225e-03, // 81
            -7.358629955283671e-02, // 82
             1.170736583986732e+01, // 83
            -6.554783773385704e+00, // 84
            -2.504951398288827e+00, // 85
             2.223089287969986e-03, // 86
             9.747169777602457e-01, // 87
             8.812738475450660e-01, // 88
            -1.187659299601829e-01, // 89
             2.453451563944304e+00, // 90
             1.612997372228605e+00, // 91
             4.975108518716281e+00, // 92
             8.061015219199991e+00, // 93
            -1.243976189142013e-01, // 94
             4.717529505812305e+00, // 95
            -5.689459794021250e-02, // 96
            -5.252383918852110e-02, // 97
            -3.333788252981905e-02, // 98
            -2.746601444331877e+00, // 99
             1.712333690926624e+00, // 100
             6.747731047370071e+00, // 101
            -8.531781997783510e+00, // 102
            -1.066964735890106e-01, // 103
            -5.259521629190689e+00, // 104
             5.430902091408301e-02, // 105
             3.055776467090027e+00, // 106
            -2.499297779337763e+00, // 107
            -2.262287479479043e+00, // 108
            -4.849919009931579e+00, // 109
             1.117880520555618e-01, // 110
             1.572187051667365e+01, // 111
            -8.087416060126944e-02, // 112
             5.573787807973442e+00, // 113
            -5.706089919530725e-02, // 114
             6.783210139216116e-02, // 115
             4.673759229288594e+00, // 116
            -3.527624862675605e-02, // 117
             4.007150142467246e+00, // 118
             1.849748062125000e+00, // 119
             2.512166883015098e-03, // 120
             5.411693582867585e-02, // 121
            -4.520035742962915e-01, // 122
             6.001212984825992e+00, // 123
            -1.098582142640980e+01, // 124
             1.608721498115200e+00, // 125
            -7.527075808020294e+00, // 126
             8.456205209473650e-02, // 127
            -1.808769165251939e+00, // 128
             2.807718146539471e-01, // 129
             2.542277813914775e+00, // 130
            -2.989121324095351e-01, // 131
             1.949879343546627e+00, // 132
            -8.577815582588929e+00, // 133
            -1.733149347230874e+00, // 134
             4.722360264532841e-02, // 135
            -2.454996165134943e+00, // 136
             2.508917880282612e+00, // 137
            -1.982608299510560e-01, // 138
             1.874757800477777e-02, // 139
            -5.430178549571291e+00, // 140
            -4.933174935982692e+00, // 141
             2.115141844541924e-03, // 142
            -2.593439814227163e+00, // 143
             3.559750569361740e-01, // 144
            -1.891920862693904e+00, // 145
             7.183188720465240e-02, // 146
            -1.550538913257653e-01, // 147
             1.339192430723058e+00, // 148
            -9.409161247512263e-01, // 149
            -5.764707091502964e+00, // 150
            -4.099002442472462e-02, // 151
             2.272993239724152e-01, // 152
             6.877130638030310e-02, // 153
             1.148720932678486e+00, // 154
             3.780474510004141e+00, // 155
             1.765193319484586e-02, // 156
             4.261427449092903e+00, // 157
             1.318014011670042e+00, // 158
             4.921431710958223e-02, // 159
             2.192895540823795e-01, // 160
             1.935191835127353e+00, // 161
            -2.096662172085008e+00, // 162
             2.060313727080563e+00, // 163
             3.669741538260437e+00, // 164
             8.250850176571214e-02, // 165
             1.948738986155918e+00, // 166
             1.039699004388605e+01, // 167
            -1.201622870405440e-01, // 168
             3.156881895564944e+00, // 169
            -8.928868076517160e-02, // 170
             9.878316045470340e-01, // 171
             7.813297932993519e-01, // 172
            -3.373735016965507e+00, // 173
            -4.161892178431281e+00, // 174
             4.046030913747598e+00, // 175
             1.523916239004319e+00, // 176
            -6.367501843460051e+00, // 177
             2.125755807254886e+00, // 178
            -2.159323529942674e+00, // 179
            -5.730826043067860e+00, // 180
            -4.299605991815664e-02, // 181
            -3.367155256539849e-01, // 182
             3.671842834019631e-01, // 183
             2.464226590004468e+00, // 184
             3.063929835496449e-01, // 185
             5.268271315016895e+00, // 186
             7.644144162877521e+00, // 187
            -1.123857257983671e+00, // 188
             1.140296167084711e+01, // 189
             8.809807564975848e+00, // 190
             2.171240824358824e-02, // 191
             2.440701825319023e-01, // 192
             6.732581417161682e-01, // 193
            -2.710770829133488e+00, // 194
             3.148950883301212e+00, // 195
             1.691108983209215e-02, // 196
            -1.519497347942965e+00, // 197
             1.166937518946979e-01, // 198
             3.854364531290923e-02, // 199
            -9.071404443243404e-01, // 200
             3.627902835256771e-03, // 201
             8.657301303572879e-02, // 202
             3.399472676718787e-01, // 203
            -2.552916482230969e-01, // 204
             7.253254049385321e+00, // 205
            -1.116732795218704e-01, // 206
            -9.595718861835111e+00, // 207
             3.382767069236742e-02, // 208
             1.687830384201882e+00, // 209
             3.659900119183443e+00, // 210
             1.431000809186601e-01, // 211
            -8.308906351417338e-02, // 212
            -5.856339366349485e+00, // 213
             2.629187868506676e+00, // 214
            -6.427401565962681e+00, // 215
            -1.317084239227138e+00, // 216
             1.220911348859394e+00, // 217
            -1.773398954455080e+00, // 218
             3.502180882166567e+00, // 219
            -2.155162412646093e-01, // 220
             9.962844411111504e-02, // 221
            -1.458341566514101e-01, // 222
            -2.872083098323456e+00, // 223
            -1.346997890982987e+01, // 224
            -8.094696832898932e-03, // 225
            -2.736935078989263e+00, // 226
            -9.745365553611158e-04, // 227
            -2.017558318180047e-01, // 228
            -6.889435668901804e+00, // 229
             2.774766443015775e-01, // 230
             1.669402454943881e-01, // 231
             3.688088196124693e-02, // 232
             2.647852048115966e+00, // 233
             4.266419141461747e-02, // 234
            -9.212399086978847e+00, // 235
             2.309860054909894e+00, // 236
             1.319958688253828e-01, // 237
             7.164936886488375e-02, // 238
            -3.442995530118020e+00, // 239
            -3.230839090178099e-01, // 240
            -4.466905130073747e+00, // 241
            -1.087689237812563e+00, // 242
            -1.235572993332932e+00, // 243
             5.488745858846029e+00, // 244
             2.121299741082140e-01, // 245
            -3.038895083124385e-01, // 246
             3.132070821617974e+00, // 247
            -3.896406381011970e+00, // 248
            -1.148013894724111e+01, // 249
             2.036708858047252e-02, // 250
            -1.111395534194044e-01, // 251
            -6.495893639837876e+00, // 252
            -2.517693813639191e+00, // 253
             2.629174506856737e-02, // 254
            -4.678400219176606e+00, // 255
             1.100301978767469e+00, // 256
            -1.970620683834761e-01, // 257
            -6.743432593561449e-02, // 258
            -1.001123194878663e+00, // 259
             8.460194474658680e+00, // 260
             1.214064610370679e-01, // 261
            -5.045477267349117e+00, // 262
             4.847163145339555e-02, // 263
             3.998899683956227e+00, // 264
            -2.745498647413954e-01}; // 265
        m_k_x_inter_A_A_0 =  1.244277387934733e+00; // A^(-1))
        m_d_x_inter_A_A_0 =  5.701501264744776e+00; // A^(-1))
        m_k_x_intra_A_B_1 =  1.487372640271976e-03; // A^(-1))
        m_d_x_intra_A_B_1 =  2.991930120243535e+00; // A^(-1))
        m_k_x_inter_A_B_0 =  2.205457840065838e+00; // A^(-1))
        m_d_x_inter_A_B_0 =  1.862092832490158e+00; // A^(-1))
        m_k_x_intra_B_B_1 =  7.572254854896845e-04; // A^(-1))
        m_d_x_intra_B_B_1 =  5.111788710064341e-03; // A^(-1))
        m_k_x_inter_B_B_0 =  6.073417761499551e-01; // A^(-1))
        m_d_x_inter_B_B_0 =  1.348545171945606e+00; // A^(-1))
        m_ri =  5.500000000000000e+00; // A
        m_ro =  6.500000000000000e+00; // A

    } // end if mon1 == "h2o" and mon2 == "h2o" and mon3 == "h2o" and mon4 == "h2o"
    // =====>> END SECTION CONSTRUCTOR <<=====
}

//----------------------------------------------------------------------------//

double mbnrg_A1B2_A1B2_A1B2_A1B2_deg4_v1::f_switch(const double r, double& g)
{
    if (r > m_ro) {
        g = 0.0;
        return 0.0;
    } else if (r > m_ri) {
        const double t1 = M_PI/(m_ro - m_ri);
        const double x = (r - m_ri)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

//----------------------------------------------------------------------------//

 double mbnrg_A1B2_A1B2_A1B2_A1B2_deg4_v1::eval(const double *xyz1, const double *xyz2, const double *xyz3, const double *xyz4, const size_t n) {
    std::vector<double> energies(n,0.0);
    std::vector<double> energies_sw(n,0.0);

    std::vector<double> xyz(36);
    double sw = 0.0;
    polynomial my_poly;

    for (size_t j = 0; j < n; j++) {
        const double *mon1 = xyz1 + 9*j;
        const double *mon2 = xyz2 + 9*j;
        const double *mon3 = xyz3 + 9*j;
        const double *mon4 = xyz4 + 9*j;


        const double d12[3] =
                         {mon1[0] - mon2[0],
                          mon1[1] - mon2[1],
                          mon1[2] - mon2[2]};

        const double d12rsq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
        const double d12r = std::sqrt(d12rsq);

        const double d13[3] =
                         {mon1[0] - mon3[0],
                          mon1[1] - mon3[1],
                          mon1[2] - mon3[2]};

        const double d13rsq = d13[0]*d13[0] + d13[1]*d13[1] + d13[2]*d13[2];
        const double d13r = std::sqrt(d13rsq);

        const double d14[3] =
                         {mon1[0] - mon4[0],
                          mon1[1] - mon4[1],
                          mon1[2] - mon4[2]};

        const double d14rsq = d14[0]*d14[0] + d14[1]*d14[1] + d14[2]*d14[2];
        const double d14r = std::sqrt(d14rsq);

        const double d23[3] =
                         {mon2[0] - mon3[0],
                          mon2[1] - mon3[1],
                          mon2[2] - mon3[2]};

        const double d23rsq = d23[0]*d23[0] + d23[1]*d23[1] + d23[2]*d23[2];
        const double d23r = std::sqrt(d23rsq);

        const double d24[3] =
                         {mon2[0] - mon4[0],
                          mon2[1] - mon4[1],
                          mon2[2] - mon4[2]};

        const double d24rsq = d24[0]*d24[0] + d24[1]*d24[1] + d24[2]*d24[2];
        const double d24r = std::sqrt(d24rsq);

        const double d34[3] =
                         {mon3[0] - mon4[0],
                          mon3[1] - mon4[1],
                          mon3[2] - mon4[2]};

        const double d34rsq = d34[0]*d34[0] + d34[1]*d34[1] + d34[2]*d34[2];
        const double d34r = std::sqrt(d34rsq);

        if (true  && d12r > m_ro  && d13r > m_ro  && d14r > m_ro  && d23r > m_ro  && d24r > m_ro  && d34r > m_ro ) {
             continue;
        }

        std::copy(mon1, mon1 + 9, xyz.begin() + 0);

        std::copy(mon2, mon2 + 9, xyz.begin() + 9);

        std::copy(mon3, mon3 + 9, xyz.begin() + 18);

        std::copy(mon4, mon4 + 9, xyz.begin() + 27);


        const double* coords_A_1_a = xyz.data() + 0;

        const double* coords_B_1_a = xyz.data() + 3;

        const double* coords_B_2_a = xyz.data() + 6;

        const double* coords_A_2_b = xyz.data() + 9;

        const double* coords_B_3_b = xyz.data() + 12;

        const double* coords_B_4_b = xyz.data() + 15;

        const double* coords_A_3_c = xyz.data() + 18;

        const double* coords_B_5_c = xyz.data() + 21;

        const double* coords_B_6_c = xyz.data() + 24;

        const double* coords_A_4_d = xyz.data() + 27;

        const double* coords_B_7_d = xyz.data() + 30;

        const double* coords_B_8_d = xyz.data() + 33;


        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;

    
        variable vs[66];

        double xs[66];

        xs[0] = vs[0].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_1_a, coords_A_2_b);
        xs[1] = vs[1].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_1_a, coords_A_3_c);
        xs[2] = vs[2].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_1_a, coords_A_4_d);
        xs[3] = vs[3].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_1_a, coords_B_1_a);
        xs[4] = vs[4].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_1_a, coords_B_2_a);
        xs[5] = vs[5].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_3_b);
        xs[6] = vs[6].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_4_b);
        xs[7] = vs[7].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_5_c);
        xs[8] = vs[8].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_6_c);
        xs[9] = vs[9].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_7_d);
        xs[10] = vs[10].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_8_d);
        xs[11] = vs[11].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_2_b, coords_A_3_c);
        xs[12] = vs[12].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_2_b, coords_A_4_d);
        xs[13] = vs[13].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_1_a);
        xs[14] = vs[14].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_2_a);
        xs[15] = vs[15].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_2_b, coords_B_3_b);
        xs[16] = vs[16].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_2_b, coords_B_4_b);
        xs[17] = vs[17].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_5_c);
        xs[18] = vs[18].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_6_c);
        xs[19] = vs[19].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_7_d);
        xs[20] = vs[20].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_8_d);
        xs[21] = vs[21].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_3_c, coords_A_4_d);
        xs[22] = vs[22].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_1_a);
        xs[23] = vs[23].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_2_a);
        xs[24] = vs[24].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_3_b);
        xs[25] = vs[25].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_4_b);
        xs[26] = vs[26].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_3_c, coords_B_5_c);
        xs[27] = vs[27].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_3_c, coords_B_6_c);
        xs[28] = vs[28].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_7_d);
        xs[29] = vs[29].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_8_d);
        xs[30] = vs[30].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_1_a);
        xs[31] = vs[31].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_2_a);
        xs[32] = vs[32].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_3_b);
        xs[33] = vs[33].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_4_b);
        xs[34] = vs[34].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_5_c);
        xs[35] = vs[35].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_6_c);
        xs[36] = vs[36].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_4_d, coords_B_7_d);
        xs[37] = vs[37].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_4_d, coords_B_8_d);
        xs[38] = vs[38].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_1_a, coords_B_2_a);
        xs[39] = vs[39].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_3_b);
        xs[40] = vs[40].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_4_b);
        xs[41] = vs[41].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_5_c);
        xs[42] = vs[42].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_6_c);
        xs[43] = vs[43].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_7_d);
        xs[44] = vs[44].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_8_d);
        xs[45] = vs[45].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_3_b);
        xs[46] = vs[46].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_4_b);
        xs[47] = vs[47].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_5_c);
        xs[48] = vs[48].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_6_c);
        xs[49] = vs[49].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_7_d);
        xs[50] = vs[50].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_8_d);
        xs[51] = vs[51].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_3_b, coords_B_4_b);
        xs[52] = vs[52].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_5_c);
        xs[53] = vs[53].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_6_c);
        xs[54] = vs[54].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_7_d);
        xs[55] = vs[55].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_8_d);
        xs[56] = vs[56].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_5_c);
        xs[57] = vs[57].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_6_c);
        xs[58] = vs[58].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_7_d);
        xs[59] = vs[59].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_8_d);
        xs[60] = vs[60].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_5_c, coords_B_6_c);
        xs[61] = vs[61].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_5_c, coords_B_7_d);
        xs[62] = vs[62].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_5_c, coords_B_8_d);
        xs[63] = vs[63].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_6_c, coords_B_7_d);
        xs[64] = vs[64].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_6_c, coords_B_8_d);
        xs[65] = vs[65].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_7_d, coords_B_8_d);

        double gsw12 = 0.0;
        double sw12 = f_switch(d12r, gsw12);
        double gsw13 = 0.0;
        double sw13 = f_switch(d13r, gsw13);
        double gsw14 = 0.0;
        double sw14 = f_switch(d14r, gsw14);
        double gsw23 = 0.0;
        double sw23 = f_switch(d23r, gsw23);
        double gsw24 = 0.0;
        double sw24 = f_switch(d24r, gsw24);
        double gsw34 = 0.0;
        double sw34 = f_switch(d34r, gsw34);

        sw = sw12*sw13*sw14*sw23*sw24*sw34; // cutoff based on largest d-OO
        // sw = 3*sw12*sw13*sw14*sw23*sw24*sw34 - sw13*sw14*sw23*sw24*sw34 - sw12*sw14*sw23*sw24*sw34 - sw12*sw13*sw23*sw24*sw34 - sw12*sw13*sw14*sw24*sw34 - sw12*sw13*sw14*sw23*sw34 - sw12*sw13*sw14*sw23*sw24 + sw12*sw13*sw14 + sw12*sw23*sw24 + sw13*sw23*sw34 + sw14*sw24*sw34; // cutoff based on center-3-neighbor criteria

        energies[j] = my_poly.eval(xs,coefficients.data());
        energies_sw[j] = energies[j]*sw;

    }

    double energy = 0.0;
    for (size_t i = 0; i < n; i++)
        energy += energies_sw[i];

    return energy;

}

//----------------------------------------------------------------------------//

double mbnrg_A1B2_A1B2_A1B2_A1B2_deg4_v1::eval(const double *xyz1, const double *xyz2, const double *xyz3, const double *xyz4, double *grad1, double *grad2, double *grad3, double *grad4 , const size_t n, std::vector<double> *virial) {
    std::vector<double> energies(n,0.0);
    std::vector<double> energies_sw(n,0.0);

    std::vector<double> xyz(36);
    double sw = 0.0;
    polynomial my_poly;

    for (size_t j = 0; j < n; j++) {
        const double *mon1 = xyz1 + 9*j;
        const double *mon2 = xyz2 + 9*j;
        const double *mon3 = xyz3 + 9*j;
        const double *mon4 = xyz4 + 9*j;


        const double d12[3] =
                         {mon1[0] - mon2[0],
                          mon1[1] - mon2[1],
                          mon1[2] - mon2[2]};

        const double d12rsq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
        const double d12r = std::sqrt(d12rsq);

        const double d13[3] =
                         {mon1[0] - mon3[0],
                          mon1[1] - mon3[1],
                          mon1[2] - mon3[2]};

        const double d13rsq = d13[0]*d13[0] + d13[1]*d13[1] + d13[2]*d13[2];
        const double d13r = std::sqrt(d13rsq);

        const double d14[3] =
                         {mon1[0] - mon4[0],
                          mon1[1] - mon4[1],
                          mon1[2] - mon4[2]};

        const double d14rsq = d14[0]*d14[0] + d14[1]*d14[1] + d14[2]*d14[2];
        const double d14r = std::sqrt(d14rsq);

        const double d23[3] =
                         {mon2[0] - mon3[0],
                          mon2[1] - mon3[1],
                          mon2[2] - mon3[2]};

        const double d23rsq = d23[0]*d23[0] + d23[1]*d23[1] + d23[2]*d23[2];
        const double d23r = std::sqrt(d23rsq);

        const double d24[3] =
                         {mon2[0] - mon4[0],
                          mon2[1] - mon4[1],
                          mon2[2] - mon4[2]};

        const double d24rsq = d24[0]*d24[0] + d24[1]*d24[1] + d24[2]*d24[2];
        const double d24r = std::sqrt(d24rsq);

        const double d34[3] =
                         {mon3[0] - mon4[0],
                          mon3[1] - mon4[1],
                          mon3[2] - mon4[2]};

        const double d34rsq = d34[0]*d34[0] + d34[1]*d34[1] + d34[2]*d34[2];
        const double d34r = std::sqrt(d34rsq);

        if (true  && d12r > m_ro  && d13r > m_ro  && d14r > m_ro  && d23r > m_ro  && d24r > m_ro  && d34r > m_ro ) {
             continue;
        }

        std::vector<double> gradients(36,0.0);

        std::copy(mon1, mon1 + 9, xyz.begin() + 0);

        std::copy(mon2, mon2 + 9, xyz.begin() + 9);

        std::copy(mon3, mon3 + 9, xyz.begin() + 18);

        std::copy(mon4, mon4 + 9, xyz.begin() + 27);
        const double* coords_A_1_a = xyz.data() + 0;

        const double* coords_B_1_a = xyz.data() + 3;

        const double* coords_B_2_a = xyz.data() + 6;

        const double* coords_A_2_b = xyz.data() + 9;

        const double* coords_B_3_b = xyz.data() + 12;

        const double* coords_B_4_b = xyz.data() + 15;

        const double* coords_A_3_c = xyz.data() + 18;

        const double* coords_B_5_c = xyz.data() + 21;

        const double* coords_B_6_c = xyz.data() + 24;

        const double* coords_A_4_d = xyz.data() + 27;

        const double* coords_B_7_d = xyz.data() + 30;

        const double* coords_B_8_d = xyz.data() + 33;


        double* coords_A_1_a_g = gradients.data() + 0;

        double* coords_B_1_a_g = gradients.data() + 3;

        double* coords_B_2_a_g = gradients.data() + 6;

        double* coords_A_2_b_g = gradients.data() + 9;

        double* coords_B_3_b_g = gradients.data() + 12;

        double* coords_B_4_b_g = gradients.data() + 15;

        double* coords_A_3_c_g = gradients.data() + 18;

        double* coords_B_5_c_g = gradients.data() + 21;

        double* coords_B_6_c_g = gradients.data() + 24;

        double* coords_A_4_d_g = gradients.data() + 27;

        double* coords_B_7_d_g = gradients.data() + 30;

        double* coords_B_8_d_g = gradients.data() + 33;



        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;

    
        variable vs[66];

        double xs[66];


        double gxs[66];

        xs[0] = vs[0].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_1_a, coords_A_2_b);
        xs[1] = vs[1].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_1_a, coords_A_3_c);
        xs[2] = vs[2].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_1_a, coords_A_4_d);
        xs[3] = vs[3].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_1_a, coords_B_1_a);
        xs[4] = vs[4].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_1_a, coords_B_2_a);
        xs[5] = vs[5].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_3_b);
        xs[6] = vs[6].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_4_b);
        xs[7] = vs[7].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_5_c);
        xs[8] = vs[8].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_6_c);
        xs[9] = vs[9].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_7_d);
        xs[10] = vs[10].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_1_a, coords_B_8_d);
        xs[11] = vs[11].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_2_b, coords_A_3_c);
        xs[12] = vs[12].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_2_b, coords_A_4_d);
        xs[13] = vs[13].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_1_a);
        xs[14] = vs[14].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_2_a);
        xs[15] = vs[15].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_2_b, coords_B_3_b);
        xs[16] = vs[16].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_2_b, coords_B_4_b);
        xs[17] = vs[17].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_5_c);
        xs[18] = vs[18].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_6_c);
        xs[19] = vs[19].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_7_d);
        xs[20] = vs[20].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_2_b, coords_B_8_d);
        xs[21] = vs[21].v_exp0(m_d_x_inter_A_A_0, m_k_x_inter_A_A_0, coords_A_3_c, coords_A_4_d);
        xs[22] = vs[22].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_1_a);
        xs[23] = vs[23].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_2_a);
        xs[24] = vs[24].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_3_b);
        xs[25] = vs[25].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_4_b);
        xs[26] = vs[26].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_3_c, coords_B_5_c);
        xs[27] = vs[27].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_3_c, coords_B_6_c);
        xs[28] = vs[28].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_7_d);
        xs[29] = vs[29].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_3_c, coords_B_8_d);
        xs[30] = vs[30].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_1_a);
        xs[31] = vs[31].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_2_a);
        xs[32] = vs[32].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_3_b);
        xs[33] = vs[33].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_4_b);
        xs[34] = vs[34].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_5_c);
        xs[35] = vs[35].v_exp0(m_d_x_inter_A_B_0, m_k_x_inter_A_B_0, coords_A_4_d, coords_B_6_c);
        xs[36] = vs[36].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_4_d, coords_B_7_d);
        xs[37] = vs[37].v_exp0(m_d_x_intra_A_B_1, m_k_x_intra_A_B_1, coords_A_4_d, coords_B_8_d);
        xs[38] = vs[38].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_1_a, coords_B_2_a);
        xs[39] = vs[39].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_3_b);
        xs[40] = vs[40].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_4_b);
        xs[41] = vs[41].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_5_c);
        xs[42] = vs[42].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_6_c);
        xs[43] = vs[43].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_7_d);
        xs[44] = vs[44].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_1_a, coords_B_8_d);
        xs[45] = vs[45].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_3_b);
        xs[46] = vs[46].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_4_b);
        xs[47] = vs[47].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_5_c);
        xs[48] = vs[48].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_6_c);
        xs[49] = vs[49].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_7_d);
        xs[50] = vs[50].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_2_a, coords_B_8_d);
        xs[51] = vs[51].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_3_b, coords_B_4_b);
        xs[52] = vs[52].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_5_c);
        xs[53] = vs[53].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_6_c);
        xs[54] = vs[54].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_7_d);
        xs[55] = vs[55].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_3_b, coords_B_8_d);
        xs[56] = vs[56].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_5_c);
        xs[57] = vs[57].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_6_c);
        xs[58] = vs[58].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_7_d);
        xs[59] = vs[59].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_4_b, coords_B_8_d);
        xs[60] = vs[60].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_5_c, coords_B_6_c);
        xs[61] = vs[61].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_5_c, coords_B_7_d);
        xs[62] = vs[62].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_5_c, coords_B_8_d);
        xs[63] = vs[63].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_6_c, coords_B_7_d);
        xs[64] = vs[64].v_exp0(m_d_x_inter_B_B_0, m_k_x_inter_B_B_0, coords_B_6_c, coords_B_8_d);
        xs[65] = vs[65].v_exp0(m_d_x_intra_B_B_1, m_k_x_intra_B_B_1, coords_B_7_d, coords_B_8_d);

        double gsw12 = 0.0;
        double sw12 = f_switch(d12r, gsw12);
        double gsw13 = 0.0;
        double sw13 = f_switch(d13r, gsw13);
        double gsw14 = 0.0;
        double sw14 = f_switch(d14r, gsw14);
        double gsw23 = 0.0;
        double sw23 = f_switch(d23r, gsw23);
        double gsw24 = 0.0;
        double sw24 = f_switch(d24r, gsw24);
        double gsw34 = 0.0;
        double sw34 = f_switch(d34r, gsw34);

        sw = sw12*sw13*sw14*sw23*sw24*sw34; // cutoff based on largest d-OO
        // sw = 3*sw12*sw13*sw14*sw23*sw24*sw34 - sw13*sw14*sw23*sw24*sw34 - sw12*sw14*sw23*sw24*sw34 - sw12*sw13*sw23*sw24*sw34 - sw12*sw13*sw14*sw24*sw34 - sw12*sw13*sw14*sw23*sw34 - sw12*sw13*sw14*sw23*sw24 + sw12*sw13*sw14 + sw12*sw23*sw24 + sw13*sw23*sw34 + sw14*sw24*sw34; // cutoff based on center-3-neighbor criteria

        energies[j] = my_poly.eval(xs,coefficients.data(),gxs);
        energies_sw[j] = energies[j]*sw;

        for (size_t i = 0; i < 66; i++) {
            gxs[i] *= sw;
        }

        vs[0].grads(gxs[0], coords_A_1_a_g, coords_A_2_b_g, coords_A_1_a, coords_A_2_b);
        vs[1].grads(gxs[1], coords_A_1_a_g, coords_A_3_c_g, coords_A_1_a, coords_A_3_c);
        vs[2].grads(gxs[2], coords_A_1_a_g, coords_A_4_d_g, coords_A_1_a, coords_A_4_d);
        vs[3].grads(gxs[3], coords_A_1_a_g, coords_B_1_a_g, coords_A_1_a, coords_B_1_a);
        vs[4].grads(gxs[4], coords_A_1_a_g, coords_B_2_a_g, coords_A_1_a, coords_B_2_a);
        vs[5].grads(gxs[5], coords_A_1_a_g, coords_B_3_b_g, coords_A_1_a, coords_B_3_b);
        vs[6].grads(gxs[6], coords_A_1_a_g, coords_B_4_b_g, coords_A_1_a, coords_B_4_b);
        vs[7].grads(gxs[7], coords_A_1_a_g, coords_B_5_c_g, coords_A_1_a, coords_B_5_c);
        vs[8].grads(gxs[8], coords_A_1_a_g, coords_B_6_c_g, coords_A_1_a, coords_B_6_c);
        vs[9].grads(gxs[9], coords_A_1_a_g, coords_B_7_d_g, coords_A_1_a, coords_B_7_d);
        vs[10].grads(gxs[10], coords_A_1_a_g, coords_B_8_d_g, coords_A_1_a, coords_B_8_d);
        vs[11].grads(gxs[11], coords_A_2_b_g, coords_A_3_c_g, coords_A_2_b, coords_A_3_c);
        vs[12].grads(gxs[12], coords_A_2_b_g, coords_A_4_d_g, coords_A_2_b, coords_A_4_d);
        vs[13].grads(gxs[13], coords_A_2_b_g, coords_B_1_a_g, coords_A_2_b, coords_B_1_a);
        vs[14].grads(gxs[14], coords_A_2_b_g, coords_B_2_a_g, coords_A_2_b, coords_B_2_a);
        vs[15].grads(gxs[15], coords_A_2_b_g, coords_B_3_b_g, coords_A_2_b, coords_B_3_b);
        vs[16].grads(gxs[16], coords_A_2_b_g, coords_B_4_b_g, coords_A_2_b, coords_B_4_b);
        vs[17].grads(gxs[17], coords_A_2_b_g, coords_B_5_c_g, coords_A_2_b, coords_B_5_c);
        vs[18].grads(gxs[18], coords_A_2_b_g, coords_B_6_c_g, coords_A_2_b, coords_B_6_c);
        vs[19].grads(gxs[19], coords_A_2_b_g, coords_B_7_d_g, coords_A_2_b, coords_B_7_d);
        vs[20].grads(gxs[20], coords_A_2_b_g, coords_B_8_d_g, coords_A_2_b, coords_B_8_d);
        vs[21].grads(gxs[21], coords_A_3_c_g, coords_A_4_d_g, coords_A_3_c, coords_A_4_d);
        vs[22].grads(gxs[22], coords_A_3_c_g, coords_B_1_a_g, coords_A_3_c, coords_B_1_a);
        vs[23].grads(gxs[23], coords_A_3_c_g, coords_B_2_a_g, coords_A_3_c, coords_B_2_a);
        vs[24].grads(gxs[24], coords_A_3_c_g, coords_B_3_b_g, coords_A_3_c, coords_B_3_b);
        vs[25].grads(gxs[25], coords_A_3_c_g, coords_B_4_b_g, coords_A_3_c, coords_B_4_b);
        vs[26].grads(gxs[26], coords_A_3_c_g, coords_B_5_c_g, coords_A_3_c, coords_B_5_c);
        vs[27].grads(gxs[27], coords_A_3_c_g, coords_B_6_c_g, coords_A_3_c, coords_B_6_c);
        vs[28].grads(gxs[28], coords_A_3_c_g, coords_B_7_d_g, coords_A_3_c, coords_B_7_d);
        vs[29].grads(gxs[29], coords_A_3_c_g, coords_B_8_d_g, coords_A_3_c, coords_B_8_d);
        vs[30].grads(gxs[30], coords_A_4_d_g, coords_B_1_a_g, coords_A_4_d, coords_B_1_a);
        vs[31].grads(gxs[31], coords_A_4_d_g, coords_B_2_a_g, coords_A_4_d, coords_B_2_a);
        vs[32].grads(gxs[32], coords_A_4_d_g, coords_B_3_b_g, coords_A_4_d, coords_B_3_b);
        vs[33].grads(gxs[33], coords_A_4_d_g, coords_B_4_b_g, coords_A_4_d, coords_B_4_b);
        vs[34].grads(gxs[34], coords_A_4_d_g, coords_B_5_c_g, coords_A_4_d, coords_B_5_c);
        vs[35].grads(gxs[35], coords_A_4_d_g, coords_B_6_c_g, coords_A_4_d, coords_B_6_c);
        vs[36].grads(gxs[36], coords_A_4_d_g, coords_B_7_d_g, coords_A_4_d, coords_B_7_d);
        vs[37].grads(gxs[37], coords_A_4_d_g, coords_B_8_d_g, coords_A_4_d, coords_B_8_d);
        vs[38].grads(gxs[38], coords_B_1_a_g, coords_B_2_a_g, coords_B_1_a, coords_B_2_a);
        vs[39].grads(gxs[39], coords_B_1_a_g, coords_B_3_b_g, coords_B_1_a, coords_B_3_b);
        vs[40].grads(gxs[40], coords_B_1_a_g, coords_B_4_b_g, coords_B_1_a, coords_B_4_b);
        vs[41].grads(gxs[41], coords_B_1_a_g, coords_B_5_c_g, coords_B_1_a, coords_B_5_c);
        vs[42].grads(gxs[42], coords_B_1_a_g, coords_B_6_c_g, coords_B_1_a, coords_B_6_c);
        vs[43].grads(gxs[43], coords_B_1_a_g, coords_B_7_d_g, coords_B_1_a, coords_B_7_d);
        vs[44].grads(gxs[44], coords_B_1_a_g, coords_B_8_d_g, coords_B_1_a, coords_B_8_d);
        vs[45].grads(gxs[45], coords_B_2_a_g, coords_B_3_b_g, coords_B_2_a, coords_B_3_b);
        vs[46].grads(gxs[46], coords_B_2_a_g, coords_B_4_b_g, coords_B_2_a, coords_B_4_b);
        vs[47].grads(gxs[47], coords_B_2_a_g, coords_B_5_c_g, coords_B_2_a, coords_B_5_c);
        vs[48].grads(gxs[48], coords_B_2_a_g, coords_B_6_c_g, coords_B_2_a, coords_B_6_c);
        vs[49].grads(gxs[49], coords_B_2_a_g, coords_B_7_d_g, coords_B_2_a, coords_B_7_d);
        vs[50].grads(gxs[50], coords_B_2_a_g, coords_B_8_d_g, coords_B_2_a, coords_B_8_d);
        vs[51].grads(gxs[51], coords_B_3_b_g, coords_B_4_b_g, coords_B_3_b, coords_B_4_b);
        vs[52].grads(gxs[52], coords_B_3_b_g, coords_B_5_c_g, coords_B_3_b, coords_B_5_c);
        vs[53].grads(gxs[53], coords_B_3_b_g, coords_B_6_c_g, coords_B_3_b, coords_B_6_c);
        vs[54].grads(gxs[54], coords_B_3_b_g, coords_B_7_d_g, coords_B_3_b, coords_B_7_d);
        vs[55].grads(gxs[55], coords_B_3_b_g, coords_B_8_d_g, coords_B_3_b, coords_B_8_d);
        vs[56].grads(gxs[56], coords_B_4_b_g, coords_B_5_c_g, coords_B_4_b, coords_B_5_c);
        vs[57].grads(gxs[57], coords_B_4_b_g, coords_B_6_c_g, coords_B_4_b, coords_B_6_c);
        vs[58].grads(gxs[58], coords_B_4_b_g, coords_B_7_d_g, coords_B_4_b, coords_B_7_d);
        vs[59].grads(gxs[59], coords_B_4_b_g, coords_B_8_d_g, coords_B_4_b, coords_B_8_d);
        vs[60].grads(gxs[60], coords_B_5_c_g, coords_B_6_c_g, coords_B_5_c, coords_B_6_c);
        vs[61].grads(gxs[61], coords_B_5_c_g, coords_B_7_d_g, coords_B_5_c, coords_B_7_d);
        vs[62].grads(gxs[62], coords_B_5_c_g, coords_B_8_d_g, coords_B_5_c, coords_B_8_d);
        vs[63].grads(gxs[63], coords_B_6_c_g, coords_B_7_d_g, coords_B_6_c, coords_B_7_d);
        vs[64].grads(gxs[64], coords_B_6_c_g, coords_B_8_d_g, coords_B_6_c, coords_B_8_d);
        vs[65].grads(gxs[65], coords_B_7_d_g, coords_B_8_d_g, coords_B_7_d, coords_B_8_d);
        // cutoff based on largest d-OO
        gsw12 *= (sw13*sw14*sw23*sw24*sw34)*energies[j]/d12r;
        gsw13 *= (sw12*sw14*sw23*sw24*sw34)*energies[j]/d13r;
        gsw14 *= (sw12*sw13*sw23*sw24*sw34)*energies[j]/d14r;
        gsw23 *= (sw12*sw13*sw14*sw24*sw34)*energies[j]/d23r;
        gsw24 *= (sw12*sw13*sw14*sw23*sw34)*energies[j]/d24r;
        gsw34 *= (sw12*sw13*sw14*sw23*sw24)*energies[j]/d34r;

        // cutoff based on center-3-neighbor criteria
        //gsw12 *= (3*sw13*sw14*sw23*sw24*sw34 - sw14*sw23*sw24*sw34 - sw13*sw23*sw24*sw34 - sw13*sw14*sw24*sw34 - sw13*sw14*sw23*sw34 - sw13*sw14*sw23*sw24 + sw13*sw14 + sw23*sw24)*energies[j]/d12r;
        //gsw13 *= (3*sw12*sw14*sw23*sw24*sw34 - sw14*sw23*sw24*sw34 - sw12*sw23*sw24*sw34 - sw12*sw14*sw24*sw34 - sw12*sw14*sw23*sw34 - sw12*sw14*sw23*sw24 + sw12*sw14 + sw23*sw34)*energies[j]/d13r;
        //gsw14 *= (3*sw12*sw13*sw23*sw24*sw34 - sw13*sw23*sw24*sw34 - sw12*sw23*sw24*sw34 - sw12*sw13*sw24*sw34 - sw12*sw13*sw23*sw34 - sw12*sw13*sw23*sw24 + sw12*sw13 + sw24*sw34)*energies[j]/d14r;
        //gsw23 *= (3*sw12*sw13*sw14*sw24*sw34 - sw13*sw14*sw24*sw34 - sw12*sw14*sw24*sw34 - sw12*sw13*sw24*sw34 - sw12*sw13*sw14*sw34 - sw12*sw13*sw14*sw24 + sw12*sw24 + sw13*sw34)*energies[j]/d23r;
        //gsw24 *= (3*sw12*sw13*sw14*sw23*sw34 - sw13*sw14*sw23*sw34 - sw12*sw14*sw23*sw34 - sw12*sw13*sw23*sw34 - sw12*sw13*sw14*sw34 - sw12*sw13*sw14*sw23 + sw12*sw23 + sw14*sw34)*energies[j]/d24r;
        //gsw34 *= (3*sw12*sw13*sw14*sw23*sw24 - sw13*sw14*sw23*sw24 - sw12*sw14*sw23*sw24 - sw12*sw13*sw23*sw24 - sw12*sw13*sw14*sw24 - sw12*sw13*sw14*sw23 + sw13*sw23 + sw14*sw24)*energies[j]/d34r;


        for (size_t i = 0; i < 3; i++) {
            gradients[0 + i] += 0.0 + (gsw12*d12[i])+ (gsw13*d13[i])+ (gsw14*d14[i]);
            gradients[9 + i] += 0.0 - (gsw12*d12[i])+ (gsw23*d23[i])+ (gsw24*d24[i]);
            gradients[18 + i] += 0.0 - (gsw13*d13[i])- (gsw23*d23[i])+ (gsw34*d34[i]);
            gradients[27 + i] += 0.0 - (gsw14*d14[i])- (gsw24*d24[i])- (gsw34*d34[i]);
        }


        for (size_t i = 0; i < 9; i++) {
            grad1[i + j*9] += gradients[0 + i];
        }


        for (size_t i = 0; i < 9; i++) {
            grad2[i + j*9] += gradients[9 + i];
        }


        for (size_t i = 0; i < 9; i++) {
            grad3[i + j*9] += gradients[18 + i];
        }


        for (size_t i = 0; i < 9; i++) {
            grad4[i + j*9] += gradients[27 + i];
        }

        
        if (virial != 0) {
        
            (*virial)[0] += -coords_A_1_a[0]*coords_A_1_a_g[0]
                        -coords_B_1_a[0]*coords_B_1_a_g[0]
                        -coords_B_2_a[0]*coords_B_2_a_g[0]
                        -coords_A_2_b[0]*coords_A_2_b_g[0]
                        -coords_B_3_b[0]*coords_B_3_b_g[0]
                        -coords_B_4_b[0]*coords_B_4_b_g[0]
                        -coords_A_3_c[0]*coords_A_3_c_g[0]
                        -coords_B_5_c[0]*coords_B_5_c_g[0]
                        -coords_B_6_c[0]*coords_B_6_c_g[0]
                        -coords_A_4_d[0]*coords_A_4_d_g[0]
                        -coords_B_7_d[0]*coords_B_7_d_g[0]
                        -coords_B_8_d[0]*coords_B_8_d_g[0];

        
            (*virial)[1] += -coords_A_1_a[0]*coords_A_1_a_g[1]
                        -coords_B_1_a[0]*coords_B_1_a_g[1]
                        -coords_B_2_a[0]*coords_B_2_a_g[1]
                        -coords_A_2_b[0]*coords_A_2_b_g[1]
                        -coords_B_3_b[0]*coords_B_3_b_g[1]
                        -coords_B_4_b[0]*coords_B_4_b_g[1]
                        -coords_A_3_c[0]*coords_A_3_c_g[1]
                        -coords_B_5_c[0]*coords_B_5_c_g[1]
                        -coords_B_6_c[0]*coords_B_6_c_g[1]
                        -coords_A_4_d[0]*coords_A_4_d_g[1]
                        -coords_B_7_d[0]*coords_B_7_d_g[1]
                        -coords_B_8_d[0]*coords_B_8_d_g[1];

        
            (*virial)[2] += -coords_A_1_a[0]*coords_A_1_a_g[2]
                        -coords_B_1_a[0]*coords_B_1_a_g[2]
                        -coords_B_2_a[0]*coords_B_2_a_g[2]
                        -coords_A_2_b[0]*coords_A_2_b_g[2]
                        -coords_B_3_b[0]*coords_B_3_b_g[2]
                        -coords_B_4_b[0]*coords_B_4_b_g[2]
                        -coords_A_3_c[0]*coords_A_3_c_g[2]
                        -coords_B_5_c[0]*coords_B_5_c_g[2]
                        -coords_B_6_c[0]*coords_B_6_c_g[2]
                        -coords_A_4_d[0]*coords_A_4_d_g[2]
                        -coords_B_7_d[0]*coords_B_7_d_g[2]
                        -coords_B_8_d[0]*coords_B_8_d_g[2];

        
            (*virial)[4] += -coords_A_1_a[1]*coords_A_1_a_g[1]
                        -coords_B_1_a[1]*coords_B_1_a_g[1]
                        -coords_B_2_a[1]*coords_B_2_a_g[1]
                        -coords_A_2_b[1]*coords_A_2_b_g[1]
                        -coords_B_3_b[1]*coords_B_3_b_g[1]
                        -coords_B_4_b[1]*coords_B_4_b_g[1]
                        -coords_A_3_c[1]*coords_A_3_c_g[1]
                        -coords_B_5_c[1]*coords_B_5_c_g[1]
                        -coords_B_6_c[1]*coords_B_6_c_g[1]
                        -coords_A_4_d[1]*coords_A_4_d_g[1]
                        -coords_B_7_d[1]*coords_B_7_d_g[1]
                        -coords_B_8_d[1]*coords_B_8_d_g[1];

        
            (*virial)[5] += -coords_A_1_a[1]*coords_A_1_a_g[2]
                        -coords_B_1_a[1]*coords_B_1_a_g[2]
                        -coords_B_2_a[1]*coords_B_2_a_g[2]
                        -coords_A_2_b[1]*coords_A_2_b_g[2]
                        -coords_B_3_b[1]*coords_B_3_b_g[2]
                        -coords_B_4_b[1]*coords_B_4_b_g[2]
                        -coords_A_3_c[1]*coords_A_3_c_g[2]
                        -coords_B_5_c[1]*coords_B_5_c_g[2]
                        -coords_B_6_c[1]*coords_B_6_c_g[2]
                        -coords_A_4_d[1]*coords_A_4_d_g[2]
                        -coords_B_7_d[1]*coords_B_7_d_g[2]
                        -coords_B_8_d[1]*coords_B_8_d_g[2];

        
            (*virial)[8] += -coords_A_1_a[2]*coords_A_1_a_g[2]
                        -coords_B_1_a[2]*coords_B_1_a_g[2]
                        -coords_B_2_a[2]*coords_B_2_a_g[2]
                        -coords_A_2_b[2]*coords_A_2_b_g[2]
                        -coords_B_3_b[2]*coords_B_3_b_g[2]
                        -coords_B_4_b[2]*coords_B_4_b_g[2]
                        -coords_A_3_c[2]*coords_A_3_c_g[2]
                        -coords_B_5_c[2]*coords_B_5_c_g[2]
                        -coords_B_6_c[2]*coords_B_6_c_g[2]
                        -coords_A_4_d[2]*coords_A_4_d_g[2]
                        -coords_B_7_d[2]*coords_B_7_d_g[2]
                        -coords_B_8_d[2]*coords_B_8_d_g[2];

            (*virial)[3] = (*virial)[1];
            (*virial)[6] = (*virial)[2];
            (*virial)[7] = (*virial)[5];


        }


    }

    double energy = 0.0;
    for (size_t i = 0; i < n; i++)
        energy += energies_sw[i];

    return energy;

}

//----------------------------------------------------------------------------//
} // namespace mbnrg_A1B2_A1B2_A1B2_A1B2_deg4

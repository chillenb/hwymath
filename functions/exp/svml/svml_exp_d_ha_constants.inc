#ifndef SVML_TABLE  /* storage class of the tables (C++ users may define it as static constexpr) */
#define SVML_TABLE static const
#endif
/* vector path (__svml_dexp_ha_data_internal_avx512) */
SVML_TABLE double vTh[16] = {  /* 2^(j/16), rounded */
    0x1.0000000000000p+0,
    0x1.0b5586cf9890fp+0,
    0x1.172b83c7d517bp+0,
    0x1.2387a6e756238p+0,
    0x1.306fe0a31b715p+0,
    0x1.3dea64c123422p+0,
    0x1.4bfdad5362a27p+0,
    0x1.5ab07dd485429p+0,
    0x1.6a09e667f3bcdp+0,
    0x1.7a11473eb0187p+0,
    0x1.8ace5422aa0dbp+0,
    0x1.9c49182a3f090p+0,
    0x1.ae89f995ad3adp+0,
    0x1.c199bdd85529cp+0,
    0x1.d5818dcfba487p+0,
    0x1.ea4afa2a490dap+0,
};
SVML_TABLE double vTlr[16] = {  /* 2^(j/16) = Th[j] * (1 + Tlr[j]) */
    0x0.0p+0,
    0x1.79aa65d837b6dp-54,
    -0x1.01b15eaa59348p-55,
    0x1.68efde3a8a894p-54,
    0x1.34d754db0abb6p-55,
    0x1.59f48a72a4c6dp-55,
    0x1.690cebb7aafb0p-56,
    0x1.063e1e21c5409p-54,
    -0x1.3b3efbf5e2228p-54,
    -0x1.b32dcb94da51dp-56,
    0x1.db72fc1f0eab4p-55,
    0x1.1affc2b91ce27p-56,
    0x1.c1a7792cb3387p-55,
    0x1.36eae30af0cb3p-56,
    0x1.4a385a63d07a7p-56,
    -0x1.ff7128fd391f0p-55,
};
#define vLog2e       0x1.71547652b82fep+0
#define vShifter     0x1.8000000003ff0p+48
#define vLn2Hi       0x1.62e42fefa39efp-1
#define vLn2Lo       0x1.abc9e3b39803fp-56
#define vC6          0x1.7411836940c04p-10
#define vC5          0x1.1101cbbc265c0p-7
#define vC4          0x1.55557242d68fep-5
#define vC3          0x1.5555553939732p-3
#define vC2          0x1.000000000d008p-1
#define vC1          0x1.fffffffffff70p-1
#define vThreshold   0x1.61da04cbafe44p+9
#define vRMaskBits   0xbfffffffffffffffULL
#define vAbsMaskBits 0x7fffffffffffffffULL

/* scalar callout (_imldExpHATab) */
SVML_TABLE double sTab[128] = {  /* {hi, lo} of 2^(j/64), j = 0..63 */
    0x1.0000000000000p+0, 0x0.0p+0,
    0x1.02c9a3e778061p+0, -0x1.160139cd8dc5dp-56,
    0x1.059b0d3158574p+0, 0x1.cd2523567f613p-55,
    0x1.0874518759bc8p+0, 0x1.0f74e61e6c861p-57,
    0x1.0b5586cf9890fp+0, 0x1.79aa65d837b6dp-54,
    0x1.0e3ec32d3d1a2p+0, 0x1.ebe3d702f9cd1p-60,
    0x1.11301d0125b51p+0, -0x1.556522a2fbd0ep-54,
    0x1.1429aaea92de0p+0, -0x1.1c923b9d5f416p-54,
    0x1.172b83c7d517bp+0, -0x1.01b15eaa59348p-55,
    0x1.1a35beb6fcb75p+0, 0x1.b898c3f1353bfp-55,
    0x1.1d4873168b9aap+0, 0x1.aecf73e3a2f60p-54,
    0x1.2063b88628cd6p+0, 0x1.a6f4144a6c38dp-55,
    0x1.2387a6e756238p+0, 0x1.68efde3a8a894p-54,
    0x1.26b4565e27cddp+0, 0x1.0472b981fe7f2p-55,
    0x1.29e9df51fdee1p+0, 0x1.2f7e16d09ab31p-55,
    0x1.2d285a6e4030bp+0, 0x1.b3782720c0ab4p-55,
    0x1.306fe0a31b715p+0, 0x1.34d754db0abb6p-55,
    0x1.33c08b26416ffp+0, 0x1.fdd395dd3f84ap-55,
    0x1.371a7373aa9cbp+0, -0x1.24aedcc4b5068p-54,
    0x1.3a7db34e59ff7p+0, -0x1.1d1e83e9436d2p-56,
    0x1.3dea64c123422p+0, 0x1.59f48a72a4c6dp-55,
    0x1.4160a21f72e2ap+0, -0x1.8a78f4817895bp-58,
    0x1.44e086061892dp+0, 0x1.363ed60c2ac11p-59,
    0x1.486a2b5c13cd0p+0, 0x1.ecce1daa10379p-57,
    0x1.4bfdad5362a27p+0, 0x1.690cebb7aafb0p-56,
    0x1.4f9b2769d2ca7p+0, -0x1.f94340071a38ep-55,
    0x1.5342b569d4f82p+0, -0x1.8dec6bd0f385fp-56,
    0x1.56f4736b527dap+0, 0x1.3350518fdd78ep-54,
    0x1.5ab07dd485429p+0, 0x1.063e1e21c5409p-54,
    0x1.5e76f15ad2148p+0, 0x1.432e62b64c035p-54,
    0x1.6247eb03a5585p+0, -0x1.c33c53bef4da8p-55,
    0x1.6623882552225p+0, -0x1.3cedd78565858p-54,
    0x1.6a09e667f3bcdp+0, -0x1.3b3efbf5e2228p-54,
    0x1.6dfb23c651a2fp+0, -0x1.367efb86da9eep-57,
    0x1.71f75e8ec5f74p+0, -0x1.81f647e5a3ecfp-56,
    0x1.75feb564267c9p+0, -0x1.619321e55e68ap-55,
    0x1.7a11473eb0187p+0, -0x1.b32dcb94da51dp-56,
    0x1.7e2f336cf4e62p+0, 0x1.5ebe1abd66c55p-57,
    0x1.82589994cce13p+0, -0x1.369b6f13b3734p-54,
    0x1.868d99b4492edp+0, -0x1.4d450d872576ep-54,
    0x1.8ace5422aa0dbp+0, 0x1.db72fc1f0eab4p-55,
    0x1.8f1ae99157736p+0, 0x1.bf68359f35f44p-56,
    0x1.93737b0cdc5e5p+0, -0x1.da9b88b6c1e29p-58,
    0x1.97d829fde4e50p+0, -0x1.2434322f4f9aap-54,
    0x1.9c49182a3f090p+0, 0x1.1affc2b91ce27p-56,
    0x1.a0c667b5de565p+0, -0x1.7c50422622263p-55,
    0x1.a5503b23e255dp+0, -0x1.1bbd1d3bcbb15p-54,
    0x1.a9e6b5579fdbfp+0, 0x1.469846e735ab3p-55,
    0x1.ae89f995ad3adp+0, 0x1.c1a7792cb3387p-55,
    0x1.b33a2b84f15fbp+0, -0x1.5c3d956dcaebap-58,
    0x1.b7f76f2fb5e47p+0, -0x1.8d6f438ad9334p-57,
    0x1.bcc1e904bc1d2p+0, 0x1.4ffd70a5fddcdp-56,
    0x1.c199bdd85529cp+0, 0x1.36eae30af0cb3p-56,
    0x1.c67f12e57d14bp+0, 0x1.4e08fd10959acp-55,
    0x1.cb720dcef9069p+0, 0x1.76b2c6c921968p-57,
    0x1.d072d4a07897cp+0, -0x1.fad5d3ffffa6fp-55,
    0x1.d5818dcfba487p+0, 0x1.4a385a63d07a7p-56,
    0x1.da9e603db3285p+0, 0x1.e5a50d5c192acp-55,
    0x1.dfc97337b9b5fp+0, -0x1.2d52107b43e1fp-55,
    0x1.e502ee78b3ff6p+0, 0x1.4b604603a88d3p-56,
    0x1.ea4afa2a490dap+0, -0x1.ff7128fd391f0p-55,
    0x1.efa1bee615a27p+0, 0x1.ec3bc41aa2008p-55,
    0x1.f50765b6e4540p+0, 0x1.a64a931d185eep-55,
    0x1.fa7c1819e90d8p+0, 0x1.7893b4d91cd9dp-56,
};
#define s64OverLn2   0x1.71547652b82fep+6
#define sShifter     0x1.8000000000000p+52
#define sC2          0x1.0000000000000p-1
#define sC3          0x1.55555555548f8p-3
#define sC4          0x1.5555555558fccp-5
#define sC5          0x1.111123aaf20d3p-7
#define sC6          0x1.6c16a1c2a3ffdp-10
#define sOverflowX   0x1.62e42fefa39efp+9
#define sDenormX     -0x1.6232bdd7abcd2p+9
#define sUnderflowX  -0x1.74910d52d3051p+9
#define sLn2o64Hi    0x1.62e42fefa0000p-7
#define sLn2o64Lo    0x1.cf79abc9e3b3ap-46
#define sTiny        0x1.0000000000001p-1022
#define sHuge        0x1.fffffffffffffp+1023
#define sZero        0x0.0p+0
#define sOne         0x1.0000000000000p+0
#define sTwo         0x1.0000000000000p+1
#define sTwoM60      0x1.0000000000000p-60
#define sSplitter    0x1.8000000000000p+32

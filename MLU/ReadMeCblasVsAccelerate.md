# Summary of findings

1. Results using Apple Accelerate framework agree with those using GSL CBLAS.
2. Apple Accelerate framework gives ~3% performance improvement over GSL CBLAS.
3. Previous version was correctly using Apple Accelerate.

Relative differences between Apple Accelerate framework and GSL CBLAS appear in:

- central replica ~1e-5
- other replicas ~ 1e-3

I.e. results are compatible to first 5 digits.
Given solver precision set to 7 digits, this seems quite good. 


# Correctness testing

After changing to use CBLAS interface directly (12 Apr 2024), compared results against prior '`golden`' version of code which:

- used GSL's high level BLAS implementation
- linked with Apple Accelerate framework
- was run on 13 Feb 2024.

Tests were performed using both Apple Accelerate framework and GSL CBLAS.

## Test results

### GSL CBLAS vs Golden (prior) results

Comparison command:

    h5diff -p 1e-6 -rn 10 -c Cont.{gslcblas,golden}/linear/Simul/renormE3-CZ34/F3_K_Ds.corr_f0_fplus.g5P_g5W.model.h5; echo $?

Results

    dataset: </model/Summary> and </model/Summary>
    size:           [2x216]           [2x216]
    position        Summary         Summary         difference      relative       
    ------------------------------------------------------------------------
    [ 0 29 ]          0.709774075     0.7097750452    9.702654662e-07 1.367006066e-06
    [ 0 30 ]          -0.03085269777  -0.03086094262  8.244852662e-06 0.0002672327951
    [ 0 30 ]          0.2251301748    0.2251305074    3.326359659e-07 1.477527241e-06
    [ 0 30 ]          0.3261056545    0.3261051517    5.028088268e-07 1.541858658e-06
    [ 0 30 ]          0.3906623032    0.3906634446    1.141398485e-06 2.921701111e-06
    [ 0 31 ]          -0.2083200794   -0.2083197858   2.936196843e-07 1.409464153e-06
    [ 0 31 ]          -0.144820188    -0.1448210991   9.111370115e-07 6.291505515e-06
    [ 0 46 ]          0.01457666076   0.01454349765   3.316311013e-05 0.002275082796
    [ 0 46 ]          0.2206262635    0.2206265669    3.033032697e-07 1.374737825e-06
    [ 0 46 ]          0.3050880381    0.3050876527    3.854609612e-07 1.263441738e-06
    [ 0 46 ]          0.3602888849    0.3602896036    7.187421802e-07 1.994905229e-06
    11 differences found
    dataset: </model/data_C> and </model/data_C>
    size:           [216]           [216]
    position        data_C          data_C          difference      relative       
    ------------------------------------------------------------------------
    [ 30 ]          0.3261056545    0.3261051517    5.028088268e-07 1.541858658e-06
    [ 31 ]          -0.2083200794   -0.2083197858   2.936196843e-07 1.409464153e-06
    [ 46 ]          0.3050880381    0.3050876527    3.854609612e-07 1.263441738e-06
    [ 47 ]          -0.1505596136   -0.150559401    2.125867755e-07 1.411977425e-06
    4 differences found
    dataset: </model/data_S> and </model/data_S>
    size:           [10000x216]           [10000x216]
    position        data_S          data_S          difference      relative       
    ------------------------------------------------------------------------
    [ 0 30 ]          0.3773719079    0.377374073     2.165121649e-06 5.737368374e-06
    [ 0 31 ]          -0.122279615    -0.1222807469   1.131910585e-06 9.256739848e-06
    [ 0 46 ]          0.3284664676    0.3284679993    1.531708786e-06 4.663212038e-06
    [ 0 47 ]          -0.06748072153  -0.06748159128  8.697496609e-07 1.288886131e-05
    [ 0 49 ]          0.9438220371    0.94382065      1.387102138e-06 1.469664919e-06
    [ 0 50 ]          -1.216784049    -1.216782207    1.841788462e-06 1.5136527e-06 
    [ 0 51 ]          0.4514879545    0.4514871437    8.107657355e-07 1.795763824e-06
    [ 1 30 ]          0.1912562733    0.191255435     8.383815978e-07 4.383550841e-06
    [ 1 31 ]          -0.1589258479   -0.1589248354   1.012494791e-06 6.37086292e-06
    [ 1 46 ]          0.2227919997    0.2227913962    6.034569489e-07 2.708611394e-06
    10 differences found
    1

### Apple Accelerate vs Golden (prior) results

Comparison command:

    h5diff -c Cont.{Accelerate,golden}/linear/Simul/renormE3-CZ34/F3_K_Ds.corr_f0_fplus.g5P_g5W.model.h5; echo $?

Results

    0

### GSL CBLAS vs Apple Accelerate

Comparison command:

    h5diff -p 1e-6 -rn 10 -c Cont.{gslcblas,Accelerate}/linear/Simul/renormE3-CZ34/F3_K_Ds.corr_f0_fplus.g5P_g5W.model.h5; echo $?

Results

    dataset: </model/Summary> and </model/Summary>
    size:           [2x216]           [2x216]
    position        Summary         Summary         difference      relative       
    ------------------------------------------------------------------------
    [ 0 29 ]          0.709774075     0.7097750452    9.702654662e-07 1.367006066e-06
    [ 0 30 ]          -0.03085269777  -0.03086094262  8.244852662e-06 0.0002672327951
    [ 0 30 ]          0.2251301748    0.2251305074    3.326359659e-07 1.477527241e-06
    [ 0 30 ]          0.3261056545    0.3261051517    5.028088268e-07 1.541858658e-06
    [ 0 30 ]          0.3906623032    0.3906634446    1.141398485e-06 2.921701111e-06
    [ 0 31 ]          -0.2083200794   -0.2083197858   2.936196843e-07 1.409464153e-06
    [ 0 31 ]          -0.144820188    -0.1448210991   9.111370115e-07 6.291505515e-06
    [ 0 46 ]          0.01457666076   0.01454349765   3.316311013e-05 0.002275082796
    [ 0 46 ]          0.2206262635    0.2206265669    3.033032697e-07 1.374737825e-06
    [ 0 46 ]          0.3050880381    0.3050876527    3.854609612e-07 1.263441738e-06
    [ 0 46 ]          0.3602888849    0.3602896036    7.187421802e-07 1.994905229e-06
    11 differences found
    dataset: </model/data_C> and </model/data_C>
    size:           [216]           [216]
    position        data_C          data_C          difference      relative       
    ------------------------------------------------------------------------
    [ 30 ]          0.3261056545    0.3261051517    5.028088268e-07 1.541858658e-06
    [ 31 ]          -0.2083200794   -0.2083197858   2.936196843e-07 1.409464153e-06
    [ 46 ]          0.3050880381    0.3050876527    3.854609612e-07 1.263441738e-06
    [ 47 ]          -0.1505596136   -0.150559401    2.125867755e-07 1.411977425e-06
    4 differences found
    dataset: </model/data_S> and </model/data_S>
    size:           [10000x216]           [10000x216]
    position        data_S          data_S          difference      relative       
    ------------------------------------------------------------------------
    [ 0 30 ]          0.3773719079    0.377374073     2.165121649e-06 5.737368374e-06
    [ 0 31 ]          -0.122279615    -0.1222807469   1.131910585e-06 9.256739848e-06
    [ 0 46 ]          0.3284664676    0.3284679993    1.531708786e-06 4.663212038e-06
    [ 0 47 ]          -0.06748072153  -0.06748159128  8.697496609e-07 1.288886131e-05
    [ 0 49 ]          0.9438220371    0.94382065      1.387102138e-06 1.469664919e-06
    [ 0 50 ]          -1.216784049    -1.216782207    1.841788462e-06 1.5136527e-06 
    [ 0 51 ]          0.4514879545    0.4514871437    8.107657355e-07 1.795763824e-06
    [ 1 30 ]          0.1912562733    0.191255435     8.383815978e-07 4.383550841e-06
    [ 1 31 ]          -0.1589258479   -0.1589248354   1.012494791e-06 6.37086292e-06
    [ 1 46 ]          0.2227919997    0.2227913962    6.034569489e-07 2.708611394e-06
    10 differences found
    1

# Benchmarking `ContVariations.sh`

Benchmark 12 Apr 2024 on:

- Apple Mac Studio
- Chip M1 Ultra
- 20 CPU
- 48 GPU

## GSL CBLAS

That is, link with `-l gslcblas`

    StudioM:NoSync mike$ time ContVariations.sh
    Making fit variations cubic shrink linear

    real  1m34.646s
    user  26m34.836s
    sys   1m8.696s
    StudioM:NoSync mike$ time Make= ContCompare.sh

    real  0m34.636s
    user  0m27.039s
    sys   0m5.270s

## Apple Accelerate

That is, link with `-framework Accelerate`

    StudioM:NoSync mike$ time ContVariations.sh
    Making fit variations cubic shrink linear

    real  1m31.307s
    user  25m45.026s
    sys   1m8.523s
    StudioM:NoSync mike$ time Make= ContCompare.sh

    real  0m34.834s
    user  0m27.118s
    sys   0m5.362s

## Conclusion

Apple Accelerate gives ~3% performance improvement.

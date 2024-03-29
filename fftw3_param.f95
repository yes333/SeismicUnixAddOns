! FFTW3 constant and related global parameters

module fftw3_constants

      INTEGER(4), PARAMETER ::  FFTW_R2HC=0
      INTEGER(4), PARAMETER ::  FFTW_HC2R=1
      INTEGER(4), PARAMETER ::  FFTW_DHT=2
      INTEGER(4), PARAMETER ::  FFTW_REDFT00=3
      INTEGER(4), PARAMETER ::  FFTW_REDFT01=4
      INTEGER(4), PARAMETER ::  FFTW_REDFT10=5
      INTEGER(4), PARAMETER ::  FFTW_REDFT11=6
      INTEGER(4), PARAMETER ::  FFTW_RODFT00=7
      INTEGER(4), PARAMETER ::  FFTW_RODFT01=8
      INTEGER(4), PARAMETER ::  FFTW_RODFT10=9
      INTEGER(4), PARAMETER ::  FFTW_RODFT11=10
      INTEGER(4), PARAMETER ::  FFTW_FORWARD=-1
      INTEGER(4), PARAMETER ::  FFTW_BACKWARD=+1
      INTEGER(4), PARAMETER ::  FFTW_MEASURE=0
      INTEGER(4), PARAMETER ::  FFTW_DESTROY_INPUT=1
      INTEGER(4), PARAMETER ::  FFTW_UNALIGNED=2
      INTEGER(4), PARAMETER ::  FFTW_CONSERVE_MEMORY=4
      INTEGER(4), PARAMETER ::  FFTW_EXHAUSTIVE=8
      INTEGER(4), PARAMETER ::  FFTW_PRESERVE_INPUT=16
      INTEGER(4), PARAMETER ::  FFTW_PATIENT=32
      INTEGER(4), PARAMETER ::  FFTW_ESTIMATE=64
      INTEGER(4), PARAMETER ::  FFTW_ESTIMATE_PATIENT=128
      INTEGER(4), PARAMETER ::  FFTW_BELIEVE_PCOST=256
      INTEGER(4), PARAMETER ::  FFTW_NO_DFT_R2HC=512
      INTEGER(4), PARAMETER ::  FFTW_NO_NONTHREADED=1024
      INTEGER(4), PARAMETER ::  FFTW_NO_BUFFERING=2048
      INTEGER(4), PARAMETER ::  FFTW_NO_INDIRECT_OP=4096
      INTEGER(4), PARAMETER ::  FFTW_ALLOW_LARGE_GENERIC=8192
      INTEGER(4), PARAMETER ::  FFTW_NO_RANK_SPLITS=16384
      INTEGER(4), PARAMETER ::  FFTW_NO_VRANK_SPLITS=32768
      INTEGER(4), PARAMETER ::  FFTW_NO_VRECURSE=65536
      INTEGER(4), PARAMETER ::  FFTW_NO_SIMD=131072
      INTEGER(4), PARAMETER ::  FFTW_NO_SLOW=262144
      INTEGER(4), PARAMETER ::  FFTW_NO_FIXED_RADIX_LARGE_N=524288
      INTEGER(4), PARAMETER ::  FFTW_ALLOW_PRUNING=1048576
      INTEGER(4), PARAMETER ::  FFTW_WISDOM_ONLY=2097152

end module fftw3_constants

module fftw3_global

    INTEGER(8), SAVE    :: fftw_plan_fw = 0
    INTEGER(8), SAVE    :: fftw_plan_bw = 0

    contains

    subroutine destroy_plan()
        if ( fftw_plan_fw /= 0 ) call sfftw_destroy_plan(fftw_plan_fw)
        if ( fftw_plan_bw /= 0 ) call sfftw_destroy_plan(fftw_plan_bw)
    end subroutine destroy_plan
    
end module fftw3_global

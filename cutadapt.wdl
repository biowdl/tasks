task cutadapt {
    String? format
    Int? cores = 1
    String? adapter3end
    String? adapter5end
    String? adapterAnywhere
    Float? errorRate
    Boolean? noIndels
    Int? count
    Int? minLength
    Boolean? matchReadWildcards
    Boolean? noMatchAdapterWildcards
    Boolean? noTrim
    String? maskAdapter
    Int? cutLength

}
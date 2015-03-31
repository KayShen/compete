ds <- function(M){
    n= nrow(M)
    ds_p = DS(M)
    n_ds_p = NormDS(ds_p)
    ds_d = Dyadic(M)
    n_ds_d = NormDS(ds_d)
    r = order(n_ds_p)[n:1]
    r2 = order(n_ds_d)[n:1]
    return(list(rank_p = r, rank_d = r2))
}

proportion <- function(M){
    M_temp = M + t(M)
    M_temp[M_temp == 0] = 0.01
    M_p = M/M_temp
    return(M_p)
}

DS <- function(M){
    M_p = proportion(M)
    w = rowSums(M_p)
    l = colSums(M_p)
    M_p_weighted = t(w*t(M_p))
    M_p_weighted2 = l*M_p
    w2 = rowSums(M_p_weighted)
    l2 = colSums(M_p_weighted2)
    result = w + w2 - l - l2
    return(result)
}

NormDS <- function(DS){
    n = length(DS)
    result = (DS + n*(n-1)/2)/n
    return(result)
}

Dyadic <- function(M){
    M_p = proportion(M)
    pp = M + t(M)
    pp2 = 1/(pp)
    pp = 1/(1+pp)
    pp[pp2==Inf] = 0
    D = M_p - ((M_p-0.5)*pp)
    M_p = D
    w = rowSums(M_p)
    l = colSums(M_p)
    M_p_weighted = t(w*t(M_p))
    M_p_weighted2 = l*M_p
    w2 = rowSums(M_p_weighted)
    l2 = colSums(M_p_weighted2)
    result = w + w2 - l - l2
    return(result)
}


# Inigo Martincorena - 2019
# EM algorithm to identify clusters of variants from VAF using a beta-binomial model.
#
# verbose=T outputs more information on each iteration
# min_vafdist>0 would not accept adding a new cluster if it is closer than this VAF distance from another cluster
#
# Example (simulation)
# n = 50; p_vec = c(0.21,0.79); m_vec = c(0.1,0.5); cov_vec = rep(80,n);
# x = rbinom(n=n, size=cov_vec, prob=sample(x=m_vec, size=n, prob=p_vec, replace=T)); counts_table = cbind(x,cov_vec);
# clusters = vafclusters_mixbbin(counts_table, num_runs=50, bb_rho=1e-4)

vafclusters_mixbbin = function(counts_table, bb_rho=1e-4, max_clusters=min(nrow(counts_table),10), num_runs=100, verbose=F, min_vafdist=0) {
    
    # Subfunction mixture of binomials
    mixbbinEM = function(counts_table, k) {
        
        bestLL = -Inf
        
        for (run in 1:num_runs) {
            
            maxiter = 100; mindist = 0.0001
            m_est = sort(runif(k,0,0.5)) # Random initial guess
            count = 0; m_est_prev = rep(-1,k); # Initialise
            
            while (count<maxiter & sum(abs(m_est-m_est_prev))>mindist) {
                count = count+1; m_est_prev = m_est
                # 1. Expectation
                probs = array(NA,dim=c(dim(counts_table)[1],k))
                for (h in 1:k) {
                    probs[,h] = VGAM::dbetabinom(x=counts_table[,1], size=counts_table[,2], prob=m_est[h], rho=bb_rho)
                }
                probs = pmax(probs,1e-20)
                probs = probs/apply(probs,1,sum)
                # 2. Maximisation (updating the estimates of the means)
                for (h in 1:k) {
                    m_est[h] = sum(counts_table[,1]*probs[,h])/sum(counts_table[,2]*probs[,h])
                }
            }
            
            # Log-likelihood
            weights = apply(probs,2,sum)/sum(probs)
            loglik = array(NA,dim=c(nrow(counts_table),k))
            for (h in 1:k) {
                loglik[,h] = log(weights[h])+VGAM::dbetabinom(x=counts_table[,1], size=counts_table[,2], prob=m_est[h], rho=bb_rho, log=T)
            }
            LL = sum(apply(loglik,1,function(x_row) log(sum(exp(x_row - max(x_row)))) + max(x_row) ))
            
            # Saving if best solution so far
            if (LL>bestLL) {
                bestLL = LL
                params = list(m_est=m_est, weights=weights, LL=LL, probs=probs)
            }
        }
        return(params)
    }
    
    # Finding the optimal number of clusters using an iterative LRT
    LL0 = -Inf
    for (k in 1:max_clusters) {
        em = mixbbinEM(counts_table, k)
        pval = 1-pchisq(2*(em$LL-LL0),2)
        if (verbose) {
            print(sprintf("%0.0f clusters: LL=%0.3f, pval=%0.3g", k, em$LL, pval)) # Printing the progress
        }
        if (pval<0.01 & !any(diff(sort(em$m_est))<min_vafdist)) {
            best_EM = em; LL0 = em$LL
        } else {
            break
        }
    }
    return(best_EM)
}

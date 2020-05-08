scale_01 = function(x) {
  (x - min(x)) / (max(x) - min(x))
}

scale_01_mat = function(m) {
  (m - min(c(m))) / (max(c(m)) - min(c(m)))
}

calc_geom_mean = function(x) {
  x = na.omit(x)
  prod(x) ^ (1 / length(x))
}

calc_cv = function(x) {
  sd(x) / abs(mean(x))
}

calc_mean_sd = function(x, name) {
  d = data.frame(mean(x,na.rm=T), sd(x,na.rm=T))
  colnames(d) = paste(name, c("mean", "sd"), sep="_")
  d
}

calc_kl = function(x,y) {
  sum(x * log2(x/y))
}

calc_jsd = function(x,y) {
  x = x[x>0]
  y = y[y>0]
  m = (x + y) / 2
  sqrt(.5 * calc_kl(x,m) + .5 * calc_kl(y,m))
}

calc_eu = function(x, y) {
  sqrt(sum((x - y)^2))
}

sem = function(x) {
  x = na.omit(x)
  sd(x)/sqrt(length(x))
}

mean_na = function(x) {
  x = x[is.finite(x)]
  mean(x, na.rm=T)
}

median_na = function(x) {
  x = x[is.finite(x)]
  median(x, na.rm=T)
}

iqr = function(x) {
  qs = quantile(na.omit(x), probs = c(.25, .75))
  diff(qs)
}

sd_na = function(x) {
  x = x[is.finite(x)]
  sd(x, na.rm=T)
}

sem_na = function(x) {
  x = na.omit(x)
  sd(x)/length(x)
}

#' https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
#' bootstrap validation 
weighted_sem = function(x, w) {
  na_ind = is.na(x)
  x = x[!na_ind]
  w = w[!na_ind]
  w1 = (w - min(w)) / (max(w) - min(w))
  ws = sum(w1)
  n = length(x)
  std_var = (n / ((n-1) * ws^2)) * sum(w1^2 * (x - mean(x))^2)
  sqrt(std_var)
}

var_na = function(x) {
  x = x[is.finite(x)]
  var(x, na.rm=T)
}
q75 = function(x) {
  x = na.omit(x)
  quantile(x, probs=c(.75))
}

zscore_clean = function(d) {
  d = na.omit(d)
  (d - mean(d)) / sd(d)
}
zscore_my = function(d, drel=NULL, default_sd=1, baseline=NULL) {
  if (is.null(drel) & !is.null(baseline)) {
    drel = d[baseline]
  }
  if (length(drel) < 2) {
    sd_value = default_sd
  } else {
    sd_value = sd(drel, na.rm=T)    
  }
  mean_value = mean(drel, na.rm=T)
  (d - mean_value) / sd_value
}

frac_to_perc = function(x) {
  x * 100
}

perc_change = function(d, drel=NULL, baseline=NULL) {
  if (is.null(drel) & !is.null(baseline)) {
    drel = d[baseline]
  }
  baseline_mean =  mean(drel, na.rm=T)
  100 * (d - baseline_mean) / baseline_mean
}

calc_diff = function(d, drel=NULL, baseline=NULL) {
  if (is.null(drel) & !is.null(baseline)) {
    drel = d[baseline]
  }
  baseline_mean =  mean(drel, na.rm=T)
  d - baseline_mean
}

#' Tissue specificty score
#'
#' @param x 
#'
#' @return scalar
#' @export
#' @references http://www.sciencedirect.com/science/article/pii/S0092867410000796
#' @examples
tissue_specificity = function(x, scale=F) {
  if (scale)
    x = (x - (min(x))) / (max(x) - (min(x)))
  x = x[x>0]
  x_norm = x /sum(x)
  #sum(x_norm * log(x_norm / (1/length(x))))
  sum(x_norm * log(x_norm / mean(x)))
}

calc_gene_specificity = function(x, scale=F) {
  if (scale)
    x = (x - (min(x))) / (max(x) - (min(x)))
  #x = x[x>0]
  x_norm = x /sum(x)
  x_norm * log(x_norm / mean(x_norm))
}

# calc_tau = function(x) {
#   x_hat = x / max(x)
#   tau = 
# }

calc_specificity = function(mat) {
   t(apply(mat, 1, function(x) x / mean(x)))
}

specificity_correlate = function(mat_a, mat_b, method="spearman") {
  
  mat_a_spec = calc_specificity(mat_a)
  mat_b_spec = na.omit(calc_specificity(mat_b))
  #mat_a_spec = t(apply(mat_a, 1, function(x)  calc_gene_specificity(x)))
  #mat_b_spec = na.omit(t(apply(mat_b, 1, function(x) calc_gene_specificity(x))))
  #mat_a_spec = t(apply(mat_a, 1, scale, scale=F))
  #mat_b_spec = na.omit(t(apply(mat_b, 1, scale, scale=F)))
  #colnames(mat_a_spec) = colnames(mat_a)
  #colnames(mat_b_spec) = colnames(mat_b)
  #mat_a_spec = log(mat_a_spec)
  #mat_b_spec = log(mat_b_spec)
  mat_a_spec = mat_a_spec[match(rownames(mat_b_spec), rownames(mat_a_spec)),]
  obj_cor = cor(mat_a_spec, mat_b_spec, method=method)
  obj_cor
}

cor_sig_matrix = function(mat_a, mat_b, method="pearson", adjust_method="BH", symmetric=T) {
  obj_cor_sig = matrix(NA, 
                       nrow = ncol(mat_a),
                       ncol=ncol(mat_b), 
                       dimnames = list(colnames(mat_a), 
                                       colnames(mat_b)))
  
  if (!symmetric) {
    i_range = 1:ncol(mat_a)
    j_range = 1:ncol(mat_b)
  }
  if (symmetric) {
    i_range = 1:(ncol(mat_a)-1)
    j_range = (i+1):ncol(mat_b)
  }
  for (i in i_range ) {
    for (j in j_range) {
      obj_cor_sig[i,j] = cor.test(mat_a[,i], mat_b[,j], method=method, exact=F)$p.value
    }
  }
  obj_cor_sig_df = na.omit(melt(obj_cor_sig))
  obj_cor_sig_df = obj_cor_sig_df %>%
    mutate(value_adj = p.adjust(value, method=adjust_method))
  obj_cor_sig_mat = acast(obj_cor_sig_df, Var1~Var2, value.var="value_adj", fill=1)
  obj_cor_sig_mat
}





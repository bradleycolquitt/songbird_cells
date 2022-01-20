region_correlate = function(mat_a, mat_b, method="spearman") {
  mat_a_spec = t(apply(mat_a, 1, function(x) x / mean(x)))
  mat_b_spec = na.omit(t(apply(mat_b, 1, function(x) x / mean(x))))
  mat_a_spec = mat_a_spec[match(rownames(mat_b_spec), rownames(mat_a_spec)),]
  obj_cor = cor(mat_a_spec, mat_b_spec, method=method)
  obj_cor
}

region_calc_summary = function(obj_cor_df,
                               taxon_md_df,
                               cor_quantile_filter = .99, 
                               value_n_quantile_filter = .95,
                               grouping_factor_compare = "pallial",
                               grouping_factor_data = "position2",
                               grouping_factor_region = "acronym"
                               ) {
  
  grouping_factor_compare_sym = as.symbol(grouping_factor_compare)
  grouping_factor_data_sym = as.symbol(grouping_factor_data)
  grouping_factor_region_sym = as.symbol(grouping_factor_region)
  grouping1 = vars(!!grouping_factor_compare_sym, !!grouping_factor_data_sym)

  cor_quantile_filter_sym = as.symbol(sprintf("value_%s", sub("0\\.", "q", as.character(cor_quantile_filter) )))
  
  obj_cor_df_filt_sum = obj_cor_df %>%
    filter(value > !!cor_quantile_filter_sym) %>% 
    #left_join(finch_samp) %>%
    left_join(taxon_md_df %>% select(!!grouping_factor_compare_sym, !!grouping_factor_region_sym, acronym)) %>% 
    group_by_at(grouping1)
  if (grouping_factor_region!="acronym") {
    obj_cor_df_filt_sum = obj_cor_df_filt_sum %>%
      group_by(!!grouping_factor_region_sym, add = T)
  }
  obj_cor_df_filt_sum = obj_cor_df_filt_sum %>%
    summarize(value_mean = mean(value),
      value_n_n = n()) %>%
    ungroup() %>% 
    arrange(!!grouping_factor_data_sym) %>% 
    mutate(!!grouping_factor_compare_sym := factor(!!grouping_factor_compare_sym, levels=unique(as.data.frame(taxon_md_df)[,grouping_factor_compare]))) %>%
    mutate(!!grouping_factor_data_sym := factor(!!grouping_factor_data_sym, levels=unique(obj_cor_df[,grouping_factor_data]))) %>%
    #select(!!grouping_factor_compare_sym, !!grouping_factor_data_sym, value_n_n, value_mean) %>% 
    complete(!!grouping_factor_compare_sym, !!grouping_factor_data_sym, fill=list(value_n_n=0, value_mean=0)) %>%
    na.omit()

  nrep= 100 
  obj_cor_shuf_sum = map(1:nrep, function(i) {
    
    shuf_df = obj_cor_df %>% 
      ungroup() %>% 
      distinct(!!grouping_factor_region_sym) %>%
      mutate(tmp = sample(!!grouping_factor_region_sym))
    obj_cor_shuf_sum = obj_cor_df %>% 
      ungroup() %>%
      left_join(shuf_df) %>%
      select(-!!grouping_factor_region_sym) %>% 
      rename(!!grouping_factor_region_sym := tmp) %>% 
      ungroup() %>% 
      #left_join(finch_samp) %>%
      left_join(taxon_md_df %>% select(!!grouping_factor_compare_sym, acronym, !!grouping_factor_region_sym)) %>%
      filter(value > !!cor_quantile_filter_sym) %>% 
      group_by_at(grouping1) %>%
      summarize(value_mean = mean(value),
        value_n_n = n()) %>%
      ungroup() %>% 
      arrange(!!grouping_factor_data_sym) %>% 
      mutate(!!grouping_factor_compare_sym := factor(!!grouping_factor_compare_sym, levels=unique(as.data.frame(taxon_md_df)[,grouping_factor_compare]))) %>%
      mutate(!!grouping_factor_data_sym := factor(!!grouping_factor_data_sym, levels=unique(obj_cor_df[,grouping_factor_data]))) %>%
      complete(!!grouping_factor_compare_sym, !!grouping_factor_data_sym, fill=list(value_n_n=0, value_mean=0)) %>% 
      na.omit() %>%
      mutate(rep = i)

    obj_cor_shuf_sum
    
  }) %>% bind_rows()
  
  obj_cor_shuf_sum_stat = obj_cor_shuf_sum %>%
    group_by(!!grouping_factor_data_sym, !!grouping_factor_compare_sym) %>%
    summarize(value_mean_shuf = mean(value_mean),
      value_n_n_q95 = quantile(value_n_n, value_n_quantile_filter),
      value_n_n_med = median(value_n_n)) 
  
  obj_cor_df_filt_sum3 = obj_cor_df_filt_sum %>% left_join(obj_cor_shuf_sum_stat) %>%
    mutate(value_n_n_fc_shuf = value_n_n / value_n_n_med )
  
  obj_cor_df_filt_sum3
}
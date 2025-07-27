library(Biostrings)
library(ggplot2)
source("scripts/handle_EVO2_output_v1.R") # needed for initialize_df in evo2_analysis_functions_v1.R
source("scripts/evo2_analysis_functions_v1.R") # importing initialize_df 

all_variants <- readDNAStringSet("input/double_genes.fasta") 
seq_names <- names(all_variants)
"seq_names: double_83_ref double_83_TTT double_random"

# don't reuse has specific file name
get_dfs <- function(){
  dfs <- list()
  for (name in seq_names){
    npy_fn <- paste0("input/input_", name, "_logits.npy")
    dfs[[name]] <- initialize_df(npy_fn, name, all_variants)
  }
  return(dfs)
}

all_gene_dfs <- get_dfs()
# Print summary of each data frame
for (name in names(all_gene_dfs)) {
  cat("\nFirst 6 rows of", name, ":\n")
  print(head(all_gene_dfs[[name]]))
}

# returns vector with first - second half of both variants
get_2h_lls <- function(all_gene_dfs){
  # initialize first half vec
  firsth_ll_vec <- numeric(length(all_gene_dfs))
  # firsth_ll_vec
  names(firsth_ll_vec) <- names(all_gene_dfs)
  names(firsth_ll_vec)
  # initialize second half vec
  secondh_ll_vec <- numeric(length(all_gene_dfs))
  # secondh_ll_vec
  names(secondh_ll_vec) <- paste0("secondh_", names(all_gene_dfs))
  names(secondh_ll_vec)
  
  for (i in 1:length(all_gene_dfs)) { 
    firsth_ll_vec[i] <- sum(all_gene_dfs[[i]]$log_likelihood[1:2625])
    secondh_ll_vec[i] <- sum(all_gene_dfs[[i]]$log_likelihood[(1+2625):(2625+2625)])
  }
  return(secondh_ll_vec)
}

# Function to create difference data frame
get_diff_df <- function(df1, df2) {
  # Create a new data frame with the difference in log likelihoods
  diff_df <- df1
  diff_df$log_likelihood <- df1$log_likelihood - df2$log_likelihood
  return(diff_df)
}

dfs <- get_dfs()
get_2h_lls(dfs)

##################################
### Making pdfs ###################
##################################
# Create difference data frame
ref_diff_df <- get_diff_df(dfs$double_83_ref, dfs$double_random)
mut_diff_df <- get_diff_df(dfs$double_83_TTT, dfs$double_random)

#diff plots
save_as_pdf(plot_log_likelihood("double_ref-random", ref_diff_df, index_rows = 1:5250, 
            highlighted = 2626:2628), "double_ref-random.pdf", out_dir = "figures", width = 12, height = 5)
save_as_pdf(plot_log_likelihood("double_mut-random", mut_diff_df, index_rows = 1:5250, 
            highlighted = 2626:2628), "double_mut-random.pdf", out_dir = "figures", width = 12, height = 5)

# plotting mutation of second part of gene
save_as_pdf(plot_log_likelihood("83_TTT_2nd", dfs$double_83_TTT, index_rows = (2625+200):(2625+300), 
            highlighted = (2625+247):(2625+249)), "83_TTT_2nd.pdf", out_dir = "figures", width = 12, height = 5)
save_as_pdf(plot_log_likelihood("83_ref_2nd", dfs$double_83_ref, index_rows = (2625+200):(2625+300), 
            highlighted = (2625+247):(2625+249)), "83_ref_2nd.pdf", out_dir = "figures", width = 12, height = 5)

# plotting transversion of second part of gene
save_as_pdf(plot_log_likelihood("trans_83_TTT", dfs$double_83_TTT, index_rows = (2620):(2600+100),
          highlighted = 2626:2628), "trans_83_TTT.pdf", out_dir = "figures", width = 12, height = 5)
save_as_pdf(plot_log_likelihood("trans_83_S1", dfs$double_83_ref, index_rows = (2620):(2600+100), 
            highlighted = 2626:2628), "trans_83_S1.pdf", out_dir = "figures", width = 12, height = 5)
save_as_pdf(plot_log_likelihood("trans_random", dfs$double_random, index_rows = (2620):(2600+100), 
            highlighted = 2626:2628), "trans_random.pdf", out_dir = "figures", width = 12, height = 5)        

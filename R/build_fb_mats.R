build_fb_mats <- function(epimodel) {
          
          # get indices for tpm products corresponding to inter-observation time tpms
          .interobs_tpm_inds <- which(epimodel$config_mat[,"ID"] == 0)
          
          # calculate the emission probability-weighted initial distribution
          pi_0 <- normalize(epimodel$.initdist * epimodel$.emission_mat[,1, drop = FALSE])
          
          # compute the first fb matrix

          epimodel$.fb_mats[[1]] <- normalize((pi_0 %*% t(epimodel$.emission_mat[, 2, drop = FALSE])) * epimodel$.tpm_products[[1]])
          
          for(m in 2:length(epimodel$.fb_mats)) {
                    
                    epimodel$.fb_mats[[m]] <- normalize((colSums(epimodel$.fb_mats[[m-1]]) %*% t(epimodel$.emission_mat[, m+1, drop = FALSE])) * epimodel$.tpm_products[[.interobs_tpm_inds[m]]])
                    
          }
          
}
#!/usr/bin/env python

"""
PRS combination weight for MIX.

"""

import numpy as np

def compute_regression_variable(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, n_prs_beta, indep_approx):
    beta_ldblk_beta = {}
    beta_beta_hat = {}

    # Initialize sum variables
    total_beta_ldblk_beta = None
    total_beta_beta_hat = None

    for chrom in range(1, 23):
        beta_std_dict[chrom] = beta_std_dict[chrom].reshape(-1, 1)  # Ensure correct shape

        beta_ldblk_beta[chrom] = np.zeros((n_prs_beta, n_prs_beta))
        beta_beta_hat[chrom] = np.zeros((n_prs_beta, 1))

        n_blk = len(ld_blk[chrom])
        mm = 0  # Start index for LD blocks

        for kk in range(n_blk):  # Iterate over LD blocks
            if blk_size[chrom][kk] == 0:
                continue  # Skip empty blocks
            
            # Get block index range
            idx_blk = np.arange(mm, mm + blk_size[chrom][kk])  # Convert range to NumPy array
            
            # Compute regression variables
            if indep_approx == 'TRUE':
                beta_ldblk_beta[chrom] += prs_beta_matrix[chrom][idx_blk, :].T @ prs_beta_matrix[chrom][idx_blk, :]
            else:
                beta_ldblk_beta[chrom] += prs_beta_matrix[chrom][idx_blk, :].T @ ld_blk[chrom][kk] @ prs_beta_matrix[chrom][idx_blk, :]

            beta_beta_hat[chrom] += prs_beta_matrix[chrom][idx_blk, :].T @ beta_std_dict[chrom][idx_blk, :]

            # Move to the next block
            mm += blk_size[chrom][kk]

        # Aggregate results across chromosomes
        if total_beta_ldblk_beta is None:
            total_beta_ldblk_beta = beta_ldblk_beta[chrom]
            total_beta_beta_hat = beta_beta_hat[chrom]
        else:
            total_beta_ldblk_beta += beta_ldblk_beta[chrom]
            total_beta_beta_hat += beta_beta_hat[chrom]

    return total_beta_ldblk_beta, total_beta_beta_hat


def linear_regression(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, n_prs_beta, indep_approx, non_negative_weights, pop, out_dir, out_name):
    print('... MIX linear regression weight calculating ...')

    beta_ldblk_beta, beta_beta_hat = compute_regression_variable(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, n_prs_beta, indep_approx)

    weights = np.linalg.solve(beta_ldblk_beta, beta_beta_hat)

    if non_negative_weights == 'TRUE':
        weights = np.maximum(weights,0)
        weights = weights / weights.sum()
        weight_file = out_dir + '/' + '%s_%s_linear_weights_approx%s_non_negative.txt' % (out_name, pop, indep_approx)
    else:
        weight_file = out_dir + '/' + '%s_%s_linear_weights_approx%s.txt' % (out_name, pop, indep_approx)
        
    with open(weight_file, 'w') as ff:
        ff.write("\t".join(f"{w:.6e}" for w in weights.flatten()))
            
    print(f"Weights saved to: {weight_file}")


def compute_aic_bic(ld_blk, blk_size, beta_std_dict, prs_beta_matrix_subset, n_selected_prs_beta, n_gwas, indep_approx):
    # Compute OLS estimates using the feature subset
    beta_ldblk_beta, beta_beta_hat = compute_regression_variable(ld_blk, blk_size, beta_std_dict, prs_beta_matrix_subset, n_selected_prs_beta, indep_approx)
    weights = np.linalg.solve(beta_ldblk_beta, beta_beta_hat)

    # Compute AIC: AIC = 2d + n log(RSS/n)
    aic = 2 * n_selected_prs_beta + n_gwas * np.log(1 - 2 * weights.T @ beta_beta_hat + weights.T @ beta_ldblk_beta @ weights)

    # Compute BIC: BIC = 2d + n log(RSS/n)
    bic = n_selected_prs_beta * np.log(n_gwas) + n_gwas * np.log(1 - 2 * weights.T @ beta_beta_hat + weights.T @ beta_ldblk_beta @ weights)

    return weights, aic, bic


def linear_regression_with_forward_selection(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, n_prs_beta, n_gwas, indep_approx, selection_criterion, max_features, non_negative_weights, pop, out_dir, out_name):
    print('... MIX linear regression with forward selection weight calculating ...')

    remaining_features = list(range(n_prs_beta))  # Indices of all predictors
    selected_features = []  # Best feature subset
    current_evaluation = float('inf')  # Start with a large AIC/BIC value

    while remaining_features:
        evaluation_scores = []

        prs_beta_matrix_subset = {}

        for feature in remaining_features:

            for chrom in range(1, 23):

                prs_beta_matrix_subset[chrom] = np.column_stack([prs_beta_matrix[chrom][:, f] for f in (selected_features + [feature])])
            
            n_selected_prs_beta = len(selected_features) + 1

            # Compute AIC and BIC
            weights, aic_value, bic_value = compute_aic_bic(ld_blk, blk_size, beta_std_dict, prs_beta_matrix_subset, n_selected_prs_beta, n_gwas, indep_approx)

            if selection_criterion == "AIC":
                evaluation_value = aic_value
            elif selection_criterion == "BIC":
                evaluation_value = bic_value

            # require positive or not
            if non_negative_weights == 'TRUE':
                if np.all(weights >= 0):
                    evaluation_scores.append((evaluation_value, feature))
            else:
                evaluation_scores.append((evaluation_value, feature))

        if not evaluation_scores:
            break

        # Select feature with the lowest AIC/BIC among valid candidates
        evaluation_scores.sort()
        best_evaluation, best_feature = evaluation_scores[0]

        # If AIC/BIC improves, update selected features
        if best_evaluation < current_evaluation:
            current_evaluation = best_evaluation
            selected_features.append(best_feature)
            remaining_features.remove(best_feature)
        else:
            break  # Stop if no improvement

        # Stop if max features are reached
        if max_features and len(selected_features) >= max_features:
            break

    # Compute final OLS with selected features
    prs_beta_matrix_final = {}

    if selected_features:

        for chrom in range(1, 23):

            prs_beta_matrix_final[chrom] = np.column_stack([prs_beta_matrix[chrom][:, f] for f in selected_features])

        n_selected_prs_beta = len(selected_features)
        beta_final, _, _ = compute_aic_bic(ld_blk, blk_size, beta_std_dict, prs_beta_matrix_final, n_selected_prs_beta, n_gwas, indep_approx)

    else:
        beta_final = np.array([])  # No features selected

    # Assign regression coefficients to selected features
    weights = np.zeros(n_prs_beta)
    weights[selected_features] = beta_final.flatten()

    # Save weights to file
    if non_negative_weights == 'TRUE':
        weights = weights / weights.sum()
        weight_file = out_dir + '/' + '%s_%s_linear_weights_with_forward_selection_%s_maxfeatures%s_approx%s_non_negative.txt' % (out_name, pop, selection_criterion, max_features, indep_approx)
    else:
        weight_file = out_dir + '/' + '%s_%s_linear_weights_with_forward_selection_%s_maxfeatures%s_approx%s.txt' % (out_name, pop, selection_criterion, max_features, indep_approx)

    with open(weight_file, 'w') as ff:
        ff.write("\t".join(f"{w:.6e}" for w in weights.flatten()))

    print(f"Weights saved to: {weight_file}")


def compute_non_negative_weights(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, indep_approx, active_set, prev_weights, tol):
    weights_active_full = prev_weights.copy()
    local_active_set = active_set.copy()

    while True:
        # Build the subproblem for the *current* local_active_set
        n_selected_prs_beta = np.sum(local_active_set)

        # Subset the matrix
        prs_beta_matrix_subset = {}
        for chrom in range(1, 23):
            sub = prs_beta_matrix[chrom][:, local_active_set]
            # Ensure it's 2D
            if sub.ndim == 1:
                sub = sub.reshape(-1, 1)
            prs_beta_matrix_subset[chrom] = sub
            
        # Solve the unconstrained least-squares subproblem
        beta_ldblk_beta, beta_beta_hat = compute_regression_variable(ld_blk, blk_size, beta_std_dict, prs_beta_matrix_subset, n_selected_prs_beta, indep_approx)
        weights_active_sub = np.linalg.solve(beta_ldblk_beta, beta_beta_hat)

        # Check if any coefficients are negative
        negative_mask = (weights_active_sub < -tol)
        if not np.any(negative_mask):
            # No negativity => clamp small negatives => done
            weights_active_sub = np.clip(weights_active_sub, 0, None)
            # Embed in the full solution
            weights_active_full[local_active_set] = weights_active_sub
            weights_active_full[~local_active_set] = 0.0
            return weights_active_full, local_active_set
        
        # Partial line search: remove negative coefficients
        neg_idx = np.flatnonzero(negative_mask)
        prev_sol_sub = weights_active_full[local_active_set]

        alpha_vals = prev_sol_sub[neg_idx] / (prev_sol_sub[neg_idx] - weights_active_sub[neg_idx])
        alpha_min = np.min(alpha_vals)
        alpha_min_idx = np.flatnonzero(np.isclose(alpha_vals, alpha_min, atol=tol))

        stepped_sol_sub = prev_sol_sub + alpha_min*(weights_active_sub - prev_sol_sub)
        stepped_sol_sub[neg_idx[alpha_min_idx]] = 0.0

        # Update weights and remove them from local_active_set
        weights_active_full[local_active_set] = stepped_sol_sub
        weights_active_full[~local_active_set] = 0.0

        local_inds_now = np.flatnonzero(local_active_set)
        zeroed_global_inds = local_inds_now[neg_idx[alpha_min_idx]]
        local_active_set[zeroed_global_inds] = False


def non_negative_linear_regression(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, n_prs_beta, n_gwas, indep_approx, pop, out_dir, out_name, tol=1e-10, max_iter=20):
    print('... MIX non-negative least squares regression calculating ...')

    # Initialize weights, active set
    weights = np.zeros((n_prs_beta,1))
    active_set = np.zeros(n_prs_beta, dtype=bool)
    
    # Obtain linear_weights for further comparison
    beta_ldblk_beta, beta_beta_hat = compute_regression_variable(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, n_prs_beta, indep_approx)
    linear_weights = np.linalg.solve(beta_ldblk_beta, beta_beta_hat)
    linear_weights = np.maximum(linear_weights,0)
    linear_weights_normalize = linear_weights / linear_weights.sum()
    RMSE_linear = 1 - 2 * linear_weights_normalize.T @ beta_beta_hat + linear_weights_normalize.T @ beta_ldblk_beta @ linear_weights_normalize

    for iter_num in range(max_iter):
        # Compute gradient
        gradient = n_gwas * (beta_beta_hat - beta_ldblk_beta @ weights)
        
        # Find most correlated inactive feature (largest positive gradient)
        inactive_gradient = gradient.flatten() * (~active_set)  # Only consider inactive features
        if np.all(inactive_gradient < -tol):  # Stop if no positive gradients remain
            break
        idx = np.argmax(inactive_gradient) 
        active_set[idx] = True

        # Compute new weights
        new_weights, new_active_set = compute_non_negative_weights(ld_blk, blk_size, beta_std_dict, prs_beta_matrix, indep_approx, active_set, weights, tol)

        # Update weights and activate_set for the next iteration
        weights = new_weights.copy()
        active_set = new_active_set.copy()

    # Decide simple linear regression or NNLS
    weights_normalize = weights / weights.sum()
    RMSE_NNLS = 1 - 2 * weights_normalize.T @ beta_beta_hat + weights_normalize.T @ beta_ldblk_beta @ weights_normalize

    if RMSE_linear.item() < RMSE_NNLS.item():
        final_weights = linear_weights
    else:
        final_weights = weights

    # Normalize the final weights and save to file
    final_weights = final_weights / final_weights.sum()  # Normalize final weights
    weight_file = out_dir + '/' + '%s_%s_non_negative_linear_weights_approx%s.txt' % (out_name, pop, indep_approx)

    with open(weight_file, 'w') as ff:
        ff.write("\t".join(f"{w:.6e}" for w in final_weights.flatten()))

    print(f"Weights saved to: {weight_file}")

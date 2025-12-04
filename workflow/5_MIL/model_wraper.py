import torch
from data import MyDataset,DataLoader

def train(model, loader_pos, loader_neg, device, optimizer, epochs=100, top_k=6, margin=4.0):
    """
    Train the model using top_k instance selection for both positive and negative bags.
    
    For each positive bag and for each negative bag, the function:
      - Computes model scores.
      - Selects the top_k scores from the bag.
      - Pairs the top_k scores (assuming they are sorted in descending order) and computes 
        a hinge loss: loss = max(0, (negative_score - positive_score + margin))
      - Sums the loss over the k pairs and (optionally) adds a small L1 regularization.
    """
    import torch
    all_losses = []
    lambda_reg = 0.09#001 # best 0.05

    for e in range(epochs):
        total_loss = 0.0
        model.train()

        # Loop over positive bags
        for pbag in loader_pos:
            pbag = pbag.float().to(device)
            # pbag[0] is assumed to be the tensor of instances in the bag
            p_scores = model(pbag[0])
            # Select top_k positive scores (largest scores)
            p_top, _ = torch.topk(p_scores.ravel(), k=top_k, largest=False)
            pos_loss = torch.tensor(0.0, device=device)
            
            # Loop over negative bags
            for nbag in loader_neg:
                nbag = nbag.float().to(device)
                n_scores = model(nbag[0])
                # Select top_k negative scores (largest scores)
                n_top, _ = torch.topk(n_scores.ravel(), k=top_k, largest=True)
                
                # Compute the hinge loss for each paired top-k instance.
                # Here we assume that the ordering of p_top and n_top is meaningful,
                # i.e. the i-th highest instance in the positive bag is paired with
                # the i-th highest instance in the negative bag.
                # The loss for each pair is: max(0, n_top[i] - p_top[i] + margin)
                losses = torch.clamp(n_top - p_top + margin, min=0)
                loss_value = losses.sum()  # or losses.mean() if you prefer averaging
                
                # Optionally add L1 regularization (applied on weights only)
                l1_reg = torch.tensor(0.0, device=device)
                for name, param in model.named_parameters():
                    if 'weight' in name:
                        l1_reg += torch.norm(param, 1)
                
                pos_loss += loss_value + lambda_reg * l1_reg
            
            optimizer.zero_grad()
            pos_loss.backward() 
            optimizer.step()
            total_loss += pos_loss.item()

        all_losses.append(total_loss)
        print(f'Epoch {e+1}/{epochs}, Loss: {total_loss:.8f}')

# def train(model, loader_pos,loader_neg,device,optimizer,epochs = 100):
#     all_losses = []
#     lambda_reg = 1e-5
#     for e in range(epochs):
#         total_loss = 0.0
#         model.train()

#         # Positive bag loop
#         for pbag in loader_pos:
#             pbag = pbag.float().to(device)
#             p_scores = model(pbag[0])
#             #max_p = torch.median(p_scores)

#             #UPDATE1
#             max_p = torch.median(p_scores).unsqueeze(0)
#             pos_loss = torch.tensor(0.0, device=device)
#             #Negative bags
#             for nbag in loader_neg:
#                 nbag = nbag.float().to(device)
#                 n_scores = model(nbag[0])
#                 max_n = torch.median(n_scores).unsqueeze(0)
#                 z = torch.zeros(1, device=device)
#                 loss = torch.max(z, (max_n - max_p)+2)

#                 l2_reg = torch.tensor(0.0, device=device)
#                 for name,param in model.named_parameters():
#                     if 'weight' in name:
#                         l2_reg += torch.norm(param,1)  # L2 regularization term (sum of squared weights)

#                 # Add regularization to pos_loss
#                 #pos_loss += lambda_reg * l2_reg
#                 pos_loss += loss[0]
#             optimizer.zero_grad()
#             pos_loss.backward() 
#             optimizer.step()
#             total_loss += pos_loss.item()

#         all_losses.append(total_loss)
#         print(f'Epoch {e+1}/{epochs}, Loss: {total_loss:.8f}')


def test(model, data_dict, wsi_pred_dict, itr_dict, fidx=0, test_indices=None, top_k = 5):
    """
    Test the model using the data stored in the dictionary and update the provided dictionaries directly.
    
    Args:
    - model: The trained model to be evaluated.
    - data_dict: Dictionary containing test data such as bags, sample IDs, labels, and region IDs.
    - wsi_pred_dict: Dictionary to store whole-slide image (WSI) level predictions.
    - itr_dict: Dictionary to store individual region (IT region) level predictions.
    - fidx: Fold index for cross-validation (default: 0).

    Returns:
    - None (The results are updated directly into wsi_pred_dict and itr_dict).
    """
    
    # Extract required data from the dictionary
    if test_indices is not None:
        bags_ts = data_dict['bags'][test_indices]
        samp_idx = data_dict['samples_id'][test_indices]
        samp_lbls = data_dict['labels'][test_indices]
        its_idx = data_dict['samples_its_id'][test_indices]
    else:
        bags_ts = data_dict['bags']
        samp_idx = data_dict['samples_id']
        samp_lbls = data_dict['labels']
        its_idx = data_dict['samples_its_id']

    # Testing loop
    model.eval()  # Put the model in evaluation mode
    test_dataset = MyDataset(bags_ts)
    loader_ts = DataLoader(test_dataset, batch_size=1)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    with torch.no_grad():  # No need to calculate gradients during testing
        for idx, tsbag in enumerate(loader_ts):
            tsbag = tsbag.float().to(device)
            scores = model(tsbag[0])

            effects = tsbag[0]*model.out.weight
            effects = effects.cpu().numpy()
            
            # Store the max score for predictions (sample level)
            sample_id = samp_idx[idx]
            top_scores, _ = torch.topk(scores.ravel(), k=top_k, largest=True)
            max_score = float(torch.mean(top_scores))
            
            # Update WSI-level predictions
            if sample_id not in wsi_pred_dict:
                wsi_pred_dict[sample_id] = {
                    'scores': [],
                    'labels': [],
                    'fold_idx': []
                }
            wsi_pred_dict[sample_id]['scores'] = max_score
            wsi_pred_dict[sample_id]['labels'] = samp_lbls[idx]
            wsi_pred_dict[sample_id]['fold_idx'] = fidx
            # Update IT region level predictions
            scores = scores.cpu()
            for it_idx, it_id in enumerate(its_idx[idx]):
                if it_id not in itr_dict:
                    itr_dict[it_id] = {
                        'scores': [],
                        'effects': [],
                        'labels': [],
                        'fold_idx': []
                    }
                itr_dict[it_id]['scores'] = float(scores[it_idx])
                itr_dict[it_id]['labels'] = samp_lbls[idx]
                itr_dict[it_id]['fold_idx'] = fidx
                itr_dict[it_id]['effects'] = effects[it_idx]
    

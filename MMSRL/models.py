import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
from sklearn.cluster import KMeans
import torch.optim as optim
from random import shuffle
import pandas as pd
import numpy as np
import scanpy as sc
from functools import partial
from MMSRL.layers import GraphConvolution,SelfAttention,MLP,MGCN,decoder
import sys

# Define a function to compute the Symmetric Contrastive Embedding (SCE) loss
def sce_loss(x, y, alpha=3):
    # Apply L2 normalization to inputs x and y
    x = F.normalize(x, p=2, dim=-1)
    y = F.normalize(y, p=2, dim=-1)

    # Compute the SCE loss: (1 - x·y)^alpha
    # The previously commented-out loss formulation used dot product and Euclidean distance
    loss = (1 - (x * y).sum(dim=-1)).pow_(alpha)

    # Return the mean loss
    loss = loss.mean()
    return loss


# Define the MMSRL model class, which inherits from torch.nn.Module
class MMSRL(nn.Module):
    def __init__(self, nfeatX, nfeatI, hidden_dims):
        super(MMSRL, self).__init__()

        # Initialize various layers, including MGCN, SelfAttention, fully connected layer, and decoder
        self.mgcn = MGCN(nfeatX, nfeatI, hidden_dims)
        self.attlayer1 = SelfAttention(dropout=0.1)  # Attention layer for scRNA-seq data
        self.attlayer2 = SelfAttention(dropout=0.1)  # Attention layer for scATAC data
        self.fc = nn.Linear(hidden_dims[1] * 2, hidden_dims[1])  # Fully connected layer
        self.mlp = MLP(hidden_dims[1], dropout_rate=0.1)  # MLP layer
        self.ZINB = decoder(hidden_dims[1], nfeatX)  # Decoder for ZINB modeling

        # Initialize mask token for random masking
        self.enc_mask_token = nn.Parameter(torch.zeros(1, nfeatX))  # Input features are 3000-dimensional node features; initialized as zero vector
        self._mask_rate = 0.8  # Masking ratio
        self.criterion = self.setup_loss_fn(loss_fn='sce')  # Set loss function (default: SCE loss)


    # Set up the loss function
    def setup_loss_fn(self, loss_fn, alpha_l=3):
        if loss_fn == "mse":
            criterion = nn.MSELoss()  # Mean Squared Error loss
        elif loss_fn == "sce":
            criterion = partial(sce_loss, alpha=3)  # Symmetric Contrastive Embedding (SCE) loss
        else:
            raise NotImplementedError  # Raise error if the loss function is not supported
        return criterion

    # Define the forward propagation process
    def forward(self, x, i, a):
        # Apply masking to input data and add noise
        a, x, (mask_nodes, keep_nodes) = self.encoding_mask_noise(a, x, self._mask_rate)

        emb_x, emb_i = self.mgcn(x, i, a)

        E_G = emb_x
        E_H = emb_i

        ## Attention for omics-specific information in scRNA-seq
        att_weights_x, att_emb_x = self.attlayer1(emb_x, emb_x, emb_x)
        Z_G = att_emb_x

        ## Attention for omics-specific information in scATAC
        att_weights_i, att_emb_i = self.attlayer2(emb_i, emb_i, emb_i)
        Z_H = att_emb_i

        q_x, q_i = self.mlp(emb_x, emb_i)

        # cl_loss = crossview_contrastive_Loss(q_x, q_i)

        # Capture consistency information
        emb_con = torch.cat([q_x, q_i], dim=1)
        z_xi = self.fc(emb_con)

        Z_F = z_xi

        z_I = 20 * att_emb_x + 1 * att_emb_i + 10 * z_xi
        Z = z_I

        [defeat_f, pi, disp, mean] = self.ZINB(z_I)


        # Construct reconstruction loss
        recon = defeat_f.clone()  #
        x_init = x[mask_nodes]   # Original data at masked positions
        x_rec = recon[mask_nodes]   # Reconstructed data at masked positions

        # Compute loss
        loss = self.criterion(x_rec, x_init)  #
        # loss = self.weights_mask_loss * loss

        gcn_G = emb_x
        gcn_H = emb_i

        return z_I, q_x, q_i, pi, disp, mean, loss, E_G, E_H, Z_G, Z_H, Z_F, Z
    

    # Randomly mask some nodes in the input data; masked node values are replaced by the mask token.
    # mask_rate=0.3 means 30% of the nodes will be masked.
    def encoding_mask_noise(self, adj, x, mask_rate=0.3):
        num_nodes = adj.shape[0]  # Get the number of nodes; adj.shape[0] equals the total number of nodes in the graph
        perm = torch.randperm(num_nodes, device=x.device)  # Randomly permute node indices

        # Random masking
        num_mask_nodes = int(mask_rate * num_nodes)  # Compute the number of nodes to mask
        mask_nodes = perm[:num_mask_nodes]  # Select the first 'num_mask_nodes' from the shuffled indices as masked nodes
        keep_nodes = perm[num_mask_nodes:]  # Remaining nodes (to be kept)

        out_x = x.clone()  # Deep copy the input features to avoid modifying the original data
        token_nodes = mask_nodes  # Masked nodes
        # out_x[mask_nodes] = 0.0  # Previously, masked node values were set to zero

        # enc_mask_token is a learnable parameter whose value is updated during training
        out_x[token_nodes] += self.enc_mask_token  # Add mask token to masked node features
        use_adj = adj.clone()  # Clone the adjacency matrix

        return use_adj, out_x, (mask_nodes, keep_nodes)



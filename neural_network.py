from itertools import *

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class Highway(nn.Module):

    def __init__(self, size, num_layers, f):
        super(Highway, self).__init__()

        self.num_layers = num_layers
        self.nonlinear = nn.ModuleList([nn.Linear(size, size) for _ in range(num_layers)])
        self.linear = nn.ModuleList([nn.Linear(size, size) for _ in range(num_layers)])
        self.gate = nn.ModuleList([nn.Linear(size, size) for _ in range(num_layers)])
        self.f = f

    def forward(self, x):
        for layer in range(self.num_layers):
            gate = F.sigmoid(self.gate[layer](x))
            nonlinear = self.f(self.nonlinear[layer](x))
            linear = self.linear[layer](x)
            x = gate * nonlinear + (1 - gate) * linear
        return x


class mvml_survival_Analysis(nn.Module):

    def __init__(self):
        nn.Module.__init__(self)
        self.miRNA_encoding = nn.Sequential(nn.Linear(1881, 80),
                                            nn.Tanh(),
                                            nn.BatchNorm1d(80),
                                            Highway(80, 5, f=F.relu),
                                            )

        self.RNASeq_encoding = nn.Sequential(nn.Linear(17720, 80),
                                             nn.Tanh(),
                                             nn.BatchNorm1d(80),
                                             Highway(80, 5, f=F.relu),
                                             )

        self.miRNA_BN = nn.BatchNorm1d(80)
        self.RNASeq_BN = nn.BatchNorm1d(80)
        self.survival_weight = nn.Linear(80, 1, bias=False)

    def get_RNASeq_feature(self, gene):
        gene = gene.view(gene.shape[0], -1)
        # gene = F.tanh(self.fcg(gene))
        # gene = self.bn1_fcg(gene)
        # gene = self.fcg_highway(gene)
        # # gene = F.dropout(gene, 0.3)
        # # gene =F.sigmoid(self.bn2_fcg(gene))
        # return gene
        gene_encoding = F.tanh(self.RNASeq_BN(self.RNASeq_encoding(gene)))
        return gene_encoding

    def get_miRNA_feature(self, miRNA):
        miRNA = miRNA.view(miRNA.shape[0], -1)
        # microRNA = F.tanh(self.fcm(microRNA))
        # microRNA = self.bn1_fcm(microRNA)
        # microRNA = self.fcm_highway(microRNA)
        # # microRNA  = F.dropout(microRNA, 0.3)
        # # microRNA_feature =F.sigmoid(self.bn2_fcm(microRNA))
        # return microRNA
        miRNA_encoding = F.tanh(self.miRNA_BN(self.miRNA_encoding(miRNA)))
        return miRNA_encoding

    def get_survival_result(self, gene, miRNA):
        RNASeq_feature = self.get_RNASeq_feature(gene)
        miRNA_feature = self.get_miRNA_feature(miRNA)
        return self.survival_weight((RNASeq_feature + miRNA_feature) / 2)

    def get_survival_feature(self, gene, miRNA):
        RNASeq_feature = self.get_RNASeq_feature(gene)
        miRNA_feature = self.get_miRNA_feature(miRNA)
        return (RNASeq_feature + miRNA_feature) / 2

    def get_cox_weight(self):
        return self.survival_weight.weight

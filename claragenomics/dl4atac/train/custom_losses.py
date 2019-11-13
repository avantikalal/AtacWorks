#
# Copyright (c) 2019, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#

import torch
import torch.nn as nn
import torch.nn.functional as F
from claragenomics.dl4atac.train.metrics import CorrCoef

class PearsonLoss(nn.Module):

    def __init__(self):
        super(PearsonLoss, self).__init__()

    def forward(self, input, targets):
        r = CorrCoef()(input, targets)
        r_loss = 1 - r
        return r_loss

class FocalBCELoss(nn.Module):

        def __init__(self):
            super(FocalBCELoss, self).__init__()
        
        def forward(self, input, targets):
            self.gamma = 2
            logpt = -F.binary_cross_entropy(input, targets, reduction='none')
            #print(logpt)
            pt = torch.exp(logpt)
            #print(pt)
            loss_func = -((1-pt)**self.gamma) * logpt
            #print(loss_func)
            loss_func = loss_func.mean()
            #print(loss_func)
            return loss_func

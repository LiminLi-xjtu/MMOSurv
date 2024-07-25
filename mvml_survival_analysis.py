import sys

sys.path.append("../..")
import os
import numpy as np
# %matplotlib inline
import torch
import random
from torch.autograd import Variable
from operator import add
import time
import argparse
import json
from model.neural_network import mvml_survival_Analysis


def do_base_learning(mvml_model, x_RNASeq_batch, x_miRNA_batch, R_matrix_batch, ystatus_batch,
                     TRAIN_INNER_LR, TRAIN_N_INNER, REG_SCALE):
    new_mvml_model = mvml_survival_Analysis()
    if torch.cuda.is_available():
        new_mvml_model = new_mvml_model.cuda()

    new_mvml_model.load_state_dict(mvml_model.state_dict())  # copy? looks okay
    new_mvml_model.train()
    # inner_optimizer = torch.optim.Adam([{'params': new_mvml_model.parameters()}], lr=TRAIN_INNER_LR, weight_decay=REG_SCALE)
    inner_optimizer = torch.optim.SGD(new_mvml_model.parameters(), lr=TRAIN_INNER_LR, weight_decay=REG_SCALE)
    x_RNASeq_batch = torch.tensor(x_RNASeq_batch, dtype=torch.float)
    x_miRNA_batch = torch.tensor(x_miRNA_batch, dtype=torch.float)
    R_matrix_batch = torch.tensor(R_matrix_batch, dtype=torch.float)
    ystatus_batch = torch.tensor(ystatus_batch, dtype=torch.float)
    real_batch_size = x_RNASeq_batch.shape[0]
    if torch.cuda.is_available():
        x_RNASeq_batch = x_RNASeq_batch.cuda()
        x_miRNA_batch = x_miRNA_batch.cuda()
        R_matrix_batch = R_matrix_batch.cuda()
        ystatus_batch = ystatus_batch.cuda()
    for _ in range(TRAIN_N_INNER):
        j = np.random.choice(range(real_batch_size), real_batch_size)
        feature_RNASeq = new_mvml_model.get_RNASeq_feature(x_RNASeq_batch)
        feature_miRNA = new_mvml_model.get_miRNA_feature(x_miRNA_batch)

        norm_gene = torch.sqrt(torch.sum(torch.mul(feature_RNASeq, feature_RNASeq), dim=1))
        norm_miRNA = torch.sqrt(torch.sum(torch.mul(feature_miRNA, feature_miRNA), dim=1))
        same_sim = torch.sum(torch.mul(feature_RNASeq, feature_miRNA), 1) / torch.mul(norm_gene, norm_miRNA)
        diff_sim = torch.sum(torch.mul(feature_RNASeq, feature_miRNA[j,]), 1) / torch.mul(norm_gene, norm_miRNA[j])
        cos_loss = torch.mean(torch.max(torch.tensor(0.0).cuda(), 0.5 - (same_sim - diff_sim)))

        theta_mv = new_mvml_model.get_survival_result(x_RNASeq_batch, x_miRNA_batch)
        exp_theta_mv = torch.reshape(torch.exp(theta_mv), [x_RNASeq_batch.shape[0]])
        theta_mv = torch.reshape(theta_mv, [x_RNASeq_batch.shape[0]])

        loss_mvml = -torch.mean(torch.mul((theta_mv - torch.log(torch.sum(torch.mul(exp_theta_mv,
                                                                                    R_matrix_batch), dim=1))),
                                          torch.reshape(ystatus_batch, [x_RNASeq_batch.shape[0]])))

        loss = loss_mvml + 1 * cos_loss

        inner_optimizer.zero_grad()
        loss.backward()
        inner_optimizer.step()
    with open('nohup_loss' + '.log', 'a') as f:
        f.writelines("loss_mvml:" + str(loss_mvml.cpu()) + '\n')
    return new_mvml_model


def do_final_learning(mvml_model, x_RNASeq_test_support, x_miRNA_test_support, ytime_test_support,
                      ystatus_test_support, x_RNASeq_test_qeury, x_miRNA_test_qeury, ytime_test_qeury,
                      ystatus_test_qeury, TRAIN_LR, REG_SCALE, META_STEP):
    new_mvml_model = mvml_survival_Analysis()
    if torch.cuda.is_available():
        new_mvml_model = new_mvml_model.cuda()
    new_mvml_model.load_state_dict(mvml_model.state_dict())

    inner_optimizer = torch.optim.Adam(new_mvml_model.parameters(), lr=TRAIN_LR, weight_decay=REG_SCALE)
    # inner_optimizer = torch.optim.SGD(new_mvml_model.parameters(), lr=TRAIN_LR, weight_decay=REG_SCALE)
    # x_batch = Variable(torch.FloatTensor(x_batch), requires_grad=True)
    # R_matrix_batch = Variable(torch.FloatTensor(R_matrix_batch), requires_grad=True)
    # ystatus_batch = Variable(torch.FloatTensor(ystatus_batch), requires_grad=True)
    ytime_test_support = np.squeeze(ytime_test_support)
    R_matrix_batch = np.zeros([ytime_test_support.shape[0], ytime_test_support.shape[0]], dtype=int)
    for i in range(ytime_test_support.shape[0]):
        R_matrix_batch[i, :] = ytime_test_support >= ytime_test_support[i]

    x_RNASeq_batch = torch.tensor(x_RNASeq_test_support, dtype=torch.float)
    x_miRNA_batch = torch.tensor(x_miRNA_test_support, dtype=torch.float)
    R_matrix_batch = torch.tensor(R_matrix_batch, dtype=torch.float)
    ystatus_batch = torch.tensor(ystatus_test_support, dtype=torch.float)
    if torch.cuda.is_available():
        x_RNASeq_batch = x_RNASeq_batch.cuda()
        x_miRNA_batch = x_miRNA_batch.cuda()
        R_matrix_batch = R_matrix_batch.cuda()
        ystatus_batch = ystatus_batch.cuda()
    real_batch_size = x_RNASeq_batch.shape[0]

    cind_list = []
    c_index, _ = do_final_eval(new_mvml_model, x_RNASeq_test_qeury, x_miRNA_test_qeury, ytime_test_qeury,
                               ystatus_test_qeury)
    cind_list.append(c_index)

    for _ in list(range(META_STEP)):
        new_mvml_model.train()
        j = np.random.choice(range(real_batch_size), real_batch_size)
        feature_RNASeq = new_mvml_model.get_RNASeq_feature(x_RNASeq_batch)
        feature_miRNA = new_mvml_model.get_miRNA_feature(x_miRNA_batch)

        norm_gene = torch.sqrt(torch.sum(torch.mul(feature_RNASeq, feature_RNASeq), dim=1))
        norm_miRNA = torch.sqrt(torch.sum(torch.mul(feature_miRNA, feature_miRNA), dim=1))
        same_sim = torch.sum(torch.mul(feature_RNASeq, feature_miRNA), 1) / torch.mul(norm_gene, norm_miRNA)
        diff_sim = torch.sum(torch.mul(feature_RNASeq, feature_miRNA[j,]), 1) / torch.mul(norm_gene, norm_miRNA[j])
        cos_loss = torch.mean(torch.max(torch.tensor(0.0).cuda(), 0.5 - (same_sim - diff_sim)))

        theta_mv = new_mvml_model.get_survival_result(x_RNASeq_batch, x_miRNA_batch)
        exp_theta_mv = torch.reshape(torch.exp(theta_mv), [x_RNASeq_batch.shape[0]])
        theta_mv = torch.reshape(theta_mv, [x_RNASeq_batch.shape[0]])

        loss_mvml = -torch.mean(torch.mul((theta_mv - torch.log(torch.sum(torch.mul(exp_theta_mv,
                                                                                    R_matrix_batch), dim=1))),
                                          torch.reshape(ystatus_batch, [x_RNASeq_batch.shape[0]])))

        loss = loss_mvml + 1 * cos_loss

        inner_optimizer.zero_grad()
        loss.backward()
        inner_optimizer.step()

        c_index, _ = do_final_eval(new_mvml_model, x_RNASeq_test_qeury, x_miRNA_test_qeury, ytime_test_qeury,
                                   ystatus_test_qeury)
        cind_list.append(c_index)
    # return c_index
    return np.max(np.array(cind_list))


def do_final_eval(trained_model, x_RNASeq_test_query, x_miRNA_test_query, ytime_test_query, ystatus_test_query):
    trained_model.eval()
    x_RNASeq_batch = torch.FloatTensor(x_RNASeq_test_query)
    print(x_RNASeq_batch.shape)
    x_miRNA_batch = torch.FloatTensor(x_miRNA_test_query)
    if torch.cuda.is_available():
        x_RNASeq_batch = x_RNASeq_batch.cuda()
        x_miRNA_batch = x_miRNA_batch.cuda()
    pred_batch_test = trained_model.get_survival_result(x_RNASeq_batch, x_miRNA_batch)
    cind = CIndex(pred_batch_test.cpu().detach().numpy(), np.asarray(ytime_test_query), np.asarray(ystatus_test_query))

    return cind, pred_batch_test


def CIndex(pred, ytime_test, ystatus_test):
    ystatus_test = np.asarray(ystatus_test, dtype=bool)
    theta = pred

    ytime_theta = np.expand_dims(np.squeeze(ytime_test), 1)  # N*1
    ytime_theta_ex = np.expand_dims(np.squeeze(ytime_test), 0)  # 1*N
    W_ytime = ytime_theta - ytime_theta_ex  # N*N
    W_ytime = (W_ytime < 0) * np.expand_dims(np.squeeze(ystatus_test), 1)

    theta = np.expand_dims(np.squeeze(theta), 1)  # N*1
    theta_ex = np.expand_dims(np.squeeze(theta), 0)  # 1*N
    W_theta = theta - theta_ex  # N*N
    W_theta = (W_theta > 0) * np.expand_dims(np.squeeze(ystatus_test), 1)

    W = W_ytime * W_theta

    return (np.sum(W.reshape(-1, 1), 0)) / (np.sum(W_ytime.reshape(-1, 1), 0))


def meta_learn(x_RNASeq_train, x_miRNA_train, ytime_train, ystatus_train, TRAIN_ITER, TRAIN_OUTER_LR,
               TRAIN_INNER_LR, TRAIN_N_INNER, TRAIN_BATCH_N, TRAIN_SHOTS_N, x_RNASeq_target,
               x_miRNA_target, ytime_target, ystatus_target, Target_DEAD_SHOTS_N, Target_ALIVE_SHOTS_N,
               Target_TASK_NUMBER, model_parameters_path, REG_SCALE):
    ystatus_target_dead_index = np.argwhere(ystatus_target == 1)
    ystatus_target_alive_index = np.argwhere(ystatus_target == 0)

    ystatus_target_dead = ystatus_target[ystatus_target_dead_index,]
    ystatus_target_alive = ystatus_target[ystatus_target_alive_index,]
    x_RNASeq_target_dead = x_RNASeq_target[ystatus_target_dead_index,]
    x_RNASeq_target_alive = x_RNASeq_target[ystatus_target_alive_index,]
    x_miRNA_target_dead = x_miRNA_target[ystatus_target_dead_index,]
    x_miRNA_target_alive = x_miRNA_target[ystatus_target_alive_index,]
    ytime_target_dead = ytime_target[ystatus_target_dead_index,]
    ytime_target_alive = ytime_target[ystatus_target_alive_index,]
    target_test_rate = 0.2

    ystatus_train_dead_index = np.argwhere(ystatus_train == 1)
    ystatus_train_alive_index = np.argwhere(ystatus_train == 0)

    ystatus_train_dead = ystatus_train[ystatus_train_dead_index,]
    ystatus_train_alive = ystatus_train[ystatus_train_alive_index,]
    x_RNASeq_train_dead = x_RNASeq_train[ystatus_train_dead_index,]
    x_RNASeq_train_alive = x_RNASeq_train[ystatus_train_alive_index,]
    x_miRNA_train_dead = x_miRNA_train[ystatus_train_dead_index,]
    x_miRNA_train_alive = x_miRNA_train[ystatus_train_alive_index,]
    ytime_train_dead = ytime_train[ystatus_train_dead_index,]
    ytime_train_alive = ytime_train[ystatus_train_alive_index,]

    for k_num in range(20):
        torch.manual_seed(k_num)
        mvml_model = mvml_survival_Analysis()
        if torch.cuda.is_available():
            mvml_model = mvml_model.cuda()
        #  mvml_model_filepath = model_parameters_path + 'mvml_epoch' + str(0) + '_Adam' + str(TRAIN_OUTER_LR) + '.pt'
        # torch.save(mvml_model.state_dict(), mvml_model_filepath)
        optimizer = torch.optim.Adam([{'params': mvml_model.parameters()}], lr=TRAIN_OUTER_LR)

        x_RNASeq_target_train_list = []
        x_miRNA_target_train_list = []
        ystatus_target_train_list = []
        ytime_target_train_list = []

        x_RNASeq_target_test_list = []
        x_miRNA_target_test_list = []
        ystatus_target_test_list = []
        ytime_target_test_list = []

        x_RNASeq_support_list = []
        x_miRNA_support_list = []
        ystatus_support_list = []
        ytime_support_list = []

        for _ in range(Target_TASK_NUMBER):
            x_RNASeq_target_dead_index = range(x_RNASeq_target_dead.shape[0])
            ind_target_dead_test = random.sample(x_RNASeq_target_dead_index,
                                                 int(x_RNASeq_target_dead.shape[0] * target_test_rate))
            target_train_index_dead = [n for n in x_RNASeq_target_dead_index if n not in ind_target_dead_test]
            x_RNASeq_target_alive_index = range(x_RNASeq_target_alive.shape[0])
            ind_target_alive_test = random.sample(x_RNASeq_target_alive_index,
                                                  int(x_RNASeq_target_alive.shape[0] * target_test_rate))
            target_train_index_alive = [n for n in x_RNASeq_target_alive_index if n not in ind_target_alive_test]

            x_RNASeq_target_test_dead = x_RNASeq_target_dead[ind_target_dead_test,]
            x_miRNA_target_test_dead = x_miRNA_target_dead[ind_target_dead_test,]
            ystatus_target_test_dead = ystatus_target_dead[ind_target_dead_test,]
            ytime_target_test_dead = ytime_target_dead[ind_target_dead_test,]

            x_RNASeq_target_test_alive = x_RNASeq_target_alive[ind_target_alive_test,]
            x_miRNA_target_test_alive = x_miRNA_target_alive[ind_target_alive_test,]
            ystatus_target_test_alive = ystatus_target_alive[ind_target_alive_test,]
            ytime_target_test_alive = ytime_target_alive[ind_target_alive_test,]

            x_RNASeq_target_test = np.concatenate((x_RNASeq_target_test_dead, x_RNASeq_target_test_alive), axis=0)
            x_RNASeq_target_test_list.append(x_RNASeq_target_test)

            x_miRNA_target_test = np.concatenate((x_miRNA_target_test_dead, x_miRNA_target_test_alive), axis=0)
            x_miRNA_target_test_list.append(x_miRNA_target_test)

            ystatus_target_test = np.concatenate((ystatus_target_test_dead, ystatus_target_test_alive), axis=0)
            ystatus_target_test_list.append(ystatus_target_test)

            ytime_target_test = np.concatenate((ytime_target_test_dead, ytime_target_test_alive), axis=0)
            ytime_target_test_list.append(ytime_target_test)

            x_RNASeq_target_train_dead = x_RNASeq_target_dead[target_train_index_dead,]
            x_miRNA_target_train_dead = x_miRNA_target_dead[target_train_index_dead,]
            ystatus_target_train_dead = ystatus_target_dead[target_train_index_dead,]
            ytime_target_train_dead = ytime_target_dead[target_train_index_dead,]

            x_RNASeq_target_train_alive = x_RNASeq_target_alive[target_train_index_alive,]
            x_miRNA_target_train_alive = x_miRNA_target_alive[target_train_index_alive,]
            ystatus_target_train_alive = ystatus_target_alive[target_train_index_alive,]
            ytime_target_train_alive = ytime_target_alive[target_train_index_alive,]

            x_RNASeq_target_train_list.append(
                np.concatenate((x_RNASeq_target_train_dead, x_RNASeq_target_train_alive), axis=0))
            x_miRNA_target_train_list.append(
                np.concatenate((x_miRNA_target_train_dead, x_miRNA_target_train_alive), axis=0))
            ystatus_target_train_list.append(
                np.concatenate((ystatus_target_train_dead, ystatus_target_train_alive), axis=0))
            ytime_target_train_list.append(np.concatenate((ytime_target_train_dead, ytime_target_train_alive), axis=0))

            x_RNASeq_target_train_dead_index = range(x_RNASeq_target_train_dead.shape[0])
            target_ind_dead_shots = random.sample(x_RNASeq_target_train_dead_index, Target_DEAD_SHOTS_N)
            target_rest_train_index_dead = [n for n in x_RNASeq_target_train_dead_index if
                                            n not in target_ind_dead_shots]
            x_RNASeq_target_train_alive_index = range(x_RNASeq_target_train_alive.shape[0])
            target_ind_alive_shots = random.sample(x_RNASeq_target_train_alive_index, Target_ALIVE_SHOTS_N)
            target_rest_train_index_alive = [n for n in x_RNASeq_target_train_alive_index if
                                             n not in target_ind_alive_shots]

            # ind = np.random.choice(range(x_train.shape[0]), shots_n, p=(1/weight_train)/np.sum(1/weight_train))
            # print(ind)

            x_RNASeq_target_batch_dead = x_RNASeq_target_train_dead[target_ind_dead_shots,]
            x_miRNA_target_batch_dead = x_miRNA_target_train_dead[target_ind_dead_shots,]
            ystatus_target_batch_dead = ystatus_target_train_dead[target_ind_dead_shots,]
            ytime_target_batch_dead = ytime_target_train_dead[target_ind_dead_shots,]

            x_RNASeq_target_batch_alive = x_RNASeq_target_train_alive[target_ind_alive_shots,]
            x_miRNA_target_batch_alive = x_miRNA_target_train_alive[target_ind_alive_shots,]
            ystatus_target_batch_alive = ystatus_target_train_alive[target_ind_alive_shots,]
            ytime_target_batch_alive = ytime_target_train_alive[target_ind_alive_shots,]

            x_RNASeq_support = np.concatenate((x_RNASeq_target_batch_dead, x_RNASeq_target_batch_alive), axis=0)
            x_RNASeq_support_list.append(x_RNASeq_support)

            x_miRNA_support = np.concatenate((x_miRNA_target_batch_dead, x_miRNA_target_batch_alive), axis=0)
            x_miRNA_support_list.append(x_miRNA_support)

            ystatus_support = np.concatenate((ystatus_target_batch_dead, ystatus_target_batch_alive), axis=0)
            ystatus_support_list.append(ystatus_support)

            ytime_support = np.concatenate((ytime_target_batch_dead, ytime_target_batch_alive), axis=0)
            ytime_support_list.append(ytime_support)

        cind_list_iter_max = []

        for t in range(TRAIN_ITER):

            start = time.time()

            ind_train_dead_shots = np.random.choice(range(x_RNASeq_train_dead.shape[0]), int(TRAIN_SHOTS_N * 0.25))
            ind_train_alive_shots = np.random.choice(range(x_RNASeq_train_alive.shape[0]), int(TRAIN_SHOTS_N * 0.75))
            # print(ind)

            x_RNASeq_batch_train_dead = x_RNASeq_train_dead[ind_train_dead_shots,]
            x_miRNA_batch_train_dead = x_miRNA_train_dead[ind_train_dead_shots,]
            ystatus_batch_train_dead = ystatus_train_dead[ind_train_dead_shots,]
            ytime_batch_train_dead = ytime_train_dead[ind_train_dead_shots,]

            x_RNASeq_batch_train_alive = x_RNASeq_train_alive[ind_train_alive_shots,]
            x_miRNA_batch_train_alive = x_miRNA_train_alive[ind_train_alive_shots,]
            ystatus_batch_train_alive = ystatus_train_alive[ind_train_alive_shots,]
            ytime_batch_train_alive = ytime_train_alive[ind_train_alive_shots,]

            x_RNASeq_batch_train = np.concatenate((x_RNASeq_batch_train_dead, x_RNASeq_batch_train_alive), axis=0)
            x_miRNA_batch_train = np.concatenate((x_miRNA_batch_train_dead, x_miRNA_batch_train_alive), axis=0)
            ystatus_batch_train = np.concatenate((ystatus_batch_train_dead, ystatus_batch_train_alive), axis=0)
            ytime_batch_train = np.concatenate((ytime_batch_train_dead, ytime_batch_train_alive), axis=0)

            x_RNASeq_batch = np.squeeze(x_RNASeq_batch_train)
            x_miRNA_batch = np.squeeze(x_miRNA_batch_train)
            ystatus_batch = np.squeeze(ystatus_batch_train)
            ytime_batch = np.squeeze(ytime_batch_train)
            R_matrix_batch = np.zeros([ytime_batch.shape[0], ytime_batch.shape[0]], dtype=int)
            for i in range(ytime_batch.shape[0]):
                R_matrix_batch[i, :] = ytime_batch >= ytime_batch[i]
            new_mvml_model = do_base_learning(mvml_model, x_RNASeq_batch, x_miRNA_batch, R_matrix_batch, ystatus_batch,
                                              TRAIN_INNER_LR, TRAIN_N_INNER, REG_SCALE)

            diff = list()
            for p, new_p in zip(mvml_model.parameters(), new_mvml_model.parameters()):
                temp = Variable(torch.zeros(p.size()))
                temp = temp.cuda()
                temp.add_(p.data - new_p.data)
                diff.append(temp)

            for _ in range(TRAIN_BATCH_N - 1):

                ind_train_dead_shots = np.random.choice(range(x_RNASeq_train_dead.shape[0]), int(TRAIN_SHOTS_N * 0.25))
                ind_train_alive_shots = np.random.choice(range(x_RNASeq_train_alive.shape[0]),
                                                         int(TRAIN_SHOTS_N * 0.75))
                # print(ind)

                x_RNASeq_batch_train_dead = x_RNASeq_train_dead[ind_train_dead_shots,]
                x_miRNA_batch_train_dead = x_miRNA_train_dead[ind_train_dead_shots,]
                ystatus_batch_train_dead = ystatus_train_dead[ind_train_dead_shots,]
                ytime_batch_train_dead = ytime_train_dead[ind_train_dead_shots,]

                x_RNASeq_batch_train_alive = x_RNASeq_train_alive[ind_train_alive_shots,]
                x_miRNA_batch_train_alive = x_miRNA_train_alive[ind_train_alive_shots,]
                ystatus_batch_train_alive = ystatus_train_alive[ind_train_alive_shots,]
                ytime_batch_train_alive = ytime_train_alive[ind_train_alive_shots,]

                x_RNASeq_batch_train = np.concatenate((x_RNASeq_batch_train_dead, x_RNASeq_batch_train_alive), axis=0)
                x_miRNA_batch_train = np.concatenate((x_miRNA_batch_train_dead, x_miRNA_batch_train_alive), axis=0)
                ystatus_batch_train = np.concatenate((ystatus_batch_train_dead, ystatus_batch_train_alive), axis=0)
                ytime_batch_train = np.concatenate((ytime_batch_train_dead, ytime_batch_train_alive), axis=0)

                x_RNASeq_batch = np.squeeze(x_RNASeq_batch_train)
                x_miRNA_batch = np.squeeze(x_miRNA_batch_train)
                ystatus_batch = np.squeeze(ystatus_batch_train)
                ytime_batch = np.squeeze(ytime_batch_train)
                R_matrix_batch = np.zeros([ytime_batch.shape[0], ytime_batch.shape[0]], dtype=int)
                for i in range(ytime_batch.shape[0]):
                    R_matrix_batch[i, :] = ytime_batch >= ytime_batch[i]
                new_mvml_model = do_base_learning(mvml_model, x_RNASeq_batch, x_miRNA_batch, R_matrix_batch,
                                                  ystatus_batch,
                                                  TRAIN_INNER_LR, TRAIN_N_INNER, REG_SCALE)

                diff_next = list()
                for p, new_p in zip(mvml_model.parameters(), new_mvml_model.parameters()):
                    temp = Variable(torch.zeros(p.size()))
                    temp = temp.cuda()
                    temp.add_(p.data - new_p.data)
                    diff_next.append(temp)

                diff = list(map(add, diff, diff_next))

            diff_ave = [x / TRAIN_BATCH_N for x in diff]

            ind_k = 0
            for p in mvml_model.parameters():
                if p.grad is None:
                    p.grad = Variable(torch.zeros(p.size())).cuda()
                p.grad.data.add_(diff_ave[ind_k])
                ind_k = ind_k + 1
            optimizer.step()
            optimizer.zero_grad()
            end = time.time()
            with open('nohup_mvml.log', 'a') as f:
                f.writelines("第" + str(t) + "步训练开始" + '\n')

            cind_list_max = []
            for i in range(Target_TASK_NUMBER):
                CI = do_final_learning(mvml_model, x_RNASeq_support_list[i], x_miRNA_support_list[i],
                                       ytime_support_list[i],
                                       ystatus_support_list[i], x_RNASeq_target_test_list[i],
                                       x_miRNA_target_test_list[i],
                                       ytime_target_test_list[i], ystatus_target_test_list[i], TRAIN_OUTER_LR,
                                       REG_SCALE,
                                       META_STEP=5)
                cind_list_max.append(CI)

                with open('Cindex_' + str(k_num) + '.log', 'a') as f:
                    f.writelines(str(CI) + '\n')
            with open('Cindex_' + str(k_num) + '.log', 'a') as f:
                f.writelines("Meta train" + str(t) + "次迭代, Ave C_index:" + str(np.mean(np.array(cind_list_max))) + '\n')

            cind_list_iter_max.append(np.mean(np.array(cind_list_max)))

        with open('Cindex' + '.log', 'a') as f:
            f.writelines(str(np.max(np.array(cind_list_iter_max))) + '\n')


parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, default='config.json',
                    help='configuration json file')

if __name__ == '__main__':

    args = parser.parse_args()
    with open(args.config) as f:
        config = json.load(f)

    TRAIN_PATH = config['train_path']
    TRAIN_OUTER_LR = config['train_outer_lr']
    TRAIN_INNER_LR = config['train_inner_lr']
    TRAIN_SHOTS_N = config['train_shots_n']
    TRAIN_BATCH_N = config['train_batch_n']
    TRAIN_N_INNER = config['train_n_inner']
    REG_SCALE = config['reg_scale']
    TRAIN_ITER = config['train_iters']

    Target_PATH = config['target_path']
    Target_DEAD_SHOTS_N = config['target_dead_shots_n']
    Target_ALIVE_SHOTS_N = config['target_alive_shots_n']
    Target_TASK_NUMBER = config['target_task_number']

    model_parameters_path = config['model_parameter_path']

    # x_RNASeq_train = np.loadtxt(fname=config['train_RNASeq_feature'], delimiter=",", skiprows=1)
    # x_miRNA_train = np.loadtxt(fname=config['train_miRNA_feature'], delimiter=",", skiprows=1)
    # y_train = np.loadtxt(fname=config['train_time'], delimiter=",", skiprows=1)
    # ystatus_train = np.loadtxt(fname=config['train_status'], delimiter=",", skiprows=1)

    dir_res = os.listdir(TRAIN_PATH)
    Flag = True
    for i in range(len(dir_res)):
        if Flag:
            temp_train_path = TRAIN_PATH + '/' + dir_res[i]
            RNASeq_arr = np.loadtxt(fname=temp_train_path + '/RNASeq.csv', delimiter=",", skiprows=1)
            miRNA_arr = np.loadtxt(fname=temp_train_path + '/miRNA.csv', delimiter=",", skiprows=1)
            ytime_arr = np.loadtxt(fname=temp_train_path + '/ytime.csv', delimiter=",", skiprows=1)
            ystatus_arr = np.loadtxt(fname=temp_train_path + '/ystatus.csv', delimiter=",", skiprows=1)
            Flag = False
        else:
            temp_train_path = TRAIN_PATH + '/' + dir_res[i]
            temp_RNASeq_arr = np.loadtxt(fname=temp_train_path + '/RNASeq.csv', delimiter=",", skiprows=1)
            RNASeq_arr = np.concatenate((RNASeq_arr, temp_RNASeq_arr), axis=0)
            temp_miRNA_arr = np.loadtxt(fname=temp_train_path + '/miRNA.csv', delimiter=",", skiprows=1)
            miRNA_arr = np.concatenate((miRNA_arr, temp_miRNA_arr), axis=0)
            temp_ytime_arr = np.loadtxt(fname=temp_train_path + '/ytime.csv', delimiter=",", skiprows=1)
            ytime_arr = np.concatenate((ytime_arr, temp_ytime_arr), axis=0)
            temp_ystatus_arr = np.loadtxt(fname=temp_train_path + '/ystatus.csv', delimiter=",", skiprows=1)
            ystatus_arr = np.concatenate((ystatus_arr, temp_ystatus_arr), axis=0)

    x_RNASeq_train = RNASeq_arr
    x_miRNA_train = miRNA_arr
    ytime_train = ytime_arr
    ystatus_train = ystatus_arr

    x_RNASeq_target = np.loadtxt(fname=Target_PATH + "/RNASeq.csv", delimiter=",", skiprows=1)
    x_miRNA_target = np.loadtxt(fname=Target_PATH + "/miRNA.csv", delimiter=",", skiprows=1)
    ytime_target = np.loadtxt(fname=Target_PATH + "/ytime.csv", delimiter=",", skiprows=1)
    ystatus_target = np.loadtxt(fname=Target_PATH + "/ystatus.csv", delimiter=",", skiprows=1)

    RNASeq_feature = np.concatenate((x_RNASeq_train, x_RNASeq_target), axis=0)
    miRNA_feature = np.concatenate((x_miRNA_train, x_miRNA_target), axis=0)

    from sklearn import preprocessing

    standard_scaler = preprocessing.StandardScaler()
    RNASeq_feature = standard_scaler.fit_transform(RNASeq_feature)
    miRNA_feature = standard_scaler.fit_transform(miRNA_feature)
    x_RNASeq_train = RNASeq_feature[0:x_RNASeq_train.shape[0], :]
    x_RNASeq_target = RNASeq_feature[x_RNASeq_train.shape[0]:, :]
    x_miRNA_train = miRNA_feature[0:x_miRNA_train.shape[0], :]
    x_miRNA_target = miRNA_feature[x_miRNA_train.shape[0]:, ]

    random.seed(1)
    np.random.seed(1)
    meta_learn(x_RNASeq_train=x_RNASeq_train, x_miRNA_train=x_miRNA_train,
               ytime_train=ytime_train, ystatus_train=ystatus_train, TRAIN_ITER=TRAIN_ITER,
               TRAIN_OUTER_LR=TRAIN_OUTER_LR, TRAIN_INNER_LR=TRAIN_INNER_LR,
               TRAIN_N_INNER=TRAIN_N_INNER, TRAIN_BATCH_N=TRAIN_BATCH_N, TRAIN_SHOTS_N=TRAIN_SHOTS_N,
               x_RNASeq_target=x_RNASeq_target, x_miRNA_target=x_miRNA_target, ytime_target=ytime_target,
               ystatus_target=ystatus_target, Target_DEAD_SHOTS_N=Target_DEAD_SHOTS_N,
               Target_ALIVE_SHOTS_N=Target_ALIVE_SHOTS_N, Target_TASK_NUMBER=Target_TASK_NUMBER,
               model_parameters_path=model_parameters_path, REG_SCALE=REG_SCALE)

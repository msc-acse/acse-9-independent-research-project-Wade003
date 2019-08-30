from numpy import concatenate, zeros
import torch
from torch import nn
import numpy as np
from torch.autograd import Variable

LR = 0.02           # learning rate
class LSNN(nn.Module):
    def __init__(self, history_level, nodes_number,hidden_size = 30, hidden_layer = 1):
        super(LSNN, self).__init__()
        self.lstm = nn.LSTM(
            input_size = history_level,
            hidden_size = hidden_size,
            num_layers = hidden_layer,
            batch_first = True,

        )
        self.nodes_number = nodes_number
        self.hidden = (torch.autograd.Variable(torch.zeros(hidden_layer, 1, hidden_size)),torch.autograd.Variable(torch.zeros(hidden_layer, 1, hidden_size)))
        self.out = nn.Linear(hidden_size, 1)

    def forward(self,x):
        # x (batch, time_step, input_size)
        # h_state (n_layers, batch, hidden_size)
        # r_out (batch, time_step, output_size)
        r_out,self.hidden= self.lstm(x,self.hidden)
        self.hidden=(Variable(self.hidden[0]),Variable(self.hidden[1]))
        outs = []

        for node in range(self.nodes_number):
            outs.append(self.out(r_out[:, node, :]))
        return torch.stack(outs, dim=1)


def LSTM_train(dd_pod_coeffs_all, history_level, local_list, order):
    torch.manual_seed(1)    # reproducible
    print("DD_LSTM_training")
    lstmNN = []
    lstmNN_temp = []
    optimizer = []
    time_steps = dd_pod_coeffs_all.shape[1]
    nodes_number = dd_pod_coeffs_all.shape[2]
    for sub in range(len(local_list)):
        lstmNN.append(LSNN(history_level, nodes_number, hidden_size = nodes_number*len(local_list[sub])))
        optimizer.append(torch.optim.Adam(lstmNN[sub].parameters(), lr=LR))
    loss_func = nn.MSELoss()
    #loss_1 = []
    n_sub = len(local_list)
    for step in range(time_steps - history_level):
        print("percentage", float(step)/(float(time_steps) - float(history_level)))
        #prediction_list = []
        for sub in order:
            x_np = []
            y_np = [[]]
            x = []
            y = []
            for i in range(history_level):
                null_list = []
                for local_input in local_list[sub]:
                    null_list.append(dd_pod_coeffs_all[local_input][step + i])
                # null_list.append(dd_pod_coeffs_all[local_list[sub][0]][step + i])
                x_np.append(null_list)

            x = Variable(torch.from_numpy(np.array(x_np)).float())
            for local_input in local_list[sub]:
                y_np[0].append(dd_pod_coeffs_all[local_input][step + history_level])
            y = Variable(torch.from_numpy(np.array(y_np)).float())
            x = x.permute(1,2,0).view(-1, nodes_number*len(local_list[sub]), history_level)
            y = y.permute(1,2,0).view(-1, nodes_number, 1)
            prediction = lstmNN[sub](x)
            loss = loss_func(prediction, y)
            optimizer[sub].zero_grad()
            loss.backward()
            optimizer[sub].step()

        #
        # lstmNN_temp = lstmNN
        # prediction_update = []
        # if(step == 0):
        #     for sub in range(len(local_list)):
        #         prediction_update.append([])
        #         prediction_update[sub].append(dd_pod_coeffs_all[sub][step + history_level - 1])
        # else:
        #     for sub in range(len(local_list)):
        #         prediction_update.append([])
        #         prediction_update[sub].append(prediction_list[sub])
        #
        # #print("start convergence test")
        # for loc_iter in range(10):
        #     for sub in range(n_sub):
        #         x_temp = []
        #         for i in range(history_level):
        #             null_list = []
        #             if(loc_iter < history_level - 1):
        #                 if(i < history_level - 1):    # first time 0,1, second time 0
        #                     for local_input in local_list[sub]:
        #                         null_list.append(dd_pod_coeffs_all[local_input][step + i])
        #                 else:                                # first time 2, second time 1,2, third time 0,1,2 and all the same after
        #                     for local_input in local_list[sub]:
        #                         null_list.append(prediction_update[local_input][-(history_level - i)])   # -3, -2, -1
        #             else:
        #                 for local_input in local_list[sub]:
        #                     null_list.append(prediction_update[local_input][-(history_level - i)])
        #             # for local_input in local_list[sub]:
        #             #     # null_list.append(dd_pod_coeffs_all[local_input][step + i])
        #             #     null_list.append(dd_pod_coeffs_all[local_input][step + i])
        #             x_temp.append(null_list)
        #         x = Variable(torch.from_numpy(np.array(x_temp)).float())
    	#         y_temp = []
        #         y_temp.append(dd_pod_coeffs_all[sub][step + history_level])
        #         y = Variable(torch.from_numpy(np.array(y_temp)).float())
        #         #print("x", x.size(), nodes_number*len(local_list[sub]))
        #         x = x.permute(1,2,0).view(-1, nodes_number*len(local_list[sub]), history_level)
        #         y = y.permute(-1,0).view(-1, nodes_number, 1)
        #
        #         prediction = lstmNN[sub](x)
        #         prediction_update[sub].append(prediction.data.view(nodes_number).numpy())
        #         #print("sub=", sub, np.max(prediction_update[sub][-1]-prediction_update[sub][-2]))
        #         loss = loss_func(prediction, y)
        #         optimizer[sub].zero_grad()
        #         loss.backward()
        #         optimizer[sub].step()
        # lstmNN = lstmNN_temp
        # prediction_list = []
        # for sub in range(n_sub):
        #     x_np = []
        #     y_np = []
        #     x = []
        #     y = []
        #     for i in range(history_level):
        #         null_list = []
        #         for local_input in local_list[sub]:
        #             null_list.append(prediction_update[local_input][-history_level+i])
        #             #null_list.append(dd_pod_coeffs_all[local_input][step + i])
        #         # null_list.append(dd_pod_coeffs_all[local_list[sub][0]][step + i])
        #         x_np.append(null_list)
        #     x = Variable(torch.from_numpy(np.array(x_np)).float())
        #     y_np.append(dd_pod_coeffs_all[sub][step + history_level])
        #     y = Variable(torch.from_numpy(np.array(y_np)).float())
        #     x = x.permute(1,2,0).view(-1, nodes_number*len(local_list[sub]), history_level)
        #     # x = x.permute(1,2,0).view(-1, nodes_number, history_level)
        #     y = y.permute(-1,0).view(-1, nodes_number, 1)
        #
        #     # x = x.permute(1,0).view(1, nodes_number*len(local_list[sub]), history_level)
        #     # y = y.permute(1,0).view(1, nodes_number, 1)
        #     prediction = lstmNN[sub](x)
        #     prediction_list.append(prediction.data.view(nodes_number).numpy())
        #     loss = loss_func(prediction, y)
        #     optimizer[sub].zero_grad()
        #     loss.backward()
        #     optimizer[sub].step()
    return lstmNN

def LSTM_predict(dd_pod_coeffs_all, history_level, lstmNN, time_steps, local_list, order):
    print("DD_LSTM_predicting")
    nodes_number = dd_pod_coeffs_all.shape[2]
    prediction_list = []
    for sub in range(len(lstmNN)):
        prediction_list.append([])
    for step in range(time_steps - history_level):
        for sub in order:
            x_np = []
            y_np = []
            x = []
            y = []
            if(step < history_level):
                for i in range(history_level - step):
                    null_list = []
                    for local_sub in local_list[sub]:
                        null_list.append(dd_pod_coeffs_all[local_sub][step + i])
                    x_np.append(null_list)
                for i in range(step + 1):
                    if i > 0:
                        null_list = []
                        for local_sub in local_list[sub]:
                            null_list.append(prediction_list[local_sub][-i])
                        x_np.append(null_list)
            else:
                for i in range(history_level):
                    null_list = []
                    for local_sub in local_list[sub]:
                        null_list.append(prediction_list[local_sub][-history_level+i])
                    x_np.append(null_list)
            x = Variable(torch.from_numpy(np.array(x_np)).float())
            #print("x", x.size())
            x = x.permute(1,2,0).view(1, nodes_number*len(local_list[sub]), history_level)
            with torch.no_grad():
                prediction = lstmNN[sub](x)
            prediction_list[sub].append(prediction.data.view(nodes_number).numpy()[:nodes_number])
            #print("length", len(prediction_list))
    return np.array(prediction_list)

def LSTM_DD_global_train(dd_pod_coeffs_all, history_level, local_list, order):
    torch.manual_seed(1)    # reproducible
    print("DD_LSTM_training")
    optimizer = []
    time_steps = dd_pod_coeffs_all.shape[1]
    nodes_number = dd_pod_coeffs_all.shape[2]
    n_sub = len(local_list)
    lstmNN = LSNN(history_level, nodes_number * n_sub, hidden_size = nodes_number)
    optimizer=torch.optim.Adam(lstmNN.parameters(), lr=LR)
    loss_func = nn.MSELoss()
    n_sub = len(local_list)

    for step in range(time_steps - history_level):
        x_np = []
        y_np = []
        x = []
        y = []
        for i in range(history_level):
            x_np.append(dd_pod_coeffs_all[:, step + i])
        x = Variable(torch.from_numpy(np.array(x_np)).float())
        y_np.append(dd_pod_coeffs_all[:, step + history_level])
        y = Variable(torch.from_numpy(np.array(y_np)).float())
        x = x.permute(1,2,0).view(1, nodes_number * n_sub, history_level)
        y = y.permute(1,2,0).view(1, nodes_number * n_sub, 1)
        prediction = lstmNN(x)
        #prediction_list.append(prediction.data.view(nodes_number).numpy())
        loss = loss_func(prediction, y)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    return lstmNN

def LSTM_DD_global_predict(dd_pod_coeffs_all, history_level, lstmNN, time_steps, local_list, order):
    print("LSTM_predicting")
    nodes_number = dd_pod_coeffs_all.shape[2]   # for one sub domian
    prediction_list = []
    n_sub = dd_pod_coeffs_all.shape[0]
    for step in range(time_steps - history_level):
        x_np = []
        y_np = []
        x = []
        y = []
        if(step < history_level):
            for i in range(history_level - step):
                x_np.append(dd_pod_coeffs_all[:, step + i])
            for i in range(step + 1):
                if i > 0:
                    x_np.append(prediction_list[-i])
        elif(step >= history_level):
            for i in range(history_level):
                x_np.append(prediction_list[-history_level+i])
        x = Variable(torch.from_numpy(np.array(x_np)).float())
        x = x.permute(1,2,0).view(1, nodes_number * n_sub, history_level)
        with torch.no_grad():
            prediction = lstmNN(x)
        prediction_list.append(prediction.data.view(n_sub, nodes_number).numpy())
        #print("length", len(prediction_list))
    return prediction_list
